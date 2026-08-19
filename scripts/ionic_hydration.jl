# =============================================================================
#  ionic_hydration.jl — hydration of a Portland cement through its PORE SOLUTION
#
#  The clinker only dissolves — into Ca²⁺, SiO₂, AlO₂⁻, FeO₂⁻, SO₄²⁻, CO₃²⁻ and
#  H⁺ — and a Gibbs free-energy minimization decides, at every accepted step of
#  the integration, which hydrates are stable and in what amounts. No sequencing
#  rule is written anywhere: the aluminate cascade, the fate of the ettringite,
#  the pH of the pore solution are all results.
#
#  Contrast with `opc_semiadiabatic_calorimetry.jl`, which runs the same cement
#  through aggregated solid → solid reactions: there each reaction STATES its
#  products, and which one forms when has to be decided by hand.
#
#  Chemistry only. The same model coupled to a four-scale micromechanical
#  estimate of E(t) is the Applications chapter "Hydration through the pore
#  solution" of MeanFieldHomogenization.jl, which duplicates the chemistry below
#  and adds the mechanics.
#
#  Usage:
#    julia --project=scripts scripts/ionic_hydration.jl
#    or from the REPL: include("scripts/ionic_hydration.jl")
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections
using OrdinaryDiffEq
using Printf

const IONIC_CEMDATA = joinpath(pkgdir(ChemistryLab), "data", "cemdata18-thermofun.json")

# Four systems. `:opc` is the model; the three smaller ones are the ladder that
# was used to isolate what made it hard, and they are kept because reproducing an
# intermediate is the fastest way to localize a regression.
const IONIC_SYSTEMS = Dict(
    :silicates => (
        anhydrous = ["C3S", "C2S"],
        hydrates = ["Portlandite", "Jennite"],
    ),
    :aluminate => (
        anhydrous = ["C3S", "C2S", "C3A"],
        hydrates = ["Portlandite", "Jennite", "C3AH6"],
    ),
    :sulfate => (
        anhydrous = ["C3S", "C2S", "C3A", "Gp"],
        hydrates = ["Portlandite", "Jennite", "ettringite", "monosulphate12", "C3AH6"],
    ),
    :opc => (
        anhydrous = ["C3S", "C2S", "C3A", "C4AF", "Gp", "Cal"],
        hydrates = [
            "Portlandite", "Jennite", "ettringite", "monosulphate12",
            "monocarbonate", "C3AH6", "FeOOHmic",
        ],
    ),
)

"""
    IONIC_DEFAULT_SYSTEM

`:opc` — the full ordinary Portland cement: alite, belite, aluminate, ferrite,
gypsum and calcite dissolving to ions, with the assemblage left to the
thermodynamics.

Measured over 28 days: 202 accepted steps, `retcode = Success`, pore solution
holding at pH 12.58. The aluminate sequence comes out of the free-energy
minimization with no sequencing rule written anywhere — ettringite forms early,
peaks around 6 hours and converts to monosulphate once the sulfate is spent, the
AFm settling at exactly the sulfate budget.

!!! note "The reported compositions are certified"
    [`speciated_states`](@ref) passes each instant to `DualEquilibriumSolver`,
    which solves the KKT system and returns a proof: the Gibbs problem is convex,
    so stationarity of the interior species, the component balance, and
    undersaturation of every absent phase are together sufficient for **global**
    optimality. On this OPC every replayed instant is certified, with an element
    balance between 1e-11 and 1e-13 mol.

    That matters because the interior-point solver alone is not reliable here. On
    this package's own calcite reference it returns pH 6.96 where the certified
    answer is 9.90, and it misses a trace species by 147 %; it also rarely
    reports convergence at all, so its return code cannot be used to tell the two
    cases apart.
"""
const IONIC_DEFAULT_SYSTEM = :opc

# Phase families, so an assemblage can be reported by role rather than by
# CEMDATA18 symbol.
const IONIC_GROUPS = [
    "anhydrous" => ["C3S", "C2S", "C3A", "C4AF"],
    "gypsum" => "Gp",
    "calcite" => "Cal",
    "C-S-H" => "Jennite",
    "CH" => "Portlandite",
    "AFt" => "ettringite",
    "AFm" => ["monosulphate12", "monocarbonate"],
    "hydrogarnet" => "C3AH6",
    "FH3" => "FeOOHmic",
    # "water" is filled in per system by `ionic_water_group`: it must collect
    # EVERY aqueous species, not just H2O@. The liquid phase is the pore
    # solution, and at cement ionic strength the solutes are not negligible —
    # their negative partial molar volumes sum to -3e-4 of the paste volume,
    # which left as an unassigned remainder trips the volume bookkeeping.
]

"""
    ionic_water_group(cs) -> Pair

The `"water"` phase: every aqueous species of `cs`, i.e. the whole pore solution
rather than the solvent alone.
"""
ionic_water_group(cs) =
    "water" => [symbol(s) for s in cs.species if aggregate_state(s) == AS_AQUEOUS]

"""
    GEL_WATER_PER_CSH

Water held in the C-S-H gel, 1.9 mol per mole of silicon: the difference between
the C₁.₇SH₄ of a saturated gel and the 2.1 mol the `Jennite` end-member carries in
its own formula. Counting it with the C-S-H rather than with the pore solution is
the convention of Lavergne et al. (2018); pass `gel_water = 0` to
[`ionic_phase_history`](@ref) to leave it in the liquid.
"""
const GEL_WATER_PER_CSH = 4.0 - 2.1

"""
    build_ionic_system(system = IONIC_DEFAULT_SYSTEM) -> ChemicalSystem

The CEMDATA18 subset for the ionic model: the dissolving phases, the candidate
hydrates, and every aqueous species `speciation` pulls in from
`CEMDATA_PRIMARIES` — 49 of them, which is what carries the pore chemistry.
"""
function build_ionic_system(system::Symbol = IONIC_DEFAULT_SYSTEM)
    haskey(IONIC_SYSTEMS, system) || throw(
        ArgumentError(
            "system must be one of $(sort(collect(keys(IONIC_SYSTEMS)))), got :$system",
        ),
    )
    spec = IONIC_SYSTEMS[system]
    subs = build_species(IONIC_CEMDATA)
    sp = speciation(subs, vcat(spec.anhydrous, spec.hydrates); aggregate_state = [AS_AQUEOUS])
    return ChemicalSystem(sp, CEMDATA_PRIMARIES)
end

"""
    IONIC_CALIBRATION

Per-phase multipliers on the dissolution rates.

They are all 1 here, and the entry exists so that a study which needs the degrees
of hydration to follow a reference curve can scale them without touching the rate
laws. A Palandri–Kharaka parameter set would be the principled alternative, but
none is published for the clinker phases and inventing one would be a
fabrication.
"""
const IONIC_CALIBRATION = Dict(
    "C3S" => 1.0, "C2S" => 1.0, "C3A" => 1.0, "C4AF" => 1.0,
    "Gp" => 1.0, "Cal" => 1.0,
)

"""
    ionic_reactions(cs; wb, blaine, system, calibration) -> Vector{KineticReaction}

Congruent dissolution of each anhydrous phase into the primary aqueous species,
with a Parrot–Killoh rate scaled by a per-phase calibration factor.

The reactions are **balanced by ChemistryLab** from the phase and the primaries,
and come out in the expected acid-driven form, e.g.

    C₃S = 3 H₂O + 3 Ca²⁺ + SiO₂ − 6 H⁺

The negative proton coefficient is what drives the pore solution alkaline.
"""
function ionic_reactions(
        cs; wb, blaine,
        system::Symbol = IONIC_DEFAULT_SYSTEM,
        calibration = IONIC_CALIBRATION,
    )
    nmv = symbol.(cs.species)
    prim = [
        p for p in ("Ca+2", "SiO2@", "AlO2-", "FeO2-", "SO4-2", "CO3-2", "H2O@", "H+")
            if p in nmv
    ]
    α_max = powers_alpha_max(wb)

    pk_of = Dict(
        "C3S" => PK84_PARAMS_C3S, "C2S" => PK84_PARAMS_C2S,
        "C3A" => PK84_PARAMS_C3A, "C4AF" => PK84_PARAMS_C4AF,
    )

    out = KineticReaction[]
    for a in IONIC_SYSTEMS[system].anhydrous
        rxn = Reaction([cs[a]], [cs[p] for p in prim]; symbol = "$a dissolution")

        base = if haskey(pk_of, a)
            parrot_killoh_avrami(pk_of[a], a; α_max, blaine)
        else
            # Gypsum and calcite are not clinker: give them a fast first-order
            # release so sulfate and carbonate are available to the minimization
            # from the start, rather than rate-limiting it.
            k = get(calibration, a, 1.0) * 1.0e-4
            KineticFunc(
                (T, P, t, n, lna, n0) -> k * max(n[a], zero(eltype(n.data))),
                (T = 293.15u"K", P = 1.0e5u"Pa"), u"mol/s",
            )
        end

        f = get(calibration, a, 1.0)
        rate = haskey(pk_of, a) ?
            KineticFunc(
            (T, P, t, n, lna, n0) -> f * base(T, P, t, n, lna, n0),
            (T = 293.15u"K", P = 1.0e5u"Pa"), u"mol/s",
        ) : base

        rxn[:rate] = rate
        push!(out, KineticReaction(cs, rxn))
    end
    return out
end

"""
    run_ionic_hydration(; wb, clinker, gypsum, filler, blaine, tend,
                          binder_mass, calorimeter, system) -> NamedTuple

Integrate the coupled problem: dissolution kinetics on the anhydrous phases,
Gibbs minimization on everything else, once per accepted step.

`binder_mass` is the mass of binder simulated (1 kg by default, so every
extensive result is per kilogram of binder). `calorimeter`, when given, couples
the thermal balance into the same ODE.

Returns `(; cs, state0, kp, system, calorimeter, sol)`. The activity model is
`HKFActivityModel` on both halves — a cement pore solution sits at
I ≈ 0.1–0.7 mol/kg, where a dilute model is not defensible.
"""
function run_ionic_hydration(;
        wb = 0.5,
        clinker = (C3S = 0.65, C2S = 0.11, C3A = 0.11, C4AF = 0.08),
        gypsum = 0.046, filler = 0.035,
        blaine = 380.0u"m^2/kg", tend = 28 * 86400.0,
        reltol = 1.0e-7, abstol = 1.0e-10,
        system::Symbol = IONIC_DEFAULT_SYSTEM,
        binder_mass = 1.0u"kg",
        calorimeter = nothing,
    )
    cs = build_ionic_system(system)
    nmv = symbol.(cs.species)
    f_clinker = 1.0 - gypsum - filler

    mb = ustrip(us"kg", binder_mass)
    T0 = calorimeter === nothing ? 293.15u"K" : _calorimeter_T0(calorimeter)

    state0 = ChemicalState(cs; T = T0)
    for (nm, w) in pairs(clinker)
        string(nm) in nmv && set_quantity!(state0, string(nm), (mb * f_clinker * w)u"kg")
    end
    gypsum > 0 && "Gp" in nmv && set_quantity!(state0, "Gp", (mb * gypsum)u"kg")
    filler > 0 && "Cal" in nmv && set_quantity!(state0, "Cal", (mb * filler)u"kg")
    set_quantity!(state0, "H2O@", (mb * wb)u"kg")

    model = HKFActivityModel()
    krs = ionic_reactions(cs; wb, blaine, system)
    kp = if calorimeter === nothing
        KineticsProblem(
            cs, krs, state0, (0.0, tend);
            activity_model = model,
            equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()),
        )
    else
        KineticsProblem(
            cs, krs, state0, (0.0, tend);
            activity_model = model,
            equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()),
            calorimeter = calorimeter,
        )
    end
    ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = reltol, abstol = abstol)
    return (; cs, state0, kp, system, calorimeter, sol = integrate(kp, ks))
end

_calorimeter_T0(cal::IsothermalCalorimeter) = cal.T
_calorimeter_T0(cal::SemiAdiabaticCalorimeter) = cal.T0

# ── calorimetry, after Lavergne et al. (2018) §4.1 ────────────────────────────

"""
    LAVERGNE_MIX_C100

Mix proportions of the plain-cement semi-adiabatic test of Lavergne et al.
(2018), Table 11, at w/b = 0.5: 371 g of binder, 1113 g of dry sand, 196 g of
water. The sand is there to keep the temperature rise moderate, as NF EN 196-9
prescribes; it takes no part in the chemistry and enters only through its heat
capacity.
"""
const LAVERGNE_MIX_C100 = (binder = 0.371u"kg", sand = 1.113u"kg", water = 0.196u"kg")

"""
    LAVERGNE_LOSS_A, LAVERGNE_LOSS_B

Calibration of the calorimeter's heat loss, Lavergne et al. (2018) Eq. (23):
`φ(ΔT) = a ΔT + b ΔT²`, with `a = 75 J/(h·K)` and `b = 0.260 J/(h·K²)` from the
NF EN 196-9 calibration. Converted here to watts.
"""
const LAVERGNE_LOSS_A = 75.0 / 3600            # W/K
const LAVERGNE_LOSS_B = 0.260 / 3600           # W/K²

"""
    LAVERGNE_VESSEL_CP

Heat capacity of the calorimeter vessel, **380 J/K**.

The paper prints "about 380 kJ/K", and that cannot be the figure its own results
correspond to. Its Table 11 mix holds 371 g of binder, which releases of order
420 J/g, so about 156 kJ in total; against 380 kJ/K the temperature would rise by
0.4 K, where the test reports tens of kelvin. The rest of the setup is consistent
with joules: sand and water alone contribute roughly 1.6 kJ/K, so a vessel of
380 J/K puts the total near 2.1 kJ/K and the adiabatic rise near 75 K, which is
the order the measurements show. It is read as 380 J/K here, and this note is
deliberate: the alternative is to change a published number in silence.
"""
const LAVERGNE_VESSEL_CP = 380.0               # J/K

const _SAND_CP_PER_KG = Ref{Float64}(NaN)

"""
    sand_heat_capacity(mass) -> Float64

Heat capacity [J/K] of `mass` of quartz sand, from the `Qtz` entry of CEMDATA18
rather than from a remembered figure.
"""
function sand_heat_capacity(mass)
    if isnan(_SAND_CP_PER_KG[])
        q = first(s for s in build_species(IONIC_CEMDATA) if symbol(s) == "Qtz")
        cp = ustrip(us"J/K/mol", q[:Cp⁰](T = 293.15, P = 1.0e5, unit = true))
        _SAND_CP_PER_KG[] = cp / ustrip(us"kg/mol", q[:M])
    end
    return _SAND_CP_PER_KG[] * ustrip(us"kg", mass)
end

"""
    lavergne_semiadiabatic(; mix = LAVERGNE_MIX_C100, T0 = 293.15u"K", T_env = T0)

The NF EN 196-9 device of Lavergne et al. (2018), as a `SemiAdiabaticCalorimeter`.

`Cp` here is what the ODE does NOT compute for itself: the vessel and the inert
sand. The paste's own heat capacity is `Σᵢ nᵢ Cp⁰ᵢ(T)`, which `ChemistryLab` adds
at every step from the database, so it must not be counted twice.
"""
function lavergne_semiadiabatic(;
        mix = LAVERGNE_MIX_C100, T0 = 293.15u"K", T_env = T0,
    )
    Cp_fixed = LAVERGNE_VESSEL_CP + sand_heat_capacity(mix.sand)
    return SemiAdiabaticCalorimeter(;
        Cp = Cp_fixed * u"J/K",
        heat_loss = ΔT -> LAVERGNE_LOSS_A * ΔT + LAVERGNE_LOSS_B * ΔT^2,
        T_env = T_env,
        T0 = T0,
    )
end

"""
    langavant_temperature(t, qdot_per_g, states; mix = LAVERGNE_MIX_C100, T_env = 293.15)

Temperature history of the NF EN 196-9 cell, integrated from a heat rate measured
at `T_env`.

`qdot_per_g` is in W per gram of binder — what [`heat_release`](@ref) returns,
divided by the binder mass of the run. The total heat capacity is the vessel, the
sand, and the paste's own `Σᵢ nᵢ Cp⁰ᵢ(T)` read from `states` and scaled to the
mass of binder in `mix`.

!!! warning "The Arrhenius feedback is not included"
    The rate handed in was computed at `T_env`. The temperature reached in the
    cell accelerates the reactions — Parrot–Killoh carries activation energies of
    42, 21, 54 and 32 kJ/mol for C₃S, C₂S, C₃A and C₄AF — so the true peak comes
    earlier and higher. Closing that loop needs the heat source inside the ODE,
    which under partial equilibrium requires differentiating the equilibrium map;
    `KineticsProblem` refuses that combination with a warning rather than
    returning a number it cannot support.
"""
function langavant_temperature(
        t, qdot_per_g, states;
        mix = LAVERGNE_MIX_C100, T_env = 293.15,
    )
    m_binder_g = ustrip(us"kg", mix.binder) * 1000
    C_fixed = LAVERGNE_VESSEL_CP + sand_heat_capacity(mix.sand)
    φ(ΔT) = LAVERGNE_LOSS_A * ΔT + LAVERGNE_LOSS_B * ΔT^2

    T = fill(T_env, length(t))
    for i in 2:length(t)
        C_paste = ustrip(us"J/K", heat_capacity(states[i])) * m_binder_g / 1000
        C_tot = C_fixed + C_paste
        q = qdot_per_g[i] * m_binder_g                     # W
        T[i] = T[i - 1] + (t[i] - t[i - 1]) * (q - φ(T[i - 1] - T_env)) / C_tot
    end
    return T
end

# ── post-processing ───────────────────────────────────────────────────────────

"""
    ionic_phase_history(run, times; gel_water = GEL_WATER_PER_CSH, states = nothing)
        -> (times, Vector{Dict{String,Float64}}, Vector{Float64}, Vector{NamedTuple})

Volume fractions of the phase families of [`IONIC_GROUPS`](@ref) at each instant,
the pore-solution pH, and the sealed-curing porosity.

Each instant is re-speciated by [`speciated_states`](@ref), which certifies it
against the KKT conditions, so the fractions reported are those of a composition
proved optimal.

The porosity is the two-argument form, referred to the **fresh** specimen: a
sealed paste keeps the volume it was cast with, and the empty porosity left by
the Le Chatelier contraction is not a species at all. Dividing by the current,
shrunken volume instead understates it by some five points on a w/c = 0.5 paste
at 28 days.

Pass `states` to reuse a replay already in hand: certifying forty instants is the
dominant cost of the whole calculation, and the calorimetry needs the very same
states, so computing them twice doubles the run for nothing.
"""
function ionic_phase_history(run, times; gel_water = GEL_WATER_PER_CSH, states = nothing)
    V_ref = volume(run.state0).total
    fs = Vector{Dict{String, Float64}}(undef, length(times))
    phs = Vector{Float64}(undef, length(times))
    pors = Vector{NamedTuple}(undef, length(times))
    nmv = symbol.(run.cs.species)
    groups = vcat(
        [
            k => v for (k, v) in IONIC_GROUPS
                if any(x -> x in nmv, v isa AbstractString ? [v] : v)
        ],
        [ionic_water_group(run.cs)],
    )
    sts = states === nothing ? speciated_states(run.sol, run.kp; times = times) : states
    for (i, st) in enumerate(sts)
        f = volume_fractions(st, groups; reference = run.state0)
        d = Dict{String, Float64}(k => v for (k, v) in f)

        if gel_water > 0
            n_csh = ustrip(us"mol", moles(st, "Jennite"))
            V_H2O = ustrip(
                run.cs["H2O@"][:V⁰](
                    T = temperature(st), P = pressure(st); unit = true,
                ) / V_ref,
            )
            f_gel = min(n_csh * gel_water * V_H2O, get(d, "water", 0.0))
            d["C-S-H"] = get(d, "C-S-H", 0.0) + f_gel
            d["water"] = get(d, "water", 0.0) - f_gel
        end
        fs[i] = d
        phs[i] = something(pH(st), NaN)
        pors[i] = porosity(st, run.state0)
    end
    return times, fs, phs, pors
end

# ── driver ────────────────────────────────────────────────────────────────────

function main()
    clinker = (C3S = 0.65, C2S = 0.11, C3A = 0.11, C4AF = 0.08)
    tend = 28 * 86400.0
    times = 10 .^ range(log10(0.05 * 86400), log10(tend); length = 40)

    println("\nCEM I 52.5 N of Lavergne et al. (2018) Table 9, w/b = 0.50, 28 days.")
    println("Run twice: with 3.5 % limestone filler, and without it.\n")

    run_cal = run_ionic_hydration(;
        wb = 0.5, clinker, gypsum = 0.046, filler = 0.035, tend,
    )
    run_nol = run_ionic_hydration(;
        wb = 0.5, clinker, gypsum = 0.046, filler = 0.0, tend,
    )
    @printf "with limestone: %d accepted steps, retcode = %s\n" length(run_cal.sol.t) run_cal.sol.retcode
    @printf "no limestone  : %d accepted steps, retcode = %s\n\n" length(run_nol.sol.t) run_nol.sol.retcode

    # The certified replay is the expensive part: do it once and hand it to both
    # the phase history and the calorimetry.
    st_c = speciated_states(run_cal.sol, run_cal.kp; times = times)
    st_n = speciated_states(run_nol.sol, run_nol.kp; times = times)

    _, fr_c, pH_c, po_c = ionic_phase_history(run_cal, times; states = st_c)
    _, fr_n, _, _ = ionic_phase_history(run_nol, times; states = st_n)
    t_cal, Q_c, qd_c = heat_release(run_cal.sol, run_cal.kp; times, states = st_c)
    _, Q_n, qd_n = heat_release(run_nol.sol, run_nol.kp; times, states = st_n)
    g = 1000.0                                    # the runs simulate 1 kg of binder

    println("── the aluminate sequence, decided by the free energy alone ──")
    for (lbl, f) in ("with limestone" => fr_c, "no limestone" => fr_n)
        aft = [get(fi, "AFt", 0.0) for fi in f]
        j = argmax(aft)
        @printf "%-15s  AFt peak %.4f at %5.2f d   AFt(28 d) %.4f   AFm(28 d) %.4f\n" lbl aft[j] (
            times[j] / 86400
        ) aft[end] get(f[end], "AFm", 0.0)
    end

    @printf "\n── pore solution ──\npH from %.2f to %.2f\n" first(pH_c) last(pH_c)

    println("\n── isothermal calorimetry at 20 °C ──")
    @printf "monotone: with limestone %s, without %s\n" all(diff(Q_c) .>= -1.0e-9) all(
        diff(Q_n) .>= -1.0e-9
    )
    for (lbl, Q, qd) in (("with limestone", Q_c, qd_c), ("no limestone", Q_n, qd_n))
        j = argmax(qd)
        @printf "%-15s  Q: %5.1f (1 d) %5.1f (7 d) %5.1f (28 d) J/g   peak %.2f mW/g at %.2f h\n" lbl (
            Q[argmin(abs.(t_cal .- 86400))] / g
        ) (Q[argmin(abs.(t_cal .- 7 * 86400))] / g) (Q[end] / g) (qd[j] / g * 1000) (
            t_cal[j] / 3600
        )
    end

    println("\n── semi-adiabatic cell, NF EN 196-9 ──")
    m_binder_g = ustrip(us"kg", LAVERGNE_MIX_C100.binder) * 1000
    for (lbl, qd, Q, st) in (("with limestone", qd_c, Q_c, st_c), ("no limestone", qd_n, Q_n, st_n))
        T = langavant_temperature(t_cal, qd ./ g, st)
        j = argmax(T)
        C_tot = LAVERGNE_VESSEL_CP + sand_heat_capacity(LAVERGNE_MIX_C100.sand) +
            ustrip(us"J/K", heat_capacity(st[end])) * m_binder_g / 1000
        @printf "%-15s  ΔT max %.1f K at %.1f h    adiabatic ΔT(28 d) %.1f K\n" lbl (
            T[j] - 293.15
        ) (t_cal[j] / 3600) (Q[end] / g * m_binder_g / C_tot)
    end

    println("\n── sealed-curing porosity, referred to the fresh paste ──")
    @printf "%6s  %8s  %8s  %8s\n" "t [d]" "φ" "S" "void"
    for t in (0.25, 1.0, 7.0, 28.0)
        i = argmin(abs.(times ./ 86400 .- t))
        @printf "%6.2f  %8.4f  %8.4f  %8.4f\n" (times[i] / 86400) po_c[i].total (
            po_c[i].liquid / po_c[i].total
        ) get(fr_c[i], "void", 0.0)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main()
end
