# =============================================================================
# opc_semiadiabatic_calorimetry.jl
#
# Simulation of the heat release of a CEM I paste in a semi-adiabatic
# calorimeter using the method of Lavergne et al. (2018).
#
# Hydration kinetics: Parrot & Killoh (1984) model
#   via ChemistryLab.parrot_killoh (KineticFunc),
#   temperature: Arrhenius correction per phase,
#   available water: maximum hydration limit α_max (Powers 1948).
#
# Calorimeter: semi-adiabatic, quadratic heat loss
#   φ(ΔT) = a·ΔT + b·ΔT²
#   (ChemistryLab.SemiAdiabaticCalorimeter)
#
# References:
#   Lavergne F., Ben Fraj A., Bayane I., Barthelemy J.-F. (2018).
#   Estimating the mechanical properties of hydrating blended cementitious
#   materials: an investigation based on micromechanics.
#   Cem. Concr. Res. 104, 37-60.
#   https://doi.org/10.1016/j.cemconres.2017.10.018
#
#   Parrot L.J. & Killoh D.C. (1984).
#   Prediction of cement hydration.
#   Br. Ceram. Proc. 35, 41-53.
#
# Usage:
#   julia --project=scripts scripts/opc_semiadiabatic_calorimetry.jl
#   or from the REPL:  include("scripts/opc_semiadiabatic_calorimetry.jl")
#   or, in VS Code, just run the file: it activates `scripts/` itself, and every
#   data file it reads is named with `datapath`, so no working directory matters.

import Pkg
Pkg.activate(@__DIR__; io = devnull)

using ChemistryLab
using OrdinaryDiffEq
using DynamicQuantities
using OrderedCollections
using Printf

# ── 1. ChemicalSystem from the CEMDATA18 database ───────────────────────────

const DATA_FILE = datapath("cemdata18-thermofun.json")

substances = build_species(DATA_FILE)

input_species = split(
    "C3S C2S C3A C4AF " *
        "Portlandite Jennite ettringite monosulphate12 C3AH6 C3FH6 " *
        "H2O@",
)

species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])
cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

@info "Chemical system: $(length(cs.species)) species, $(length(cs.reactions)) reactions"

# ── 2. Cement paste composition ──────────────────────────────────────────────
#
# Typical CEM I 52.5, from Lavergne et al. (2018)
# Mass fractions in the clinker

const PHASE_MASS_FRAC = (C3S = 0.619, C2S = 0.165, C3A = 0.08, C4AF = 0.087)

# Water/cement ratio
const WC = 0.4

# Maximum degree of hydration (Powers 1948 law: α_max ≤ w/c / 0.42)
const ALPHA_MAX = min(1.0, WC / 0.42)

state0 = ChemicalState(cs)
for (name, frac) in pairs(PHASE_MASS_FRAC)
    set_quantity!(state0, string(name), frac * u"kg")
end
set_quantity!(state0, "H2O@", WC * u"kg")

# ── 3. Parrot & Killoh kinetic models ────────────────────────────────────────
#
# Parameters K₁, N₁, K₂, N₂, K₃, N₃, B, Ea — Parrot & Killoh (1984)
# Ea — Schindler & Folliard (2005) / van Breugel (1991) [J/mol]

const PK_SMOOTHED_C3S = (
    K₁ = 1.5u"1/d", N₁ = 3.3, K₂ = 0.018u"1/d", N₂ = 2.5,
    K₃ = 0.0024u"1/d", N₃ = 4.0, B = 0.5,
    Ea = 41_570.0u"J/mol", T_ref = 293.15u"K",
)
const PK_SMOOTHED_C2S = (
    K₁ = 0.95u"1/d", N₁ = 0.5, K₂ = 0.0005u"1/d", N₂ = 2.5,
    K₃ = 0.0024u"1/d", N₃ = 4.0, B = 0.2,
    Ea = 43_670.0u"J/mol", T_ref = 293.15u"K",
)
const PK_SMOOTHED_C3A = (
    K₁ = 0.082u"1/d", N₁ = 0.87, K₂ = 0.00024u"1/d", N₂ = 2.0,
    K₃ = 0.0024u"1/d", N₃ = 4.0, B = 0.04,
    Ea = 54_040.0u"J/mol", T_ref = 293.15u"K",
)
const PK_SMOOTHED_C4AF = (
    K₁ = 0.165u"1/d", N₁ = 3.7, K₂ = 0.0015u"1/d", N₂ = 2.5,
    K₃ = 0.0024u"1/d", N₃ = 4.0, B = 0.5,
    Ea = 34_420.0u"J/mol", T_ref = 293.15u"K",
)

pk_C3S = parrot_killoh(PK_SMOOTHED_C3S, "C3S"; α_max = ALPHA_MAX)
pk_C2S = parrot_killoh(PK_SMOOTHED_C2S, "C2S"; α_max = ALPHA_MAX)
pk_C3A = parrot_killoh(PK_SMOOTHED_C3A, "C3A"; α_max = ALPHA_MAX)
pk_C4AF = parrot_killoh(PK_SMOOTHED_C4AF, "C4AF"; α_max = ALPHA_MAX)

# ── 4. Kinetic reaction list ─────────────────────────────────────────────────
#
# Balanced hydration reactions (CEMDATA18).
# Jennite in CEMDATA18 = (SiO₂)(CaO)_{5/3}(H₂O)_{21/10}, Ca:Si = 5/3.
# ΔᵣH⁰ is computed automatically from species ΔₐH⁰ (requires balanced reactions).
#
#   C₃S  + 103/30 H₂O  →  Jennite  + 4/3 Portlandite       (ΔᵣH⁰ ≈ −124 kJ/mol)
#   C₂S  +  73/30 H₂O  →  Jennite  + 1/3 Portlandite       (ΔᵣH⁰ ≈  −48 kJ/mol)
#   C₃A  +  6     H₂O  →  C₃AH₆                            (ΔᵣH⁰ ≈ −261 kJ/mol)
#   C₄AF + 2 Portlandite + 10 H₂O → C₃AH₆ + C₃FH₆         (ΔᵣH⁰ ≈ −147 kJ/mol)

sp(name) = cs[name]

rxn_C3S = Reaction(
    OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 103 / 30),
    OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 4 / 3);
    symbol = "C₃S hydration",
)
rxn_C3S[:rate] = pk_C3S

rxn_C2S = Reaction(
    OrderedDict(sp("C2S") => 1.0, sp("H2O@") => 73 / 30),
    OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 1 / 3);
    symbol = "C₂S hydration",
)
rxn_C2S[:rate] = pk_C2S

rxn_C3A = Reaction(
    OrderedDict(sp("C3A") => 1.0, sp("H2O@") => 6.0),
    OrderedDict(sp("C3AH6") => 1.0);
    symbol = "C₃A hydration",
)
rxn_C3A[:rate] = pk_C3A

rxn_C4AF = Reaction(
    OrderedDict(sp("C4AF") => 1.0, sp("Portlandite") => 2.0, sp("H2O@") => 10.0),
    OrderedDict(sp("C3AH6") => 1.0, sp("C3FH6") => 1.0);
    symbol = "C₄AF hydration",
)
rxn_C4AF[:rate] = pk_C4AF

kinetic_reactions = [rxn_C3S, rxn_C2S, rxn_C3A, rxn_C4AF]

# ── 6. Semi-adiabatic calorimeter ────────────────────────────────────────────
#
# Langavant-type device (standard NF EN 196-9).
# Quadratic heat loss φ(ΔT) = a·ΔT + b·ΔT² (Lavergne et al. 2018).
# Heat capacity Cp [J/K] for 1 kg cement + WC kg water + flask:
#   cp_cement ≈ 800 J/(kg·K), cp_water = 4186 J/(kg·K), cp_flask ≈ 900 J/(kg·K)

const T0_K = 20.0 + 273.15    # initial temperature [K]
const T_ENV_K = 20.0 + 273.15  # ambient temperature [K]
const CP = 1.0 * 800.0 + WC * 4186.0 + 1.0 * 900.0   # ≈ 3449 J/K

const A_CAL = 0.3   # linear coefficient  [W/K]
const B_CAL = 0.003  # quadratic coefficient [W/K²]

cal = SemiAdiabaticCalorimeter(;
    Cp = CP * u"J/K",
    T_env = T_ENV_K * u"K",
    heat_loss = ΔT -> A_CAL * ΔT + B_CAL * ΔT^2,
    T0 = T0_K * u"K",
)

# ── 7. Kinetics problem ─────────────────────────────────────────────────────

const TSPAN = (0.0, 7.0 * 86400.0)

kp = KineticsProblem(
    cs,
    kinetic_reactions,
    state0,
    TSPAN;
    calorimeter = cal,
    equilibrium_solver = nothing,
)

# ── 8. Integration ───────────────────────────────────────────────────────────

@info "Solving (7 days, Rodas5P)..."
ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-9)
sol = integrate(kp, ks)
@info "Done: $(length(sol.t)) accepted steps."

# ── 9. Post-processing ──────────────────────────────────────────────────────

t_h = sol.t ./ 3600.0   # seconds → hours

t_T, T_K_vec = temperature_profile(sol, cal)
t_Q, Q_J_vec = cumulative_heat(sol, cal)
T_°C_vec = T_K_vec .- 273.15
Q_kJ_vec = Q_J_vec ./ 1000.0

_, qdot_W = heat_flow(sol, cal)   # [W]

# Degree of hydration per phase
n0_kin = [sol.prob.p.n_initial_full[i] for i in kp.idx_kinetic]
n_kin = [[u[i] for u in sol.u] for i in eachindex(n0_kin)]

# idx_kinetic lists the kinetic species in order of appearance.
# The position of each clinker phase depends on cs.species; we retrieve
# it by comparing with kp.idx_kinetic.
function phase_alpha(cs, kp, sol, n0_kin, n_kin, name)
    sp_idx = findfirst(sp -> ChemistryLab.symbol(sp) == name, cs.species)
    pos = findfirst(==(sp_idx), kp.idx_kinetic)
    isnothing(pos) && return fill(NaN, length(sol.t))
    return 1.0 .- n_kin[pos] ./ n0_kin[pos]
end

α_C3S_vec = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C3S")
α_C2S_vec = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C2S")
α_C3A_vec = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C3A")
α_C4AF_vec = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C4AF")

# Mass-weighted mean degree of hydration
α_mean = (
    PHASE_MASS_FRAC.C3S .* α_C3S_vec .+
        PHASE_MASS_FRAC.C2S .* α_C2S_vec .+
        PHASE_MASS_FRAC.C3A .* α_C3A_vec .+
        PHASE_MASS_FRAC.C4AF .* α_C4AF_vec
)

# ── Summary ──────────────────────────────────────────────────────────────────
println()
println("╔══════════════════════════════════════════════╗")
println("║    Semi-adiabatic calorimetry results         ║")
println("║    CEM I — Lavergne et al. (2018)             ║")
println("╠══════════════════════════════════════════════╣")
println("║  Duration    = 7 days")
@printf "║  T final     = %.2f °C\n" T_°C_vec[end]
@printf "║  ΔT max      = %.2f °C\n" maximum(T_°C_vec) - (T0_K - 273.15)
@printf "║  α(C3S)      = %.4f\n"   α_C3S_vec[end]
@printf "║  α(C2S)      = %.4f\n"   α_C2S_vec[end]
@printf "║  α(C3A)      = %.4f\n"   α_C3A_vec[end]
@printf "║  α(C4AF)     = %.4f\n"   α_C4AF_vec[end]
@printf "║  α mean      = %.4f\n"   α_mean[end]
@printf "║  Q total     = %.1f kJ/kg_cement\n" Q_kJ_vec[end]
println("╚══════════════════════════════════════════════╝")

# ── 10. Plots ────────────────────────────────────────────────────────────────

using Plots
gr()

p1 = plot(
    t_T ./ 3600, T_°C_vec;
    xlabel = "Time [h]", ylabel = "Temperature [°C]",
    title = "Temperature profile — semi-adiabatic calorimeter",
    label = "T(t)", lw = 2, color = :red,
)
hline!(p1, [T0_K - 273.15]; linestyle = :dash, color = :gray, label = "T₀ = T_env")

p2 = plot(
    t_h, [α_C3S_vec α_C2S_vec α_C3A_vec α_C4AF_vec α_mean];
    xlabel = "Time [h]", ylabel = "Degree of hydration α",
    title = "Clinker phase hydration (Parrot-Killoh model)",
    label = ["C₃S" "C₂S" "C₃A" "C₄AF" "ᾱ mean"],
    lw = 2, ls = [:solid :dash :dot :dashdot :solid],
)
hline!(p2, [ALPHA_MAX]; linestyle = :dash, color = :black, label = "α_max")

p3 = plot(
    t_T ./ 3600, qdot_W ./ 1000;
    xlabel = "Time [h]", ylabel = "Heat flow [kW/kg_cement]",
    title = "Instantaneous heat flow",
    label = "q̇(t)", lw = 2, color = :orange,
)

p4 = plot(
    t_Q ./ 3600, Q_kJ_vec;
    xlabel = "Time [h]", ylabel = "Q [kJ/kg_cement]",
    title = "Cumulative heat",
    label = "Q(t)", lw = 2, color = :purple,
)

display(
    plot(
        p1, p2, p3, p4; layout = (2, 2), size = (1200, 800),
        plot_title = "Lavergne et al. (2018) — CEM I, w/c = $WC",
    )
)
