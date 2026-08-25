# =============================================================================
# cement_clinker_kinetics.jl
#
# Kinetic simulation of OPC clinker hydration via ChemistryLab:
#   - ChemicalSystem built from the CEMDATA18 database
#   - KineticsProblem with parrot_killoh (KineticFunc) for the 4 clinker phases
#   - Heat tracking via SemiAdiabaticCalorimeter
#
# This script demonstrates the full workflow:
#   1. ChemicalSystem  → 2. ChemicalState  → 3. list of annotated Reactions
#   → 4. KineticsProblem → 5. integrate → 6. post-processing
#
# Kinetics: Parrot & Killoh (1984), Arrhenius correction from Schindler &
#   Folliard (2005).
# Semi-adiabatic calorimeter with quadratic heat losses (Lavergne et al. 2018).
#
# Usage:
#   julia --project=scripts scripts/cement_clinker_kinetics.jl
#   or from the REPL:  include("scripts/cement_clinker_kinetics.jl")
#   or, in VS Code, just run the file: it activates `scripts/` itself, and every
#   data file it reads is named with `datapath`, so no working directory matters.

import Pkg
Pkg.activate(@__DIR__; io = devnull)

using ChemistryLab
using OrdinaryDiffEq
using DynamicQuantities
using OrderedCollections
using Printf

# ── 1. ChemicalSystem from the CEMDATA18 database ────────────────────────────
#
# We select:
#   - the 4 clinker phases (kinetic species)
#   - the main hydration products
#   - liquid water (solvent)
# CEMDATA_PRIMARIES provides the independent components of the system.

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

# ── 2. Initial state ─────────────────────────────────────────────────────────
#
# Typical CEM I 52.5 R composition (Lavergne et al. 2018):
# mass fractions in the clinker. We work with 1 kg of cement.

const WC = 0.4         # water/cement ratio [-]
const COMPOSITION = (            # clinker phase mass fractions [kg/kg cement]
    C3S = 0.619,
    C2S = 0.165,
    C3A = 0.08,
    C4AF = 0.087,
)

state0 = ChemicalState(cs)

for (name, frac) in pairs(COMPOSITION)
    set_quantity!(state0, string(name), frac * u"kg")
end
set_quantity!(state0, "H2O@", WC * u"kg")

# ── 3. Parrot & Killoh kinetic models ─────────────────────────────────────────
#
# Maximum degree of hydration according to Powers (1948): α_max ≤ w/c / 0.42
# We create models with the computed α_max rather than the default value 1.0.

const α_max = min(1.0, WC / 0.42)

pk_C3S = parrot_killoh(PK_PARAMS_C3S, "C3S"; α_max)
pk_C2S = parrot_killoh(PK_PARAMS_C2S, "C2S"; α_max)
pk_C3A = parrot_killoh(PK_PARAMS_C3A, "C3A"; α_max)
pk_C4AF = parrot_killoh(PK_PARAMS_C4AF, "C4AF"; α_max)

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
#
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

# ── 5. Kinetics problem + calorimeter ────────────────────────────────────────
#
# tspan: 7 days in seconds
# equilibrium_solver = nothing → no aqueous re-speciation (sufficient for PK)

const TSPAN = (0.0, 7.0 * 86400.0)

# Semi-adiabatic calorimeter parameters (Lavergne et al. 2018)
#   Cp [J/K]: 1 kg cement + WC kg water + Dewar flask
#   quadratic losses: φ(ΔT) = a·ΔT + b·ΔT²
cal = SemiAdiabaticCalorimeter(;
    Cp = (1.0 * 800.0 + WC * 4186.0 + 1.0 * 900.0) * u"J/K",   # ≈ 3449 J/K
    T_env = 293.15u"K",
    heat_loss = ΔT -> 0.3 * ΔT + 0.003 * ΔT^2,
    T0 = 293.15u"K",
)

kp = KineticsProblem(
    cs,
    kinetic_reactions,
    state0,
    TSPAN;
    calorimeter = cal,
    equilibrium_solver = nothing,
)

# ── 6. Integration ───────────────────────────────────────────────────────────

@info "Integration in progress (7 days, Rodas5P)..."
ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-9)
sol = integrate(kp, ks)
@info "Done: $(length(sol.t)) accepted steps."

# ── 7. Post-processing ──────────────────────────────────────────────────────

t_h = sol.t ./ 3600.0   # [s] → [h]

# Temperature and heat from the calorimeter
t_T, T_K_vec = temperature_profile(sol, cal)
t_Q, Q_J_vec = cumulative_heat(sol, cal)
T_°C_vec = T_K_vec .- 273.15
Q_kJ_vec = Q_J_vec ./ 1000.0

# Degree of hydration per phase
n0_kin = [sol.prob.p.n_initial_full[i] for i in kp.idx_kinetic]
n_kin = [[u[i] for u in sol.u] for i in eachindex(n0_kin)]

function phase_alpha(cs, kp, sol, n0_kin, n_kin, name)
    sp_idx = findfirst(sp -> ChemistryLab.symbol(sp) == name, cs.species)
    pos = findfirst(==(sp_idx), kp.idx_kinetic)
    isnothing(pos) && return fill(NaN, length(sol.t))
    return 1.0 .- n_kin[pos] ./ n0_kin[pos]
end

α_C3S = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C3S")
α_C2S = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C2S")
α_C3A = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C3A")
α_C4AF = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C4AF")

# Mass-weighted mean degree of hydration
w = COMPOSITION
α_mean = (
    w.C3S .* α_C3S .+ w.C2S .* α_C2S .+
        w.C3A .* α_C3A .+ w.C4AF .* α_C4AF
) ./
    (w.C3S + w.C2S + w.C3A + w.C4AF)

# Summary
println()
println("╔══════════════════════════════════════════════════╗")
println("║   Clinker kinetics — KineticsProblem              ║")
println("║   CEM I, w/c = $WC, T₀ = 20 °C, 7 days            ║")
println("╠══════════════════════════════════════════════════╣")
@printf "║  ΔT max      = %6.2f °C\n"  maximum(T_°C_vec) - 20.0
@printf "║  α(C3S)      = %6.4f\n"    α_C3S[end]
@printf "║  α(C2S)      = %6.4f\n"    α_C2S[end]
@printf "║  α(C3A)      = %6.4f\n"    α_C3A[end]
@printf "║  α(C4AF)     = %6.4f\n"    α_C4AF[end]
@printf "║  ᾱ mean       = %6.4f\n"    α_mean[end]
@printf "║  Q total     = %6.1f kJ/kg\n" Q_kJ_vec[end]
println("╚══════════════════════════════════════════════════╝")

# ── 8. Plots ─────────────────────────────────────────────────────────────────

using Plots
gr()

p1 = plot(
    t_T ./ 3600, T_°C_vec;
    xlabel = "Time [h]", ylabel = "T [°C]",
    title = "Temperature (semi-adiabatic calorimeter)",
    label = "T(t)", lw = 2, color = :red,
)
hline!(p1, [20.0]; linestyle = :dash, color = :gray, label = "T₀ = T_env")

p2 = plot(
    t_h, [α_C3S α_C2S α_C3A α_C4AF α_mean];
    xlabel = "Time [h]", ylabel = "Degree of hydration α",
    title = "Clinker phase hydration",
    label = ["C₃S" "C₂S" "C₃A" "C₄AF" "ᾱ mean"],
    lw = 2, ls = [:solid :dash :dot :dashdot :solid],
)
hline!(p2, [α_max]; linestyle = :dash, color = :black, label = "α_max (Powers)")

p3 = plot(
    t_Q ./ 3600, Q_kJ_vec;
    xlabel = "Time [h]", ylabel = "Q [kJ/kg cement]",
    title = "Cumulative heat",
    label = "Q(t)", lw = 2, color = :purple,
)

display(
    plot(
        p1, p2, p3; layout = (1, 3), top_margin = 7Plots.mm, left_margin = 8Plots.mm, bottom_margin = 8Plots.mm, size = (1400, 420),
        plot_title = "ChemistryLab — KineticsProblem cement — CEM I w/c=$WC",
    )
)
