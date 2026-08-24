# =============================================================================
# clinker_hydrate_equilibration.jl
#
# Compute the equilibrium mineral assemblage of a hydrated OPC clinker paste
# (C3S, C2S, C3A, C4AF + gypsum + water) by Gibbs energy minimization.
# Compares the result against a Reaktoro reference solution.
#
# Solver: OptimaSolver (default) or Ipopt — selected via USE_OPTIMA flag.
#
# Usage:
#   julia --project=scripts scripts/clinker_hydrate_equilibration.jl
#   or from the REPL:  include("scripts/clinker_hydrate_equilibration.jl")
#   or, in VS Code, just run the file: it activates `scripts/` itself, and every
#   data file it reads is named with `datapath`, so no working directory matters.

import Pkg
Pkg.activate(@__DIR__; io = devnull)

using ChemistryLab
using Optimization, OptimizationIpopt
using OptimaSolver
using DynamicQuantities
using LinearAlgebra

# ── Solver selection ─────────────────────────────────────────────────────────

const USE_OPTIMA = true

# ── Species and chemical system ──────────────────────────────────────────────

substances = build_species(datapath("cemdata18-thermofun.json"))
input_species = split("C3S C2S C3A C4AF Gp Anh Portlandite Jennite H2O@")
species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])

cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

# ── Initial state ────────────────────────────────────────────────────────────
# Typical OPC clinker composition (mass fractions), w/c = 0.4

state = ChemicalState(cs)
compo = [
    "C3S" => 67.8 / 100,
    "C2S" => 16.6 / 100,
    "C3A" => 4 / 100,
    "C4AF" => 7.2 / 100,
    "Gp" => 2.8 / 100,
]
for (name, frac) in compo
    set_quantity!(state, name, frac * u"kg")
end
w = 0.4 * sum(last.(compo))
set_quantity!(state, "H2O@", w * u"kg")
rescale!(state, 1.0u"kg")
# H⁺ and OH⁻ auto-seeded at neutral pH when water was added

# ── Equilibrate ──────────────────────────────────────────────────────────────

if USE_OPTIMA
    state_eq = equilibrate(
        state,
        OptimaOptimizer(tol = 1.0e-10, verbose = true);
        variable_space = Val(:linear),
    )
else
    state_eq = equilibrate(state; variable_space = Val(:linear), verbose = 5)
end

display(state_eq)

# ── Gibbs energy comparison ──────────────────────────────────────────────────

μ = build_potentials(cs, DiluteSolutionModel())
p = ChemistryLab._build_params(state; ϵ = 1.0e-16)
Gini = μ(ustrip.(state.n), p) ⋅ ustrip.(state.n)
Gfin = μ(ustrip.(state_eq.n), p) ⋅ ustrip.(state_eq.n)

println("G initial  = ", round(Gini, digits = 6))
println("G final    = ", round(Gfin, digits = 6))
println("ΔG         = ", round(Gfin - Gini, digits = 6))

# ── Reference solution (Reaktoro) ────────────────────────────────────────────

values = [
    6.3267e+0, 5.2875e-6, 4.7599e-11, 6.0853e-10, 1.0e-16, 2.4934e-10,
    1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16,
    3.1801e-1, 1.0e-16, 6.5027e-13, 1.4563e-6, 4.3002e-1, 2.1376e-1,
    1.0e-16, 1.0e-16, 6.9034e-8, 6.8628e-6, 7.9351e-3, 1.1895e-15,
    6.9034e-8, 1.337e-3, 1.7474e-16, 1.0e-16, 8.0774e-16, 1.7135e-4,
    1.0e-16, 1.0e-16, 1.0e-16, 2.5682e-7, 1.0e-16, 1.0e-16, 1.0e-16,
    8.1095e-4, 5.8068e-9, 2.8552e+0, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16,
    3.534e+0, 1.1724e-1, 4.806e-16,
]
state_rkt = copy(state)
state_rkt.n .= values * u"mol"
ChemistryLab._update_derived!(state_rkt)
Grkt = μ(ustrip.(state_rkt.n), p) ⋅ ustrip.(state_rkt.n)

println("G (Reaktoro ref) = ", round(Grkt, digits = 6))
println("ΔG vs Reaktoro   = ", round(Gfin - Grkt, digits = 6))
