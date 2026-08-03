# =============================================================================
# slag_clinker_equilibration.jl
#
# Equilibrium mineral assemblage for a hydrated OPC clinker paste with slag-
# related secondary phases (hydrotalcite, AFm, C3AFS). Extended species list
# compared to clinker_hydrate_equilibration.jl.
#
# Usage:
#   julia --project scripts/slag_clinker_equilibration.jl
# =============================================================================

using Pkg
Pkg.activate(@__DIR__)

using ChemistryLab
using DynamicQuantities
using Optimization, OptimizationIpopt
using OptimaSolver
using LinearAlgebra

# ── Species and chemical system ──────────────────────────────────────────────

substances = build_species("data/cemdata18-thermofun.json")

clinker = "C3S C2S C3A C4AF Gp Anh"
hydrates = "Portlandite Jennite H2O@ ettringite monosulphate12 C3AH6 C3FH6 C4FH13"
slag = "hydrotalcite C3AFS0.84H4.32 C3AS0.84H4.32 C4AH13"

species = speciation(
    substances,
    split(clinker * " " * hydrates * " " * slag);
    aggregate_state = [AS_AQUEOUS],
)

cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

# ── Initial state ────────────────────────────────────────────────────────────
# Clinker + gypsum + MgSO4 additive, w/c = 0.4, normalized to 1 kg total

state = ChemicalState(cs)
compo = [
    "C3S" => 67.8 / 100,
    "C2S" => 16.6 / 100,
    "C3A" => 4 / 100,
    "C4AF" => 7.2 / 100,
    "Gp" => 2.8 / 100,
    "MgSO4@" => 1.6 / 100,
]
c = sum(last.(compo))
wc = 0.4
w = wc * c
mtot = c + w
for (name, frac) in compo
    set_quantity!(state, name, frac / mtot * u"kg")
end
set_quantity!(state, "H2O@", w / mtot * u"kg")
# H⁺ and OH⁻ auto-seeded at neutral pH when water was added

# ── Equilibrate ──────────────────────────────────────────────────────────────

state_eq = equilibrate(state)

display(state_eq)

# ── Gibbs energy check ───────────────────────────────────────────────────────

μ = build_potentials(cs, DiluteSolutionModel())
p = ChemistryLab._build_params(state; ϵ = 1.0e-16)
Gini = μ(ustrip.(state.n), p) ⋅ ustrip.(state.n)
Gfin = μ(ustrip.(state_eq.n), p) ⋅ ustrip.(state_eq.n)

println("G initial = ", round(Gini, digits = 6))
println("G final   = ", round(Gfin, digits = 6))
println("ΔG        = ", round(Gfin - Gini, digits = 6))
