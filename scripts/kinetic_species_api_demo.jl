# =============================================================================
# demo_kinetic_species.jl
#
# Demonstration of the new kinetic_species API of ChemicalSystem with
# full ODE integration (Leal et al. 2017 formulation).
#
# Workflow:
#   1. ChemicalSystem with kinetic_species → auto-generated reactions
#   2. ChemicalState (initial OPC composition)
#   3. KineticsProblem with SemiAdiabaticCalorimeter
#   4. ODE integration via OrdinaryDiffEq (Rodas5P)
#   5. Post-processing: mole and temperature tracking
#
# Usage:
#   julia --project=scripts scripts/kinetic_species_api_demo.jl
#   or from the REPL:  include("scripts/kinetic_species_api_demo.jl")
#   or, in VS Code, just run the file: it activates `scripts/` itself, and every
#   data file it reads is named with `datapath`, so no working directory matters.

import Pkg
Pkg.activate(@__DIR__; io = devnull)

using ChemistryLab
using OrdinaryDiffEq
using DynamicQuantities
using Printf

# ── 1. Load CEMDATA18 and select species ─────────────────────────────────────

const DATA_FILE = datapath("cemdata18-thermofun.json")

substances = build_species(DATA_FILE)

input_species = split(
    "C3S C2S C3A C4AF " *
        "Portlandite Jennite ettringite monosulphate12 C3AH6 C3FH6 " *
        "H2O@",
)

species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])

# ── 2. Parrot & Killoh kinetic models ────────────────────────────────────────

const WC = 0.4
const α_max = min(1.0, WC / 0.42)

# Blaine fineness of the binder, used by the canonical rate law through
# `blaine_factor`. 380 m²/kg is an ordinary CEM I.
const BLAINE = 380.0u"m^2/kg"

pk_C3S = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S"; α_max, blaine = BLAINE)
pk_C2S = parrot_killoh_avrami(PK84_PARAMS_C2S, "C2S"; α_max, blaine = BLAINE)
pk_C3A = parrot_killoh_avrami(PK84_PARAMS_C3A, "C3A"; α_max, blaine = BLAINE)
pk_C4AF = parrot_killoh_avrami(PK84_PARAMS_C4AF, "C4AF"; α_max, blaine = BLAINE)

# ── 3. ChemicalSystem with kinetic_species ───────────────────────────────────

cs = ChemicalSystem(
    species, CEMDATA_PRIMARIES;
    kinetic_species = Dict(
        "C3S" => pk_C3S,
        "C2S" => pk_C2S,
        "C3A" => pk_C3A,
        "C4AF" => pk_C4AF,
    ),
)

println("System: $(length(cs.species)) species, $(length(cs.SM.primaries)) primaries")
println("Kinetic species: ", join(symbol.(kinetic_species(cs)), ", "))
println("Kinetic reactions: $(length(cs.reactions))")

# ── 4. Initial state ─────────────────────────────────────────────────────────

const COMPOSITION = (C3S = 0.619, C2S = 0.165, C3A = 0.08, C4AF = 0.087)

state0 = ChemicalState(cs)
for (name, frac) in pairs(COMPOSITION)
    set_quantity!(state0, string(name), frac * u"kg")
end
set_quantity!(state0, "H2O@", WC * u"kg")

# ── 5. KineticsProblem ───────────────────────────────────────────────────────

kp = KineticsProblem(cs, state0, (0.0, 7 * 86400.0))

println("\nODE state dimension: ", length(build_u0(kp)))
println("  nₖ = ", length(kp.idx_kinetic), " kinetic species")

# ── 6. Integration ───────────────────────────────────────────────────────────

println("\nODE integration (7 days)...")
sol = integrate(kp)
println("  $(length(sol.t)) time steps, t_final = $(@sprintf("%.1f", sol.t[end] / 86400)) days")

# ── 7. Results ───────────────────────────────────────────────────────────────

println("\nKinetic species moles (initial → final):")
u0 = build_u0(kp)
u_end = sol.u[end]
for (j, idx) in enumerate(kp.idx_kinetic)
    sp = cs.species[idx]
    n0 = u0[j]
    nf = u_end[j]
    α = 1.0 - nf / n0
    @printf(
        "  %-6s : %.4f → %.4f mol  (α = %.1f%%)\n",
        symbol(sp), n0, nf, 100α
    )
end

println("\nChecks:")
println("  max |A·N| = ", maximum(abs, cs.SM.A * cs.SM.N))
println("  Integration completed successfully ✓")
