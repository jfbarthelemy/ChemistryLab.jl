# =============================================================================
# aqueous_equilibration_basic.jl
#
# Minimal example: equilibrate a dilute NaOH solution (1 mmol in 1 kg H₂O).
# Demonstrates the simplest ChemicalSystem → ChemicalState → equilibrate flow.
#
# Usage:
#   julia --project=scripts scripts/aqueous_equilibration_basic.jl
#   or from the REPL:  include("scripts/aqueous_equilibration_basic.jl")
#   or, in VS Code, just run the file: it activates `scripts/` itself, and every
#   data file it reads is named with `datapath`, so no working directory matters.

import Pkg
Pkg.activate(@__DIR__; io = devnull)

using ChemistryLab
using Optimization, OptimizationIpopt
using OptimaSolver
using DynamicQuantities

# ── Species and chemical system ──────────────────────────────────────────────

substances = build_species(datapath("slop98-inorganic-thermofun.json"))
dict = Dict(symbol(s) => s for s in substances)
species = [dict[s] for s in split("H2O@ Na+ NaOH@ H+ OH-")]

cs = ChemicalSystem(species)

# ── Initial state ────────────────────────────────────────────────────────────

state = ChemicalState(cs)
set_quantity!(state, "NaOH@", 1.0e-3u"mol")
set_quantity!(state, "H2O@", 1.0u"kg")
# H⁺ and OH⁻ auto-seeded at neutral pH when water was added

# ── Equilibrate ──────────────────────────────────────────────────────────────

state_eq = equilibrate(state)

display(state_eq)
