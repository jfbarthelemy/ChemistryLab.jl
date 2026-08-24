# =============================================================================
# titration_maleic_acid.jl
#
# Titration curve of maleic acid (diprotic, 0.1 M, 100 mL) by NaOH (2 M).
# Two equivalence points at ~5 mL and ~10 mL; pKa₁ = 1.92, pKa₂ = 6.27.
#
# Solver: OptimaSolver (default) or Ipopt — selected via USE_OPTIMA flag.
#
# Usage:
#   julia --project=scripts scripts/titration_maleic_acid.jl
#   or from the REPL:  include("scripts/titration_maleic_acid.jl")
#   or, in VS Code, just run the file: it activates `scripts/` itself, and every
#   data file it reads is named with `datapath`, so no working directory matters.

import Pkg
Pkg.activate(@__DIR__; io = devnull)

using ChemistryLab
using Optimization, OptimizationIpopt
using OptimaSolver
using DynamicQuantities
using ProgressMeter

# ── Solver selection ─────────────────────────────────────────���───────────────

const USE_OPTIMA = true

# ── Species and chemical system ──────────────────────────────────────────────

substances_inorg = build_species(datapath("slop98-inorganic-thermofun.json"))
substances_org = build_species(datapath("slop98-organic-thermofun.json"))

dict_all_species = merge(
    Dict(symbol(s) => s for s in substances_inorg),
    Dict(symbol(s) => s for s in substances_org),
)
species = [dict_all_species[s] for s in split("H2O@ Na+ NaOH@ H+ OH- MalH2@ MalH- Mal-2")]

cs = ChemicalSystem(species, ["H2O@", "H+", "Mal-2", "Na+", "Zz"])

# ── Equilibrium solver ───────────────────────────────────────────────────────

solver = if USE_OPTIMA
    EquilibriumSolver(
        cs, DiluteSolutionModel(),
        OptimaOptimizer(tol = 1.0e-12, warm_start = true);
        variable_space = Val(:linear),
    )
else
    EquilibriumSolver(
        cs, DiluteSolutionModel(),
        IpoptOptimizer(mu_strategy = "adaptive");
        variable_space = Val(:linear),
        abstol = 1.0e-8, reltol = 1.0e-8, maxiters = 100, verbose = 0,
    )
end

# ── Titration parameters ─────────────────────────────────────────────────────

const V_acid = 100.0e-3    # acid solution volume [L]
const c_acid = 0.1          # maleic acid concentration [mol/L]
const c_base = 2.0          # NaOH concentration [mol/L]
const n_H2A = V_acid * c_acid   # total moles of H₂A

const V_eq1 = n_H2A / c_base * 1.0e3       # first equivalence point [mL]
const V_eq2 = 2 * n_H2A / c_base * 1.0e3   # second equivalence point [mL]

# ── Titration loop ───────────────────────────────────────────────────────────

volumes_NaOH = range(0, 15; length = 201)   # [mL]
pH_vals = Float64[]

s = ChemicalState(cs)
@showprogress for V_mL in volumes_NaOH
    V_NaOH = V_mL * 1.0e-3        # [L]
    n_NaOH = c_base * V_NaOH      # [mol] NaOH added
    V_total = V_acid + V_NaOH     # total volume [L]

    set_quantity!(s, "MalH2@", n_H2A * u"mol")
    set_quantity!(s, "NaOH@", n_NaOH * u"mol")
    set_quantity!(s, "H2O@", V_total * u"kg")
    set_neutral_pH!(s)

    s_eq = solve(solver, s)
    push!(pH_vals, pH(s_eq))
end

# ── Results ──────────────────────────────────────────────────────────────────

println("pH at V = 0 mL   (pure acid)         : ", round(pH_vals[1], digits = 2))
println("pH at V = 2.5 mL (½ PE₁, ≈ pKa₁)   : ", round(pH_vals[6], digits = 2))
println("pH at V = 5 mL   (PE₁)              : ", round(pH_vals[11], digits = 2))
println("pH at V = 7.5 mL (½ PE₂, ≈ pKa₂)   : ", round(pH_vals[16], digits = 2))
println("pH at V = 10 mL  (PE₂)              : ", round(pH_vals[21], digits = 2))
println("pH at V = 15 mL  (excess NaOH)       : ", round(pH_vals[26], digits = 2))

# ── Plot ─────────────────────────────────────────────────────────────────────

using Plots
gr()

pKa1 = 1.92
pKa2 = 6.27

solver_name = USE_OPTIMA ? "OptimaSolver" : "Ipopt"
p = plot(
    collect(volumes_NaOH), pH_vals;
    xlabel = "V(NaOH) [mL]", ylabel = "pH",
    label = "Titration curve ($solver_name)",
    lw = 2, marker = :circle, markersize = 3, color = :steelblue,
    title = "Titration of maleic acid (0.1 M) by NaOH (2 M)",
    ylims = (0, 14), legend = :topleft,
)
vline!(p, [V_eq1]; ls = :dash, color = :red, label = "PE₁ ($(round(V_eq1, digits = 1)) mL)")
vline!(p, [V_eq2]; ls = :dash, color = :blue, label = "PE₂ ($(round(V_eq2, digits = 1)) mL)")
hline!(p, [pKa1]; ls = :dot, color = :orange, label = "pKₐ₁ = $pKa1")
hline!(p, [pKa2]; ls = :dot, color = :green, label = "pKₐ₂ = $pKa2")
display(p)
