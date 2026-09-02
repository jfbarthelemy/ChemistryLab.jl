# =============================================================================
# titration_malonic_acid.jl
#
# Titration curve of malonic acid (diprotic, 0.1 M, 100 mL) by NaOH (2 M),
# with both pKa values DERIVED from the database rather than tabulated.
#
# The species are malonate, not maleate. SLOP98 names them `MALONIC-ACID,AQ`,
# `H-MALONATE,AQ` and `MALONATE,AQ`, and their formula carries three carbons:
# malonic acid is HOOC-CH2-COOH (propanedioic acid, C3H4O4), while maleic acid
# is the cis-butenedioic acid C4H4O4. This script was called "maleic" and
# printed maleic acid's tabulated pKa (1.92 / 6.27) over a curve computed from
# malonate data, so the two dotted lines crossed the curve at neither
# half-equivalence point. Neither maleic nor fumaric acid is in SLOP98.
#
# Solver: OptimaSolver (default) or Ipopt — selected via USE_OPTIMA flag.
#
# Usage:
#   julia --project=scripts scripts/titration_malonic_acid.jl
#   or from the REPL:  include("scripts/titration_malonic_acid.jl")
#   or, in VS Code, just run the file: it activates `scripts/` itself, and every
#   data file it reads is named with `datapath`, so no working directory matters.

import Pkg
Pkg.activate(@__DIR__; io = devnull)

using ChemistryLab
using Optimization, OptimizationIpopt
using OptimaSolver
using DynamicQuantities
using ProgressMeter

# ── Solver selection ─────────────────────────────────────────────────────────

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

# ── Dissociation constants, from the database ────────────────────────────────
#
# Nothing here is tabulated. Both constants come from the standard Gibbs
# energies of formation the database ships, exactly as in
# `titration_acetic_acid.jl` and `co2_carbonate_system.jl`. The literature
# values are printed alongside for comparison, never used in the calculation.

R = ustrip(Constants.R)
T = 298.15   # K
RT = R * T

sp = Dict(symbol(s) => s for s in cs.species)
_pKa(num, den) = -log10(
    exp(-(sum(sp[x].ΔₐG⁰(T = T) for x in num) - sp[den].ΔₐG⁰(T = T)) / RT),
)

# MalH2@ ⇌ MalH⁻ + H⁺   and   MalH⁻ ⇌ Mal²⁻ + H⁺
pKa1 = _pKa(("MalH-", "H+"), "MalH2@")
pKa2 = _pKa(("Mal-2", "H+"), "MalH-")

println("pKa1 (MalH2@ / MalH⁻) = ", round(pKa1, digits = 2), "   (lit. malonic 2.83)")
println("pKa2 (MalH⁻  / Mal²⁻) = ", round(pKa2, digits = 2), "   (lit. malonic 5.69)")

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

const V_acid = 100.0e-3     # acid solution volume [L]
const c_acid = 0.1          # malonic acid concentration [mol/L]
const c_base = 2.0          # NaOH concentration [mol/L]
const n_H2A = V_acid * c_acid   # total moles of H₂A

const V_eq1 = n_H2A / c_base * 1.0e3       # first equivalence point [mL]
const V_eq2 = 2 * n_H2A / c_base * 1.0e3   # second equivalence point [mL]

# ── Titration loop ───────────────────────────────────────────────────────────

volumes_NaOH = collect(range(0, 15; length = 201))   # [mL], step 0.075 mL
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
#
# Indices are looked up from the volume, never hard-coded: the grid has 201
# points over 15 mL, so a step is 0.075 mL. An earlier version carried indices
# from a 61-point grid and mislabeled all six lines — `pH_vals[6]` was announced
# as 2.5 mL when it is 0.375 mL.

_at(V) = argmin(abs.(volumes_NaOH .- V))

for (V, what) in (
        (0.0, "pure acid"),
        (V_eq1 / 2, "½ PE₁, ≈ pKa₁"),
        (V_eq1, "PE₁"),
        (V_eq1 + V_eq1 / 2, "½ PE₂, ≈ pKa₂"),
        (V_eq2, "PE₂"),
        (15.0, "excess NaOH"),
    )
    i = _at(V)
    println(
        "pH at V = ", lpad(round(volumes_NaOH[i], digits = 3), 6), " mL  ",
        rpad("($what)", 18), " : ", round(pH_vals[i], digits = 2),
    )
end

# ── Plot ─────────────────────────────────────────────────────────────────────

using Plots
gr()

solver_name = USE_OPTIMA ? "OptimaSolver" : "Ipopt"
p = plot(
    volumes_NaOH, pH_vals;
    xlabel = "V(NaOH) [mL]", ylabel = "pH",
    label = "Titration curve ($solver_name)",
    lw = 2, marker = :circle, markersize = 3, color = :steelblue,
    title = "Titration of malonic acid (0.1 M) by NaOH (2 M)",
    ylims = (0, 14), legend = :topleft,
)
vline!(p, [V_eq1]; ls = :dash, color = :red, label = "PE₁ ($(round(V_eq1, digits = 1)) mL)")
vline!(p, [V_eq2]; ls = :dash, color = :blue, label = "PE₂ ($(round(V_eq2, digits = 1)) mL)")
hline!(
    p, [pKa1]; ls = :dot, color = :orange,
    label = "pKₐ₁ = $(round(pKa1, digits = 2)) (from the database)",
)
hline!(
    p, [pKa2]; ls = :dot, color = :green,
    label = "pKₐ₂ = $(round(pKa2, digits = 2)) (from the database)",
)
display(p)
