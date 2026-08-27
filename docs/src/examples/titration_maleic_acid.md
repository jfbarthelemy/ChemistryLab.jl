# [Titration of Maleic Acid by NaOH](@id sec-titration-maleic)

This example simulates the **potentiometric titration** of maleic acid (H₂A, a diprotic weak acid)
by sodium hydroxide (NaOH, a strong base), both at 0.1 M.
Maleic acid has two well-separated dissociation constants (pKₐ₁ = 1.92, pKₐ₂ = 6.27),
producing two distinct inflection points on the titration curve.

---

## System setup

The species are distributed across two SLOP98 databases:
the inorganic database provides H₂O, H⁺, OH⁻ and Na⁺;
the organic database provides maleic acid and its conjugate bases.

| Symbol | Species | Phase |
|:-------|:--------|:------|
| `MalH2@` | H₂A — maleic acid | aqueous solute |
| `MalH-`  | HA⁻ — hydrogen maleate | aqueous solute |
| `Mal-2`  | A²⁻ — maleate | aqueous solute |
| `Na+`    | Na⁺ | aqueous solute |
| `H+`     | H⁺ | aqueous solute |
| `OH-`    | OH⁻ | aqueous solute |
| `H2O@`   | H₂O | aqueous solvent |

```@setup titration_setup
using ChemistryLab
using DynamicQuantities

substances_inorg = build_species(datapath("slop98-inorganic-thermofun.json"))
substances_org   = build_species(datapath("slop98-organic-thermofun.json"))

dict_all_species = merge(Dict(symbol(s) => s for s in substances_inorg), Dict(symbol(s) => s for s in substances_org))
species = [dict_all_species[s] for s in split("H2O@ Na+ NaOH@ H+ OH- MalH2@ MalH- Mal-2")]

cs = ChemicalSystem(species, ["H2O@", "H+", "Mal-2", "Na+", "Zz"])
```

```julia
using ChemistryLab
using DynamicQuantities

substances_inorg = build_species(datapath("slop98-inorganic-thermofun.json"))
substances_org   = build_species(datapath("slop98-organic-thermofun.json"))

dict_all_species = merge(Dict(symbol(s) => s for s in substances_inorg), Dict(symbol(s) => s for s in substances_org))
species = [dict_all_species[s] for s in split("H2O@ Na+ NaOH@ H+ OH- MalH2@ MalH- Mal-2")]

cs = ChemicalSystem(species, ["H2O@", "H+", "Mal-2", "Na+", "Zz"])
```

Build the [`EquilibriumSolver`](@ref) once — it is reused for each of the 66 titration points:

```@setup titration_setup
using Optimization, OptimizationIpopt

solver = EquilibriumSolver(
    cs,
    DiluteSolutionModel(),
    IpoptOptimizer(
        mu_strategy = "adaptive",
    );
    variable_space = Val(:linear),
    abstol  = 1e-8,
    reltol  = 1e-8,
    maxiters = 100,
    verbose = false,
)
```

```julia
using Optimization, OptimizationIpopt

solver = EquilibriumSolver(
    cs,
    DiluteSolutionModel(),
    IpoptOptimizer(
        mu_strategy = "adaptive",
    );
    variable_space = Val(:linear),
    abstol  = 1e-8,
    reltol  = 1e-8,
    maxiters = 100,
    verbose = false,
)
```

---

## Running the titration

At each titration point the total composition is reset from scratch and the equilibrium is recomputed.
The conservation constraint (total moles of each element) is automatically enforced by the solver.

```@example titration_setup
V_acid = 100e-3   # volume of acid solution, L
c_acid = 0.1     # maleic acid concentration, mol/L
c_base = 2     # NaOH concentration, mol/L
n_H2A  = V_acid * c_acid   # total moles of H₂A = 2.5 mmol

ρ_water = 1.   # kg/L

volumes_NaOH = range(0, 15; length = 101)   # mL
pH_vals = Float64[]

s = ChemicalState(cs)
for V_mL in volumes_NaOH
    V_NaOH  = V_mL * 1e-3           # L
    n_NaOH    = c_base * V_NaOH        # mol of NaOH (= mol of Na⁺ added)
    V_total = V_acid + V_NaOH        # total volume, L

    set_quantity!(s, "MalH2@", n_H2A   * u"mol")
    set_quantity!(s, "NaOH@", n_NaOH * u"mol")
    set_quantity!(s, "H2O@",   ρ_water * V_total * u"kg")

    V_liq = volume(s).liquid
    set_quantity!(s, "H+",  1e-7u"mol/L" * V_liq)   # pH-neutral seed
    set_quantity!(s, "OH-", 1e-7u"mol/L" * V_liq)

    s_eq = solve(solver, s)
    push!(pH_vals, pH(s_eq))
end

println("pH at V = 0 mL (pure acid)        : ", round(pH_vals[1],  digits = 2))
println("pH at V = 15 mL  (excess NaOH)    : ", round(pH_vals[100], digits = 2))

pprint(cs.SM)
```

---

## Titration curve

```julia
using Plots

pKa1 = 1.92
pKa2 = 6.27
V_eq1 = n_H2A / c_base * 1e3    # first equivalence point  = 25 mL
V_eq2 = 2 * n_H2A / c_base * 1e3  # second equivalence point = 50 mL

p = plot(
    collect(volumes_NaOH), pH_vals;
    xlabel     = "V(NaOH) (mL)",
    ylabel     = "pH",
    label      = "Titration curve",
    linewidth  = 2,
    marker     = :circle,
    markersize = 3,
    color      = :steelblue,
    title      = "Titration of maleic acid (0.1 M) by NaOH (2 M)",
    ylims      = (0, 14),
    legend     = :topleft,
)
vline!(p, [V_eq1]; linestyle = :dash, color = :red,    label = "PE₁ ($(round(V_eq1,digits = 1)) mL)")
vline!(p, [V_eq2]; linestyle = :dash, color = :blue,   label = "PE₂ ($(round(V_eq2,digits = 1)) mL)")
hline!(p, [pKa1];  linestyle = :dot,  color = :orange, label = "pKₐ₁ = $pKa1")
hline!(p, [pKa2];  linestyle = :dot,  color = :green,  label = "pKₐ₂ = $pKa2")
```

![Maleitric titration curve](../assets/maleic_titration.png)

---

## Analysis

The titration curve shows five characteristic zones:

| Zone | V(NaOH) | Dominant species | pH |
|:-----|:--------|:-----------------|:---|
| Initial state | 0 mL | H₂A | Low, controlled by pKₐ₁ |
| First buffer | 0–5 mL | H₂A / HA⁻ | ≈ pKₐ₁ = 1.92 at V = 2.5 mL |
| First equivalence point (PE₁) | 5 mL | HA⁻ | First inflection |
| Second buffer | 5–10 mL | HA⁻ / A²⁻ | ≈ pKₐ₂ = 6.27 at V = 7.5 mL |
| Second equivalence point (PE₂) | 10 mL | A²⁻ | Second inflection |
| Excess base | > 10 mL | A²⁻ + OH⁻ | Controlled by excess NaOH |

- **V = 0 mL** — The pH is low, determined mainly by the first dissociation (pKₐ₁ = 1.92).
- **V = 5 mL (PE₁)** — The first proton is fully neutralized; the dominant species transitions from H₂A to HA⁻.
- **V = 2.5 mL (half-equivalence 1)** — pH ≈ pKₐ₁ = 1.92 (Henderson–Hasselbalch condition).
- **5 mL < V < 10 mL** — The HA⁻/A²⁻ couple acts as a buffer; at V = 7.5 mL, pH ≈ pKₐ₂ = 6.27.
- **V = 10 mL (PE₂)** — The second proton is fully neutralized; the dominant species is A²⁻.
- **V > 10 mL** — pH rises steeply, controlled by the concentration of free OH⁻ from excess NaOH.

!!! note "Δ pKₐ and resolution"
    The two dissociation constants of maleic acid are well separated (Δ pKₐ ≈ 4.35).
    This large gap produces two clearly resolved inflection points, making it an ideal
    model compound for potentiometric titration analysis.
