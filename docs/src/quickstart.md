# Getting started

This quickstart shows a few common, minimal examples to get you productive with ChemistryLab. It demonstrates loading species from a database, building reactions and solving a thermodynamic equilibrium problem.

## Simplified example

Let us start with a minimal example in which we compute the thermodynamic properties of a reaction. As a first illustration, we consider the equilibrium of calcite in water. This equilibrium can be written as:

$\ce{CaCO3 <=> Ca^2+ + CO3^2-}$

It is possible to calculate the thermodynamic properties of the reaction, in particular the solubility constant of the reaction ($\ln K$) which is related to the Gibbs free energy of the reaction ($\Delta_r G^°$). This solubility constant is a function of temperature and the calculation is performed at a reference temperature of 298 K and at a pressure of 1 Atm using the following equation:

$\Delta_r G^° = - RT\ln K$

where $\Delta_r G^°$ is deduced from the Gibbs energies of formation ($\Delta_f {G_i}^°$) of the other chemical species involved in the reaction:

$\Delta_r G^° = \sum_i \nu_i \Delta_f {G_i}^°$

------------------------

To do this, we load the species from one of the databases integrated into ChemistryLab, filter those relevant to the calcite–water system, then build the stoichiometric matrix and derive the reactions.

In this example, the database is [cemdata](https://www.empa.ch/web/s308/thermodynamic-data). The `.json` file is included in ChemistryLab but is a copy of a file which can be found in [ThermoHub](https://github.com/thermohub).

`build_species` reads the database file and returns a `Vector{Species}` with compiled thermodynamic functions. `speciation` then filters this list to the species whose atomic composition is a subset of the seed atoms (here Ca, C, H and O from `Cal`, `H2O@` and `CO2`):

```@example from_scratch
using ChemistryLab
using DynamicQuantities # for unit management

all_species = build_species(datapath("cemdata18-thermofun.json"))
species_calcite = speciation(all_species, split("Cal H2O@ CO2");
                             aggregate_state=[AS_AQUEOUS],
                             exclude_species=split("H2@ O2@ CH4@"))
dict_species_calcite = Dict(symbol(s) => s for s in species_calcite)
```

During species creation, ChemistryLab calculates the molar mass of the species. It also constructs thermodynamic functions (heat capacity, entropy, enthalpy, and Gibbs free energy of formation) as a function of temperature. The evolution of thermodynamic properties as a function of temperature, such as heat capacity, can thus be easily plotted.

```@example from_scratch
dict_species_calcite["Cal"]
```

```@example from_scratch
using Plots

p1 = plot(xlabel="Temperature [°C]", ylabel="Cp⁰ [J/mol/K]", title="Heat capacity of calcite \nas a function of temperature")
plot!(p1, θ -> dict_species_calcite["Cal"].Cp⁰(T = θ*ua"degC"), 0:0.1:100, label="Cp⁰")
```

Obtaining stoichiometric matrices requires the choice of a species-independent basis.

```@example from_scratch
primaries = [dict_species_calcite[s] for s in split("H2O@ H+ CO3-2 Ca+2")]
SM = StoichMatrix(collect(values(dict_species_calcite)), primaries)
pprint(SM)
```

These stoichiometric matrices thus allow us to write the chemical reactions at work.

```@example from_scratch
list_reactions = reactions(SM)
dict_reactions_calcite = Dict(r.symbol => r for r in list_reactions)
```

Again, when constructing the reactions, the thermodynamic properties of the reactions as a function of temperature are deduced. It is thus possible to see, for example, the expression for the solubility product of calcite for the reaction under study and to plot its evolution.

```@example from_scratch
dict_reactions_calcite["Cal"].logK⁰
```

```@example from_scratch
p2 = plot(xlabel="Temperature [°C]", ylabel="pKs", title="Solubility product (pKs) of calcite \nas a function of temperature")
plot!(p2, θ -> dict_reactions_calcite["Cal"].logK⁰(T = θ*ua"degC"), 0:0.1:100, label="pKs")
```

## Equilibrium solving

The previous section computed thermodynamic properties of reactions analytically. ChemistryLab can go further and solve the full **thermodynamic equilibrium** — that is, find the species amounts that minimize the Gibbs free energy of the system given initial conditions.

Three objects are needed:

| Object | Role |
|--------|------|
| [`ChemicalSystem`](@ref) | Immutable description of the system: species list, primary species, stoichiometric matrices, and derived index maps. Built once and reused. |
| [`ChemicalState`](@ref) | Mutable thermodynamic state: amounts `n` (mol), temperature `T` and pressure `P`. Modified in-place before and after solving. |
| `equilibrate` | Convenience function that wraps a `ChemicalSystem` + `ChemicalState` into an optimization problem and solves it. Returns a new equilibrated `ChemicalState`. |

```@example from_scratch
using Optimization, OptimizationIpopt
using DynamicQuantities

# ChemicalSystem: declare the species and which four are the independent basis
primaries_eq = [dict_species_calcite[s] for s in split("H2O@ H+ CO3-2 Ca+2")]
cs = ChemicalSystem(collect(values(dict_species_calcite)), primaries_eq)

# ChemicalState: set the initial amounts
state = ChemicalState(cs)
set_quantity!(state, "Cal",  1e-3u"mol")   # 1 mmol of calcite
set_quantity!(state, "H2O@", 1.0u"kg")     # 1 kg of water

# Seed H⁺ and OH⁻ at pH 4 (trace amounts to help convergence)
V = volume(state)
set_quantity!(state, "H+",  1e-4u"mol/L" * V.liquid)
set_quantity!(state, "OH-", 1e-10u"mol/L" * V.liquid)

# Solve: find the Gibbs-energy minimum
state_eq = equilibrate(state)
```

```@example from_scratch
println("pH = ", round(pH(state_eq), digits = 2))
```

Derived quantities such as pH, pOH, phase volumes and individual species amounts are all accessible on the returned `ChemicalState`. For a detailed description of the solver options, activity models, and temperature sweeps, see the equilibrium tutorial (`Tutorial → Chemical Equilibrium`).

## Notes and next steps

- The `Formula`, `Species`, `Reaction` and `StoichMatrix` APIs are intentionally small and composable — explore the `docs/src/` pages for detailed examples.
- For equilibrium calculations, see `docs/src/man/equilibrium.md` and the worked examples `co2_carbonate_system` and `cement_carbonation`.
- For cement-specific workflows, use `CemSpecies` and the `databases` utilities to convert between oxide- and atom-based representations.

Now try the `quickstart` examples interactively in the REPL and then follow the next pages of the tutorial for deeper coverage.
