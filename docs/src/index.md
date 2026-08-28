```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "ChemistryLab.jl"
  text: "Chemistry you can script"
  tagline: Formulas, species, reactions, thermodynamic databases and equilibrium — for aqueous geochemistry and cement chemistry.
  image:
    src: /logo.png
    alt: ChemistryLab
  actions:
    - theme: brand
      text: Get started
      link: /quickstart
    - theme: alt
      text: Tutorials
      link: /tutorials/formula_manipulation
    - theme: alt
      text: Examples
      link: /examples/co2_carbonate_system
    - theme: alt
      text: View on GitHub
      link: https://github.com/MicroPoroChemoMechanics/ChemistryLab.jl

features:
  - icon: 🚀
    title: Getting started
    details: Install the package, compute the solubility constant of a reaction, then solve a first equilibrium.
    link: /quickstart
  - icon: 🎓
    title: Tutorials
    details: The library topic by topic — formulas, species, databases, stoichiometric matrices, reactions, equilibrium, kinetics, coupling.
    link: /tutorials/formula_manipulation
  - icon: 🧪
    title: Examples
    details: Complete problems solved end to end — carbonate speciation, cement carbonation, acid-base titrations, clinker dissolution and hydration.
    link: /examples/co2_carbonate_system
  - icon: 🗄️
    title: Databases
    details: Read, merge and write ThermoFun and Cemdata sources, and pull species and reactions out of them programmatically.
    link: /tutorials/databases
  - icon: ⚖️
    title: Equilibrium and kinetics
    details: Gibbs-energy minimization with dilute-solution, HKF and Davies activity models, and time-resolved dissolution-precipitation.
    link: /tutorials/equilibrium
  - icon: 📖
    title: API reference
    details: Every exported function, grouped by topic, from parsing tools to thermodynamic models.
    link: /api/thermo_functions
---
```

## What it does

`ChemistryLab` handles chemical formulas, species and reactions as first-class
objects, reads thermodynamic data from ThermoFun and Cemdata, and solves
equilibrium by Gibbs-energy minimization. Kinetics and chemo-mechanical coupling
build on the same objects.

It is written for work that has to be reproducible and scripted: aqueous
geochemistry, cement chemistry, and any problem where speciation, a database and
a solver have to be driven from code rather than from a dialog box.

## A first calculation

Calcite in water. The species come from a thermodynamic database, four of them
are declared as the independent basis, and the equilibrium state follows from a
Gibbs-energy minimization:

```@example home
using ChemistryLab, DynamicQuantities
using OptimaSolver          # the default equilibrium back-end

species = speciation(build_species(datapath("cemdata18-thermofun.json")),
                     split("Cal H2O@ CO2");
                     aggregate_state = [AS_AQUEOUS],
                     exclude_species = split("H2@ O2@ CH4@"))
byname = Dict(symbol(s) => s for s in species)

system = ChemicalSystem(collect(values(byname)),
                        [byname[s] for s in split("H2O@ H+ CO3-2 Ca+2")])

state = ChemicalState(system)
set_quantity!(state, "Cal", 1e-3u"mol")     # 1 mmol of calcite
set_quantity!(state, "H2O@", 1.0u"kg")      # in 1 kg of water
V = volume(state)
set_quantity!(state, "H+",  1e-4u"mol/L" * V.liquid)
set_quantity!(state, "OH-", 1e-10u"mol/L" * V.liquid)

equilibrated = equilibrate(state)
nothing # hide
```

The state carries everything derived from it — pH, phase volumes, individual
amounts:

```@example home
using Printf
@printf("pH                = %.2f\n", pH(equilibrated))
@printf("dissolved Ca(2+)  = %.3e mol\n", ustrip(moles(equilibrated, "Ca+2")))
@printf("remaining calcite = %.3e mol\n", ustrip(moles(equilibrated, "Cal")))
```

[Getting started](@ref) takes the same problem more slowly, and shows how the
solubility constant is obtained analytically before any solver is involved.
