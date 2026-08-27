```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "ChemistryLab.jl"
  text: "Chemistry you can script"
  tagline: Formulas, species, reactions, thermodynamic databases and equilibrium — programmatic and precise, for geochemistry and cement chemistry.
  image:
    src: /logo.png
    alt: ChemistryLab
  actions:
    - theme: brand
      text: Get started
      link: /quickstart
    - theme: alt
      text: Tutorial
      link: /man/formula_manipulation
    - theme: alt
      text: API
      link: /api/thermo_functions
    - theme: alt
      text: View on GitHub
      link: https://github.com/MicroPoroChemoMechanics/ChemistryLab.jl

features:
  - icon: 🧬
    title: Formulas and species
    details: Parse and manipulate chemical formulas — charges, subscripts, rational stoichiometry — then build Species and cement-specific CemSpecies and attach properties to them.
    link: /man/formula_manipulation
  - icon: ⚗️
    title: Reactions and stoichiometry
    details: Canonical stoichiometric matrices, converted to Reaction objects that can be combined, simplified and transformed.
    link: /man/reactions
  - icon: 🗄️
    title: Thermodynamic databases
    details: Read, merge and write ThermoFun and Cemdata sources, and pull species and reactions out of them programmatically.
    link: /man/databases
  - icon: ⚖️
    title: Equilibrium
    details: Gibbs-energy minimization with built-in activity models — dilute solution, HKF, Davies — and an extension API for your own.
    link: /man/equilibrium
  - icon: 📊
    title: Worked examples
    details: Runnable end to end — carbonate speciation, cement carbonation, acid-base titrations, clinker dissolution, Bogue calculation.
    link: /examples/co2_carbonate_system
  - icon: 📖
    title: API reference
    details: Every exported function, grouped by topic, from parsing tools to thermodynamic models.
    link: /api/thermo_functions
---
```

## What it does

`ChemistryLab` is a Julia toolkit for parsing and manipulating chemical
formulas, creating species, building stoichiometric matrices and reactions,
importing thermodynamic data from ThermoFun and Cemdata, and solving
equilibrium problems by Gibbs-energy minimization.

It is written for scientists, engineers and scientific programmers who want
**programmatic** control over chemical species and thermodynamic workflows:
chemists, geochemists, civil-materials researchers working on cement chemistry,
and anyone who needs precise, scriptable handling of formulas, reactions and
databases.

Temperature-dependent sweeps and speciation diagrams come with it, and the
whole surface is designed to be driven from a script rather than a dialog box.

## A first look

```julia
using ChemistryLab

Species("CaCO3")            # a species from its formula
Formula("SO4-2")            # charges and subscripts parsed
```

Continue with [Getting started](@ref) for installation, a worked
speciation problem and a first equilibrium solve.
