```@meta
CurrentModule = ChemistryLab
```

# ChemistryLab

Welcome to the ChemistryLab documentation. This set of pages guides you through the core features of ChemistryLab — a Julia toolkit for parsing and manipulating chemical formulas, creating species, building stoichiometric matrices and reactions, and importing thermodynamic data (ThermoFun / Cemdata).

## Audience

This documentation is intended for scientists, engineers and scientific programmers who want programmatic control over chemical species and simple thermodynamic workflows in Julia. It suits chemists, geochemists, civil-materials researchers (cement chemistry), and anyone who needs precise, scriptable handling of formulas, reactions and databases.

## Installation

To install ChemistryLab.jl, use the Julia package manager:

- From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add ChemistryLab
```

- Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("ChemistryLab")
```

## Citation

If you use ChemistryLab in your work, please cite the following:

```bibtex
@software{chemistrylab_jl,
  author       = {Barthélémy, Jean-François and
                  Soive, Anthony},
  title        = {ChemistryLab.jl: Numerical laboratory for
                   computational chemistry},
  doi          = {10.5281/zenodo.17756074},
  url          = {https://doi.org/10.5281/zenodo.17756074},
}
```

## Features

- Parse and manipulate chemical formulas (supporting charges, subscripts/superscripts and rational stoichiometry).
- Create and inspect `Species` and `CemSpecies` objects and attach thermodynamic or auxiliary properties.
- Build canonical stoichiometric matrices and convert them to `Reaction` objects.
- Read, merge and write ThermoFun / Cemdata sources and extract species and reactions programmatically.
- Combine, simplify and transform reactions using the `Reaction` API.
- Solve thermodynamic equilibrium problems by Gibbs energy minimization (`equilibrate`, `EquilibriumSolver`).
- Built-in activity models (`DiluteSolutionModel`, `HKFActivityModel`, `DaviesActivityModel`) and extension API for custom models.
- Temperature-dependent thermodynamic sweeps and speciation diagrams.

## Documentation structure

- `introduction` — this page: goals, audience, and how to run the examples.
- `quickstart` — minimal examples to get you started quickly.
- `tutorial` — basic examples to get you started more deeply in the library.
  - `formula_manipulation` — parsing, conversions (Phreeqc ⇄ Unicode), formatting and arithmetic on `Formula` objects.
  - `species` — constructing `Species` and `CemSpecies`, managing properties and basic queries.
  - `cement_species` — cement-specific species and oxide nomenclature.
  - `stoich_matrices` — building canonical stoichiometric matrices and converting them to reaction objects.
  - `reactions` — parsing, combining and simplifying chemical reactions.
  - `databases` — reading ThermoFun JSON and Cemdata `.dat` files, merging reactions and exporting results.
  - `equilibrium` — setting up and solving chemical equilibrium problems (`ChemicalSystem`, `ChemicalState`, `equilibrate`, activity models, temperature sweeps).
  - `advanced` — advanced patterns: formula arithmetic, reaction algebra, database operations.
- `examples/` — runnable worked examples (also built to HTML under `docs/build/examples/`):
  - `co2_carbonate_system` — CO₂ dissolution in water, pKₐ extraction and carbonate speciation diagram.
  - `cement_carbonation` — progressive carbonation of hydrated cement paste (pH, portlandite, calcite).
  - `titration_acetic_acid`, `titration_maleic_acid` — acid–base titration curves.
  - `simplified_clinker_dissolution` — clinker phase dissolution.
  - `bogue_calculation` — Bogue calculation for cement clinker.

## Try it locally

The repository contains a `docs` folder with Documenter sources and a `make.jl` script. To build the documentation locally from the project root, activate the docs environment and run the make script:

```powershell
# from the project root (Windows PowerShell)
julia --project=docs docs/make.jl
```

Alternatively, ensure the package environment is prepared and run Documenter through the package project:

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate(); include(\"docs/make.jl\")"
```

After the build completes, open `docs/build/index.html` in a browser to read the rendered tutorial and follow the example pages.

## Quick tips

- In the REPL try small calls like `using ChemistryLab; Species("CaCO3")` and `Formula("SO4-2")` to explore parsing behavior interactively.
- Start with `docs/src/examples/example_stoich_matrix.md` to see a concise, runnable example converting a stoichiometric matrix into reactions.
- For equilibrium calculations, see `docs/src/man/equilibrium.md` for the minimal workflow, then explore the `co2_carbonate_system` and `cement_carbonation` examples.
- If you plan to work with ThermoFun/Cemdata sources, run the examples in `docs/src/man/databases.md` after placing the required `.json`/`.dat` data files in the `data/` directory.

## Next steps

You can see the `examples` section for more advanced runnable examples and small worked problems, including CO₂ dissolution, carbonate speciation, cement carbonation, and clinker dissolution.

Happy exploring — this tutorial aims to be practical and runnable, so please tell me which example you want expanded into a fully reproducible script.
