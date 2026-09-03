# Kinetics API

Chemical kinetics module: rate models, kinetic reactions, ODE problem setup, solvers, and calorimetry.

## Rate constants and models

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/rate_models.jl"]
```

## Kinetic reactions and surface area

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/kinetics_reactions.jl"]
```

## Kinetics problem

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/kinetics_problems.jl"]
```

## Kinetics solver

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/kinetics_solver.jl"]
```

## Calorimetry

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/calorimetry.jl"]
```

## Post-processing a solution

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/kinetics_postprocessing.jl"]
```

## The implicit kinetic step

One fully implicit problem per step, with the reaction extents as unknowns of the
same Gibbs minimization. See
[Which route: two ways to advance in time](@ref sec-kinetics-routes).

```@autodocs
Modules = [ChemistryLab]
Pages   = ["kinetics/implicit_step.jl"]
```
