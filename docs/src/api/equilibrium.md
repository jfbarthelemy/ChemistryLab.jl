# Equilibrium

```@index
Pages = ["equilibrium.md"]
```

## Activity models

```@autodocs
Modules = [ChemistryLab]
Pages = ["equilibrium/activities.jl"]
```

## Solid solutions

```@autodocs
Modules = [ChemistryLab]
Pages = ["chemical_structs/solid_solutions.jl"]
```

## Certified equilibrium

Newton on the KKT system in element-potential space, and the certificate that
proves a composition optimal. See
[Proving that an answer is the answer](@ref) for the derivation, and
[`equilibrate_certified`](@ref) for the route `equilibrate` takes by default.

```@autodocs
Modules = [ChemistryLab]
Pages = ["equilibrium/dual_solver.jl", "equilibrium/certified.jl"]
```

## Constraints

What is held fixed while the Gibbs energy is minimized. See
[Constraints other than fixed T and P](@ref sec-equilibrium-constraints).

```@autodocs
Modules = [ChemistryLab]
Pages = ["equilibrium/constraints.jl"]
```

## Problem and solver

```@autodocs
Modules = [ChemistryLab]
Pages = ["equilibrium/equilibrium_problems.jl", "equilibrium/equilibrium_solver.jl"]
```
