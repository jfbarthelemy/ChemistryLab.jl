# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using SciMLBase

# ── Default solver registry ───────────────────────────────────────────────────

"""
    _DEFAULT_KINETICS_SOLVER_FACTORY

`Ref` holding a zero-argument factory `() -> KineticsSolver` for the default
ODE solver used by `integrate(kp)` (no explicit solver argument).

Priority:
- `KineticsOrdinaryDiffEqExt.__init__` registers `Rodas5P()` on load.
- Only one registration; extensions override it as needed.

Users can override with:
```julia
ChemistryLab._DEFAULT_KINETICS_SOLVER_FACTORY[] = () -> KineticsSolver(ode_solver=MyAlg())
```
"""
const _DEFAULT_KINETICS_SOLVER_FACTORY = Ref{Union{Nothing, Function}}(nothing)

# ── KineticsSolver ────────────────────────────────────────────────────────────

"""
    struct KineticsSolver{ODE_S, ES}

Bundles the ODE algorithm and optional equilibrium solver used by [`integrate`](@ref).
Construct once, reuse across multiple [`KineticsProblem`](@ref) instances.

# Fields

  - `ode_solver`: any `SciMLBase.AbstractODEAlgorithm` (e.g. `Rodas5P()` from
    `OrdinaryDiffEq`), or `:auto` to let `OrdinaryDiffEq` choose — its default
    polyalgorithm detects stiffness at run time and switches. `nothing` keeps
    `Rodas5P()`, which stays the default because on the problems here `:auto`
    lands on the same stiff method and returns the same answers: switching would
    change nothing while removing a reproducible choice. Use `nothing` before
    `OrdinaryDiffEq` is loaded; an error will be raised at solve time.
  - `equilibrium_solver`: optional [`EquilibriumSolver`](@ref) to re-equilibrate
    aqueous speciation at each ODE evaluation. When `nothing`, the kinetic minerals
    evolve without re-speciation (faster, less accurate).
  - `kwargs`: keyword arguments forwarded to `DifferentialEquations.solve`
    (e.g. `reltol`, `abstol`, `saveat`, `maxiters`).

# Examples

```julia
using OrdinaryDiffEq          # activates KineticsOrdinaryDiffEqExt
using Optimization, OptimizationIpopt   # needed for equilibrium_solver

es = EquilibriumSolver(cs, HKFActivityModel(), IpoptOptimizer())
ks = KineticsSolver(; ode_solver=Rodas5P(), equilibrium_solver=es,
                     reltol=1e-8, abstol=1e-10)
sol = integrate(kp, ks)
```
"""
struct KineticsSolver{ODE_S, ES}
    ode_solver::ODE_S
    equilibrium_solver::ES
    kwargs::Base.Pairs
end

"""
    KineticsSolver(; ode_solver=nothing, equilibrium_solver=nothing, kwargs...) -> KineticsSolver

Construct a [`KineticsSolver`](@ref).

`kwargs` are forwarded to `DifferentialEquations.solve` (e.g. `reltol=1e-8`,
`abstol=1e-10`, `saveat=0:60:3600`).
"""
function KineticsSolver(;
        ode_solver = nothing,
        equilibrium_solver = nothing,
        kwargs...,
    )
    return KineticsSolver{typeof(ode_solver), typeof(equilibrium_solver)}(
        ode_solver, equilibrium_solver, kwargs
    )
end

"""
    integrate(kp::KineticsProblem; kwargs...) -> ODESolution

Shortcut that uses the default solver registered by [`_DEFAULT_KINETICS_SOLVER_FACTORY`](@ref).

Requires `OrdinaryDiffEq` to be loaded (which sets the default to `Rodas5P()`).

```julia
using OrdinaryDiffEq
sol = integrate(kp)              # uses Rodas5P() by default
sol = integrate(kp; reltol=1e-6) # forward kwargs to the ODE solver
```
"""
function integrate(kp::KineticsProblem; kwargs...)
    factory = _DEFAULT_KINETICS_SOLVER_FACTORY[]
    isnothing(factory) && error(
        "integrate requires OrdinaryDiffEq.jl to be loaded.\n" *
            "Add `using OrdinaryDiffEq` before calling this function.",
    )
    ks = factory()
    return integrate(kp, ks; kwargs...)
end
