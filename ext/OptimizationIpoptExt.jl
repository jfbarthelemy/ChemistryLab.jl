# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

module OptimizationIpoptExt

using ChemistryLab
import ChemistryLab:
    EquilibriumProblem,
    EquilibriumSolver,
    ChemicalState,
    AbstractActivityModel,
    DiluteSolutionModel,
    _build_params,
    _build_n0,
    _solution_transform,
    _update_derived!
using Optimization
using OptimizationIpopt
using SciMLBase
using LinearAlgebra: dot, mul!
using DynamicQuantities
using ForwardDiff

# ── OptimizationProblem conversions ──────────────────────────────────────────

"""
    SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:linear}; kwargs...)

Convert an `EquilibriumProblem` to an `OptimizationProblem` in linear variable space.
The objective minimizes G = n⋅μ(n) subject to the mass conservation constraint A*n = b.
"""
function SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:linear}; kwargs...)
    Gibbs_energy(x, p) = dot(x, ep.μ(x, p))
    conservation_constraints(res, x, _) = mul!(res, ep.A, x) .-= ep.b

    optf = OptimizationFunction(
        Gibbs_energy,
        Optimization.AutoForwardDiff();
        cons = conservation_constraints,
    )

    return OptimizationProblem(
        optf,
        ep.u0,
        ep.p;
        lb = ep.lb,
        ub = ep.ub,
        lcons = zeros(size(ep.A, 1)),
        ucons = zeros(size(ep.A, 1)),
        kwargs...,
    )
end

"""
    SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:log}; kwargs...)

Convert an `EquilibriumProblem` to an `OptimizationProblem` in log variable space.
Solves for x = log(n), which automatically enforces positivity and is more robust
for systems spanning many orders of magnitude.
"""
function SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:log}; kwargs...)
    Gibbs_energy_log(x, p) = (n = exp.(x); dot(n, ep.μ(n, p)))
    conservation_constraints_log(res, x, _) = (n = exp.(x); mul!(res, ep.A, n); res .-= ep.b)

    optf = OptimizationFunction(
        Gibbs_energy_log,
        Optimization.AutoForwardDiff();
        cons = conservation_constraints_log,
    )

    return OptimizationProblem(
        optf,
        log.(ep.u0),
        ep.p;
        lb = log.(ep.lb),
        ub = log.(ep.ub),
        lcons = zeros(size(ep.A, 1)),
        ucons = zeros(size(ep.A, 1)),
        kwargs...,
    )
end

# ── solve(EquilibriumProblem, solver) ─────────────────────────────────────────

"""
    SciMLBase.solve(ep::EquilibriumProblem, solver; variable_space=Val(:linear), kwargs...)

Solve an `EquilibriumProblem`. The solution is transformed back to physical (mole) space.
"""
function SciMLBase.solve(
        ep::EquilibriumProblem,
        solver;
        variable_space = Val(:linear),
        kwargs...,
    )
    opt_prob = SciMLBase.OptimizationProblem(ep, variable_space)
    sol = ChemistryLab._check_converged(
        SciMLBase.solve(opt_prob, solver; kwargs...),
        "equilibrium solve",
    )
    transform = _solution_transform(variable_space)
    sol.u .= transform.(sol.u)
    return sol
end

# ── solve(EquilibriumSolver, ChemicalState) ───────────────────────────────────

"""
    SciMLBase.solve(esolver::EquilibriumSolver, state::ChemicalState; ϵ=1e-16) -> ChemicalState

Solve a chemical equilibrium problem from an initial `ChemicalState`.

The conservation matrix `A` is taken from `state.system.SM.A`.
The initial mole vector `n0` and thermodynamic parameters `ΔₐG⁰/RT`
are extracted from `state` at its current `T` and `P`.

Returns a new `ChemicalState` with equilibrium mole amounts,
sharing the same `ChemicalSystem` as the input.

# Arguments

  - `esolver`: the `EquilibriumSolver` to use.
  - `state`: initial state — defines `T`, `P`, and initial composition.
  - `ϵ`: regularization floor (default: `1e-16`).

# Examples
```julia
solver = EquilibriumSolver(cs, DiluteSolutionModel(), IpoptOptimizer();
                           variable_space=Val(:log), abstol=1e-10)
state0 = ChemicalState(cs, n0; T=298.15u"K", P=1u"bar")
state_eq = solve(solver, state0)
```
"""
function SciMLBase.solve(
        esolver::EquilibriumSolver,
        state::ChemicalState;
        ϵ::Float64 = 1.0e-16,
        b = nothing,
    )
    # A composition carrying dual numbers takes the implicit-function route:
    # primal solve, then sensitivities from the optimality conditions. No solver
    # is asked to iterate on dual numbers.
    if eltype(state.n) <: DynamicQuantities.AbstractQuantity{<:ForwardDiff.Dual}
        return ChemistryLab._solve_dual(esolver, state, ϵ; b = b)
    end

    n0 = max.(_build_n0(state), ϵ)
    p = _build_params(state; ϵ = ϵ)

    # `b` given explicitly is Leal's φ(b): minimize G subject to A n = b, with
    # `state` supplying only the starting guess and the T, P conditions. The
    # element totals then come from the caller — the ODE state of a kinetics
    # run — instead of being derived from a composition that may not carry them.
    prob = isnothing(b) ?
        EquilibriumProblem(state.system.SM.A, esolver.μ, n0; p = p) :
        EquilibriumProblem(state.system.SM.A, esolver.μ, n0; b = collect(b), p = p)
    sol = SciMLBase.solve(prob, esolver.solver; variable_space = esolver.variable_space, esolver.kwargs...)

    state_eq = copy(state)
    for (i, nᵢ) in enumerate(sol.u)
        state_eq.n[i] = max(nᵢ, ϵ) * u"mol"
    end
    _update_derived!(state_eq)

    return state_eq
end

# ── Default Ipopt solver factory ──────────────────────────────────────────────

function _default_ipopt_solver()
    return IpoptOptimizer(
        acceptable_tol = 1.0e-12,
        dual_inf_tol = 1.0e-12,
        acceptable_iter = 1000,
        constr_viol_tol = 1.0e-12,
        warm_start_init_point = "no",
    )
end

# ── __init__: register default solver (low priority) ─────────────────────────

function __init__()
    ChemistryLab.register_solver_factory!(_default_ipopt_solver)
    return if isnothing(ChemistryLab._DEFAULT_SOLVER_FACTORY[])
        ChemistryLab._DEFAULT_SOLVER_FACTORY[] = _default_ipopt_solver
    end
end

end # module OptimizationIpoptExt
