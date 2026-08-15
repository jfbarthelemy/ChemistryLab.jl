# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

module OptimaSolverExt

using ChemistryLab
import ChemistryLab:
    EquilibriumProblem,
    EquilibriumSolver,
    ChemicalState,
    _build_params,
    _build_n0,
    _solution_transform,
    _update_derived!
using OptimaSolver: OptimaOptimizer
using SciMLBase
using LinearAlgebra: dot, mul!
using DynamicQuantities
using ForwardDiff

# ── OptimizationProblem helpers (NoAD — OptimaOptimizer handles gradients) ────

# The element-conservation matrix and vector are passed through the parameters.
# Without them `OptimaSolver` falls back to rebuilding `A` by finite differences
# on the constraint function, which caps the achievable feasibility at ~1e-6
# whatever tolerance is requested — `A` is known exactly, so it is handed over.
# Only in the linear parameterization: in log space the constraint is
# A·exp(x) = b, which is not linear in the optimization variables, so handing
# over `A` would be wrong there.
_with_constraints(q::NamedTuple, A, b) = merge(q, (A = A, b = b))

function _build_optima_opt_prob(ep::EquilibriumProblem, μ, ::Val{:linear})
    f_gibbs(x, q) = dot(x, μ(x, q))
    cons!(res, x, _) = mul!(res, ep.A, x) .-= ep.b
    optf = SciMLBase.OptimizationFunction{true}(f_gibbs; cons = cons!)
    return SciMLBase.OptimizationProblem(
        optf, ep.u0, _with_constraints(ep.p, ep.A, ep.b);
        lb = ep.lb, ub = ep.ub,
        lcons = zeros(size(ep.A, 1)),
        ucons = zeros(size(ep.A, 1)),
    )
end

function _build_optima_opt_prob(ep::EquilibriumProblem, μ, ::Val{:log})
    f_gibbs(x, q) = (n = exp.(x); dot(n, μ(n, q)))
    cons!(res, x, _) = (n = exp.(x); mul!(res, ep.A, n); res .-= ep.b)
    optf = SciMLBase.OptimizationFunction{true}(f_gibbs; cons = cons!)
    return SciMLBase.OptimizationProblem(
        optf, log.(ep.u0), ep.p;
        lb = log.(ep.lb), ub = log.(ep.ub),
        lcons = zeros(size(ep.A, 1)),
        ucons = zeros(size(ep.A, 1)),
    )
end

# ── solve(EquilibriumSolver{OptimaOptimizer}, ChemicalState) ──────────────────

"""
    SciMLBase.solve(esolver::EquilibriumSolver{F,<:OptimaOptimizer,V},
                   state::ChemicalState; ϵ=1e-16) -> ChemicalState

Solve a chemical equilibrium problem using an `OptimaOptimizer` solver.
Loaded automatically when `using OptimaSolver` is active.
"""
function SciMLBase.solve(
        esolver::EquilibriumSolver{F, <:OptimaOptimizer, V},
        state::ChemicalState;
        ϵ::Float64 = 1.0e-16,
        b = nothing,
    ) where {F, V}
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
    opt_prob = _build_optima_opt_prob(prob, esolver.μ, esolver.variable_space)

    sol = ChemistryLab._check_converged(
        SciMLBase.solve(opt_prob, esolver.solver; esolver.kwargs...),
        "equilibrium solve",
    )
    transform = _solution_transform(esolver.variable_space)

    state_eq = copy(state)
    for (i, nᵢ) in enumerate(sol.u)
        state_eq.n[i] = max(transform(nᵢ), ϵ) * u"mol"
    end
    _update_derived!(state_eq)

    return state_eq
end

# ── __init__: register default solver (high priority — always overrides) ──────
#
# Measured on the cement equilibria this package targets, once the exact
# conservation matrix is handed over (see `_with_constraints` above):
# OptimaSolver reaches a 4e-14 element balance and is 3 to 26 times faster than
# Ipopt. Pass a solver explicitly to override the choice.

_default_optima_solver() = OptimaOptimizer()

function __init__()
    return ChemistryLab._DEFAULT_SOLVER_FACTORY[] = _default_optima_solver
end

end # module OptimaSolverExt
