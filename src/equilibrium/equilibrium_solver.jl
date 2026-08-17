# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using ForwardDiff
using LinearAlgebra: pinv
using OrderedCollections
using SciMLBase

# ── EquilibriumSolver ─────────────────────────────────────────────────────────

"""
    struct EquilibriumSolver{F<:Function, S, V<:Val}

Encapsulates all fixed ingredients of a chemical equilibrium calculation:
the potential function, the SciML solver, and the variable space.

Construct once, call repeatedly with different `ChemicalState` inputs.

# Fields

  - `μ`: chemical potential closure `μ(n, p) -> Vector{Float64}`.
  - `solver`: any Optimization.jl-compatible solver (e.g. `IpoptOptimizer()`).
  - `variable_space`: variable space — `Val(:linear)` or `Val(:log)`.
  - `kwargs`: solver keyword arguments forwarded to `solve`.

# Examples
```julia
julia> cs = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("H+";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
           Species("OH-"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
       ]);

julia> solver = EquilibriumSolver(cs, DiluteSolutionModel(), IpoptOptimizer());

julia> solver isa EquilibriumSolver
true
```
"""
struct EquilibriumSolver{F <: Function, S, V <: Val, M <: AbstractActivityModel}
    μ::F                    # potential closure — built once from cs and model
    solver::S               # SciML-compatible solver
    variable_space::V       # Val(:linear) or Val(:log)
    kwargs::Base.Pairs      # forwarded to solve
    model::M                # the activity model μ was built from
end

"""
    EquilibriumSolver(cs, model, solver; variable_space=Val(:linear), kwargs...)

Construct an `EquilibriumSolver` from a `ChemicalSystem`, an activity model,
and a SciML solver.

The potential function is built once at construction time from `cs` and `model`.
Repeated calls to `solve` with different `ChemicalState` inputs reuse it.

# Arguments

  - `cs`: the `ChemicalSystem` defining species and conservation matrix.
  - `model`: an `AbstractActivityModel` (e.g. `DiluteSolutionModel()`).
  - `solver`: any Optimization.jl solver.
  - `variable_space`: `Val(:linear)` (default) or `Val(:log)`.
  - `kwargs...`: forwarded to the underlying `solve` call (tolerances, verbosity...).
"""
function EquilibriumSolver(
        cs::ChemicalSystem,
        model::AbstractActivityModel,
        solver::S;
        variable_space::V = Val(:linear),
        kwargs...,
    ) where {S, V <: Val}
    μ = build_potentials(cs, model)     # built once — captures indices and constants
    # `model` is kept alongside `μ`, which is opaque once compiled. A coupled
    # kinetics run has to rebuild the solver for the equilibrium SUB-system, and
    # without this it had no way to know which activity model to rebuild with —
    # it silently fell back to the problem's own default and ran a dilute solve
    # on a pore solution the user had asked to treat with HKF.
    return EquilibriumSolver{typeof(μ), S, V, typeof(model)}(
        μ, solver, variable_space, kwargs, model
    )
end

"""
    activity_model(solver::EquilibriumSolver) -> AbstractActivityModel

The activity model the solver's potential function was built from.
"""
activity_model(solver::EquilibriumSolver) = solver.model

# ── Internal helper: build p from ChemicalState ───────────────────────────────

"""
    _build_params(state::ChemicalState; ϵ=1e-16) -> NamedTuple

Extract dimensionless parameters from a `ChemicalState`.
`ΔₐG⁰overT` is evaluated at the current `T` and `P` of the state.
Units are stripped — compatible with ForwardDiff dual numbers.

# Returned fields

  - `ΔₐG⁰overT`: vector of standard Gibbs energies divided by RT (dimensionless).
  - `T`: temperature in K (plain number, Dual-safe).
  - `P`: pressure in Pa (plain number, Dual-safe).
  - `ϵ`: regularization floor (default `1e-16`).

`T` and `P` are included so that temperature-dependent activity models
(e.g. [`HKFActivityModel`](@ref) with `temperature_dependent=true`) can
recompute their parameters inside the potential closure.
"""
function _build_params(state::ChemicalState; ϵ::Float64 = 1.0e-16)
    T = temperature(state)
    P = pressure(state)
    R = Constants.R
    RT = R * T                  # keeps units — division below strips them

    # ustrip without forced Float64 conversion — preserves Dual if T is Dual
    ΔₐG⁰overT = [
        ustrip(s[:ΔₐG⁰](T = T, P = P; unit = true) / RT)
            for s in state.system.species
    ]

    T_K = ustrip(us"K", T)   # Quantity{Dual} → Dual, Float64 → Float64
    P_Pa = ustrip(us"Pa", P)

    return (ΔₐG⁰overT = ΔₐG⁰overT, T = T_K, P = P_Pa, ϵ = ϵ)
end

"""
    _build_n0(state::ChemicalState) -> Vector

Extract the dimensionless mole vector from a `ChemicalState`.
Type is inferred from the state — compatible with ForwardDiff dual numbers.
"""
function _build_n0(state::ChemicalState)
    return ustrip.(us"mol", state.n)    # no Float64 cast — type follows eltype(state.n)
end

# ── Default solver factory ────────────────────────────────────────────────────

"""
    _DEFAULT_SOLVER_FACTORY

Internal `Ref{Union{Nothing, Function}}` — populated by extension `__init__`
functions to register a default solver factory.

  - `OptimizationIpoptExt.__init__`: registers only if nothing is set (low priority).
  - `OptimaSolverExt.__init__`: always overrides (high priority).

Result: OptimaSolver wins whenever loaded, regardless of load order.
"""
const _DEFAULT_SOLVER_FACTORY = Ref{Union{Nothing, Function}}(nothing)

# ── Differentiating an equilibrium ────────────────────────────────────────────

"""
    _primal(state::ChemicalState) -> ChemicalState

The same state with every dual number replaced by its value.
"""
function _primal(state::ChemicalState)
    n = [_plain(ustrip(us"mol", nᵢ)) * u"mol" for nᵢ in state.n]
    T = _plain(ustrip(us"K", temperature(state))) * u"K"
    P = _plain(ustrip(us"Pa", pressure(state))) * u"Pa"
    return ChemicalState(state.system, n; T = T, P = P)
end

"""
    _equilibrium_sensitivity(A, H, gθ, bdot, nstar; maxpin = 8) -> Vector

Sensitivity of an equilibrium composition, from the optimality conditions.

Equilibrium is the Gibbs minimization of [Leal2017](@cite), the problem Reaktoro
solves: `min G(n)` subject to `A n = b`, `n >= 0`, whose first-order conditions
are `grad G(n) - A' y - z = 0`, `A n = b`, `n_i z_i = 0`, with `y` the element
potentials and `z >= 0` the stability multipliers.

Differentiating them gives, on the FREE set (species present, `z_i = 0`),

```
    | H   A' | | ndot |   | -dgradG |
    | A   0  | | ydot | = |   bdot  |
```

with `ndot = 0` on the complement and `H = grad^2 G` at the solution. One
factorization serves every partial derivative, and the answer is exact — no
finite difference, no step size.

The complementarity block is not optional, and dropping it fails loudly rather
than subtly: on calcite + CO2 in water **with a gas phase declared**, the
unreduced system puts the whole perturbation into the absent gas species —
`n(CO2,g) = 5.8e-9`, held at its bound — returning `ndot = e_CO2(g)`, which
satisfies `A ndot = bdot` to 4e-16 and means nothing.

No back-end returns `z`, so the active set is recovered here: a species that is
negligible on the scale of the system yet takes a leading share of the response
is pinned, and the system re-solved. Each pass pins at least one species, so the
loop terminates.

`H` is singular by construction, and correctly so — a pure phase has unit
activity, hence a zero row. The saddle-point form handles that; any method
inverting `H` does not.
"""
function _equilibrium_sensitivity(A, H, gθ, bdot, nstar; maxpin::Int = 8)
    ns = length(nstar)
    m = size(A, 1)
    scale = maximum(abs, nstar)
    pinned = falses(ns)
    ndot = zeros(ns)

    for _ in 0:maxpin
        free = findall(!, pinned)
        K = [H[free, free] A[:, free]'; A[:, free] zeros(m, m)]
        rhs = vcat(-gθ[free], bdot)
        sol = try
            K \ rhs
        catch
            pinv(Matrix(K)) * rhs
        end
        fill!(ndot, 0.0)
        ndot[free] .= sol[1:length(free)]

        big = maximum(abs, ndot)
        offenders = [
            i for i in free
                if nstar[i] < 1.0e-6 * scale && abs(ndot[i]) > 0.1 * big
        ]
        isempty(offenders) && return ndot
        pinned[offenders] .= true
    end
    return ndot
end

"""
    STRICT_CONVERGENCE

Whether a non-converged equilibrium solve raises (`true`) or warns (`false`,
the default).

The default is *not* strict, and deliberately so: the back-end's convergence
flag is unreliable in both directions on these problems. It reports `MaxIters`
on points that are numerically excellent — pure water comes back flagged while
giving `[H⁺]/[OH⁻] = 1.000003` — and reports success on points that are not the
minimum. Raising on the flag alone would reject good answers and would still
miss the bad ones, so it is offered as an opt-in for callers who want the
strictest possible reading.
"""
const STRICT_CONVERGENCE = Ref(false)

"""
    NONCONVERGED :: Ref{Int}

Running count of equilibrium solves that returned a non-success retcode.

With `STRICT_CONVERGENCE[] = false` — the default — a non-converged solve is a
`@warn` at `maxlog = 1` and its result is used anyway. Over the thousands of
steps of a coupled kinetics run that is one warning for an arbitrary number of
bad speciations, and it never reached the failure count reported by
[`integrate`](@ref), which only saw solves that actually *threw*.

Reset it with `ChemistryLab.NONCONVERGED[] = 0` before a run and read it after;
[`integrate`](@ref) does exactly that and reports the total.
"""
const NONCONVERGED = Ref(0)

"""
    EXACT_HESSIAN

Whether to hand the back-end the exact Gibbs Hessian diagonal `∂μ/∂n` computed
by `ForwardDiff`, instead of letting it approximate one. Default `false`.

It is off by default only because it currently trades one defect for another,
and both are measured:

| | pure water | calcite + CO₂ |
|:--|--:|--:|
| `false` (back-end approximation) | `[H⁺]/[OH⁻] = 3.78` ✗ | worst ×19.8 |
| `true` (exact `∂μ/∂n`) | `[H⁺]/[OH⁻] = 1.000003` ✓ | worst ×3751 ✗ |

Turn it on for an **aqueous-only** system, where it makes the water
autoprotolysis come out right. Leave it off when a pure phase is present: the
exact curvature of such a phase is zero, the interior-point iteration then
stalls essentially at its starting point, and the whole speciation is wrong.

That stall is the open problem. It is not a matter of iterating longer (the
answer is identical at 300 and at 200 000 iterations) nor of the near-singular
Hessian entry (capping the inverse curvature over five orders of magnitude
changes nothing, because the solve stops before the barrier has decayed). The
C++ Optima this back-end is ported from offers a `Nullspace` linear solver in
addition to the `Rangespace` one implemented here — the latter carries the
warning that it suits diagonal Hessians only, and it is the one that inverts
`H`. Porting the nullspace path, which never inverts `H`, is the identified
next step.
"""
const EXACT_HESSIAN = Ref(false)

"""
    NULLSPACE_STEP

Whether the back-end computes its Newton step by the nullspace method rather
than the Schur complement. Default `true`.

The Schur complement forms `S = A H⁻¹ Aᵀ`, so it needs `H` invertible — and a
pure phase has unit activity, hence exactly zero curvature. The nullspace method
writes `dn = dnₚ + Z dz` with `Z` a basis of `null(A)` and solves
`(Zᵀ H Z) dz = −Zᵀ(ex + H dnₚ)`, in which `H` appears only as a product. This is
the route the C++ Optima takes by default, and its `Rangespace` counterpart —
the Schur complement — is documented there as suitable for invertible diagonal
Hessians only.

It is what makes the water autoprotolysis come out right: pure water gives
`[H⁺]/[OH⁻] = 1.0` and `pKw = 13.9994`, against 3.78 and 13.9897 through the
Schur complement, with no change to mixed solid/aqueous systems.
"""
const NULLSPACE_STEP = Ref(true)

"""
    _check_converged(sol, what) -> sol

Return `sol`, raising or warning according to [`STRICT_CONVERGENCE`](@ref) when
the optimizer did not converge. Before this existed, neither extension looked at
the return code and a non-converged iterate was written into the state as though
it were the equilibrium.
"""
function _check_converged(sol, what::AbstractString)
    SciMLBase.successful_retcode(sol) && return sol
    NONCONVERGED[] += 1
    if !STRICT_CONVERGENCE[]
        @warn "$what returned `$(sol.retcode)`; the composition may not be an \
               equilibrium. Set `ChemistryLab.STRICT_CONVERGENCE[] = true` to \
               raise instead." maxlog = 1
        return sol
    end
    throw(
        ErrorException(
            """
            $what did not converge: the optimizer returned `$(sol.retcode)`.

            Raise `max_iter`, loosen `tol`, or start from a better-conditioned \
            composition.
            """
        ),
    )
end

"""
    _solve_dual(esolver, state, ϵ) -> ChemicalState

Equilibrium of a composition carrying dual numbers.

No optimization solver is asked to iterate on dual numbers — most cannot, and
Ipopt never will, being a C library. The equilibrium is solved once at the
primal values and the sensitivities come from `_equilibrium_sensitivity`, the
implicit-function-theorem route on the optimality conditions.

Called from the back-end `solve` methods, which dispatch on the solver type;
making this a method of `solve` dispatching on the *state* would be ambiguous
with them.
"""
function _solve_dual(
        esolver::EquilibriumSolver,
        state::ChemicalState{C, S, Q, R},
        ϵ::Float64;
        b = nothing,
    ) where {C, S, Q, R <: ForwardDiff.Dual}

    state_v = _primal(state)
    eq_v = SciMLBase.solve(esolver, state_v; ϵ = ϵ, b = isnothing(b) ? nothing : _plain.(b))
    nstar = Float64[ustrip(us"mol", nᵢ) for nᵢ in eq_v.n]

    A = Float64.(state.system.SM.A)
    p_v = _build_params(state_v; ϵ = ϵ)
    H = ForwardDiff.jacobian(n -> esolver.μ(n, p_v), nstar)

    # The parameter enters through the potentials and through the element
    # amounts; both partial derivatives are read off the dual parts.
    p_d = _build_params(state; ϵ = ϵ)
    μ_d = esolver.μ(nstar, p_d)
    n0_d = _build_n0(state)

    npart = ForwardDiff.npartials(R)
    ns = length(nstar)
    ndot = Matrix{Float64}(undef, ns, npart)
    for k in 1:npart
        gθ = Float64[ForwardDiff.partials(μᵢ, k) for μᵢ in μ_d]
        bdot = isnothing(b) ?
            A * Float64[ForwardDiff.partials(nᵢ, k) for nᵢ in n0_d] :
            Float64[ForwardDiff.partials(bᵢ, k) for bᵢ in b]
        @views ndot[:, k] .= _equilibrium_sensitivity(A, H, gθ, bdot, nstar)
    end

    Tag = ForwardDiff.tagtype(R)
    n_dual = [
        ForwardDiff.Dual{Tag}(nstar[i], ForwardDiff.Partials(ntuple(k -> ndot[i, k], npart)))
            for i in 1:ns
    ]
    return ChemicalState(
        state.system, n_dual .* u"mol";
        T = temperature(state), P = pressure(state),
    )
end

# ── equilibrate ──────────────────────────────────────────────────────────────

"""
    equilibrate(state::ChemicalState, solver; model=..., variable_space=..., ϵ=...) -> ChemicalState
    equilibrate(state::ChemicalState; kwargs...) -> ChemicalState

Compute the chemical equilibrium state by minimizing the Gibbs free energy.

**Two-argument form** (solver explicit, always available once an extension is loaded):

```julia
using Optimization, OptimizationIpopt
state_eq = equilibrate(state, IpoptOptimizer())

using OptimaSolver
state_eq = equilibrate(state, OptimaOptimizer())
```

**One-argument form** (uses the default solver of the loaded extension):

```julia
state_eq = equilibrate(state)
```

When both extensions are loaded, `OptimaSolverExt` always takes priority.

# Arguments

  - `state`: initial `ChemicalState` — defines the system, T, P, and composition.
  - `solver`: any SciML-compatible solver (e.g. `IpoptOptimizer()`, `OptimaOptimizer()`).
  - `model`: activity model (default: `DiluteSolutionModel()`).
  - `variable_space`: `Val(:linear)` (default) or `Val(:log)`.
  - `ϵ`: regularization floor for mole amounts (default: `1e-16`).
  - `kwargs...`: forwarded to the underlying solver.
"""
function equilibrate(
        state::ChemicalState,
        solver;
        model::AbstractActivityModel = DiluteSolutionModel(),
        variable_space::Val = Val(:linear),
        ϵ::Float64 = 1.0e-16,
        kwargs...,
    )
    esolver = EquilibriumSolver(
        state.system, model, solver;
        variable_space = variable_space,
        kwargs...,
    )
    return SciMLBase.solve(esolver, state; ϵ = ϵ)
end

function equilibrate(state::ChemicalState; kwargs...)
    f = _DEFAULT_SOLVER_FACTORY[]
    if isnothing(f)
        error(
            "equilibrate without explicit solver requires an extension.\n" *
                "Add `using Optimization, OptimizationIpopt` or `using OptimaSolver`, " *
                "or call `equilibrate(state, solver; ...)` explicitly.",
        )
    end
    return equilibrate(state, f(); kwargs...)
end
