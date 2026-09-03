# ── certified.jl ──────────────────────────────────────────────────────────────
#
# An equilibrium that comes with a proof, or says it has none.
#
# The back ends do not agree, and the disagreement is not small. Measured on the
# same problems, with the element balance judged row by row against each row's
# own budget:
#
#   | case                        | interior point | dual Newton |
#   |-----------------------------|----------------|-------------|
#   | calcite in water, 10-40 °C  | 3.0e-2         | 3.0e-12     |
#   | calcite + 1 mmol CO2        | 1.0e-3         | 6.3e-13     |
#   | calcite + 50 mmol CO2       | 1.3e-16        | 8.5e-14     |
#   | pure water                  | 1.3e-16        | 1.3e-16     |
#   | CEM I paste, w/c 0.45, 0.60 | 1.5e-14        | 2.0e-14     |
#   | CEM I paste, w/c 0.30       | 7.2e-16        | 3.3e-10, not certified |
#
# The interior point is wrong in the second digit of the charge balance on a
# calcite solution, because the fraction-to-boundary rule caps its step and the
# residual stops moving; the dual Newton fails to admit a supersaturated phase on
# a low-water cement. Neither is the right default on its own.
#
# What settles it is that the problem is convex whenever the mixing terms are, so
# the KKT conditions are sufficient and `optimality_certificate` DECIDES rather
# than ranks. Given a decision procedure, solving by every available route and
# keeping a proved answer is exact, not heuristic.

"""
    equilibrate_certified(state; model, ϵ, b, verbose) -> (state, certificate)

Equilibrium composition together with a proof of its global optimality, obtained
by solving from every registered back end and keeping the answer
[`optimality_certificate`](@ref) proves optimal.

`certificate.optimal == true` is a **proof**, valid because the Gibbs
minimization is convex when the mixing terms are — ideal mixing and any activity
model whose excess Gibbs energy is convex in the amounts. It is not a proof for a
model that is not, and none of the activity models that ship here have been shown
to violate it; `HKFActivityModel`, `DaviesActivityModel` and the Redlich–Kister
solid solutions are used within their stated ranges.

When no route yields a proof, the answer with the smallest KKT error is returned,
its certificate says so, and a warning names the residual. That is the honest
outcome, and it is not the same thing as a failure: on a low-water cement the
returned composition satisfies the element balance to 1e-15 and has a
supersaturated phase left out, which the certificate reports as
`worst_supersaturation > 0`.

Requires `OptimaSolver` (the dual Newton lives there). Systems without an aqueous
phase, or without `H2O@`, cannot use the dual route; for those, this falls back to
the plain [`equilibrate`](@ref) and returns `nothing` as the certificate.

# Example

```julia
using ChemistryLab, OptimaSolver
eq, cert = equilibrate_certified(state)
cert.optimal          # true — proved globally optimal
cert.balance          # element balance residual
cert.worst_supersaturation   # negative: every absent phase undersaturated
```
"""
function equilibrate_certified(
        state::ChemicalState;
        model::AbstractActivityModel = DiluteSolutionModel(),
        b = nothing,
        ϵ::Float64 = 1.0e-16,
        verbose::Bool = false,
        constraint::EquilibriumConstraint = FixedTP(),
        parameters::Union{Nothing, Base.RefValue} = nothing,
        kwargs...,
    )
    if !_DUAL_AVAILABLE[]
        error(
            "equilibrate_certified needs `OptimaSolver`: the KKT solver and the " *
                "certificate live there. Add `using OptimaSolver` — and " *
                "optionally `using Optimization, OptimizationIpopt` as a second " *
                "starting route, since neither back end certifies every case.",
        )
    end

    if !_dual_applicable(state.system)
        constraint isa FixedTP || throw(
            ArgumentError(
                "a constraint other than `FixedTP` needs the dual route, which " *
                    "requires an aqueous phase and `H2O@` among the species: the " *
                    "prescribed property is an unknown of that solver's system, " *
                    "and the interior-point back ends have nowhere to put it.",
            )
        )
        # No dual route: return the plain answer and say there is no proof.
        eq = equilibrate(state; model = model, ϵ = ϵ, certify = false, kwargs...)
        return (eq, nothing)
    end

    # Duals take the implicit-function route, dispatched on the state's element
    # type rather than tested for. See `_certified_primal_then_derivative`.
    dual_route = _certified_dual_route(
        _amount_number_type(state), state, model, b, ϵ, verbose, constraint,
        parameters, kwargs,
    )
    dual_route === nothing || return dual_route

    des = DualEquilibriumSolver(state.system, model; verbose = verbose)

    # `b` is fixed ONCE, from the state as given. Letting each start define its
    # own would pose a different problem for each: a start that violates the
    # balance — the interior point does, by 3e-6 mol on this class of problem —
    # shifts the component totals by exactly its own infeasibility, and the dual
    # solve then certifies the answer to the shifted problem. Measured, that gave
    # two "certified" compositions 0.2 % apart on dissolved calcium, which on a
    # convex problem with one minimum can only mean two different problems.
    bfix = b === nothing ?
        des.A * Float64[ustrip(us"mol", x) for x in state.n] :
        Float64.(collect(b))

    starts = ChemicalState[]
    for f in _SOLVER_FACTORIES
        try
            esolver = EquilibriumSolver(state.system, model, f(); kwargs...)
            push!(starts, SciMLBase.solve(esolver, state; ϵ = ϵ, b = bfix))
        catch err
            verbose && @info "start rejected" backend = f err
        end
    end
    # The state as given is a legitimate start too, and the only one available if
    # every back end threw.
    push!(starts, state)

    eq, cert = solve_certified(
        des, starts; b = bfix, ϵ = ϵ, constraint = constraint, parameters = parameters,
    )
    if !cert.optimal
        @warn "no route produced a certifiable equilibrium; returning the answer with the smallest KKT error — audit it with `optimality_certificate`" stationarity = cert.stationarity balance = cert.balance worst_supersaturation = cert.worst_supersaturation maxlog = 1
    end
    return (eq, cert)
end

"""
    _dual_applicable(system) -> Bool

Whether [`DualEquilibriumSolver`](@ref) can be built for `system`: it needs an
aqueous phase, and `H2O@` among the species, because it parameterizes the interior
variables by the solvent's chemical potential.
"""
function _dual_applicable(system::ChemicalSystem)
    isempty(system.idx_aqueous) && return false
    return haskey(system.dict_species, "H2O@")
end


"""
    _certified_dual_route(state, model, b, ϵ, verbose, constraint, parameters, kwargs)

`nothing` for a real-valued composition; the certified answer with its derivative
attached for one carrying `ForwardDiff.Dual` amounts.

Dispatched on the element type, positionally, so the choice is the type system's
and neither path pays for the other. Two paths are needed for a mathematical
reason, not for want of a generic element type: making the component totals
generic would let duals flow into the solver, and what came back would be the
derivative of the **algorithm** — an active set decided by sign tests, a line
search with branches, an iteration count that varies with the data — rather than
the derivative of the **solution**. The map `b ↦ n*(b)` is smooth only piecewise,
and on each piece the implicit function theorem gives its derivative at the
solution with the active set frozen. That is how `OptimaSolver` computes its own
`Sensitivity`, and how Optima does upstream.
"""
_amount_number_type(state::ChemicalState) = _number_type(eltype(state.n))
_number_type(::Type{<:DynamicQuantities.AbstractQuantity{T}}) where {T} = T
_number_type(::Type{T}) where {T <: Real} = T

_certified_dual_route(::Type{<:Real}, state, model, b, ϵ, verbose, constraint, parameters, kwargs) =
    nothing

function _certified_dual_route(
        ::Type{<:ForwardDiff.Dual}, state, model, b, ϵ, verbose, constraint,
        parameters, kwargs,
    )
    eq_v, cert = equilibrate_certified(
        _primal(state); model = model, ϵ = ϵ, verbose = verbose,
        constraint = constraint, parameters = parameters,
        b = b === nothing ? nothing : _plain.(b), kwargs...,
    )
    nstar = Float64[ustrip(us"mol", x) for x in eq_v.n]
    μ = build_potentials(state.system, model)
    return (_attach_sensitivity(state, nstar, μ, ϵ; b = b), cert)
end
