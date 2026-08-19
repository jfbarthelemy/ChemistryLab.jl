# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── dual_solver.jl ────────────────────────────────────────────────────────────
#
# Chemical equilibrium by exact Newton on the KKT system, certified.
#
# THE ALGORITHM IS NOT HERE. It lives in `OptimaSolver.dual_newton_solve`, where
# it belongs: the active set on the bound-constrained variables, the degeneracy
# test on the rows of `A`, the two-level Newton and the certificate are all
# statements about a convex program, not about chemistry. What this file supplies
# is the four things that ARE chemistry:
#
#   1. the conservation matrix `A` and the reference potentials `g = Δ_a G⁰/RT`;
#   2. the activity model, as the callback `h` with `∇f = g + h(n)`;
#   3. which species are strictly positive at any solution — the aqueous ones,
#      parameterized by `ln n` — and which may vanish — the pure phases;
#   4. which species plays the reference role: the solvent, whose activity is a
#      mole fraction and therefore bounded above, so that its stationarity cannot
#      be inverted and must be carried by the outer system.
#
# The mathematics is documented in the manual under
# "Proving that an answer is the answer", and the general form in OptimaSolver's
# theory page.

using LinearAlgebra

"""
    DualEquilibriumSolver(system, model; tol, maxit, max_active_updates, si_tol, verbose)

Equilibrium by Newton on the KKT system, in element-potential space.

Where [`EquilibriumSolver`](@ref) minimizes `G` by an interior-point method over
all species and stops on `MaxIters` for a cement, this solves the stationarity
conditions directly and returns a result that
[`optimality_certificate`](@ref) can **prove** optimal — the Gibbs problem being
convex, the KKT conditions are sufficient.

Writing `u = −Aᵀy` for the element potentials, an aqueous species obeys the
mass-action law `aᵢ = exp(uᵢ − gᵢ)` and a pure phase is present exactly when
`uᵢ = gᵢ`, absent when undersaturated: the classical phase-stability criterion.

Starting from the answer of an [`EquilibriumSolver`](@ref) is the intended use —
the interior-point method reaches a neighborhood, this reaches the conditions.

See also: [`optimality_certificate`](@ref), [`speciated_states`](@ref).
"""
struct DualEquilibriumSolver{L, M <: AbstractActivityModel}
    system::ChemicalSystem
    lna::L
    model::M
    idx_aq::Vector{Int}
    idx_pure::Vector{Int}
    j_solvent::Int
    ss_groups::Vector{Vector{Int}}   # end-members of each declared solid solution
    A::Matrix{Float64}
    opts::NamedTuple
end

function DualEquilibriumSolver(
        system::ChemicalSystem,
        model::AbstractActivityModel = DiluteSolutionModel();
        tol::Float64 = 1.0e-10,
        maxit::Int = 200,
        max_active_updates::Int = 200,
        si_tol::Float64 = 1.0e-8,
        verbose::Bool = false,
    )
    idx_aq = [i for (i, s) in enumerate(system.species) if aggregate_state(s) == AS_AQUEOUS]
    idx_pure = [i for (i, s) in enumerate(system.species) if aggregate_state(s) != AS_AQUEOUS]
    isempty(idx_aq) && throw(
        ArgumentError(
            "DualEquilibriumSolver needs an aqueous phase: the element potentials are " *
                "determined by the mass-action laws of the aqueous species."
        )
    )
    jw = something(findfirst(i -> symbol(system.species[i]) == "H2O@", idx_aq), 0)
    jw == 0 && throw(
        ArgumentError(
            "DualEquilibriumSolver needs `H2O@` among the species: the solvent's " *
                "activity is a mole fraction, hence bounded above, so its stationarity " *
                "cannot be inverted and is carried by the outer system instead."
        )
    )
    # A solid solution is a MIXING phase, exactly like the aqueous one: its
    # end-members have composition-dependent activities, so none of them is ever
    # exactly absent while the phase exists, and the active set belongs at the
    # level of the phase. They are therefore removed from the bound-constrained
    # set, where they would have been treated as pure phases.
    ss = system.solid_solutions
    ss_groups = Vector{Int}[]
    if ss !== nothing
        byname = Dict(symbol(sp) => i for (i, sp) in enumerate(system.species))
        for phase in ss
            idx = [byname[symbol(em)] for em in phase.end_members if haskey(byname, symbol(em))]
            length(idx) == length(phase.end_members) && push!(ss_groups, idx)
        end
    end
    in_ss = Set(vcat(ss_groups...))
    idx_pure = [i for i in idx_pure if !(i in in_ss)]

    return DualEquilibriumSolver(
        system, activity_model(system, model), model,
        idx_aq, idx_pure, jw, ss_groups, Float64.(system.SM.A),
        (; tol, maxit, max_active_updates, si_tol, verbose),
    )
end

activity_model(des::DualEquilibriumSolver) = des.model

"""
    _dual_problem(des, p, n0) -> DualNewtonProblem

Package the chemistry as the convex program `OptimaSolver` solves. Built per
solve because the reference potentials `Δ_a G⁰/RT` depend on temperature and
pressure.
"""
function _dual_problem(des::DualEquilibriumSolver, p, n0)
    # One mixing phase for the aqueous solution — always present, the solvent as
    # its reference — and one more per declared solid solution, whose presence
    # the tangent-plane test decides.
    #
    # The reference is the member whose stationarity the OUTER system carries
    # instead of inverting, so it has to be the one the phase is mostly made of.
    # For the aqueous phase that is the solvent, and not by convention: water is
    # the only species there whose activity is a mole fraction, bounded above, so
    # its law cannot be inverted at all. Inside a solid solution every member is
    # a mole fraction and any of them could serve — but the choice is not
    # indifferent. The outer unknown is `ln x_ref`, and picking a member that
    # happens to be minor sends that unknown towards −∞ and the Jacobian with it.
    # So the reference is chosen by magnitude, from the composition the caller
    # supplied, which is what the solvent already is for the aqueous phase.
    phases = [
        (
            members = des.idx_aq, j_ref = des.j_solvent,
            always_present = true, mole_fraction = false,
        ),
    ]
    for grp in des.ss_groups
        j_ref = argmax(@view n0[grp])
        push!(
            phases,
            (members = grp, j_ref = j_ref, always_present = false, mole_fraction = true),
        )
    end
    return _optima_dual_problem(
        des.A, Float64.(p.ΔₐG⁰overT), des.lna, phases, des.idx_pure, p,
    )
end

"""
    SciMLBase.solve(des::DualEquilibriumSolver, state; b = nothing, ϵ = 1e-16)
        -> ChemicalState

Equilibrium composition. `state` supplies the temperature, the pressure and the
starting guess; `b` the component totals, defaulting to those of `state`.

The amounts are returned as computed, without a floor: an amount of `e⁻³⁰⁰` IS
the mass-action answer for a species that is not there, and raising it to `ϵ`
falsifies its activity by hundreds of `RT` units — which is enough to make
[`optimality_certificate`](@ref) report a residual of 74 on a composition solved
to `5e-12`.
"""
function SciMLBase.solve(
        des::DualEquilibriumSolver,
        state::ChemicalState;
        b = nothing,
        ϵ::Float64 = 1.0e-16,
    )
    p = _build_params(state; ϵ = ϵ)
    n0 = Float64[ustrip(us"mol", x) for x in state.n]
    bv = b === nothing ? des.A * n0 : Float64.(collect(b))

    res = _optima_dual_solve(_dual_problem(des, p, n0), bv, n0, des.opts)

    res.converged || begin
        NONCONVERGED[] += 1
        @warn """the dual equilibrium solve did not certify optimality; audit it with \
        `optimality_certificate`.""" maxlog = 1
    end

    return ChemicalState(
        des.system, [nᵢ * u"mol" for nᵢ in res.x];
        T = temperature(state), P = pressure(state),
    )
end

"""
    optimality_certificate(des, state; b = nothing, ϵ = 1e-16, floor = 1e-25)
        -> (; stationarity, balance, worst_supersaturation, n_interior,
             n_absent_component, optimal)

Check the KKT conditions at a composition, independently of how it was obtained.

For a convex problem these conditions are sufficient, so `optimal = true` is a
proof of **global** optimality. Use it to audit any solver — including
[`EquilibriumSolver`](@ref), whose interior-point iteration reports `MaxIters` on
a cement equilibrium and cannot say whether the point it returns is the answer.

The three quantities are the stationarity of the interior species, the component
balance, and the worst saturation index among absent phases (negative when every
one of them is undersaturated, as optimality requires).
"""
function optimality_certificate(
        des::DualEquilibriumSolver, state::ChemicalState;
        b = nothing, ϵ::Float64 = 1.0e-16, floor::Float64 = 1.0e-25,
    )
    p = _build_params(state; ϵ = ϵ)
    n = Float64[ustrip(us"mol", x) for x in state.n]
    bv = b === nothing ? des.A * n : Float64.(collect(b))

    c = _optima_kkt_certificate(
        _dual_problem(des, p, n), n, bv, floor, des.opts.tol, des.opts.si_tol,
    )
    return (;
        stationarity = c.stationarity, balance = c.feasibility,
        worst_supersaturation = c.worst_violation, n_interior = c.n_interior,
        n_absent_component = c.n_forced_zero, optimal = c.optimal,
    )
end

# ── hooks filled in by the OptimaSolver extension ────────────────────────────
#
# The algorithm is `OptimaSolver`'s, and that package is a weak dependency here,
# so the three entry points are indirected through functions the extension
# overrides. Called without it loaded, they say what is missing.

_optima_dual_problem(args...) = _need_optima()
_optima_dual_solve(args...) = _need_optima()
_optima_kkt_certificate(args...) = _need_optima()

_need_optima() = error(
    "DualEquilibriumSolver needs OptimaSolver ≥ 0.3: the KKT solver and its " *
        "certificate live there. Add `using OptimaSolver`."
)
