# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
#
# Post-processing of a kinetics ODESolution: reaction extents, reconstructed
# chemical states, degrees of hydration.
#
# The ODE state carries the kinetic species moles and the EXTENTS OF REACTION ξ
# (`dξ/dt = r`), plus element amounts and temperature where applicable. The moles
# of every other species are not stored — the parameter tuple keeps a single
# `n_full` buffer, mutated in place, which after a run holds the last composition
# only — but they are recovered exactly from `n = n(0) + νᵀ ξ`, which is what
# this file does.

using OrderedCollections

# ── reaction_extents ─────────────────────────────────────────────────────────

"""
    reaction_extents(sol, kp::KineticsProblem; times = sol.t) -> Matrix{Float64}

Extent of reaction ``ξ_j(t) = ∫_0^t r_j\\,dt`` for every kinetic reaction of
`kp`, along the solution `sol`.

Returns a `length(times) × M` matrix, `M` being the number of kinetic reactions,
in mol.

The extents are **carried by the ODE state** — the integrator advances
`dξ/dt = r` alongside the kinetic moles — so this is a read of the solution, not
a quadrature. It is exact to the solver's own tolerance at any instant, and
costs nothing.

# Arguments

  - `sol`: the `ODESolution` returned by [`integrate`](@ref).
  - `kp`: the [`KineticsProblem`](@ref) that produced it.
  - `times`: instants at which the extents are wanted. Any instants within
    `tspan`, in any order; the dense output is interpolated.

See also: [`state_at`](@ref), [`degrees_of_hydration`](@ref),
[`extent_residual`](@ref).
"""
function reaction_extents(sol, kp::KineticsProblem; times = sol.t)
    p = sol.prob.p
    M = p.n_rxn_state
    off = p.n_be + p.n_nk
    ts = collect(float.(times))
    ξ = zeros(Float64, length(ts), M)
    for (i, t) in enumerate(ts)
        u = sol(t)
        for j in 1:M
            ξ[i, j] = u[off + j]
        end
    end
    return ξ
end

"""
    extent_residual(sol, kp::KineticsProblem, ξ; times = sol.t) -> Float64

Largest absolute discrepancy, in mol, between the kinetic-species moles carried
by the integrator and those implied by the extents `ξ`, that is
``\\max_{i,t} |n_k(t) - n_k(0) - (ν_k^{\\mathsf T} ξ(t))_i|``.

The two are redundant by construction — `nₖ` integrates `νₖᵀ r` while `ξ`
integrates `r` — so this measures how far the integrator has let them drift, and
should sit at the solver's tolerance. A large value means the solve is too loose,
not that the post-processing is inaccurate.
"""
function extent_residual(sol, kp::KineticsProblem, ξ::AbstractMatrix; times = sol.t)
    p = sol.prob.p
    nk0 = p.n_initial_full[kp.idx_kinetic]
    worst = 0.0
    for (i, t) in enumerate(times)
        u = sol(t)
        nk = u[(p.n_be + 1):(p.n_be + p.n_nk)]
        pred = nk0 .+ kp.νk' * ξ[i, :]
        worst = max(worst, maximum(abs, nk .- pred))
    end
    return worst
end

# ── state_at ─────────────────────────────────────────────────────────────────

"""
    state_at(sol, kp::KineticsProblem, t; ξ = nothing) -> ChemicalState

Reconstruct the full [`ChemicalState`](@ref) of the system at time `t` from a
kinetics solution.

Every species — kinetic and non-kinetic alike — is obtained from the same
stoichiometric replay ``n = n_0 + ν^{\\mathsf T} ξ(t)``, so the returned state
satisfies element conservation exactly. Temperature is taken from the solution
when a [`SemiAdiabaticCalorimeter`](@ref) is in use, and from the problem
otherwise.

The extents are read from the ODE state, so this is exact to the solver's
tolerance. Pass a precomputed extent matrix as `ξ` — any matrix whose **last
row** holds ξ(t) — to avoid re-interpolating the solution when sweeping many
instants.

!!! note "Re-speciation is not replayed"
    When the run used an equilibrium solver, the aqueous partition was
    re-speciated at every accepted step and that redistribution is *not*
    recoverable from the stoichiometry alone. The returned state is then the
    purely kinetic reconstruction; call [`equilibrate`](@ref) on it to recover
    the speciated composition.

# Examples

```julia
ξ = reaction_extents(sol, kp)                      # once
states = [state_at(sol, kp, t; ξ = ξ[i:i, :]) for (i, t) in enumerate(sol.t)]
```

See also: [`reaction_extents`](@ref), [`degrees_of_hydration`](@ref),
[`volume_fractions`](@ref).
"""
function state_at(sol, kp::KineticsProblem, t::Real; ξ = nothing)
    ξ_t = if ξ === nothing
        vec(reaction_extents(sol, kp; times = [float(t)]))
    else
        vec(ξ[end, :])
    end

    p = sol.prob.p
    n_new = p.n_initial_full .+ kp.ν' * ξ_t
    @. n_new = max(n_new, 0.0)

    T_new = kp.calorimeter isa SemiAdiabaticCalorimeter ? sol(t)[end] : p.T
    return ChemicalState(
        kp.system;
        T = T_new * u"K",
        P = p.P * u"Pa",
        n = [nᵢ * u"mol" for nᵢ in n_new],
    )
end

# ── degrees_of_hydration ─────────────────────────────────────────────────────

"""
    degrees_of_hydration(sol, kp::KineticsProblem; times = sol.t) -> OrderedDict{String, Vector{Float64}}

Degree of reaction ``α_i(t) = 1 - n_i(t)/n_i(0)`` of every kinetic species,
keyed by species symbol.

Species with a zero initial amount are omitted rather than returned as `NaN`:
a phase that was never present has no degree of reaction.

# Examples

```julia
α = degrees_of_hydration(sol, kp)
α["C3S"][end]                       # final degree of hydration of alite
```

See also: [`mean_degree_of_hydration`](@ref), [`state_at`](@ref).
"""
function degrees_of_hydration(sol, kp::KineticsProblem; times = sol.t)
    p = sol.prob.p
    out = OrderedDict{String, Vector{Float64}}()
    us = [sol(t) for t in times]
    for (j, idx) in enumerate(kp.idx_kinetic)
        n0 = p.n_initial_full[idx]
        n0 > 0 || continue
        sym = symbol(kp.system.species[idx])
        out[sym] = [1.0 - u[p.n_be + j] / n0 for u in us]
    end
    return out
end

"""
    mean_degree_of_hydration(sol, kp::KineticsProblem; times = sol.t, weights = :mass) -> Vector{Float64}

Degree of reaction of the binder as a whole, averaged over the kinetic species.

`weights` selects the averaging:

  - `:mass` (default) — weighted by the initial **mass** of each species, the
    convention used when a single ᾱ is quoted for a cement.
  - `:mole` — weighted by initial moles.

Only species with a non-zero initial amount take part, consistently with
[`degrees_of_hydration`](@ref).
"""
function mean_degree_of_hydration(
        sol, kp::KineticsProblem; times = sol.t, weights::Symbol = :mass
    )
    weights in (:mass, :mole) ||
        throw(ArgumentError("weights must be :mass or :mole, got :$weights"))
    p = sol.prob.p
    α = degrees_of_hydration(sol, kp; times = times)
    w = Float64[]
    for idx in kp.idx_kinetic
        n0 = p.n_initial_full[idx]
        n0 > 0 || continue
        sp = kp.system.species[idx]
        if weights === :mass
            haskey(properties(sp), :M) || throw(
                ArgumentError(
                    "species \"$(symbol(sp))\" has no molar mass :M, so a " *
                        "mass-weighted mean is undefined; pass weights = :mole"
                )
            )
            push!(w, n0 * ustrip(us"kg/mol", sp[:M]))
        else
            push!(w, n0)
        end
    end
    total = sum(w)
    out = zeros(Float64, length(times))
    for (k, αᵢ) in enumerate(values(α))
        @. out += w[k] * αᵢ
    end
    return out ./ total
end

# ── speciated_states ─────────────────────────────────────────────────────────

"""
    speciated_states(sol, kp::KineticsProblem; times = sol.t) -> Vector{ChemicalState}

The **speciated** compositions along a solution: the kinetic species read from
the ODE state, and the equilibrium partition recovered by re-solving `φ(bₑ)` at
each instant with the element totals the run carried.

This is what [`state_at`](@ref) deliberately does not do. That function returns
the purely kinetic reconstruction `n(0) + νᵀξ`, because the redistribution
performed by the equilibrium solve is not recoverable from the stoichiometry —
for a cement that reconstruction is meaningless, putting every dissolved element
in solution with not one hydrate. Nor is the composition left in the solver's
own buffers usable: that buffer is rewritten at every right-hand-side
evaluation, Jacobian differences and rejected steps included, so it is not the
accepted composition at any particular time.

Without an equilibrium solver on `kp` there is nothing to replay, and this
returns [`state_at`](@ref) for each instant.

# The instants must be ascending

Each solve is warm-started from the previous one, and the guess is first capped
at the element budget and projected back into `{Aₑn = bₑ, n ≥ 0}`. Both matter,
and neither is optional:

  - the warm start is the equilibrium of the *previous* `bₑ`, so once an element
    has been spent — the sulfate of an ordinary Portland cement, after the gypsum
    is gone — it demands more of it than now exists and the interior-point solve
    begins outside its own feasible set;
  - solved instead from a cold guess, `φ(bₑ)` does not converge at all for a
    cement: the H⁺ component of `bₑ` reaches −14 mol, and a solve started from
    pure water returns no hydrates and a pore solution at pH 6 while the run
    itself computed 2.2 mol of C-S-H.

# Examples

```julia
states = speciated_states(sol, kp; times = [1, 7, 28] .* 86400.0)
pH(states[end])                      # pore solution at 28 days
volume_fractions(states[end], groups; reference = state0)
```

See also: [`state_at`](@ref), [`reaction_extents`](@ref), [`volume_fractions`](@ref).
"""
function speciated_states(sol, kp::KineticsProblem; times = sol.t)
    p = sol.prob.p

    # No equilibrium partition: the kinetic reconstruction IS the composition.
    p.n_be == 0 && return ChemicalState[state_at(sol, kp, t) for t in times]

    issorted(times) || throw(
        ArgumentError(
            "`times` must be ascending: each speciation warm-starts from the previous one."
        )
    )

    # NOTE the absent `ϵ`. `build_kinetics_params` defaults it to 1e-30, which as
    # a lower bound is fourteen orders below the solve's own 1e-16, and an
    # interior-point method pinned that close to zero does not recover: replayed
    # with 1e-30 a full OPC came out at pH 14.7 with 0.44 mol of ettringite
    # against 0.27 mol of sulfate, and with the default at pH 12.6 and the
    # sulfate budget closing exactly.
    sub = _equilibrium_subsystem(kp.system, kp.idx_equilibrium)
    # A FRESH back-end instance, not the one the run used. `p.eq_solver.solver`
    # is the very object that drove thousands of solves during the integration,
    # and replaying through it returned pH 14.2 with 0.31 mol of ettringite and
    # no AFm, where a clean instance of the same type and settings gives pH 12.58
    # with the sulfate entirely in AFm — the same trajectory, the same `bₑ`, the
    # same guess. The back-end therefore carries state across solves; until that
    # is fixed upstream, a replay must not inherit it.
    es = EquilibriumSolver(
        sub, activity_model(p.eq_solver), typeof(p.eq_solver.solver)();
        variable_space = p.eq_solver.variable_space,
    )
    guess = Float64[max(x, _EQ_GUESS_FLOOR) for x in p.n_eq_init]
    n_sp = length(kp.system.species)

    # Walk exactly the instants asked for, and never before the first of them.
    #
    # Refining the walk with the integrator's own accepted steps looks safer and
    # is not: `sol.t` begins at t = 0, where `bₑ` is mixing water and nothing
    # else. The equilibrium there is degenerate, and warm-starting the next solve
    # from it poisons the whole chain — replayed that way a full OPC came out at
    # pH 14.2 with 0.31 mol of ettringite and no AFm at all, against pH 12.58 and
    # the sulfate entirely in AFm when the walk starts at the first requested
    # instant.
    #
    # The consequence is that the FIRST requested instant matters: ask for one
    # late enough that the paste has reacted, early enough that the jump from the
    # cast composition is not too large. A first point within the first hours is
    # what the cement cases here use.
    out = ChemicalState[]
    for t in times
        u = sol(t)
        be = collect(@view u[1:(p.n_be)])

        _budget_clip!(guess, p.Ae, be)
        _restore_feasibility!(guess, p.Ae, be)

        eq = SciMLBase.solve(
            es,
            ChemicalState(sub, guess .* u"mol"; T = p.T_q[], P = p.P_q[]);
            b = be,
        )
        guess = Float64[max(ustrip(us"mol", x), _EQ_GUESS_FLOOR) for x in eq.n]

        n = zeros(Float64, n_sp)
        for (j, idx) in enumerate(kp.idx_equilibrium)
            n[idx] = ustrip(us"mol", eq.n[j])
        end
        for (j, idx) in enumerate(kp.idx_kinetic)
            n[idx] = max(u[p.n_be + j], 0.0)
        end

        push!(
            out, ChemicalState(
                kp.system; T = p.T * u"K", P = p.P * u"Pa",
                n = [nᵢ * u"mol" for nᵢ in n],
            )
        )
    end
    return out
end
