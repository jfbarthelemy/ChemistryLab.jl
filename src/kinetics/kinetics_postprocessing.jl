# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
#
# Post-processing of a kinetics ODESolution: reaction extents, reconstructed
# chemical states, degrees of hydration.
#
# The ODE state vector `u` carries only the KINETIC species moles (plus element
# amounts and temperature when applicable). The moles of every other species are
# never stored: the integrator's parameter tuple keeps a single `n_full` buffer
# that is mutated in place, so after a run it holds the LAST composition only.
# Recovering the composition at an arbitrary time therefore requires replaying
# the stoichiometry, which is what this file does.

using OrderedCollections

# ── reaction_extents ─────────────────────────────────────────────────────────

"""
    reaction_extents(sol, kp::KineticsProblem; times = sol.t, nsub = 4) -> Matrix{Float64}

Extent of reaction ``ξ_j(t) = ∫_0^t r_j\\,dt`` for every kinetic reaction of
`kp`, evaluated along the solution `sol`.

Returns a `length(times) × M` matrix, `M` being the number of kinetic reactions,
in mol.

The rates are not stored by the integrator, so the extents are recovered by
quadrature on the solver's dense output. The integration always runs on the
**union of `times` and the accepted steps** `sol.t`, each resulting interval
being split into `nsub` sub-intervals and integrated by composite Simpson's
rule; the requested rows are then extracted. Asking for a coarse output grid
therefore costs accuracy nothing — the early transient of a hydration run is
resolved by the solver's own steps, and a uniform user grid would step straight
over it.

The accuracy floor is set by the dense output, not by `nsub`: beyond a handful
of sub-intervals the residual reported by [`extent_residual`](@ref) stops
falling and tracks the integrator's `reltol` instead. Tighten the solve, not
`nsub`, if it is too large.

# Arguments

  - `sol`: the `ODESolution` returned by [`integrate`](@ref).
  - `kp`: the [`KineticsProblem`](@ref) that produced it.
  - `times`: instants at which the extents are wanted, ascending. Defaults to the
    accepted steps. They need **not** start at `kp.tspan[1]`: the quadrature
    always runs from the beginning of the problem, so `times = [t]` returns the
    extent accumulated over the whole run up to `t`.
  - `nsub`: sub-intervals per output interval for the quadrature.

See also: [`state_at`](@ref), [`degrees_of_hydration`](@ref),
[`extent_residual`](@ref).
"""
function reaction_extents(sol, kp::KineticsProblem; times = sol.t, nsub::Int = 4)
    nsub >= 1 || throw(ArgumentError("nsub must be ≥ 1, got $nsub"))
    isodd(nsub) && (nsub += 1)   # Simpson needs an even number of sub-intervals
    p = sol.prob.p
    M = length(kp.kinetic_reactions)
    ts = collect(float.(times))
    issorted(ts) || throw(ArgumentError("`times` must be ascending"))

    # Integrate on the union of the requested grid and the solver's accepted
    # steps: a coarse user grid would otherwise step over the early transient,
    # where nearly all of the extent accumulates.
    #
    # The quadrature ALWAYS starts at `kp.tspan[1]`, never at `times[1]`. An
    # extent is measured from the start of the run, so asking for a single late
    # instant must integrate everything before it — returning zero there would be
    # silently wrong rather than an error.
    t_start = float(kp.tspan[1])
    grid = sort!(unique!(vcat(t_start, ts, [t for t in sol.t if t_start <= t <= ts[end]])))

    ξ_grid = zeros(Float64, length(grid), M)
    rates = zeros(Float64, M)
    acc = zeros(Float64, M)

    for i in 2:lastindex(grid)
        t0, t1 = grid[i - 1], grid[i]
        h = (t1 - t0) / nsub
        fill!(acc, 0.0)
        for k in 0:nsub
            # Simpson weights 1, 4, 2, 4, …, 4, 1
            w = (k == 0 || k == nsub) ? 1.0 : (isodd(k) ? 4.0 : 2.0)
            _rates_at!(rates, sol, p, kp, t0 + k * h)
            @. acc += w * rates
        end
        @. ξ_grid[i, :] = ξ_grid[i - 1, :] + acc * h / 3
    end

    # Extract the requested rows (`grid` is sorted and contains every `ts`).
    rows = [searchsortedfirst(grid, t) for t in ts]
    return ξ_grid[rows, :]
end

# Evaluate every kinetic reaction rate at time `t` along the solution.
function _rates_at!(rates::Vector{Float64}, sol, p, kp::KineticsProblem, t::Real)
    u = sol(t)
    n_full = copy(p.n_initial_full)
    nk = @view u[(p.n_be + 1):(p.n_be + p.n_nk)]
    for (j, idx) in enumerate(p.idx_kinetic)
        n_full[idx] = max(nk[j], p.ϵ)
    end
    T_curr = kp.calorimeter isa SemiAdiabaticCalorimeter ? u[end] : p.T
    n_sv = StateView(n_full, p.species_index)
    lna_sv = StateView(p.lna_fn(n_full, p), p.species_index)
    n0_sv = StateView(p.n_initial_full, p.species_index)
    for (j, kr) in enumerate(p.kin_rxns)
        rates[j] = kr.rate_fn(T_curr, p.P, t, n_sv, lna_sv, n0_sv)
    end
    return rates
end

"""
    extent_residual(sol, kp::KineticsProblem, ξ; times = sol.t) -> Float64

Largest absolute discrepancy, in mol, between the kinetic-species moles carried
by the integrator and those implied by the extents `ξ`, that is
``\\max_{i,t} |n_k(t) - n_k(0) - (ν_k^{\\mathsf T} ξ(t))_i|``.

This is the quadrature error of [`reaction_extents`](@ref), and it is the
quantity to watch when deciding whether `nsub` is large enough: the
reconstruction of [`state_at`](@ref) is mass-balanced by construction, so this
residual is the only place the quadrature can show up.
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
    state_at(sol, kp::KineticsProblem, t; ξ = nothing, nsub = 4) -> ChemicalState

Reconstruct the full [`ChemicalState`](@ref) of the system at time `t` from a
kinetics solution.

Every species — kinetic and non-kinetic alike — is obtained from the same
stoichiometric replay ``n = n_0 + ν^{\\mathsf T} ξ(t)``, so the returned state
satisfies element conservation exactly. Temperature is taken from the solution
when a [`SemiAdiabaticCalorimeter`](@ref) is in use, and from the problem
otherwise.

Pass a precomputed extent matrix as `ξ` (from [`reaction_extents`](@ref) with
`times = [kp.tspan[1], t]`, or any matrix whose **last row** holds ξ(t)) to
avoid re-integrating the rates when sweeping many instants; see the example
below.

!!! note "Re-speciation is not replayed"
    When the run used an equilibrium solver, the aqueous partition was
    re-speciated at every accepted step and that redistribution is *not*
    recoverable from the stoichiometry alone. The returned state is then the
    purely kinetic reconstruction; call [`equilibrate`](@ref) on it to recover
    the speciated composition.

# Examples

```julia
ξ = reaction_extents(sol, kp)                      # once
states = [state_at(sol, kp, t; ξ = ξ[1:i, :]) for (i, t) in enumerate(sol.t)]
```

See also: [`reaction_extents`](@ref), [`degrees_of_hydration`](@ref),
[`volume_fractions`](@ref).
"""
function state_at(sol, kp::KineticsProblem, t::Real; ξ = nothing, nsub::Int = 4)
    ξ_t = if ξ === nothing
        reaction_extents(sol, kp; times = [kp.tspan[1], float(t)], nsub = nsub)[end, :]
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
