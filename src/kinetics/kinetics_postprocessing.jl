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
