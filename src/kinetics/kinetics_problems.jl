# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using LinearAlgebra
using SciMLBase

# ── KineticsProblem ───────────────────────────────────────────────────────────

"""
    struct KineticsProblem{CS, CAL, ES, AM}

Encapsulates a kinetics simulation following Leal et al. (2017).

The ODE state vector `u` is structured as:
  - Without re-speciation: `u = [nₖ₁, …, nₖ_K, [T]]`
  - With re-speciation:    `u = [bₑ₁, …, bₑ_C, nₖ₁, …, nₖ_K, [T]]`

where `bₑ` are the element amounts in the equilibrium partition, `nₖ` are
the moles of kinetic species, and `T` is the temperature (semi-adiabatic only).

# Fields

  - `system`: [`ChemicalSystem`](@ref).
  - `kinetic_reactions`: vector of [`KineticReaction`](@ref) objects.
  - `initial_state`: [`ChemicalState`](@ref) providing initial moles, T, P.
  - `tspan`: `(t_start, t_end)` time interval [s].
  - `calorimeter`: `nothing`, [`IsothermalCalorimeter`](@ref),
    or [`SemiAdiabaticCalorimeter`](@ref).
  - `activity_model`: [`AbstractActivityModel`](@ref) for log-activities.
  - `equilibrium_solver`: solver for re-speciation, or `nothing`.
  - `idx_kinetic`: indices of kinetic species in `system.species`.
  - `idx_equilibrium`: indices of equilibrium species.
  - `ν`: stoichiometric matrix (M × N) = `SM.N'` restricted to kinetic reactions.
  - `νe`, `νk`: partitions of `ν` for equilibrium / kinetic species.
  - `Ae`: formula matrix restricted to equilibrium species (C × Nₑ).

See also: [`integrate`](@ref), [`KineticsSolver`](@ref).
"""
struct KineticsProblem{
        CS <: ChemicalSystem,
        KR <: AbstractVector,
        CAL,
        ES,
        AM <: AbstractActivityModel,
    }
    system::CS
    kinetic_reactions::KR
    initial_state::ChemicalState
    tspan::Tuple{Float64, Float64}
    calorimeter::CAL
    activity_model::AM
    equilibrium_solver::ES
    # ── Pre-computed partitions (Leal 2017, Eq. 53) ──
    idx_kinetic::Vector{Int}
    idx_equilibrium::Vector{Int}
    ν::Matrix{Float64}          # (M × N) stoichiometric matrix of kinetic reactions
    νe::Matrix{Float64}         # (M × Nₑ) equilibrium columns
    νk::Matrix{Float64}         # (M × K)  kinetic columns
    Ae::Matrix{Float64}         # (C × Nₑ) formula matrix, equilibrium partition
end

"""
    KineticsProblem(cs, kinetic_reactions, initial_state, tspan; ...) -> KineticsProblem

Construct a [`KineticsProblem`](@ref) from an explicit list of reactions.

Each element of `kinetic_reactions` must be either a [`KineticReaction`](@ref) or
a [`Reaction`](@ref) with a `:rate` entry in its properties. Reaction objects
are automatically wrapped via
`KineticReaction(cs, rxn)`.

    KineticsProblem(cs, initial_state, tspan; ...) -> KineticsProblem

Construct from a [`ChemicalSystem`](@ref) that has `kinetic_species` declared
(reactions and rates auto-generated via the `kinetic_species` keyword).

# Arguments

  - `cs`: [`ChemicalSystem`](@ref).
  - `kinetic_reactions`: `AbstractVector` of [`KineticReaction`](@ref) or
    [`Reaction`](@ref) objects carrying a `:rate` property.
  - `initial_state`: [`ChemicalState`](@ref) providing initial moles, T, P.
  - `tspan`: `(t0, tf)` time interval. Plain `Real` → [s]; `Quantity` → converted.
  - `calorimeter`: `nothing` (no thermal coupling),
    [`IsothermalCalorimeter`](@ref), or [`SemiAdiabaticCalorimeter`](@ref).
  - `activity_model`: activity model for log-activity computation (default: dilute).
  - `equilibrium_solver`: `nothing` (no re-speciation) or an [`EquilibriumSolver`](@ref).

# Examples

```julia
# From explicit reactions
rxn = Reaction(OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 3.33),
               OrderedDict(sp("Jennite") => 0.167, sp("Portlandite") => 1.5))
rxn[:rate] = parrot_killoh(PK_PARAMS_C3S, "C3S"; α_max)
kp = KineticsProblem(cs, [rxn], state0, (0.0, 7 * 86400.0))

# From kinetic_species in ChemicalSystem
cs = ChemicalSystem(species, primaries;
    kinetic_species = Dict("C3S" => pk_C3S, "C2S" => pk_C2S))
kp = KineticsProblem(cs, state0, (0.0, 7 * 86400.0))
```
"""
function _build_kinetics_problem(
        system::ChemicalSystem,
        kin_rxns::AbstractVector{<:KineticReaction},
        initial_state::ChemicalState,
        tspan::Tuple;
        calorimeter = nothing,
        activity_model::AbstractActivityModel = DiluteSolutionModel(),
        equilibrium_solver = nothing,
    )
    n_sp = length(system.species)
    idx_kin = unique!(Int[kr.idx_mineral for kr in kin_rxns])
    idx_eq = setdiff(1:n_sp, idx_kin)
    n_rxn = length(kin_rxns)

    # Stoichiometric matrix ν (M × N) — Leal Eq. 44
    ν = zeros(Float64, n_rxn, n_sp)
    for (i, kr) in enumerate(kin_rxns)
        ν[i, :] .= kr.stoich
    end

    # Partition (Leal Eq. 53): ν = [νₑ  νₖ]
    νe = ν[:, idx_eq]
    νk = ν[:, idx_kin]

    # Formula matrix for equilibrium partition: Aₑ = CSM.A[:, idx_eq]
    # Conservation matrix of the equilibrium partition, taken from the partition
    # SUB-SYSTEM and in the basis the equilibrium solve uses (`SM.A`, with
    # respect to the primaries — not the canonical element matrix `CSM.A`).
    # Conserving the primaries is equivalent to conserving the elements, but the
    # row counts differ, and so do the parent's and the sub-system's: `bₑ` must
    # be built on exactly the matrix the solve is posed on, or every step fails
    # on a dimension mismatch.
    Ae = Float64.(_equilibrium_subsystem(system, idx_eq).SM.A)

    return KineticsProblem{
        typeof(system), typeof(kin_rxns), typeof(calorimeter),
        typeof(equilibrium_solver), typeof(activity_model),
    }(
        system,
        kin_rxns,
        initial_state,
        (Float64(safe_ustrip(us"s", tspan[1])), Float64(safe_ustrip(us"s", tspan[2]))),
        calorimeter,
        activity_model,
        equilibrium_solver,
        collect(Int, idx_kin),
        collect(Int, idx_eq),
        ν, νe, νk, Ae,
    )
end

# 4-argument form: explicit reaction list (Reaction or KineticReaction)
function KineticsProblem(
        system::ChemicalSystem,
        kinetic_reactions::AbstractVector,
        initial_state::ChemicalState,
        tspan::Tuple;
        calorimeter = nothing,
        activity_model::AbstractActivityModel = DiluteSolutionModel(),
        equilibrium_solver = nothing,
    )
    kin_rxns = [
        r isa KineticReaction ? r : KineticReaction(system, r)
            for r in kinetic_reactions
    ]
    return _build_kinetics_problem(
        system, kin_rxns, initial_state, tspan;
        calorimeter, activity_model, equilibrium_solver,
    )
end

# 3-argument form: reactions from ChemicalSystem.reactions (kinetic_species API)
function KineticsProblem(
        system::ChemicalSystem,
        initial_state::ChemicalState,
        tspan::Tuple;
        calorimeter = nothing,
        activity_model::AbstractActivityModel = DiluteSolutionModel(),
        equilibrium_solver = nothing,
    )
    isempty(system.idx_kinetic) &&
        throw(
        ArgumentError(
            "ChemicalSystem has no kinetic species. " *
                "Pass kinetic_species to the ChemicalSystem constructor, " *
                "or use the 4-argument form KineticsProblem(cs, reactions, state, tspan)."
        )
    )
    kin_rxns = [KineticReaction(system, rxn) for rxn in system.reactions]
    return _build_kinetics_problem(
        system, kin_rxns, initial_state, tspan;
        calorimeter, activity_model, equilibrium_solver,
    )
end

"""
    _with_equilibrium_solver(kp::KineticsProblem, es) -> KineticsProblem

Return `kp` with its equilibrium solver replaced by `es`, or `kp` itself when
`es` is `nothing`. Used by `integrate` so that a solver passed on the
[`KineticsSolver`](@ref) reaches the ODE; a solver already set on the problem
wins when both are given, and the mismatch is reported.
"""
function _with_equilibrium_solver(kp::KineticsProblem, es)
    isnothing(es) && return kp
    if !isnothing(kp.equilibrium_solver)
        kp.equilibrium_solver === es || @warn(
            "an equilibrium solver is set on both the KineticsProblem and the " *
                "KineticsSolver; the one on the problem is used."
        )
        return kp
    end
    return KineticsProblem{
        typeof(kp.system), typeof(kp.kinetic_reactions), typeof(kp.calorimeter),
        typeof(es), typeof(kp.activity_model),
    }(
        kp.system, kp.kinetic_reactions, kp.initial_state, kp.tspan,
        kp.calorimeter, kp.activity_model, es,
        kp.idx_kinetic, kp.idx_equilibrium, kp.ν, kp.νe, kp.νk, kp.Ae,
    )
end

# ── build_u0 ─────────────────────────────────────────────────────────────────

"""
    build_u0(kp::KineticsProblem) -> Vector{Float64}

Build the initial ODE state vector.

Structure of `u`:
  - Without re-speciation: `u = [nₖ₁, …, nₖ_K]`
  - With re-speciation:    `u = [bₑ₁, …, bₑ_C, nₖ₁, …, nₖ_K]`
  - Semi-adiabatic adds `T` at the end: `u = [..., T₀]`
"""
function build_u0(kp::KineticsProblem)
    n_mol = Float64[
        ustrip(us"mol", kp.initial_state.n[i])
            for i in eachindex(kp.system.species)
    ]
    # Kinetic species moles
    nk0 = n_mol[kp.idx_kinetic]

    u0 = if isnothing(kp.equilibrium_solver)
        copy(nk0)
    else
        # Element amounts in equilibrium partition: bₑ = Aₑ nₑ
        ne0 = n_mol[kp.idx_equilibrium]
        be0 = kp.Ae * ne0
        vcat(be0, nk0)
    end

    # Append temperature for semi-adiabatic calorimeter
    if kp.calorimeter isa SemiAdiabaticCalorimeter
        push!(u0, Float64(safe_ustrip(us"K", kp.calorimeter.T0)))
    end

    return u0
end

# ── build_kinetics_params ────────────────────────────────────────────────────

"""
    build_kinetics_params(kp::KineticsProblem; ϵ=1e-30) -> NamedTuple

Build the immutable parameter tuple `p` passed to the ODE function.

Key fields: `T`, `P`, `ϵ`, `lna_fn`, `kin_rxns`, `species_index`,
`n_initial_full`, `n_full`, `cp_fns`, `rates_buf`, index ranges
`n_be`, `n_nk`, `idx_kinetic`, `idx_equilibrium`, `νe`, `νk`, `Ae`.
"""
function build_kinetics_params(kp::KineticsProblem; ϵ::Float64 = 1.0e-30)
    state = kp.initial_state
    T_K = Float64(ustrip(us"K", temperature(state)))
    P_Pa = Float64(ustrip(us"Pa", pressure(state)))

    lna_fn = activity_model(kp.system, kp.activity_model)

    # Species name → index dict (built once, shared by all StateViews)
    species_index = Dict{String, Int}()
    for (i, sp) in enumerate(kp.system.species)
        species_index[phreeqc(formula(sp))] = i
        sym = ChemistryLab.symbol(sp)
        !isempty(sym) && (species_index[sym] = i)
    end

    n_sp = length(kp.system.species)
    n_initial_full = Float64[ustrip(us"mol", state.n[i]) for i in 1:n_sp]
    n_full = copy(n_initial_full)

    cp_fns = [haskey(sp, :Cp⁰) ? sp[:Cp⁰] : nothing for sp in kp.system.species]

    kin_rxns = kp.kinetic_reactions
    rates_buf = zeros(Float64, length(kin_rxns))

    eq_sys, eq_sub, n_eq_init = if isnothing(kp.equilibrium_solver)
        nothing, nothing, Float64[]
    else
        sys_e = _equilibrium_subsystem(kp.system, kp.idx_equilibrium)
        es = kp.equilibrium_solver
        (
            sys_e,
            EquilibriumSolver(
                sys_e, kp.activity_model, es.solver;
                variable_space = es.variable_space, es.kwargs...
            ),
            Float64[n_initial_full[i] for i in kp.idx_equilibrium],
        )
    end

    # State layout sizes
    n_be = isnothing(kp.equilibrium_solver) ? 0 : size(kp.Ae, 1)
    n_nk = length(kp.idx_kinetic)
    has_T = kp.calorimeter isa SemiAdiabaticCalorimeter

    # Calorimeter parameters (semi-adiabatic)
    cal = kp.calorimeter
    Cp_calo = cal isa SemiAdiabaticCalorimeter ? Float64(safe_ustrip(us"J/K", cal.Cp)) : 0.0
    T_env = cal isa SemiAdiabaticCalorimeter ? Float64(safe_ustrip(us"K", cal.T_env)) : T_K
    heat_loss_fn = cal isa SemiAdiabaticCalorimeter ? cal.heat_loss : identity

    return (
        T = T_K,
        P = P_Pa,
        ϵ = ϵ,
        lna_fn = lna_fn,
        kin_rxns = kin_rxns,
        species_index = species_index,
        n_initial_full = n_initial_full,
        n_full = n_full,
        cp_fns = cp_fns,
        rates_buf = rates_buf,
        # Index layout
        n_be = n_be,
        n_nk = n_nk,
        has_T = has_T,
        idx_kinetic = kp.idx_kinetic,
        idx_equilibrium = kp.idx_equilibrium,
        # Leal partitions
        νe = Float64.(kp.νe),
        νk = Float64.(kp.νk),
        Ae = Float64.(kp.Ae),
        # Calorimeter
        Cp_calo = Cp_calo,
        T_env = T_env,
        heat_loss_fn = heat_loss_fn,
        # Equilibrium — Leal et al. (2017) §5. The re-speciation φ(bₑ) is a
        # minimization over the EQUILIBRIUM PARTITION ONLY, at frozen kinetic
        # amounts. Running it over the whole system would let the kinetic
        # minerals equilibrate instantaneously, which is exactly what a kinetic
        # description exists to prevent.
        eq_system = eq_sys,
        eq_solver = eq_sub,
        n_eq_init = n_eq_init,
        n_eq_buf = similar(n_eq_init),
        xi_buf = zeros(Float64, length(kp.idx_kinetic)),
        T_q = Ref(temperature(state)),
        P_q = Ref(pressure(state)),
        # Pseudo-inverse of Aₑ, to project the previous speciation back onto the
        # element amounts the ODE carries. Built once — it is a fixed matrix.
        Ae_pinv = isnothing(kp.equilibrium_solver) ?
            zeros(Float64, 0, 0) : pinv(Float64.(kp.Ae)),
        state_ref = Ref{ChemicalState}(state),
        eq_failures = Ref(0),
    )
end

# ── Equilibrium partition sub-system ─────────────────────────────────────────

"""
    _equilibrium_subsystem(system, idx_equilibrium) -> ChemicalSystem

The chemical system restricted to the equilibrium partition, as the partitioned
formulation of [Leal2017](@cite) requires.

Its formula matrix is exactly `system.CSM.A[:, idx_equilibrium]`, in the same
species order, so the element amounts `bₑ` carried by the ODE state are handed
to it unchanged. The primaries are those of the parent system that survive the
restriction — the kinetic minerals never do, not being in the partition.
"""
function _equilibrium_subsystem(system::ChemicalSystem, idx_equilibrium)
    sub_species = system.species[idx_equilibrium]
    # The parent's primaries are the row labels of its stoichiometric matrix,
    # `system.SM.primaries` — not `idx_components`, which is a class-based index
    # and returns a different, much smaller set. Getting this wrong leaves the
    # sub-system with only a couple of conservation constraints for a dozen
    # species: the minimization is then wildly under-determined and returns an
    # assemblage with essentially no water left.
    comp_names = Set(symbol(sp) for sp in system.SM.primaries)
    prim = [sp for sp in sub_species if symbol(sp) in comp_names]
    return ChemicalSystem(sub_species, isempty(prim) ? sub_species : prim)
end

# ── respeciate! ──────────────────────────────────────────────────────────────

"""
    respeciate!(p, u) -> Bool

Solve the equilibrium sub-problem once, and write the result into the running
composition `p.n_full`. Returns `true` when a solve actually happened.

This is the second half of the operator-splitting step: the ODE advances the
kinetic minerals with the speciation held frozen, then this function
re-equilibrates the equilibrium partition under the element amounts the ODE has
just produced.

The element amounts `bₑ` carried by the state vector are the constraint of that
sub-problem (Leal et al. 2017, Eq. 54). `solve` conserves `A·n`, so what has to
be handed to it is a composition whose element totals are exactly `bₑ` — here
the previous speciation, projected onto `bₑ` through the pseudo-inverse of
`Aₑ`. Handing over `p.n_full` unchanged, as an earlier version did, discards
`bₑ` entirely and leaves the element balance to drift.
"""
function respeciate!(p, u)
    p.n_be > 0 || return false

    # φ(bₑ), Leal et al. (2017) Eq. 54: the element amounts carried by the ODE
    # state ARE the constraint of the minimization, and they are handed to the
    # solver as `b`.
    #
    # That is the whole reason the state integrates `bₑ` rather than `nₑ`. Along
    # the way an individual species may want to go negative — the generated
    # dissolution reactions are written in H⁺, and a cement paste contains no
    # acid — and it is the minimizer, not the caller, that redistributes the
    # elements over a feasible set. An earlier version reconstructed `nₑ` from
    # `bₑ` through `pinv(Aₑ)` and clamped the result at `ϵ`; the clamp destroyed
    # the balance the projection had just established, and the solve went on to
    # return amounts of 1e65.
    #
    # The composition below is a starting guess only, and does not have to carry
    # `bₑ`. It is built from the reaction extents, which come free from the
    # kinetic amounts: each reaction carries ν = −1 on its controlling mineral
    # (`_normalise_to_mineral!`), so ξⱼ = nₖⱼ(0) − nₖⱼ.
    be = collect(@view u[1:(p.n_be)])

    nk = @view u[(p.n_be + 1):(p.n_be + p.n_nk)]
    ξ = p.xi_buf
    for (j, idx) in enumerate(p.idx_kinetic)
        ξ[j] = p.n_initial_full[idx] - nk[j]
    end
    n_eq = p.n_eq_buf
    mul!(n_eq, p.νe', ξ)
    for j in eachindex(n_eq)
        n_eq[j] = max(n_eq[j] + p.n_eq_init[j], p.ϵ)
    end

    state_eq = ChemicalState(p.eq_system, n_eq .* u"mol"; T = p.T_q[], P = p.P_q[])

    # `p.eq_solver` is a prebuilt `EquilibriumSolver` over the partition — a
    # solver *object*, not a SciML algorithm — so `solve`, not `equilibrate`.
    local eq_result
    try
        eq_result = SciMLBase.solve(p.eq_solver, state_eq; ϵ = p.ϵ, b = be)
    catch err
        p.eq_failures[] += 1
        if p.eq_failures[] == 1
            @warn """re-speciation failed; the composition is left frozen for \
            this step. Later failures are counted, not reported.""" exception = err
        end
        return false
    end

    for (j, idx) in enumerate(p.idx_equilibrium)
        p.n_full[idx] = ustrip(us"mol", eq_result.n[j])
    end
    return true
end

# ── build_kinetics_ode ───────────────────────────────────────────────────────

"""
    build_kinetics_ode(kp::KineticsProblem) -> Function

Build the ODE right-hand-side `f!(du, u, p, t)` implementing Leal et al. (2017).

State layout:
  - `u[1:n_be]`         = bₑ (element amounts in equilibrium partition)
  - `u[n_be+1:n_be+n_nk]` = nₖ (moles of kinetic species)
  - `u[end]`            = T  (semi-adiabatic only)

ODE equations (Leal 2017, Eq. 66):
  - `dnₖ/dt = νₖᵀ r`
  - `dbₑ/dt = Aₑ νₑᵀ r`
  - `dT/dt  = (q̇ − φ(ΔT)) / Cp_total`  (semi-adiabatic)

where `nₑ = φ(bₑ)` is the equilibrium re-speciation constraint.
"""
function build_kinetics_ode(kp::KineticsProblem)
    function f!(du, u, p, t)
        T_elt = eltype(u)

        # ── 1. Extract state components ──────────────────────────────────
        nk = @view u[(p.n_be + 1):(p.n_be + p.n_nk)]
        T_curr = p.has_T ? u[end] : p.T

        # ── 2. Reconstruct full mole vector ──────────────────────────────
        if T_elt === Float64
            n_full = p.n_full
        else
            n_full = T_elt.(p.n_full)
        end

        # 2a. Kinetic species from nₖ
        for (j, idx) in enumerate(p.idx_kinetic)
            n_full[idx] = max(nk[j], p.ϵ)
        end

        # 2b. Equilibrium species: read the *frozen* speciation.
        #
        # The equilibrium sub-problem is NOT solved here. It is solved once per
        # accepted step, by `respeciate!` below, in an operator-splitting step.
        # Two reasons, and both matter:
        #
        #  * the right-hand side of a stiff solver is evaluated many times per
        #    step and differentiated for the Jacobian; solving an optimization
        #    problem inside it makes the cost unpredictable, and — as long as
        #    `ChemicalState` stores `Float64` moles — a `Dual` cannot even be
        #    written into the state, so the solve would be skipped exactly on
        #    the evaluations that build the Jacobian;
        #  * with the solve outside, the same speciation is seen by the residual
        #    and by its Jacobian, whatever the number type of `u`.
        #
        # `p.n_full` already carries the equilibrium partition as left by the
        # last `respeciate!`, so nothing has to be copied here.

        # ── 3. Compute log-activities ────────────────────────────────────
        lna = p.lna_fn(n_full, p)

        # ── 4. Build StateViews (O(1) named access) ─────────────────────
        n_sv = StateView(n_full, p.species_index)
        lna_sv = StateView(lna, p.species_index)
        n0_sv = StateView(p.n_initial_full, p.species_index)

        # ── 5. Evaluate kinetic rates r(T, P, t, n, lna, n₀) ────────────
        n_rxn = length(p.kin_rxns)
        rates = Vector{T_elt}(undef, n_rxn)
        for (i, kr) in enumerate(p.kin_rxns)
            rates[i] = kr.rate_fn(T_curr, p.P, t, n_sv, lna_sv, n0_sv)
            if T_elt === Float64
                p.rates_buf[i] = rates[i]
            end
        end

        # ── 6. ODE: dnₖ/dt = νₖᵀ r (Leal Eq. 56) ───────────────────────
        fill!(du, zero(T_elt))
        du_nk = p.νk' * rates
        for j in 1:(p.n_nk)
            du[p.n_be + j] = du_nk[j]
        end

        # ── 7. ODE: dbₑ/dt = Aₑ νₑᵀ r (Leal Eq. 65) ────────────────────
        if p.n_be > 0
            du_be = p.Ae * (p.νe' * rates)
            for j in 1:(p.n_be)
                du[j] = du_be[j]
            end
        end

        # ── 8. ODE: dT/dt = (q̇ − φ(ΔT)) / Cp_total (semi-adiabatic) ───
        if p.has_T
            # Heat generation: q̇ = Σᵢ rᵢ × (−ΔᵣH⁰ᵢ) [W]
            qdot = heat_rate(p.kin_rxns, rates, T_curr)

            # Total heat capacity: Cp_calo + Σᵢ nᵢ Cp°ᵢ(T)
            Cp_total = p.Cp_calo
            for (i, cp_fn) in enumerate(p.cp_fns)
                isnothing(cp_fn) && continue
                cp_i = cp_fn(; T = T_curr, unit = false)
                Cp_total = Cp_total + n_full[i] * cp_i
            end

            ΔT = T_curr - p.T_env
            du[end] = (qdot - p.heat_loss_fn(ΔT)) / Cp_total
        end

        return nothing
    end

    return f!
end
