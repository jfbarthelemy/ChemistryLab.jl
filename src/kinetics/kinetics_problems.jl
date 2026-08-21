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
  - Without re-speciation: `u = [nₖ₁, …, nₖ_K, ξ₁, …, ξ_M, [T | Q]]`
  - With re-speciation:    `u = [bₑ₁, …, bₑ_C, nₖ₁, …, nₖ_K, ξ₁, …, ξ_M, [T | Q]]`

where `bₑ` are the element amounts in the equilibrium partition, `nₖ` the moles
of kinetic species and `ξ` the extents of the kinetic reactions. The trailing
slot is present only with a calorimeter: the temperature for
[`SemiAdiabaticCalorimeter`](@ref), the accumulated heat for
[`IsothermalCalorimeter`](@ref).

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
  - Without re-speciation: `u = [nₖ₁, …, nₖ_K, ξ₁, …, ξ_M]`
  - With re-speciation:    `u = [bₑ₁, …, bₑ_C, nₖ₁, …, nₖ_K, ξ₁, …, ξ_M]`
  - Semi-adiabatic adds `T` at the end: `u = [..., T₀]`

The extents of reaction `ξ` are carried alongside the kinetic moles, integrating
`dξ/dt = r`. They are redundant with `nₖ` — the two satisfy
`nₖ = nₖ(0) + νₖᵀ ξ` — but they are what makes the **non-kinetic** amounts
available inside the residual, through `n = n(0) + νᵀ ξ`. Without them, a rate
law gating on a species that is not itself kinetic would read a value frozen at
`t = 0`. They also make [`reaction_extents`](@ref) and [`state_at`](@ref) exact
rather than quadrature-limited.

The calorimeter's slot stays last and is addressed from the end of the vector,
so it is unaffected by the presence of `ξ`.
"""
function build_u0(kp::KineticsProblem)
    n_mol = Float64[
        ustrip(us"mol", kp.initial_state.n[i])
            for i in eachindex(kp.system.species)
    ]
    # Kinetic species moles
    nk0 = n_mol[kp.idx_kinetic]
    # Extents of reaction, zero by definition at t = tspan[1]
    ξ0 = zeros(Float64, length(kp.kinetic_reactions))

    u0 = if isnothing(kp.equilibrium_solver)
        vcat(nk0, ξ0)
    else
        # Element amounts in equilibrium partition: bₑ = Aₑ nₑ
        ne0 = n_mol[kp.idx_equilibrium]
        be0 = kp.Ae * ne0
        vcat(be0, nk0, ξ0)
    end

    # One trailing slot per calorimeter, and only one: the temperature for the
    # semi-adiabatic device, the accumulated heat for the isothermal one. Both are
    # `u[end]`, and `p.has_T` / `p.has_Q` say which.
    if kp.calorimeter isa SemiAdiabaticCalorimeter
        push!(u0, Float64(safe_ustrip(us"K", kp.calorimeter.T0)))
    elseif kp.calorimeter isa IsothermalCalorimeter
        push!(u0, 0.0)
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
    h_fns = [haskey(sp, :ΔₐH⁰) ? sp[:ΔₐH⁰] : nothing for sp in kp.system.species]

    kin_rxns = kp.kinetic_reactions
    rates_buf = zeros(Float64, length(kin_rxns))

    eq_sys, eq_sub, n_eq_init = if isnothing(kp.equilibrium_solver)
        nothing, nothing, Float64[]
    else
        sys_e = _equilibrium_subsystem(kp.system, kp.idx_equilibrium)
        es = kp.equilibrium_solver
        # The sub-solver must be rebuilt because `sys_e` is a different system
        # from `kp.system`, so the compiled `μ` does not carry over. It is
        # rebuilt with the model the USER's solver was built from — not with
        # `kp.activity_model`, which defaults to `DiluteSolutionModel()` and
        # would silently downgrade an HKF run on a cement pore solution.
        es_model = activity_model(es)
        if typeof(es_model) !== typeof(kp.activity_model)
            @warn """
            The equilibrium solver and the kinetics problem carry different \
            activity models. The equilibrium sub-solve uses the solver's \
            ($(nameof(typeof(es_model)))), while the log-activities handed to the \
            rate laws use the problem's ($(nameof(typeof(kp.activity_model)))). \
            Pass `activity_model = $(nameof(typeof(es_model)))()` to \
            `KineticsProblem` to make them agree.""" maxlog = 1
        end
        (
            sys_e,
            EquilibriumSolver(
                sys_e, es_model, es.solver;
                variable_space = es.variable_space, es.kwargs...
            ),
            Float64[n_initial_full[i] for i in kp.idx_equilibrium],
        )
    end

    # State layout sizes
    n_be = isnothing(kp.equilibrium_solver) ? 0 : size(kp.Ae, 1)
    n_nk = length(kp.idx_kinetic)
    n_rxn_state = length(kp.kinetic_reactions)
    has_T = kp.calorimeter isa SemiAdiabaticCalorimeter
    has_Q = kp.calorimeter isa IsothermalCalorimeter

    # Calorimeter parameters (semi-adiabatic)
    cal = kp.calorimeter
    Cp_calo = cal isa SemiAdiabaticCalorimeter ? Float64(safe_ustrip(us"J/K", cal.Cp)) : 0.0
    T_env = cal isa SemiAdiabaticCalorimeter ? Float64(safe_ustrip(us"K", cal.T_env)) : T_K
    heat_loss_fn = cal isa SemiAdiabaticCalorimeter ? cal.heat_loss : identity

    return (
        T = T_K,
        P = P_Pa,
        ϵ = ϵ,
        has_Q = has_Q,
        lna_fn = lna_fn,
        kin_rxns = kin_rxns,
        species_index = species_index,
        n_initial_full = n_initial_full,
        n_full = n_full,
        cp_fns = cp_fns,
        h_fns = h_fns,
        rates_buf = rates_buf,
        # Index layout
        n_be = n_be,
        n_nk = n_nk,
        n_rxn_state = n_rxn_state,
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
        n_eq_buf2 = similar(n_eq_init),
        xi_buf = zeros(Float64, length(kp.idx_kinetic)),
        T_q = Ref(temperature(state)),
        P_q = Ref(pressure(state)),
        # Pseudo-inverse of Aₑ, to project the previous speciation back onto the
        # element amounts the ODE carries. Built once — it is a fixed matrix.
        Ae_pinv = isnothing(kp.equilibrium_solver) ?
            zeros(Float64, 0, 0) : pinv(Float64.(kp.Ae)),
        state_ref = Ref{ChemicalState}(state),
        eq_failures = Ref(0),
        # Worst |Aₑ n − bₑ|∞ over the run. This, not the optimizer's retcode, is
        # the criterion with physical meaning: a solve can stall short of its
        # tolerance and still satisfy the element balance to machine precision.
        eq_worst_residual = Ref(0.0),
        eq_worst_abs = Ref(0.0),
        eq_worst_abs_acc = Ref(0.0),
        on_accepted = Ref(false),
        # Set once a speciation exists, so `respeciate!` can warm-start from it.
        eq_warm = Ref(false),
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

    # Carry the solid solutions over. Dropping them — as this did until 0.8.2 —
    # is silent and total: the parent may declare CSHQ, AFm or Hydrogarnet, and
    # the partition the equilibrium is actually solved on knows nothing of them,
    # so their end-members are treated as separate pure phases and the mixing
    # entropy never enters the Gibbs energy. Measured on alite and belite with
    # the four CSHQ end-members, the run was bit-identical with and without the
    # declaration, and produced no C-S-H at all: the silicon stayed in solution
    # and the portlandite came out at 4.52 mol against 2.93 with a Jennite
    # end-member.
    #
    # A solution survives only if ALL its end-members are in the partition. One
    # whose members are split between the kinetic and equilibrium sides is not a
    # phase the equilibrium can mix, and passing it truncated would be worse than
    # dropping it.
    sub_names = Set(symbol(sp) for sp in sub_species)
    ss = system.solid_solutions
    sub_ss = if ss === nothing
        nothing
    else
        kept = [
            phase for phase in ss
                if all(em -> symbol(em) in sub_names, phase.end_members)
        ]
        isempty(kept) ? nothing : kept
    end

    return ChemicalSystem(
        sub_species, isempty(prim) ? sub_species : prim;
        solid_solutions = sub_ss,
    )
end

"""
    _EQ_GUESS_FLOOR

Lower floor applied to the starting guess of the equilibrium sub-solve, chosen
strictly above the `1e-16` lower bound that `EquilibriumProblem` imposes.

An interior-point method started on its own bound stalls short of its tolerance
and returns `MaxIters`. Loosening that tolerance is *not* the fix: measured
against Reaktoro on the calcite reference case, the worst species error grows
from 4.3 % at `tol = 1e-10` to 38 % at `1e-8` and 252 % at `1e-7`. Moving the
guess inside keeps the tight tolerance and removes the stalls.
"""
const _EQ_GUESS_FLOOR = 1.0e-10

"""
    EQ_RESIDUAL_TOL

Largest per-element balance violation, relative to that element's own total,
that still counts as a converged speciation. Above it the point is not handed
on as a warm start and the solve is retried from a budget-clipped guess.
"""
const EQ_RESIDUAL_TOL = 1.0e-8

"""
    _RETRY_ABS_TOL

Element-balance violation in moles above which a replayed speciation is solved a
second time from a guess carrying no active set. Chosen well above machine
precision and well below anything chemically meaningful.
"""
const _RETRY_ABS_TOL = 1.0e-6

"""
    _CONTINUATION_STEPS

How many bisection steps a replay may take between the last certified instant and
one it cannot certify directly. Each step halves the jump in the component
totals, so eight of them reduce it by a factor 256.
"""
const _CONTINUATION_STEPS = 8

"""
    RESTORE_MAXIT

Alternating-projection sweeps allowed when restoring the feasibility of an
in-run guess. Exposed because the right value is a trade: the projection
converges linearly and a cement can need tens of thousands of sweeps, while this
runs at every right-hand-side evaluation.
"""
const RESTORE_MAXIT = Ref(200)

# ── enthalpy tracking ────────────────────────────────────────────────────────

"""
    system_enthalpy(p, u, T) -> Float64

`Σᵢ nᵢ ΔₐH⁰ᵢ(T)` [J] over the composition of the ODE state `u`: the kinetic
minerals read from `u` itself, the equilibrium partition from `p.n_full`.

# Why the heat cannot come from the kinetic reactions here

`heat_rate` sums `rᵢ (−ΔᵣH⁰ᵢ)` over the KINETIC reactions, which is right when
those reactions produce the hydrates. Under partial equilibrium they do not: they
dissolve the anhydrous phases into ions, and the hydrates are precipitated by the
Gibbs minimization, whose heat that sum cannot see. Measured on an ordinary
Portland cement, counting only the dissolution put the semi-adiabatic temperature
rise at 207 K where the test gives some tens of kelvin.

The enthalpy of the whole system has no such blind spot. It is a state function,
so the heat released between two states at the same temperature is their
difference — reactants, ions and hydrates all counted once, with no reaction
stoichiometry to write down. This is Eq. (17)–(21) of Lavergne et al. (2018).
"""
function system_enthalpy(p, u, T)
    # The kinetic amounts are read from the ODE state, NOT from `p.n_full`.
    #
    # `respeciate!` writes only the equilibrium partition into that buffer; the
    # anhydrous amounts in it were last written by the right-hand side, which the
    # stiff solver also evaluates at Jacobian probes and rejected steps. Summing
    # the buffer as it stands therefore pairs an accepted equilibrium with some
    # neighboring point's clinker, and the resulting enthalpy is not a state of
    # the trajectory at all: the recorded heat came out NON-MONOTONE, 936 J/g at
    # one day and 631 J/g at two, which no calorimeter has ever measured.
    kin = p.idx_kinetic
    H = 0.0
    @inbounds for (i, h_fn) in enumerate(p.h_fns)
        isnothing(h_fn) && continue
        j = findfirst(==(i), kin)
        nᵢ = j === nothing ? p.n_full[i] : max(u[p.n_be + j], p.ϵ)
        H += nᵢ * h_fn(; T = T, unit = false)
    end
    return H
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
    n_eq = p.n_eq_buf

    # WARM START. `p.n_full` holds the speciation left by the previous accepted
    # step, which is an equilibrium for a nearby `bₑ` — by far the best guess
    # available. The stoichiometric reconstruction below is only the cold start,
    # used on the first call before any speciation exists.
    #
    # This matters far more than it looks. The reconstruction places every
    # dissolved element in solution with ZERO hydrates, i.e. a wildly
    # supersaturated composition, which for an aqueous-only system like the
    # calcite reference is harmless and for a cement is nearly the worst possible
    # starting point: more than half the solves failed, no hydrate ever
    # precipitated, and the pore solution came out at pH 6 instead of 12.6.
    if p.eq_warm[]
        for (j, idx) in enumerate(p.idx_equilibrium)
            n_eq[j] = max(p.n_full[idx], _EQ_GUESS_FLOOR)
        end
        # The previous speciation was an equilibrium for the PREVIOUS `bₑ`. When
        # an element has since been spent — the sulfate of an OPC once the gypsum
        # is gone — that guess demands more of it than now exists, and the solve
        # starts outside the feasible set. Clipping costs nothing when the guess
        # is already feasible, which is the ordinary case.
        _budget_clip!(n_eq, p.Ae, be)
        _restore_feasibility!(n_eq, p.Ae, be; maxit = RESTORE_MAXIT[])
        _respeciate_solve!(p, n_eq, be) && return true
        # Infeasible or stalled: fall through to the cold reconstruction rather
        # than carry the bad point into the next step.
    end

    ξ = p.xi_buf
    for (j, idx) in enumerate(p.idx_kinetic)
        ξ[j] = p.n_initial_full[idx] - nk[j]
    end
    for j in eachindex(n_eq)
        # Start from the composition the specimen was cast with — for a paste,
        # the mixing water and nothing precipitated — and let
        # `_restore_feasibility!` below carry it onto the current `bₑ`.
        #
        # The stoichiometric reconstruction `νₑᵀξ` that used to be added here
        # placed every dissolved element in solution with ZERO hydrates. For an
        # aqueous-only system that is harmless, but for a cement it is close to
        # the worst possible start: it is supersaturated in every phase at once,
        # and its H⁺ entry is strongly negative (−6 per mole of alite) so it is
        # clamped to the floor, losing the acidity that the hydroxides have to
        # balance. Started there, the back-end stops next to its own guess and
        # returned an assemblage demanding 174 % of the sulfate present.
        #
        # Floor strictly inside the box, not at `p.ϵ`: `EquilibriumProblem`
        # raises anything below 1e-16 to exactly its lower bound, and an
        # interior-point method started on its own bound stalls — the calcite
        # reference reported 6 non-converged solves out of 8 steps for that
        # reason alone.
        n_eq[j] = max(p.n_eq_init[j], _EQ_GUESS_FLOOR)
    end
    _budget_clip!(n_eq, p.Ae, be)
    _restore_feasibility!(n_eq, p.Ae, be; maxit = RESTORE_MAXIT[])

    return _respeciate_solve!(p, n_eq, be)
end

"""
    _respeciate_solve!(p, n_eq, be) -> Bool

Solve `φ(bₑ)` from the guess `n_eq`, write the result into `p.n_full`, and record
the element-balance residual. Returns `false` if the solve threw.
"""
function _respeciate_solve!(p, n_eq, be)
    ok, n_e, abs_res = _one_speciation(p, n_eq, be)
    ok || return false

    # RETRY FROM AN INDEPENDENT GUESS. The warm start carries the previous
    # active set, and where the assemblage switches it is the wrong one — an
    # interior-point method started inside it does not cross over. A guess built
    # from the cast composition carries no active set at all, so it can land in a
    # different basin; keeping whichever conserves matter better costs one extra
    # solve, and only at the steps that need it.
    #
    # This is decided on the ABSOLUTE balance, in moles, because that is the
    # quantity with a meaning: it is the matter the composition fails to account
    # for. The relative measure is for reporting.
    if abs_res > _RETRY_ABS_TOL
        cold = p.n_eq_buf2
        @inbounds for j in eachindex(cold)
            cold[j] = max(p.n_eq_init[j], _EQ_GUESS_FLOOR)
        end
        _budget_clip!(cold, p.Ae, be)
        _restore_feasibility!(cold, p.Ae, be; maxit = RESTORE_MAXIT[])
        ok2, n2, abs2 = _one_speciation(p, cold, be)
        if ok2 && abs2 < abs_res
            n_e, abs_res = n2, abs2
        end
    end

    for (j, idx) in enumerate(p.idx_equilibrium)
        p.n_full[idx] = n_e[j]
    end

    res = _row_residual(p.Ae, n_e, be)
    abs_res > p.eq_worst_abs[] && (p.eq_worst_abs[] = abs_res)

    # Separate the trajectory from the probes. The right-hand side is evaluated
    # far more often than the solution advances — Jacobian finite differences and
    # rejected steps included — and a poor speciation on a PERTURBED `bₑ` never
    # enters the result. Reporting the worst over all evaluations alarmed about
    # something that does not affect the answer: on a full OPC that figure was
    # 1.13 mol while the worst over the 201 accepted steps was 4.3e-4, with a
    # median of 4.8e-9.
    p.on_accepted[] && abs_res > p.eq_worst_abs_acc[] && (p.eq_worst_abs_acc[] = abs_res)

    # Warm-starting from an INFEASIBLE point locks the error in: the next step
    # starts where this one ended, so a single bad solve poisons every solve
    # after it. Only hand over a speciation that actually satisfies the element
    # balance; otherwise leave `eq_warm` as it was and let the caller retry.
    res <= EQ_RESIDUAL_TOL && (p.eq_warm[] = true)
    res > p.eq_worst_residual[] && (p.eq_worst_residual[] = res)
    return res <= EQ_RESIDUAL_TOL
end

"""
    _one_speciation(p, guess, be) -> (ok, n_e, abs_res)

One equilibrium solve of `φ(bₑ)` from `guess`. Returns `ok = false` only if the
solve threw; a solve that returns a poor composition still returns `true`, with
its element-balance violation in moles, so the caller can compare attempts.
"""
function _one_speciation(p, guess, be)
    state_eq = ChemicalState(p.eq_system, guess .* u"mol"; T = p.T_q[], P = p.P_q[])

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
        return false, Float64[], Inf
    end

    n_e = [ustrip(us"mol", x) for x in eq_result.n]
    return true, n_e, _abs_residual(p.Ae, n_e, be)
end


"""
    _row_residual(Ae, n_e, be) -> Float64

Largest element-balance violation `|Aₑnₑ − bₑ|` relative to that element's own
budget. Unlike a single global scale it cannot hide a small element behind a
large one, which is what let a 0.465 mol sulfur violation report as 1.4e-2.
"""
function _row_residual(Ae, n_e, be)
    r = Ae * n_e .- be
    # An element whose total is a millionth of the largest is not tracked
    # meaningfully, and judging it against its own vanishing budget turns a
    # rounding error into an alarm: a full OPC balanced to 1.4e-10 mol at 28 days
    # was reported at 3.2e-2, and the near-empty rows of the first steps saturated
    # the measure at 1.0. Floor every row at a fraction of the system scale.
    scale = 1.0e-6 * maximum(abs, be; init = 0.0)
    worst = 0.0
    for i in eachindex(r)
        flux = zero(eltype(r))
        for j in eachindex(n_e)
            flux += abs(Ae[i, j] * n_e[j])
        end
        worst = max(worst, abs(r[i]) / max(scale, abs(be[i]), flux, 1.0e-30))
    end
    return worst
end

"""
    _abs_residual(Ae, n_e, be) -> Float64

Largest element-balance violation in moles. Reported alongside the relative
measure because it is the one a chemist can judge: 1e-10 mol is machine
precision whatever the system, and 7e-2 mol is not.
"""
_abs_residual(Ae, n_e, be) = maximum(abs, Ae * n_e .- be; init = 0.0)

"""
    _budget_clip!(n_eq, Ae, be)

Clip a starting guess to the element budget: no species may exceed what the
totals `bₑ` can supply, `nⱼ ≤ minᵢ bᵢ/Aᵢⱼ` over the rows it consumes.

This does not change the feasible set — it only moves the guess into it. It is
what unblocks the OPC case: once the sulfate is spent the warm start still
carried ettringite at the aluminum budget, three times the sulfur available,
and the interior-point solve could not walk back from there.
"""
function _budget_clip!(n_eq, Ae, be)
    for j in eachindex(n_eq)
        cap = Inf
        for i in axes(Ae, 1)
            a = Ae[i, j]
            # Only rows with a POSITIVE total bound a species from above; a
            # negative total (H⁺ in a hydrating cement) is met by the hydroxides
            # and bounds nothing.
            (a > 0 && be[i] >= 0) && (cap = min(cap, be[i] / a))
        end
        isfinite(cap) && (n_eq[j] = min(n_eq[j], max(cap, _EQ_GUESS_FLOOR)))
    end
    return n_eq
end

"""
    _restore_feasibility!(n_eq, Ae, be; maxit, tol) -> n_eq

Move a starting guess into `{Aₑn = bₑ, n ≥ 0}` by alternating projection: clamp
to the box, then project onto the affine set through the 7×7 system `AₑAₑᵀ`.

The Gibbs minimization is posed with `A n = b` as a hard equality, so a guess
that violates it starts the interior-point method outside its own feasible set.
On a full OPC this was not a detail: the solve stopped on a point demanding
0.732 mol of sulfate against the 0.267 mol available — 174 % over — and no
tolerance, iteration budget, barrier setting or bound changed it, because none
of them addresses an infeasible start. A projected-gradient phase-1 proved the
set is non-empty (residual 3·10⁻¹²); this is the cheap way to land in it.

Ending on the box clamp rather than the affine step matters: the affine
projection alone leaves small negative amounts, which the barrier cannot accept.

`maxit` deserves attention. Alternating projection converges linearly, at a rate
set by the angle between the box and the affine set, and on a cement that angle
can be small: at the six-hour instant of an ordinary Portland cement — where the
iron row carries 0.013 mol across thirteen species — 200 sweeps left a residual
of 8.4e-1, 2000 left 6.7e-2, and 20 000 were needed to reach 6.7e-9.

The default is small on purpose, and measured: inside the ODE right-hand side
this runs at every evaluation, and on a full OPC the worst in-run balance is
1.1 mol at 200 sweeps against 8.5 at 2000 and 41 at 100 000. A better guess
producing a worse answer is the back-end's own unpredictability; until that is
understood the in-run budget stays where it measures best. Note that the ranking
depends on the back-end and should be re-measured if it changes.

A replay ([`speciated_states`](@ref)) runs a handful of times and buys accuracy
instead, asking for a much larger budget.
"""
function _restore_feasibility!(n_eq, Ae, be; maxit::Int = 200, tol::Float64 = 1.0e-12)
    G = cholesky(Symmetric(Ae * Ae' + 1.0e-14I))
    sc = max(1.0, maximum(abs, be; init = 1.0))
    for _ in 1:maxit
        @. n_eq = max(n_eq, _EQ_GUESS_FLOOR)
        r = Ae * n_eq .- be
        maximum(abs, r; init = 0.0) <= tol * sc && break
        n_eq .-= Ae' * (G \ r)
    end
    @. n_eq = max(n_eq, _EQ_GUESS_FLOOR)
    return n_eq
end

# ── build_kinetics_ode ───────────────────────────────────────────────────────

"""
    build_kinetics_ode(kp::KineticsProblem) -> Function

Build the ODE right-hand-side `f!(du, u, p, t)` implementing Leal et al. (2017).

State layout:
  - `u[1:n_be]`                    = bₑ (element amounts in equilibrium partition)
  - `u[n_be+1 : n_be+n_nk]`        = nₖ (moles of kinetic species)
  - `u[n_be+n_nk+1 : n_be+n_nk+M]` = ξ  (extents of the M kinetic reactions)
  - `u[end]`                       = T with a `SemiAdiabaticCalorimeter`,
                                     Q with an `IsothermalCalorimeter`,
                                     absent when there is no calorimeter

ODE equations (Leal 2017, Eq. 66):
  - `dnₖ/dt = νₖᵀ r`
  - `dbₑ/dt = Aₑ νₑᵀ r`
  - `dξ/dt  = r`
  - `dT/dt  = (q̇ − φ(ΔT)) / Cp_total`  (semi-adiabatic)
  - `dQ/dt  = q̇`                       (isothermal)

where `nₑ = φ(bₑ)` is the equilibrium re-speciation constraint.
"""
function build_kinetics_ode(kp::KineticsProblem)
    function f!(du, u, p, t)
        # Promote with `typeof(t)`, not `eltype(u)` alone. Rosenbrock methods
        # (`Rodas5P`, `Rodas4`, …) need a TIME gradient, which they obtain by
        # calling `f!` with a dual `t` and a plain `u`. A rate law that actually
        # depends on `t` — `waller`, and any user law with an explicit time
        # dependence — then returns a `Dual` that cannot be stored in a
        # `Vector{Float64}`, and the solve fails with "First call to automatic
        # differentiation for time gradient failed". Rate laws that ignore `t`,
        # like `parrot_killoh`, never exposed this.
        T_elt = promote_type(eltype(u), typeof(t))

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

        # 2a'. NON-kinetic species from the extents: n = n(0) + νᵀ ξ.
        #
        # Without this the non-kinetic amounts stay pinned to their initial
        # values for the whole run, because the buffer holding them is refreshed
        # only by `respeciate!`. A rate law gating on such a species — sulfate
        # available for ettringite, portlandite available for a pozzolanic
        # reaction — would then never see it move, and would run past depletion
        # while the kinetic mass balance stayed exact and the solver reported
        # success. Skipped when an equilibrium solver is present: there the
        # equilibrium partition is owned by `respeciate!`, which redistributes it
        # in a way the stoichiometry alone cannot reproduce.
        if p.n_be == 0
            ξ = @view u[(p.n_be + p.n_nk + 1):(p.n_be + p.n_nk + p.n_rxn_state)]
            # `νe` holds exactly the equilibrium-partition columns of ν, in the
            # order of `idx_equilibrium`.
            for (k, idx) in enumerate(p.idx_equilibrium)
                acc = zero(T_elt)
                for j in eachindex(ξ)
                    acc += p.νe[j, k] * ξ[j]
                end
                n_full[idx] = max(p.n_initial_full[idx] + acc, p.ϵ)
            end
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

        # ── 6b. ODE: dξ/dt = r, the extents of reaction ────────────────
        for j in 1:(p.n_rxn_state)
            du[p.n_be + p.n_nk + j] = rates[j]
        end

        # ── 7. ODE: dbₑ/dt = Aₑ νₑᵀ r (Leal Eq. 65) ────────────────────
        if p.n_be > 0
            du_be = p.Ae * (p.νe' * rates)
            for j in 1:(p.n_be)
                du[j] = du_be[j]
            end
        end

        # ── 8a. ODE: dQ/dt = q̇ (isothermal) ────────────────────────────
        #
        # This is the heat of the KINETIC reactions. When those reactions produce
        # the hydrates it is the heat of hydration; under partial equilibrium they
        # only dissolve the anhydrous phases into ions and the precipitation heat
        # is invisible to it, which is why `cumulative_heat` says so and
        # `heat_release` exists.
        if p.has_Q
            du[end] = heat_rate(p.kin_rxns, rates, T_curr)
        end

        # ── 8b. ODE: dT/dt = (q̇ − φ(ΔT)) / Cp_total (semi-adiabatic) ───
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
