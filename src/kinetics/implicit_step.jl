using LinearAlgebra

# ── implicit_step.jl ──────────────────────────────────────────────────────────
#
# A kinetic step as Leal et al. (2017) pose it, and as Reaktoro implements it: ONE
# fully implicit problem per step, not an operator split.
#
#     min G(n)   s.t.   [A ; Kᵀ]·n + [0 ; −I]·Δξ = [b₀ ; ξ₀]      (linear)
#                       Δξ − Δt·M·r(n) = 0,   M = KᵀK             (nonlinear)
#                       n ≥ 0
#
# The reaction extents are unknowns of the SAME system as the amounts and the
# element potentials, and the rate is evaluated at the unknown end-of-step state —
# backward Euler. Two consequences follow, and both are the point:
#
#   * the equilibrium species and the kinetic ones are solved TOGETHER, so there
#     is no frozen speciation in a right-hand side and no lag between the two;
#   * the reactivity constraint `Kᵀn − Δξ = ξ₀` is LINEAR, so it joins the
#     conservation block beside the elements and the charge. The algebraic cost of
#     kinetics is the number of REACTIONS, not the number of species.
#
# `M = KᵀK` and not the identity: the constraint defines `Δξ = Kᵀ(n − n₀)`, and
# with `n − n₀ = Kξ` that is `KᵀK ξ`. Writing `Δξ = Δt·r` instead is wrong by a
# factor `M`.
#
# This coexists with `integrate`, which puts the same chemistry in an `ODEProblem`
# and hands it to SciML. That route has adaptive time-stepping and a choice of
# stiff solvers; this one has an exact end-of-step equilibrium and no splitting
# error. Neither dominates, and they are checked against each other.

"""
    KineticStepSolver(system, model, reactions; opts...)

A solver for one fully implicit kinetic step, following Leal et al. (2017).

`reactions` are [`KineticReaction`](@ref) objects. Every species not carried by a
declared reaction is at equilibrium; every declared reaction is kinetic and must
have a rate law, which `KineticReaction` already enforces.

Requires `OptimaSolver`.

# Example

```julia
kss = KineticStepSolver(cs, DiluteSolutionModel(), [kr_calcite])
st1 = kinetic_step(kss, st0, 60.0u"s")
```

See also [`kinetic_step`](@ref), [`integrate`](@ref) for the `ODEProblem` route.
"""
struct KineticStepSolver{D, R <: AbstractVector}
    dual::D
    reactions::R
    K::Matrix{Float64}          # species × reactions, restricted to the partition
    M::Matrix{Float64}          # KᵀK (`:reactions` coupling only)
    idx_kin::Vector{Int}        # the kinetic partition, in system order
    coupling::Symbol            # :reactions or :species
end

function KineticStepSolver(
        system::ChemicalSystem, model::AbstractActivityModel, reactions;
        kinetic_species = :auto,
        coupling::Symbol = :reactions,
        kwargs...,
    )
    if coupling === :species
        throw(
            ArgumentError(
                "`coupling = :species` is not available. It was implemented and " *
                    "does not converge, so it is refused rather than shipped: " *
                    "pinning a pure phase by a linear row leaves its stationarity " *
                    "row in the system as well, made trivially satisfiable by the " *
                    "row's own multiplier, and the Jacobian degenerates — measured " *
                    "on two C3A pathways, the extents came out 4.3 times short " *
                    "with the solve reporting non-convergence. Removing those " *
                    "stationarity rows is what remains to be done.\n\n" *
                    "For a model whose hydrate assemblage IS the stoichiometry — " *
                    "the form of Lavergne et al. (2018) — use the ODE route: " *
                    "`KineticsProblem(system, reactions, state, tspan)` with " *
                    "`equilibrium_solver = nothing`. It integrates " *
                    "`dn/dt = ν' r(n)` directly, so every product is imposed, and " *
                    "it already supports several reactions per mineral.",
            )
        )
    end
    coupling === :reactions || throw(
        ArgumentError("`coupling` must be `:reactions` — one constraint per " *
                      "reaction, which is Leal's form.")
    )
    isempty(reactions) && throw(
        ArgumentError("no kinetic reaction given: with none, the step is a plain " *
                      "equilibrium — call `equilibrate` instead.")
    )
    ns = length(system.species)
    for kr in reactions
        length(kr.stoich) == ns || throw(
            DimensionMismatch("a reaction's `stoich` has $(length(kr.stoich)) " *
                              "entries for a system of $ns species.")
        )
    end
    K = _reactivity_matrix(reactions, system, kinetic_species)
    idx_kin = [i for i in 1:ns if any(!iszero, @view K[i, :])]
    if coupling === :species && length(idx_kin) > 1 && rank(K) < length(reactions)
        # already caught in `_reactivity_matrix`, kept for symmetry
    end
    des = DualEquilibriumSolver(system, model; kwargs...)
    return KineticStepSolver(
        des, collect(reactions), K, transpose(K) * K, idx_kin, coupling,
    )
end

"""
    kinetic_step(kss, state, Δt; t = 0.0, ϵ = 1e-16, parameters = nothing)
        -> ChemicalState

One implicit kinetic step of duration `Δt` from `state`.

`t` is the absolute time handed to rate laws that depend on it. Pass a `Ref` as
`parameters` to receive the reaction extents `Δξ` the step found.

`warm_start` equilibrates the starting GUESS when the system carries solid
solutions, leaving the component totals untouched. Keep it on: a mixing phase is
admitted by a tangent-plane test, and from a composition where the phase is
absent that admission fails — measured, a cold start on C₃S dissolving into a
C-S-H solution left one end-member at 2.7e-9 with a stationarity residual of 6.5.

The rate laws are evaluated at the **end-of-step** composition, which is what
makes the step implicit and stable on a stiff system. `Δt` is the caller's choice;
[`kinetic_step_adaptive`](@ref) chooses it from a local error estimate instead.
"""
function kinetic_step(
        kss::KineticStepSolver, state::ChemicalState, Δt;
        t = 0.0, ϵ::Float64 = 1.0e-16,
        warm_start::Bool = true,
        pin_minerals = :auto,
        parameters::Union{Nothing, Base.RefValue} = nothing,
        certificate::Union{Nothing, Base.RefValue} = nothing,
    )
    # Whether to hold the kinetic minerals in the active set is decided by the
    # certificate, not by taste, because neither answer works on both cases.
    #
    # Measured. Two C3A pathways: without pinning, certified at a stationarity of
    # 1.8e-12 and `Δξ` exactly `Δt·M·r`; with pinning, stationarity 7.9 and not
    # certified. C3S dissolving into a C-S-H solid solution: without pinning the
    # step runs all the way to equilibrium, `C3S = 0` and `Δξ` twenty-five times
    # too large, because the mineral is exhausted at equilibrium and the drop rule
    # removes it, after which nothing enforces its reactivity row; with pinning,
    # certified at 2.3e-13.
    #
    # The certificate DECIDES — the problem is convex and its KKT conditions are
    # sufficient — so trying both and keeping a proved answer is exact, not a
    # heuristic. Pass `true` or `false` to force one.
    if pin_minerals === :auto
        for pin in (false, true)
            cert_try = Ref{Any}(nothing)
            out = kinetic_step(
                kss, state, Δt;
                t = t, ϵ = ϵ, warm_start = warm_start, pin_minerals = pin,
                parameters = parameters, certificate = cert_try,
            )
            if cert_try[] !== nothing && cert_try[].optimal
                certificate === nothing || (certificate[] = cert_try[])
                return out
            end
            if pin        # the last attempt: hand back what there is, and say so
                certificate === nothing || (certificate[] = cert_try[])
                @warn """no implicit kinetic step was certified, with the minerals \
                pinned or free; the last answer is returned. Try a smaller Δt, and \
                audit with the certificate.""" maxlog = 1
                return out
            end
        end
    end
    Δt_s = Δt isa Number && !(Δt isa DynamicQuantities.AbstractQuantity) ?
        Float64(Δt) : ustrip(us"s", Δt)
    des = kss.dual
    p = _build_params(state; ϵ = ϵ)
    n0 = Float64[ustrip(us"mol", x) for x in state.n]

    K, M = kss.K, kss.M
    nr = size(K, 2)
    m = size(des.A, 1)
    ns = length(n0)

    A_aug, Aq, b_aug, use_M = if kss.coupling === :reactions
        # Leal's form: `[A ; Kᵀ] n + [0 ; −I] Δξ = [A n₀ ; Kᵀ n₀]`, one row per
        # reaction. The extents are functionals of the composition, and every
        # species the partition excludes is left to the minimization.
        (vcat(des.A, transpose(K)),
         vcat(zeros(m, nr), -Matrix(1.0I, nr, nr)),
         vcat(des.A * n0, transpose(K) * n0),
         true)
    else
        # One row per KINETIC SPECIES: `nᵢ − Σⱼ νᵢⱼ Δξⱼ = nᵢ(0)`. Every species in
        # the partition is then pinned by the extents, so a solid-to-solid
        # reaction imposes its products instead of proposing them — while the
        # aqueous phase still minimizes `G` under element conservation, which is
        # the coupling a plain ODE integration does not have.
        nk = length(kss.idx_kin)
        E = zeros(nk, ns)
        for (r, i) in enumerate(kss.idx_kin)
            E[r, i] = 1.0
        end
        (vcat(des.A, E),
         vcat(zeros(m, nr), -K[kss.idx_kin, :]),
         vcat(des.A * n0, n0[kss.idx_kin]),
         false)
    end

    T_K = p.T
    P_Pa = p.P
    idx = des.system.dict_species
    rates = _step_rate_closure(kss, des, T_K, P_Pa, t, n0)

    # `Δξ − Δt·M·r(n) = 0` for Leal's form, `Δξ − Δt·r(n) = 0` for the
    # species-pinned one. The `M` is not decoration: `Kᵀn − Δξ = ξ₀` defines
    # `Δξ = Kᵀ(n − n₀)`, and with `n − n₀ = Kξ` that is `KᵀK ξ`. Pinning the
    # species instead makes `Δξ` the reaction progress itself, so `M` disappears.
    cq = use_M ?
        (x, q, params) -> q .- Δt_s .* (M * rates(x, params)) :
        (x, q, params) -> q .- Δt_s .* rates(x, params)

    # The difference step is scaled to the extent the explicit rate predicts,
    # which is the only scale available before the step is taken. Below a floor
    # so a reaction starting at zero rate still has a usable step.
    r0 = rates(n0, p)
    qscale = [
        max(abs(Δt_s * (use_M ? sum(M[j, :] .* r0) : r0[j])), 1.0e-12) for j in 1:nr
    ]

    # A kinetic mineral is determined by its reactivity row, so it must not be
    # subject to an active-set sign test.
    #
    # This is not defensive: on a mineral that dissolves ENTIRELY at equilibrium
    # the drop rule removes it — its saturation index is far positive, it is on
    # its way out — and the row pinning it then cannot be satisfied. Measured on
    # C₃S dissolving into a C-S-H solution, the step ran to full equilibrium with
    # `C3S = 0` exactly and `Δξ = 0.51` against `Δt·r = 0.02`. On calcite the same
    # bug is invisible, because only a thousandth of it dissolves and the active
    # set keeps it.
    pinned = pin_minerals ? [i for i in kss.idx_kin if i in des.idx_pure] : Int[]

    prob = _optima_dual_problem(
        A_aug, Float64.(p.ΔₐG⁰overT), des.lna, _dual_phases(des, n0), des.idx_pure, p,
        nothing, nothing, cq, zeros(nr), qscale, Aq, pinned,
    )
    # The GUESS is equilibrated first; `b_aug` stays the true initial budget.
    #
    # This is not a speed-up, it is what makes a solid solution work. The dual
    # solver admits a mixing phase by a tangent-plane test inside its active-set
    # loop, and from a composition where the phase is entirely absent that
    # admission fails: measured on C₃S dissolving into a C-S-H solid solution,
    # a cold start left one end-member at 2.7e-9 with a stationarity residual of
    # 6.5, while the same problem started from an equilibrated guess certified at
    # 3.6e-13 with all four end-members present. The failure belongs to the cold
    # start, not to the kinetics — a plain equilibrium fails the same way.
    x0 = if warm_start && !isempty(des.ss_groups) && _DUAL_AVAILABLE[]
        try
            eq = first(equilibrate_certified(state; model = des.model, ϵ = ϵ))
            Float64[ustrip(us"mol", x) for x in eq.n]
        catch
            n0
        end
    else
        n0
    end

    res = _optima_dual_solve(prob, b_aug, x0, des.opts)

    res.converged || begin
        NONCONVERGED[] += 1
        @warn "the implicit kinetic step did not converge; try a smaller Δt" maxlog = 1
    end
    parameters === nothing || (parameters[] = copy(res.q))

    # The certificate must be taken on the AUGMENTED problem, not on the
    # unconstrained equilibrium. A kinetically held mineral is supersaturated by
    # construction — that is what "held back" means — so testing it against the
    # unconstrained KKT conditions reports a stationarity residual of several RT
    # and calls a correct answer wrong. The reactivity rows carry their own
    # multipliers, and those are exactly what makes the mineral's stationarity
    # satisfiable at the amount the constraint demands.
    if certificate !== nothing
        certificate[] = _optima_kkt_certificate(
            prob, res.x, b_aug, 1.0e-25, des.opts.tol, des.opts.si_tol, res.q,
        )
    end

    return ChemicalState(
        des.system, [nᵢ * u"mol" for nᵢ in res.x];
        T = temperature(state), P = pressure(state),
    )
end

"""
    _reactivity_matrix(reactions, system, kinetic_species) -> Matrix

The matrix `K` of the reactivity constraint `Kᵀn − Δξ = ξ₀`: one column per
reaction, carrying that reaction's stoichiometry restricted to the **kinetic
partition** — the species whose amounts the reactions alone decide.

The restriction is the modeling decision of the whole scheme. A rate law says how
fast a mineral reacts; it says nothing about where matter released into solution
ends up, which is the equilibrium problem's business. So an aqueous species must
be left OUT of `K`: it re-speciates, and including it makes `Kᵀn` a functional of
the products, which stops pinning the mineral. Measured on calcite dissolving at
1 µmol/s for 1000 s with the full stoichiometry in `K`, the amount came out
5.7e-5 mol below `n₀ − kΔt`, the constraint having absorbed the carbonate's
speciation into the extent. The answer solved its own equations exactly and meant
something else.

`kinetic_species = :auto` takes the partition to be **every non-aqueous
participant of a declared reaction**. That covers both styles of model:

  - **dissolution into solution** — `Cal → Ca²⁺ + CO₃²⁻`. Only `Cal` is kinetic,
    so `K = [−1]`, `M = KᵀK = 1`, `Δξ = Δt·r`, and the released matter speciates
    freely under element conservation.
  - **solid to solid** — `C₃S + H₂O → C-S-H + CH`, the form
    [Lavergne2018](@cite) uses. Every solid product is kinetic too, so the
    assemblage is **imposed by the stoichiometry** rather than found by
    minimizing `G`. That is a different model, deliberately, and it is the one to
    write when the hydrate assemblage is part of what you are prescribing.

Pass an explicit list of species (symbols or indices) to override the partition.

Several reactions may share a mineral — the whole point of a reaction-centric
scheme, and what [Lavergne2018](@cite) needs for C₃A reacting through both an
ettringite and a monosulfoaluminate pathway. Their columns must be **linearly
independent**, though: two reactions differing only in products the partition
excludes would give the same row twice, `Kᵀn` could not tell their extents apart,
and the linear system would be singular. That is checked here, with the offending
reactions named.
"""
function _reactivity_matrix(reactions, system::ChemicalSystem, kinetic_species)
    ns = length(system.species)
    syms = symbol.(system.species)
    aqueous = Set(system.idx_aqueous)

    kin = if kinetic_species === :auto
        s = Set{Int}()
        for kr in reactions, i in 1:ns
            kr.stoich[i] == 0 && continue
            i in aqueous || push!(s, i)
        end
        s
    else
        Set(
            i isa Integer ? Int(i) : begin
                j = findfirst(==(String(i)), syms)
                j === nothing && throw(
                    ArgumentError("`$i` is not a species of this system.")
                )
                j
            end for i in kinetic_species
        )
    end

    isempty(kin) && throw(
        ArgumentError("no kinetic species: every participant of every declared " *
                      "reaction is aqueous, so there is nothing whose amount the " *
                      "rate laws decide. Name the controlled species with " *
                      "`kinetic_species = [...]`.")
    )

    K = zeros(Float64, ns, length(reactions))
    for (j, kr) in enumerate(reactions), i in kin
        K[i, j] = kr.stoich[i]
    end

    for j in axes(K, 2)
        all(iszero, @view K[:, j]) && throw(
            ArgumentError("reaction $j touches no kinetic species, so its extent " *
                          "is not determined by the composition. Its solid " *
                          "participants are missing from the partition.")
        )
    end

    nr = size(K, 2)
    if nr > 1 && rank(K) < nr
        # Name the pair that collides, so the message is actionable.
        culprits = String[]
        for a in 1:nr, bb in (a + 1):nr
            rank(K[:, [a, bb]]) < 2 && push!(culprits, "$a and $bb")
        end
        throw(
            ArgumentError(
                "the reactivity constraints are linearly dependent" *
                    (isempty(culprits) ? "" : " (reactions " *
                     join(culprits, ", ") * ")") *
                    ": `Kᵀn` cannot tell those extents apart and the linear " *
                    "system would be singular. Two reactions may share a mineral " *
                    "— that is the point of a reaction-centric scheme — but they " *
                    "must differ in some species of the kinetic partition. If " *
                    "they differ only in aqueous products, add those products to " *
                    "`kinetic_species`, or model the pathway as one reaction " *
                    "with a combined rate.",
            )
        )
    end
    return K
end

"""
    _step_rate_closure(kss, des, T, P, t, n0) -> (x, params) -> Vector

The reaction rates at composition `x`, in mol/s, one per declared reaction.

Built once per step so the `StateView`s and the initial amounts are not rebuilt
at every residual evaluation — the residual is called once per Newton iteration
and once per column of the difference quotient, so `nr + m + …` times per step.
"""
function _step_rate_closure(kss::KineticStepSolver, des, T, P, t, n0)
    index = Dict(symbol(s) => i for (i, s) in enumerate(des.system.species))
    n0_sv = StateView(n0, index)
    return function (x, params)
        xv = collect(x)
        n_sv = StateView(xv, index)
        lna_sv = StateView(collect(des.lna(xv, params)), index)
        return Float64[
            kr.rate_fn(T, P, t, n_sv, lna_sv, n0_sv) for kr in kss.reactions
        ]
    end
end

# ── adaptive stepping ────────────────────────────────────────────────────────
#
# What Reaktoro does not have. Its kinetics adds a single field to the
# equilibrium options — an initial step — and the step itself is the caller's,
# with no estimate of what taking it cost. The scheme is backward Euler, so the
# error is first order in `Δt`, and a step ten times too large returns an answer
# ten times less accurate with nothing to say so.
#
# The estimate is Richardson's, on the extents: one step of `Δt` against two of
# `Δt/2`. For a first-order method the difference IS the error of the coarse step
# to leading order, and `2·(ξ_fine − ξ_coarse)` estimates the error of the fine
# one. Three implicit solves per accepted step, which is the price.

"""
    kinetic_step_adaptive(kss, state, Δt; reltol, abstol, ...) -> (state, Δt_used, Δt_next)

One kinetic step with its length chosen from a local error estimate.

Takes one step of `Δt` and two of `Δt/2`, compares the extents, and accepts the
finer pair when the estimated error is within tolerance — halving and retrying
otherwise. Returns the accepted state, the step actually taken, and a suggestion
for the next one.

The estimate is Richardson's for a first-order method: `ξ_fine − ξ_coarse` is the
error of the coarse step to leading order. The error controlled is therefore on
the **reaction extents**, which is what a kinetic step advances; the composition
follows from them through an exact equilibrium solve, so it carries no separate
error of its own.

`Δt` is the step to attempt, `max_halvings` bounds the retries, and the tolerance
is `abstol + reltol·sⱼ` per reaction, where `sⱼ` is the **amount the extent acts
on** — `Σᵢ |Kᵢⱼ| nᵢ`, the kinetic minerals of reaction `j` at the start of the
step.

Scaling by the extent instead, `reltol·|Δξⱼ|`, does not work, and the reason is
worth stating because it is the obvious thing to write. `Δξ ∝ Δt`, so that
tolerance vanishes with the step while the equilibrium solve's own noise does
not: the measured error then behaves as `noise / (reltol·Δt)` and **grows** as the
step shrinks. Measured, the controller halved to the floor without advancing and
the final answer got worse as the tolerance was tightened — 2.5e-5 at `1e-3`
against 9.2e-5 at `1e-5`. An ODE integrator's `reltol` multiplies the solution,
not the increment, for exactly this reason.

# Example

```julia
st, dt_used, dt_next = kinetic_step_adaptive(kss, st0, 3600.0u"s"; reltol = 1e-4)
```

See also [`kinetic_step`](@ref) for a step of a length you choose, and
[`integrate`](@ref) for the `ODEProblem` route, which brings SciML's own
step-size control and its choice of stiff solvers.
"""
function kinetic_step_adaptive(
        kss::KineticStepSolver, state::ChemicalState, Δt;
        t = 0.0, reltol::Float64 = 1.0e-4, abstol::Float64 = 1.0e-14,
        max_halvings::Int = 12, ϵ::Float64 = 1.0e-16,
        warm_start::Bool = true, pin_minerals = :auto,
        parameters::Union{Nothing, Base.RefValue} = nothing,
        error_estimate::Union{Nothing, Base.RefValue} = nothing,
    )
    Δt_s = Δt isa DynamicQuantities.AbstractQuantity ? ustrip(us"s", Δt) : Float64(Δt)
    Δt_s > 0 || throw(ArgumentError("`Δt` must be positive, got $Δt_s s."))

    common = (;
        ϵ = ϵ, warm_start = warm_start, pin_minerals = pin_minerals,
    )

    # The scale the tolerance is relative to: the amount each reaction acts on.
    n_now = Float64[ustrip(us"mol", x) for x in state.n]
    scales = [
        max(sum(abs(kss.K[i, j]) * n_now[i] for i in axes(kss.K, 1)), 1.0e-30)
            for j in axes(kss.K, 2)
    ]

    h = Δt_s
    for _ in 0:max_halvings
        q_c = Ref(Float64[])
        kinetic_step(kss, state, h; t = t, parameters = q_c, common...)

        q_1 = Ref(Float64[])
        mid = kinetic_step(kss, state, h / 2; t = t, parameters = q_1, common...)
        q_2 = Ref(Float64[])
        fine = kinetic_step(
            kss, mid, h / 2; t = t + h / 2, parameters = q_2, common...,
        )
        ξ_fine = q_1[] .+ q_2[]

        err = maximum(
            abs(ξ_fine[j] - q_c[][j]) / (abstol + reltol * scales[j])
                for j in eachindex(ξ_fine)
        )
        error_estimate === nothing || (error_estimate[] = err)

        if err <= 1
            parameters === nothing || (parameters[] = ξ_fine)
            # First order, so the step scales as the tolerance ratio itself,
            # capped at doubling — the usual safety factor of 0.9.
            Δt_next = h * min(2.0, max(0.2, 0.9 / max(err, 1.0e-10)))
            return (fine, h, Δt_next)
        end
        h /= 2
    end

    throw(
        ErrorException(
            "the local error on the reaction extents stayed above tolerance after " *
                "$max_halvings halvings, down to Δt = $(Δt_s / 2^max_halvings) s. " *
                "Either a rate law is discontinuous at this composition — a hard " *
                "gate on a species that has just run out will do it — or the " *
                "tolerance is below what the equilibrium solve itself delivers " *
                "(its own tolerance is $(kss.dual.opts.tol)).",
        )
    )
end
