# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# Tests for the calibration example, `scripts/hydration_calibration.jl`: the
# reading of measured calorimetry, the parameter map, the loss, and the one thing
# that actually proves the machinery works — recovering parameters from data
# generated with known ones.
#
# The coupled forward model is NOT exercised here. One solve costs tens of
# seconds and a fit costs an hour; that belongs in the script's `main()`, not in a
# test suite.

using LinearAlgebra
using Statistics

const CALIB_SCRIPT = joinpath(pkgdir(ChemistryLab), "scripts", "hydration_calibration.jl")
isdefined(Main, :read_calorimetry) || include(CALIB_SCRIPT)

@testset "read_calorimetry reads the record and its conditions" begin
    d = CEM_I_TARGET

    # The four quantities the experiment reports and a fit must not touch.
    @test d.meta.blaine ≈ 397.0
    @test d.meta.wb ≈ 0.5
    @test d.meta.T ≈ 293.15
    @test d.meta.q_before_start ≈ 12.0
    @test occursin("zenodo", d.meta.doi)
    @test occursin("CC-BY-4.0", d.meta.license)

    # The record itself.
    @test length(d) > 400
    @test issorted(d.t)
    @test d.t[1] / 3600 ≈ 1.023 atol = 1.0e-3
    @test all(>(0), d.qdot)
    @test issorted(d.Q)                       # released heat cannot decrease
    @test d.Q[end] ≈ 376.0 atol = 0.1

    # The fourth column is somebody else's fit, and must not be mistaken for data.
    @test length(d.Qref) == length(d.Q)
    @test d.Qref[end] > d.Q[end]              # their affinity model overshoots here

    # Q must be the integral of q̇. The record is a log-spaced SUBSET, so the
    # trapezoidal rule on it is only good to a few percent — enough to catch a
    # unit error or a swapped column, which is the point.
    Q_int = similar(d.Q)
    Q_int[1] = d.Q[1]
    for i in 2:lastindex(d.t)
        Q_int[i] = Q_int[i - 1] + 0.5 * (d.qdot[i] + d.qdot[i - 1]) * (d.t[i] - d.t[i - 1])
    end
    @test Q_int[end] ≈ d.Q[end] rtol = 0.05
end

@testset "resample_log subsets without inventing anything" begin
    d = CEM_I_TARGET
    r = resample_log(d, 30)

    @test length(r) <= 30
    @test r.t[1] == d.t[1]
    @test r.t[end] == d.t[end]
    @test issorted(r.t)
    @test length(r.Qref) == length(r.Q)
    # every retained value must be one of the originals: no interpolation
    @test all(in(d.t), r.t)
    @test all(in(d.Q), r.Q)

    # `tmax` bounds the range, not the number of points: thirty log-spaced targets
    # over a shorter interval still land on thirty distinct rows.
    r24 = resample_log(d, 30; tmax = 24 * 3600)
    @test r24.t[end] <= 24 * 3600
    @test r24.t[end] > 20 * 3600
    @test r24.t[1] == d.t[1]
end

@testset "the parameter map is a bijection onto the box" begin
    θ0 = prior_vector()
    @test length(θ0) == length(CALIB_SPEC)
    @test to_native(to_search(θ0)) ≈ θ0

    # The box is enforced by the transform, so no iterate can leave it however
    # far the optimizer wanders.
    for z in (-50.0, -5.0, 0.0, 5.0, 50.0)
        θ = to_native(fill(z, length(CALIB_SPEC)))
        for (p, v) in zip(CALIB_SPEC, θ)
            @test p.lo <= v <= p.hi
        end
        @test all(0 .<= box_position(θ) .<= 1)
    end

    # z = 0 is the middle of the box, in log for a rate and linearly otherwise.
    mid = to_native(zeros(length(CALIB_SPEC)))
    for (p, v) in zip(CALIB_SPEC, mid)
        expected = p.logscale ? sqrt(p.lo * p.hi) : (p.lo + p.hi) / 2
        @test v ≈ expected
    end

    @test box_position(θ0) ≈ box_position(θ0, CALIB_SPEC)
end

@testset "apply_parameters writes only what it is told to" begin
    θ0 = prior_vector()
    pk = apply_parameters(θ0)

    # At the prior the sets must be the published ones, to the unit round trip.
    for ph in ("C3S", "C2S", "C3A", "C4AF")
        for k in keys(PK84_PRIOR[ph])
            @test ustrip(pk[ph][k]) ≈ ustrip(PK84_PRIOR[ph][k]) rtol = 1.0e-12
        end
    end

    # Move one parameter; nothing else may move.
    j = findfirst(p -> p.name === :k₁_C3S, CALIB_SPEC)
    θ = copy(θ0)
    θ[j] *= 2
    pk2 = apply_parameters(θ)
    @test ustrip(us"1/s", pk2["C3S"].k₁) ≈ 2 * ustrip(us"1/s", pk["C3S"].k₁)
    @test pk2["C3S"].n₁ == pk["C3S"].n₁
    @test pk2["C2S"] == pk["C2S"]
    @test pk2["C3A"] == pk["C3A"]

    # Exponents must stay bare numbers and rates must keep a reciprocal-time unit.
    @test pk2["C3S"].n₃ isa Float64
    @test ustrip(us"1/s", pk2["C3S"].k₃) isa Float64
    @test ustrip(us"1/d", pk2["C3S"].k₃) ≈ 86400 * ustrip(us"1/s", pk2["C3S"].k₃)
end

@testset "the surrogate is self-consistent" begin
    θ0 = prior_vector()
    d = resample_log(CEM_I_TARGET, 25; tmax = 3 * 86400.0)
    run = stoichiometric_run(
        apply_parameters(θ0); wb = d.meta.wb, blaine = d.meta.blaine, tend = d.t[end],
    )

    # No amount may go negative. The first version of this surrogate sent the
    # aluminate to ettringite, which drove the gypsum below zero and the released
    # heat to 1609 J/g against a measured 376 — this assertion is what catches it.
    @test all(u -> all(>=(-1.0e-12), u), run.sol.u)
    @test SciMLBase.successful_retcode(run.sol)

    Q = forward_Q(θ0, d; mode = :surrogate)
    @test all(isfinite, Q)
    @test all(>=(-1.0e-9), diff(Q))                    # non-decreasing
    @test Q[1] ≈ 0 atol = 1.0e-8                       # referred to the first instant
    # Right order of magnitude, and low, as an imposed assemblage must be.
    @test 0.6 * d.Q[end] < Q[end] < 1.1 * d.Q[end]
end

@testset "the loss behaves like a loss" begin
    θ0 = prior_vector()
    d = resample_log(CEM_I_TARGET, 25)

    L = calorimetry_loss(θ0, d; mode = :surrogate)
    @test isfinite(L) && L > 0
    @test L < LOSS_PENALTY

    # In J/g, so it is comparable with the instrument, not a dimensionless score.
    @test L ≈ sqrt(mean(abs2, forward_Q(θ0, d; mode = :surrogate) .- d.Q))

    # The shape loss ignores the level, so a model rescaled by a constant scores
    # the same on it and worse on the absolute one.
    Ls = calorimetry_loss(θ0, d; mode = :surrogate, shape = true)
    @test isfinite(Ls) && Ls > 0
    @test Ls != L

    @test_throws ArgumentError forward_Q(θ0, d; mode = :nonsense)
end

@testset "the Jander branch does bind on alite, late" begin
    d = resample_log(CEM_I_TARGET, 40)

    # Parrott & Killoh reported no diffusion-controlled stage for C₃S, and
    # `parrot_killoh_avrami`'s docstring repeats it. Numerically that is not quite
    # what this implementation does: `α̇₂ = k₂(1-ξ)^(2/3)/(1-(1-ξ)^(1/3))` falls as
    # ξ grows, so at high degrees of hydration the Jander branch can become the
    # minimum and k₂ acquires a real, if modest, influence. Pinning the numbers
    # here so the claim in the documentation stays honest.
    spec_k2 = [CalibParameter(:k₂_C3S, "C3S", :k₂, 0.05, 0.005, 0.5, true)]
    J2 = sensitivity_matrix([0.05], d; mode = :surrogate, spec = spec_k2)
    spec_k1 = [CalibParameter(:k₁_C3S, "C3S", :k₁, 1.5, 0.375, 6.0, true)]
    J1 = sensitivity_matrix([1.5], d; mode = :surrogate, spec = spec_k1)

    @test maximum(abs, J1) > 1.0            # k₁ certainly limits
    @test maximum(abs, J2) > 0.0            # so does k₂, contrary to the docstring
    @test maximum(abs, J2) < maximum(abs, J1) / 4   # but it is the minor effect

    # And it is confined to the late part of the record: nothing at early ages.
    early = findall(<=(6 * 3600), d.t)
    @test all(iszero, J2[early, 1])
end

@testset "parameter recovery on synthetic data" begin
    # The test that would catch a sign or scaling error in the loss: generate a
    # curve from KNOWN parameters, fit it back, and require the fit to find them.
    # Two parameters, both individually identifiable — recovering six would fail
    # for a reason that is physics, not code, and §5 of the documentation page
    # measures that.
    spec = [
        CalibParameter(:k₁_C3S, "C3S", :k₁, 1.5, 1.5 / 4, 1.5 * 4, true),
        CalibParameter(:k₃_C3S, "C3S", :k₃, 1.1, 1.1 / 4, 1.1 * 4, true),
    ]
    θ_true = [1.5 * 1.35, 1.1 * 0.7]

    base = resample_log(CEM_I_TARGET, 40)
    Q_true = forward_Q(θ_true, base; mode = :surrogate, spec)
    @test all(isfinite, Q_true)

    # Deterministic perturbation at the amplitude the depositors quote for their
    # own fit (4.11 J/g), so the test cannot flake.
    noise = 4.11 .* sinpi.(2 .* (0:(length(Q_true) - 1)) ./ 7)
    synthetic = CalorimetryData(base.t, base.qdot, Q_true .+ noise, Float64[], base.meta)

    r = calibrate(
        synthetic; mode = :surrogate, θ0 = prior_vector(spec), spec, maxiters = 400,
    )

    @test r.loss < r.loss0
    @test all(0.02 .< r.box .< 0.98)              # interior, so the data chose it
    @test r.θ[1] ≈ θ_true[1] rtol = 0.15
    @test r.θ[2] ≈ θ_true[2] rtol = 0.25
end

@testset "cumulative_heat now returns a heat, and agrees with heat_release" begin
    # Regression test for a real defect: `cumulative_heat(sol, ::Isothermal…)`
    # used to return `u[end]`, which under the old state layout was the last
    # REACTION EXTENT in mol, not a heat — the isothermal calorimeter added no
    # state to the ODE and its `saved_H` fast path was unreachable. Both routes
    # to the heat are compared here, on the stoichiometric formulation where both
    # are valid.
    θ0 = prior_vector()
    pk = apply_parameters(θ0)
    tend = 3.0 * 86400.0

    cs = stoichiometric_system()
    α_max = powers_alpha_max(0.5)
    state0 = ChemicalState(cs; T = 293.15u"K")
    for (nm, w) in pairs(CALIB_CLINKER)
        set_quantity!(state0, string(nm), (0.919 * w)u"kg")
    end
    set_quantity!(state0, "Gp", CALIB_GYPSUM * u"kg")
    set_quantity!(state0, "Cal", CALIB_FILLER * u"kg")
    set_quantity!(state0, "H2O@", 0.5u"kg")

    rxn = Reaction(
        OrderedDict(cs["C3S"] => 1.0, cs["H2O@"] => 103 / 30),
        OrderedDict(cs["Jennite"] => 1.0, cs["Portlandite"] => 4 / 3);
        symbol = "C₃S hydration",
    )
    rxn[:rate] = parrot_killoh_avrami(pk["C3S"], "C3S"; α_max, blaine = 397u"m^2/kg")

    cal = IsothermalCalorimeter(293.15u"K")
    kp = KineticsProblem(
        cs, [rxn], state0, (0.0, tend); calorimeter = cal, equilibrium_solver = nothing,
    )
    sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P()))

    # The calorimeter must have added exactly one state, and it must be the heat.
    @test length(sol.u[1]) == length(
        build_u0(
            KineticsProblem(
                cs, [rxn], state0, (0.0, tend); equilibrium_solver = nothing,
            )
        )
    ) + 1

    t, Q = cumulative_heat(sol, cal)
    @test Q[1] ≈ 0 atol = 1.0e-10
    @test issorted(Q)                       # exothermic throughout
    @test Q[end] > 1.0e5                    # joules, not moles of extent

    times = 10 .^ range(log10(3600.0), log10(tend); length = 20)
    _, Q_ref, _ = heat_release(sol, kp; times)
    Qi = [Q[argmin(abs.(t .- τ))] for τ in times]
    # Same physics by two routes: ∫q̇ dt against H(t₀) − H(t). They are referred
    # to different origins, so compare the increments.
    @test (Qi[end] - Qi[1]) ≈ (Q_ref[end] - Q_ref[1]) rtol = 0.02

    tq, qdot = heat_flow(sol, cal)
    @test length(tq) == length(t)
    @test maximum(qdot) > 0
end
