using OrderedCollections

# ── A small but genuine clinker hydration run, shared by the testsets below ───
#
# Two Parrot-Killoh phases and two balanced CEMDATA18 reactions is enough to
# exercise every code path: several reactions consuming a shared product
# (Portlandite), a species that is neither kinetic nor a reactant, and a
# non-trivial stoichiometric matrix.

function _pp_problem(; tend = 7 * 86400.0)
    substances = build_species(joinpath(pkgdir(ChemistryLab), "data", "cemdata18-thermofun.json"))
    sp = speciation(
        substances, split("C3S C2S Portlandite Jennite ettringite H2O@");
        aggregate_state = [AS_AQUEOUS]
    )
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)

    state0 = ChemicalState(cs)
    set_quantity!(state0, "C3S", 0.65u"mol")
    set_quantity!(state0, "C2S", 0.11u"mol")
    set_quantity!(state0, "H2O@", 8.0u"mol")

    s(name) = cs[name]

    # CEMDATA18 stores Jennite as (SiO2)1(CaO)1.666667(H2O)2.1 — a ROUNDED 5/3,
    # not the exact fraction. Writing the classic 4/3 and 103/30 coefficients
    # (as the shipped `scripts/*_kinetics.jl` do) leaves the reactions unbalanced
    # by ~3e-6 mol per mol, which would show up here as spurious element
    # non-conservation. Derive the coefficients from the database instead.
    _ca_jen = 1.666667
    _h2o_jen = 2.1
    ν_ch3 = 3 - _ca_jen               # Portlandite from C3S
    ν_ch2 = 2 - _ca_jen               # Portlandite from C2S

    r1 = Reaction(
        OrderedDict(s("C3S") => 1.0, s("H2O@") => _h2o_jen + ν_ch3),
        OrderedDict(s("Jennite") => 1.0, s("Portlandite") => ν_ch3);
        symbol = "C3S hydration",
    )
    r1[:rate] = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S")
    r2 = Reaction(
        OrderedDict(s("C2S") => 1.0, s("H2O@") => _h2o_jen + ν_ch2),
        OrderedDict(s("Jennite") => 1.0, s("Portlandite") => ν_ch2);
        symbol = "C2S hydration",
    )
    r2[:rate] = parrot_killoh_avrami(PK84_PARAMS_C2S, "C2S")

    kp = KineticsProblem(cs, [r1, r2], state0, (0.0, tend); equilibrium_solver = nothing)
    ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-12)
    return cs, state0, kp, integrate(kp, ks)
end

@testset "reaction_extents" begin

    _, _, kp, sol = _pp_problem()

    ξ = reaction_extents(sol, kp)
    @test size(ξ) == (length(sol.t), 2)

    # ── Zero at t = 0, monotone, bounded by the initial amounts ───────────────
    @test all(ξ[1, :] .== 0.0)
    @test all(diff(ξ[:, 1]) .>= -1.0e-14)
    @test all(diff(ξ[:, 2]) .>= -1.0e-14)
    n0 = sol.prob.p.n_initial_full
    @test ξ[end, 1] <= n0[kp.idx_kinetic[1]] + 1.0e-12
    @test ξ[end, 2] <= n0[kp.idx_kinetic[2]] + 1.0e-12

    # ── Mass balance: nₖ(t) - nₖ(0) = νₖᵀ ξ(t) ────────────────────────────────
    # This is the quadrature error of `reaction_extents`, and the only place it
    # can show up: `state_at` is mass-balanced by construction.
    res = extent_residual(sol, kp, ξ)
    @test res < 1.0e-6                              # mol, on ~0.65 mol of C3S

    # Refining `nsub` does NOT keep reducing it: past a handful of sub-intervals
    # the error is set by the dense output, not by the quadrature. Both stay at
    # the solver's own accuracy, which is the claim worth testing.
    ξ_fine = reaction_extents(sol, kp; nsub = 16)
    @test extent_residual(sol, kp, ξ_fine) < 1.0e-6
    @test isapprox(ξ_fine[end, 1], ξ[end, 1]; rtol = 1.0e-8)

    # ── A coarse output grid loses no accuracy ────────────────────────────────
    # The quadrature runs on the union of `times` and `sol.t`, so a uniform user
    # grid cannot step over the early transient where the extent accumulates.
    times = collect(range(0.0, kp.tspan[2]; length = 200))
    ξ_grid = reaction_extents(sol, kp; times = times, nsub = 4)
    @test all(ξ_grid[1, :] .== 0.0)
    @test isapprox(ξ_grid[end, 1], ξ[end, 1]; rtol = 1.0e-8)

    # Even an output grid of two points recovers the same final extent
    ξ_2 = reaction_extents(sol, kp; times = [0.0, kp.tspan[2]], nsub = 4)
    @test isapprox(ξ_2[end, 1], ξ[end, 1]; rtol = 1.0e-8)

    # ── A single late instant integrates the WHOLE run up to it ───────────────
    # The extent is measured from the start of the problem, not from `times[1]`.
    # Returning zero here would be silently wrong rather than an error.
    ξ_one = reaction_extents(sol, kp; times = [kp.tspan[2]], nsub = 4)
    @test size(ξ_one) == (1, 2)
    @test isapprox(ξ_one[1, 1], ξ[end, 1]; rtol = 1.0e-8)

    i_mid = length(sol.t) ÷ 2
    ξ_mid = reaction_extents(sol, kp; times = [sol.t[i_mid]], nsub = 4)
    @test isapprox(ξ_mid[1, 1], ξ[i_mid, 1]; rtol = 1.0e-8)

    # ── Guard rails ───────────────────────────────────────────────────────────
    @test_throws ArgumentError reaction_extents(sol, kp; nsub = 0)
    @test_throws ArgumentError reaction_extents(sol, kp; times = [10.0, 0.0])

end

@testset "state_at" begin

    cs, state0, kp, sol = _pp_problem()

    # ── At t = 0 the reconstruction is the initial state ──────────────────────
    s0 = state_at(sol, kp, kp.tspan[1])
    for (i, sp) in enumerate(cs.species)
        @test isapprox(
            ustrip(us"mol", s0.n[i]), ustrip(us"mol", state0.n[i]); atol = 1.0e-12
        )
    end

    # ── `state_at` conserves exactly what the reactions conserve ──────────────
    # `n = n₀ + νᵀξ` gives, for the formula matrix A, `A(n - n₀) = (Aνᵀ)ξ` to
    # machine precision — whatever the quadrature error and whatever the quality
    # of the stoichiometry. This is the property of `state_at` itself, isolated
    # from the database.
    ξ = reaction_extents(sol, kp; nsub = 16)
    st = state_at(sol, kp, kp.tspan[2]; ξ = ξ)
    A = Float64.(cs.SM.A)
    ν = Float64.(kp.ν)

    b0 = A * [ustrip(us"mol", n) for n in state0.n]
    b1 = A * [ustrip(us"mol", n) for n in st.n]
    @test isapprox(b1 .- b0, (A * transpose(ν)) * ξ[end, :]; atol = 1.0e-12)

    # ── How well the reactions close is a property of the DATABASE ────────────
    # CEMDATA18 stores Jennite as (SiO2)1(CaO)1.666667(H2O)2.1 — a rounded 5/3 —
    # so no set of Float64 coefficients balances it better than ~1e-6 per mole
    # of reaction. Element conservation of the run inherits that floor, and it is
    # worth asserting explicitly rather than hiding inside a loose tolerance.
    closure = maximum(abs, A * transpose(ν))
    @test closure < 1.0e-5
    @test isapprox(b0, b1; atol = 10 * closure * maximum(ξ[end, :]))

    # ── The kinetic species agree with the integrator's own values ────────────
    p = sol.prob.p
    u_end = sol(kp.tspan[2])
    for (j, idx) in enumerate(kp.idx_kinetic)
        @test isapprox(
            ustrip(us"mol", st.n[idx]), u_end[p.n_be + j]; atol = 1.0e-6
        )
    end

    # ── Products appeared, reactants and water were consumed ─────────────────
    @test ustrip(us"mol", moles(st, "Jennite")) > 0
    @test ustrip(us"mol", moles(st, "Portlandite")) > 0
    @test ustrip(us"mol", moles(st, "C3S")) < ustrip(us"mol", moles(state0, "C3S"))
    @test ustrip(us"mol", moles(st, "H2O@")) < ustrip(us"mol", moles(state0, "H2O@"))

    # ── No negative amounts leak out ──────────────────────────────────────────
    @test all(ustrip(us"mol", n) >= 0 for n in st.n)

    # ── The reconstructed state feeds the volume balance ─────────────────────
    f = volume_fractions(
        st, [
            "anhydrous" => ["C3S", "C2S"], "C-S-H" => "Jennite",
            "CH" => "Portlandite", "water" => "H2O@",
        ]; reference = state0
    )
    @test isapprox(sum(values(f)), 1.0; atol = 1.0e-12)
    @test f["void"] > 0                       # sealed hydration contracts
    @test f["C-S-H"] > 0 && f["CH"] > 0

end

@testset "degrees_of_hydration" begin

    _, _, kp, sol = _pp_problem()

    α = degrees_of_hydration(sol, kp)
    @test α isa OrderedDict{String, Vector{Float64}}
    @test Set(keys(α)) == Set(["C3S", "C2S"])
    @test all(length(v) == length(sol.t) for v in values(α))

    # ── Bounds and monotonicity ───────────────────────────────────────────────
    for v in values(α)
        @test v[1] ≈ 0.0 atol = 1.0e-12
        @test all(0.0 .<= v .<= 1.0 + 1.0e-12)
        @test all(diff(v) .>= -1.0e-12)
    end

    # ── Alite hydrates faster than belite, from the first hour on ────────────
    # Not from t = 0: the Avrami branch of C3S starts from `PK_AVRAMI_SEED`,
    # which makes it slower than the pure power law governing C2S for the first
    # minutes. The crossover happens well inside the induction period and has no
    # bearing on the volume fractions.
    late = findall(>=(3600.0), sol.t)
    @test all(α["C3S"][late] .>= α["C2S"][late] .- 1.0e-12)
    @test α["C3S"][end] > α["C2S"][end]

    # Order of magnitude at 7 days: alite well advanced, belite far behind.
    # Deliberately loose — the point is coherence, not replication of a
    # particular published curve.
    @test 0.5 < α["C3S"][end] < 0.9
    @test α["C2S"][end] < 0.5

    # ── A species that was never present is omitted, not returned as NaN ─────
    @test !any(any(isnan, v) for v in values(α))

    # ── Mean degree of hydration ─────────────────────────────────────────────
    ᾱ = mean_degree_of_hydration(sol, kp)
    @test length(ᾱ) == length(sol.t)
    # An average lies between the min and the max of its components — stated
    # that way rather than as C2S ≤ ᾱ ≤ C3S, which only holds after the crossover.
    lo = min.(α["C3S"], α["C2S"])
    hi = max.(α["C3S"], α["C2S"])
    @test all(lo .- 1.0e-12 .<= ᾱ .<= hi .+ 1.0e-12)

    # Mass weighting differs from mole weighting when molar masses differ
    ᾱ_mole = mean_degree_of_hydration(sol, kp; weights = :mole)
    @test ᾱ_mole[end] != ᾱ[end]
    @test all(lo .- 1.0e-12 .<= ᾱ_mole .<= hi .+ 1.0e-12)

    @test_throws ArgumentError mean_degree_of_hydration(sol, kp; weights = :volume)

end

@testset "time-dependent rate laws with a Rosenbrock solver" begin

    # `waller` is the first shipped rate law that actually depends on `t`.
    # Rosenbrock methods need a TIME gradient, obtained by calling the RHS with a
    # dual `t` and a plain `u`; typing the rate vector from `eltype(u)` alone
    # then fails with "First call to automatic differentiation for time gradient
    # failed". Rate laws that ignore `t`, like `parrot_killoh`, never exposed it.

    substances = build_species(joinpath(pkgdir(ChemistryLab), "data", "cemdata18-thermofun.json"))
    sp = speciation(
        substances, split("Amor-Sl Portlandite Jennite H2O@");
        aggregate_state = [AS_AQUEOUS]
    )
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)

    state0 = ChemicalState(cs)
    set_quantity!(state0, "Amor-Sl", 0.2u"mol")
    set_quantity!(state0, "Portlandite", 1.0u"mol")
    set_quantity!(state0, "H2O@", 8.0u"mol")

    rxn = Reaction(
        [cs["Amor-Sl"], cs["Portlandite"], cs["H2O@"]], [cs["Jennite"]];
        symbol = "pozzolanic"
    )
    rxn[:rate] = waller(WALLER_PARAMS_SILICA_FUME, "Amor-Sl"; α_max = 0.95)

    kp = KineticsProblem(cs, [rxn], state0, (0.0, 90 * 86400.0); equilibrium_solver = nothing)

    for alg in (Rodas5P(), Rosenbrock23())
        sol = integrate(kp, KineticsSolver(; ode_solver = alg, reltol = 1.0e-6, abstol = 1.0e-10))
        @test sol.retcode == ReturnCode.Success
        α = degrees_of_hydration(sol, kp)["Amor-Sl"]
        @test all(diff(α) .>= -1.0e-10)
        @test 0.0 < α[end] <= 0.95 + 1.0e-9
    end

    # A multistep BDF solver needs no time gradient at all, and must agree
    sol_r = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-12))
    sol_t = integrate(kp, KineticsSolver(; ode_solver = FBDF(), reltol = 1.0e-8, abstol = 1.0e-12))
    @test isapprox(sol_r(kp.tspan[2])[1], sol_t(kp.tspan[2])[1]; rtol = 1.0e-4)

end
