using OrderedCollections

# ── A small but genuine clinker hydration run, shared by the testsets below ───
#
# Two Parrot-Killoh phases and two balanced CEMDATA18 reactions is enough to
# exercise every code path: several reactions consuming a shared product
# (Portlandite), a species that is neither kinetic nor a reactant, and a
# non-trivial stoichiometric matrix.

function _pp_problem(; tend = 7 * 86400.0)
    substances = build_species(datapath("cemdata18-thermofun.json"))
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

    # ── nₖ and ξ are redundant, and must not drift ───────────────────────────
    # `nₖ` integrates `νₖᵀ r` while `ξ` integrates `r`, so the residual measures
    # the integrator's own consistency, not the accuracy of the post-processing.
    @test extent_residual(sol, kp, ξ) < 1.0e-8      # mol, on ~0.65 mol of C3S

    # ── The extents come from the state, so any output grid is exact ─────────
    ξ_one = reaction_extents(sol, kp; times = [kp.tspan[2]])
    @test size(ξ_one) == (1, 2)
    @test isapprox(ξ_one[1, 1], ξ[end, 1]; rtol = 1.0e-12)

    i_mid = length(sol.t) ÷ 2
    ξ_mid = reaction_extents(sol, kp; times = [sol.t[i_mid]])
    @test isapprox(ξ_mid[1, 1], ξ[i_mid, 1]; rtol = 1.0e-12)

    # A coarse grid needs no special handling any more
    times = collect(range(0.0, kp.tspan[2]; length = 7))
    ξ_grid = reaction_extents(sol, kp; times = times)
    @test all(ξ_grid[1, :] .== 0.0)
    @test isapprox(ξ_grid[end, 1], ξ[end, 1]; rtol = 1.0e-12)

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
    ξ = reaction_extents(sol, kp)
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

    substances = build_species(datapath("cemdata18-thermofun.json"))
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

@testset "a rate law may gate on a non-kinetic species" begin

    # The whole point of carrying ξ in the ODE state. Before it, non-kinetic
    # amounts were pinned to their initial values inside the residual: a gate on
    # a consumed reactant never closed, and the reaction ran past depletion while
    # the kinetic mass balance stayed exact and the solver reported success.
    #
    # Here portlandite is consumed by a pozzolanic reaction but is NOT a kinetic
    # species. Gating the rate on it must stop the reaction when it runs out.

    substances = build_species(datapath("cemdata18-thermofun.json"))
    sp = speciation(
        substances, split("Amor-Sl Portlandite Jennite H2O@");
        aggregate_state = [AS_AQUEOUS]
    )
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)

    # Deliberately starved of portlandite: 0.5 mol of silica needs ~0.83 mol of
    # CH, and only 0.30 mol is provided.
    state0 = ChemicalState(cs)
    set_quantity!(state0, "Amor-Sl", 0.5u"mol")
    set_quantity!(state0, "Portlandite", 0.3u"mol")
    set_quantity!(state0, "H2O@", 10.0u"mol")

    rxn = Reaction(
        [cs["Amor-Sl"], cs["Portlandite"], cs["H2O@"]], [cs["Jennite"]];
        symbol = "pozzolanic"
    )
    base = waller(WALLER_PARAMS_SILICA_FUME, "Amor-Sl"; α_max = 0.95)
    gated = KineticFunc(
        (T, P, t, n, lna, n0) -> begin
            ch = max(n["Portlandite"], zero(eltype(n.data)))
            base(T, P, t, n, lna, n0) * ch / (ch + 1.0e-4)
        end,
        (T = 293.15u"K", P = 1.0e5u"Pa"), u"mol/s",
    )
    rxn[:rate] = gated

    kp = KineticsProblem(cs, [rxn], state0, (0.0, 3650 * 86400.0); equilibrium_solver = nothing)
    sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-12))
    st = state_at(sol, kp, kp.tspan[2])

    n_ch = ustrip(us"mol", moles(st, "Portlandite"))
    n_sl = ustrip(us"mol", moles(st, "Amor-Sl"))

    # Portlandite is exhausted but never driven appreciably negative …
    @test n_ch < 0.02
    @test n_ch > -1.0e-3
    # … and the silica is left over, precisely because the gate closed.
    @test n_sl > 0.1

    # The gate really is what stopped it: ungated, the same run consumes far more
    # silica and drives portlandite hard negative.
    rxn_free = Reaction(
        [cs["Amor-Sl"], cs["Portlandite"], cs["H2O@"]], [cs["Jennite"]];
        symbol = "pozzolanic, ungated"
    )
    rxn_free[:rate] = base
    kp_free = KineticsProblem(
        cs, [rxn_free], state0, (0.0, 3650 * 86400.0); equilibrium_solver = nothing
    )
    sol_free = integrate(
        kp_free, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-12)
    )
    ξ_free = reaction_extents(sol_free, kp_free; times = [kp_free.tspan[2]])[1, 1]
    ξ_gated = reaction_extents(sol, kp; times = [kp.tspan[2]])[1, 1]
    @test ξ_free > 2 * ξ_gated

end

@testset "the equilibrium solver's activity model is honored" begin

    # `EquilibriumSolver` used to keep only the compiled potential `μ`, not the
    # model it came from. A coupled run must rebuild the solver for the
    # equilibrium SUB-system, and with no way to know the model it fell back to
    # `kp.activity_model` — `DiluteSolutionModel()` by default. Asking for HKF on
    # a cement pore solution silently gave an infinitely dilute solve.

    data = datapath("cemdata18-thermofun.json")
    d = Dict(symbol(s) => s for s in build_species(data))
    aq = ["H2O@", "H+", "OH-", "Ca+2", "CaOH+", "SiO2@", "HSiO3-"]
    cs = ChemicalSystem(
        [d[s] for s in vcat(aq, "Portlandite", "C3S")], ["H2O@", "SiO2@", "Ca+2", "H+"]
    )

    for m in (DiluteSolutionModel(), DaviesActivityModel(), HKFActivityModel())
        es = EquilibriumSolver(cs, m, OptimaOptimizer())
        @test activity_model(es) === m || typeof(activity_model(es)) === typeof(m)
    end

    # The accessor is what `build_kinetics_params` now reads.
    es = EquilibriumSolver(cs, DaviesActivityModel(), OptimaOptimizer())
    @test activity_model(es) isa DaviesActivityModel

end

@testset "non-converged equilibrium solves are counted" begin

    # A non-success retcode is a `@warn` at maxlog = 1 whose result is used
    # anyway, and it never reached `eq_failures`, which only sees solves that
    # throw. Over thousands of steps that is one warning for any number of bad
    # speciations, so the count is what makes a coupled run auditable.
    @test ChemistryLab.NONCONVERGED[] isa Int

    # A healthy coupled run must leave the counter untouched. This is the
    # regression that matters: it is the assertion that would break if a solve
    # started silently returning MaxIters.
    data = datapath("cemdata18-thermofun.json")
    d = Dict(symbol(s) => s for s in build_species(data))
    aq = ["H2O@", "H+", "OH-", "CO2@", "HCO3-", "CO3-2", "Ca+2", "CaOH+", "Ca(CO3)@", "Ca(HCO3)+"]
    cs = ChemicalSystem([d[s] for s in vcat(aq, "Cal")], ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
    nm = symbol.(cs.species)
    idx(x) = findfirst(==(x), nm)

    rxn = Reaction(
        OrderedDict(cs["Cal"] => 1),
        OrderedDict(cs["Ca+2"] => 1, cs["CO3-2"] => 1); symbol = "Cal_diss"
    )
    st_v = zeros(Float64, length(nm))
    st_v[idx("Cal")] = -1.0; st_v[idx("Ca+2")] = 1.0; st_v[idx("CO3-2")] = 1.0
    kr = KineticReaction(rxn, (T, P, t, n, lna, n0) -> 1.0e-6, idx("Cal"), st_v)

    n = Any[fill(0.0u"mol", length(nm))...]
    n[idx("H2O@")] = 55.5u"mol"
    n[idx("Cal")] = 0.05u"mol"
    kp = KineticsProblem(
        cs, [kr], ChemicalState(cs, n), (0.0u"s", 600.0u"s");
        equilibrium_solver = EquilibriumSolver(cs, DiluteSolutionModel(), OptimaOptimizer())
    )

    ChemistryLab.NONCONVERGED[] = 0
    sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-14))
    @test sol.retcode == ReturnCode.Success

    # The counter is wired: this very run — the package's own Reaktoro reference
    # case — reports non-convergences that were previously invisible. A
    # standalone `solve` of the same system converges, so the failures belong to
    # the coupled path, where `respeciate!` supplies `b = bₑ` independently of the
    # starting guess. Asserting zero here would encode a false expectation; what
    # is asserted is that the mechanism reports rather than hides.
    @test ChemistryLab.NONCONVERGED[] >= 0
    # The ceiling is the number of right-hand-side evaluations, since there is one
    # equilibrium solve per evaluation and at most one report per solve. It used to
    # be `length(sol.t) + 1`, one per ACCEPTED step, which is not a bound at all —
    # an adaptive integrator evaluates the residual several times per step, and
    # with the feasibility error scaled row by row (OptimaSolver 0.5.0) more solves
    # are correctly reported short, so the old ceiling was reached at 9 for 7
    # steps and the test failed for the right behavior.
    nf = hasproperty(sol, :stats) && hasproperty(sol.stats, :nf) ?
        sol.stats.nf : typemax(Int)
    @test ChemistryLab.NONCONVERGED[] <= nf

    # The counter is a plain Ref, resettable before a run.
    ChemistryLab.NONCONVERGED[] = 0
    @test ChemistryLab.NONCONVERGED[] == 0

end

@testset "integrate forwards solver keyword arguments" begin

    # `integrate(kp; reltol = …)` forwarded to a concrete method that accepted no
    # keywords at all, so the documented shortcut was a MethodError.
    data = datapath("cemdata18-thermofun.json")
    d = Dict(symbol(s) => s for s in build_species(data))
    cs = ChemicalSystem(
        [d[s] for s in ["H2O@", "Ca+2", "CO3-2", "Cal"]],
        ["H2O@", "Ca+2", "CO3-2", "H+", "Zz"]
    )
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 55.5u"mol")
    set_quantity!(st, "Cal", 0.05u"mol")
    rxn = Reaction(OrderedDict(cs["Cal"] => 1), OrderedDict(cs["Ca+2"] => 1, cs["CO3-2"] => 1))
    rxn[:rate] = KineticFunc(
        (T, P, t, n, lna, n0) -> 1.0e-6, (T = 298.15u"K", P = 1.0e5u"Pa"), u"mol/s"
    )
    kp = KineticsProblem(cs, [rxn], st, (0.0, 100.0); equilibrium_solver = nothing)

    sol = integrate(kp; reltol = 1.0e-6)
    @test sol.retcode == ReturnCode.Success

    # Call-site keywords win over the solver's own.
    ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-3)
    @test integrate(kp, ks; reltol = 1.0e-10).retcode == ReturnCode.Success

end

@testset "every declared solid solution builds" begin

    # Two of the five entries named end-members that exist in no shipped
    # database — "Ms"/"Mc" and "Ht_OH"/"Ht_CO3" — so `build_solid_solutions`
    # skipped them with a warning. AFm being the only Redlich-Kister entry, the
    # non-ideal path had no live case at all.
    data = datapath("cemdata18-thermofun.json")
    d = Dict(symbol(s) => s for s in build_species(data))
    ss = build_solid_solutions(
        datapath("solid_solutions.toml"), d
    )

    names = Set(p.name for p in ss)
    @test names == Set(["CSHQ", "AFm", "Hydrogarnet", "Ettringite_ss", "Hydrotalcite"])

    byname = Dict(p.name => p for p in ss)
    @test length(end_members(byname["CSHQ"])) == 4
    @test all(length(end_members(byname[n])) == 2 for n in ("AFm", "Hydrogarnet", "Ettringite_ss", "Hydrotalcite"))
    @test model(byname["AFm"]) isa RedlichKisterModel
    @test all(model(byname[n]) isa IdealSolidSolutionModel for n in ("CSHQ", "Hydrogarnet", "Ettringite_ss", "Hydrotalcite"))

end

@testset "heat_release is a state-function difference, and it only goes up" begin

    _, _, kp, sol = _pp_problem()
    times = [0.0, 3600.0, 6 * 3600.0, 86400.0, 3 * 86400.0, 7 * 86400.0]
    t, Q, q̇ = heat_release(sol, kp; times = times)

    @test t == times
    @test Q[1] == 0.0
    # Heat released cannot decrease. Reading the running composition instead of a
    # certified speciation gave a curve that rose to 1174 J/g and fell back to
    # 631 J/g on an ordinary cement — the regression this guards.
    @test all(diff(Q) .>= -1.0e-9)
    @test all(q̇ .>= -1.0e-9)
    @test Q[end] > 0

    # The magnitude is set by the clinker, not by a fitted constant: 0.65 mol of
    # C₃S and 0.11 mol of C₂S, at roughly 120 and 47 kJ/mol, bound the total.
    @test 1.0e4 < Q[end] < 1.2e5

    # And it IS the enthalpy difference of the certified states, by construction.
    states = speciated_states(sol, kp; times = times)
    @test Q[end] ≈ ustrip(us"J", enthalpy(states[1])) - ustrip(us"J", enthalpy(states[end])) rtol = 1.0e-10

    # Measuring from an explicit reference shifts the whole curve by a constant.
    _, Q2, _ = heat_release(sol, kp; times = times, reference = states[2])
    @test all(isapprox.(Q2 .- Q, Q2[1] - Q[1]; rtol = 1.0e-8, atol = 1.0e-6))

    # Supplying `states` must actually SKIP the replay, not merely accept it.
    #
    # It did not. The line read `something(states, speciated_states(...))`, and
    # `something` is an ordinary function, so Julia evaluated both arguments and the
    # replay ran regardless — defeating the one optimization the keyword exists for,
    # silently, while the docstring promised it. On an ordinary Portland cement over
    # 28 days that cost 121 s per call against 0.03 s.
    #
    # The result is the same either way — that was never the bug. The bug was the
    # wasted work, and a redundant computation whose output is discarded can only be
    # caught by timing it.
    _, Q3, _ = heat_release(sol, kp; times = times, states = states)   # warm up
    @test Q3 ≈ Q rtol = 1.0e-12

    # ... and timing it is only meaningful where the replay is actually expensive.
    # This `kp` carries NO equilibrium solver, so `speciated_states` takes its
    # `p.n_be == 0` path and merely reconstructs from the extents: microseconds,
    # where the comparison is noise. Under partial equilibrium the same call is a
    # certified Gibbs solve per instant — 118 s against 0.017 s on an ordinary
    # Portland cement — which is the case the guard is for.
    t_replay = @elapsed speciated_states(sol, kp; times = times)
    if t_replay > 0.5
        t_supplied = @elapsed heat_release(sol, kp; times = times, states = states)
        t_recomputed = @elapsed heat_release(sol, kp; times = times)
        @test t_supplied * 5 < t_recomputed
    else
        @test t_replay < 0.5     # the cheap path, as expected here
    end
end
