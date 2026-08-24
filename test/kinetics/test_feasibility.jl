using LinearAlgebra

# Helpers guarding the element balance of a re-speciation. They exist because a
# full OPC coupling stopped on a composition demanding 0.732 mol of sulfate
# against the 0.267 mol present — 174 % over — while the reported residual read
# 1.4e-2, because it was normalized by the 34 mol water budget.

@testset "_row_residual scales each element by its own budget" begin

    A = Float64[1 0; 0 1]

    # A violation on a small element must not be hidden by a large one. The
    # second element is out by 0.5, scaled by the larger of its budget (1.0) and
    # the matter flowing through the row (1.5), so a third — not 0.5/100, which
    # is what normalizing by the largest total across all rows used to give.
    @test ChemistryLab._row_residual(A, [100.0, 1.5], [100.0, 1.0]) ≈ 1 / 3
    @test ChemistryLab._row_residual(A, [100.0, 1.5], [100.0, 1.0]) > 0.5 / 100

    # A row whose total is zero — electroneutrality — must not turn a rounding
    # error into a huge residual. Dividing by `b` alone reported 2.1e3 here.
    @test ChemistryLab._row_residual(A, [2.0, 1.0e-12], [2.0, 0.0]) < 1.0e-1

    # An exactly satisfied balance is exactly zero
    @test ChemistryLab._row_residual(A, [2.0, 3.0], [2.0, 3.0]) == 0.0

end

@testset "_budget_clip! caps a guess at the element budget" begin

    # One species consuming 3 units of an element that has 0.267 available: the
    # guess cannot exceed 0.089, whatever it was.
    n = [10.0]
    ChemistryLab._budget_clip!(n, Float64[3;;], [0.267])
    @test n[1] ≈ 0.267 / 3

    # A guess already within budget is untouched
    n2 = [0.01]
    ChemistryLab._budget_clip!(n2, Float64[3;;], [0.267])
    @test n2[1] ≈ 0.01

    # A negative total bounds nothing: H⁺ in a hydrating cement is met by the
    # hydroxides, and clipping on it would forbid every hydrate.
    n3 = [5.0]
    ChemistryLab._budget_clip!(n3, Float64[-2;;], [-11.8])
    @test n3[1] == 5.0

end

@testset "_restore_feasibility! lands in {A n = b, n ≥ 0}" begin

    A = Float64[1 1 0; 0 1 1]
    b = [1.0, 1.0]

    # Starting from a point that is both infeasible and negative
    n = [5.0, -3.0, 7.0]
    ChemistryLab._restore_feasibility!(n, A, b)
    @test maximum(abs, A * n .- b) < 1.0e-10
    @test all(>=(0), n)                        # ends on the box clamp, not the affine step

    # Already feasible: stays feasible, and cheaply
    n2 = [0.5, 0.5, 0.5]
    ChemistryLab._restore_feasibility!(n2, A, b)
    @test maximum(abs, A * n2 .- b) < 1.0e-10
    @test all(>=(0), n2)

    # The cement case that motivated it: ettringite carries 3 sulfate against a
    # 0.267 mol budget, so no feasible point may hold more than 0.089 of it.
    Ae = Float64[3 0; 2 2]                      # rows: sulfate, aluminum
    be = [0.267, 0.488]
    ne = [0.244, 0.0]                           # the point the coupling stalled on
    @test 3 * ne[1] > be[1]                     # infeasible as it stands
    ChemistryLab._restore_feasibility!(ne, Ae, be)
    @test maximum(abs, Ae * ne .- be) < 1.0e-10
    @test all(>=(0), ne)
    @test 3 * ne[1] <= be[1] + 1.0e-10          # sulfate no longer over-spent

end

@testset "speciated_states replays the equilibrium partition" begin

    # A small carbonate system: enough to exercise the replay without paying for
    # a cement. The point is the contract, not the chemistry.
    data = datapath("cemdata18-thermofun.json")
    subs = build_species(data)
    sp = speciation(subs, ["Cal", "Portlandite"]; aggregate_state = [AS_AQUEOUS])
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)

    st0 = ChemicalState(cs)
    set_quantity!(st0, "Cal", 0.05u"kg")
    set_quantity!(st0, "H2O@", 1.0u"kg")

    prim = [p for p in ("Ca+2", "CO3-2", "H2O@", "H+") if p in symbol.(cs.species)]
    rxn = Reaction([cs["Cal"]], [cs[p] for p in prim]; symbol = "calcite dissolution")
    rxn[:rate] = KineticFunc(
        (T, P, t, n, lna, n0) -> 1.0e-6 * max(n["Cal"], zero(eltype(n.data))),
        (T = 298.15u"K", P = 1.0e5u"Pa"), u"mol/s",
    )

    model = DiluteSolutionModel()
    kp = KineticsProblem(
        cs, [rxn], st0, (0.0, 3600.0);
        activity_model = model,
        equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()),
    )
    sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-12))

    times = [600.0, 1800.0, 3600.0]
    states = speciated_states(sol, kp; times = times)

    @test length(states) == length(times)
    @test all(s -> s isa ChemicalState, states)

    # The replayed partition must satisfy the element balance the ODE carried.
    # The FIRST instant is the loose one — it has no previous speciation to start
    # from, only the cast composition — and every one after it, warm-started,
    # lands at machine precision. Asking for an early first instant is therefore
    # not a detail of taste.
    p = sol.prob.p
    res = map(enumerate(times)) do (i, t)
        be = collect(sol(t)[1:(p.n_be)])
        ne = Float64[ustrip(us"mol", moles(states[i], symbol(cs.species[j]))) for j in kp.idx_equilibrium]
        ChemistryLab._row_residual(p.Ae, ne, be)
    end
    # Reported as a number, not as an `all`: when this fails, what matters is by
    # how much, and a bare `all` hides it.
    @test res[1] < 1.0e-2
    worst_replay = maximum(res[2:end])
    @test worst_replay < 1.0e-10

    # The kinetic species come from the ODE state, not from the equilibrium solve
    for (i, t) in enumerate(times)
        @test ustrip(us"mol", moles(states[i], "Cal")) ≈ sol(t)[p.n_be + 1] atol = 1.0e-10
    end

    # Calcite dissolves, so its amount decreases along the replay
    @test ustrip(us"mol", moles(states[1], "Cal")) > ustrip(us"mol", moles(states[end], "Cal"))

    # Descending instants are refused: each solve warm-starts from the previous
    @test_throws ArgumentError speciated_states(sol, kp; times = [3600.0, 600.0])

end

@testset "the equilibrium partition keeps its solid solutions" begin

    # `_equilibrium_subsystem` rebuilds a `ChemicalSystem` for the partition the
    # equilibrium is solved on. Until 0.8.2 it dropped `solid_solutions` on the
    # way, silently and completely: a coupled run could declare CSHQ, AFm or
    # Hydrogarnet and the solve would treat their end-members as separate pure
    # phases, the mixing entropy never entering the Gibbs energy.

    data = datapath("cemdata18-thermofun.json")
    substances = build_species(data)
    toml = datapath("solid_solutions.toml")

    members = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD"]
    sp = speciation(substances, vcat("C3S", "Portlandite", members); aggregate_state = [AS_AQUEOUS])
    byname = Dict(symbol(s) => s for s in sp)
    ss = [x for x in build_solid_solutions(toml, byname) if x.name == "CSHQ"]
    @test length(ss) == 1                       # the declaration is usable at all

    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)
    @test cs.solid_solutions !== nothing

    # The partition excludes the kinetic species, here C3S. All four CSHQ
    # end-members stay, so the solution must survive into the sub-system.
    idx_eq = [i for (i, s) in enumerate(cs.species) if symbol(s) != "C3S"]
    sub = ChemistryLab._equilibrium_subsystem(cs, idx_eq)

    @test sub.solid_solutions !== nothing
    @test length(sub.solid_solutions) == 1
    @test sub.solid_solutions[1].name == "CSHQ"
    @test Set(symbol.(sub.solid_solutions[1].end_members)) == Set(members)

    # A solution whose end-members are NOT all in the partition cannot be mixed
    # there, and must be dropped rather than passed truncated.
    idx_partial = [i for (i, s) in enumerate(cs.species) if symbol(s) != "CSHQ-JenD"]
    sub2 = ChemistryLab._equilibrium_subsystem(cs, idx_partial)
    @test sub2.solid_solutions === nothing

    # A system with no solid solutions at all still yields a partition, and
    # `nothing` rather than an empty vector.
    cs_plain = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    sub3 = ChemistryLab._equilibrium_subsystem(cs_plain, idx_eq)
    @test sub3.solid_solutions === nothing

end
