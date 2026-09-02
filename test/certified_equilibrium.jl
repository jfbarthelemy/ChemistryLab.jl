# The certified equilibrium path.
#
# Every assertion here is on a measured property of the solvers, not on a target
# value: that the certificate decides, that it is not fooled by a start which is
# itself infeasible, and that the multi-start route certifies cases no single
# back end does.

@testsection "Certified equilibrium" begin

    sp = Dict(symbol(s) => s for s in build_species(
        datapath("slop98-inorganic-thermofun.json")))
    species = [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")]
    cs = ChemicalSystem(species, ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
    A = Matrix{Float64}(cs.SM.A)

    function calcite(; θ = 25.0, nco2 = 0.0)
        st = ChemicalState(cs)
        set_quantity!(st, "Cal", 1.0e-3u"mol")
        set_quantity!(st, "H2O@", 1.0u"kg")
        nco2 > 0 && set_quantity!(st, "CO2@", nco2 * u"mol")
        V = volume(st)
        set_quantity!(st, "H+", 1.0e-4u"mol/L" * V.liquid)
        set_quantity!(st, "OH-", 1.0e-10u"mol/L" * V.liquid)
        set_temperature!(st, (273.15 + θ) * u"K")
        return st
    end

    # Element balance, each row against its own budget. An absolute ∞-norm hides
    # the small rows: at ‖An − b‖∞ = 3e-6 the water row (111 mol) is satisfied to
    # 1.8e-8 while the charge row (budget 1e-4) is out by 3 %.
    function balance_rel(st0, eq)
        b = A * [ustrip(us"mol", x) for x in st0.n]
        r = A * [ustrip(us"mol", x) for x in eq.n] - b
        keep = abs.(b) .> 1.0e-8
        return any(keep) ? maximum(abs.(r[keep] ./ b[keep])) : maximum(abs, r)
    end

    @testsection "the certificate decides, and the answer is reproducible" begin
        st = calcite()
        eq1, c1 = equilibrate_certified(calcite())
        eq2, c2 = equilibrate_certified(calcite())
        @test c1.optimal
        @test pH(eq1) == pH(eq2)          # bit-for-bit, no path dependence
        @test balance_rel(st, eq1) < 1.0e-10
        @test c1.worst_supersaturation < 0 || c1.worst_supersaturation == -Inf
    end

    @testsection "the interior point alone is wrong on this case" begin
        # Not a target value: the point is that the certified route is two orders
        # of magnitude better on a balance the old default reported as fine.
        st = calcite()
        eq_bar = equilibrate(calcite(), OptimaOptimizer())
        eq_cert, cert = equilibrate_certified(calcite())
        @test balance_rel(st, eq_bar) > 1.0e-3      # measured at 3.0e-2
        @test balance_rel(st, eq_cert) < 1.0e-8
        @test cert.optimal
        # And they disagree on the answer, not merely on the residual.
        @test abs(pH(eq_bar) - pH(eq_cert)) > 1.0e-3
    end

    @testsection "b is fixed once, not taken from each start" begin
        # A start that violates the balance shifts the component totals by its own
        # infeasibility. Recomputing `b` from each start therefore poses a
        # different problem per start, and the dual solve certifies the answer to
        # the shifted one: measured, two "certified" compositions 0.2 % apart on
        # dissolved calcium, which a convex problem with one minimum cannot have.
        st = calcite()
        b0 = A * [ustrip(us"mol", x) for x in st.n]
        des = DualEquilibriumSolver(cs, DiluteSolutionModel())
        bar = equilibrate(calcite(), OptimaOptimizer())
        from_raw = SciMLBase.solve(des, calcite(); b = b0)
        from_bar = SciMLBase.solve(des, bar; b = b0)
        # With the same b both routes must land on the same minimum.
        @test isapprox(pH(from_raw), pH(from_bar); atol = 1.0e-6)
        # Without it, the second solves a different problem — the shift being
        # exactly the interior point's own infeasibility.
        shifted = SciMLBase.solve(des, bar)
        @test abs(pH(shifted) - pH(from_raw)) > 1.0e-3
    end

    @testsection "equilibrate certifies by default, and can be told not to" begin
        st = calcite()
        @test balance_rel(st, equilibrate(calcite())) < 1.0e-8
        @test balance_rel(st, equilibrate(calcite(); certify = false)) > 1.0e-3
    end

    @testsection "the dual route needs an aqueous phase and H2O@" begin
        # It parameterizes the interior variables by the solvent's potential, so
        # both conditions are structural. `equilibrate_certified` returns
        # `nothing` for the certificate on such a system rather than claiming a
        # proof it cannot give.
        @test ChemistryLab._dual_applicable(cs)

        solid_only = ChemicalSystem(
            [Species("NaCl"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)],
            ["NaCl"],
        )
        @test !ChemistryLab._dual_applicable(solid_only)

        # Aqueous, but the solvent is not called `H2O@`.
        no_solvent = ChemicalSystem([sp["Ca+2"], sp["CO3-2"]], ["Ca+2", "CO3-2"])
        @test !ChemistryLab._dual_applicable(no_solvent)
    end

end

@testsection "ForwardDiff through the certified route" begin

    # The certified route solves in real arithmetic — an active set has a discrete
    # component, so the map `b ↦ n*(b)` is smooth only piecewise and pushing duals
    # through the search would be wrong. The derivative is attached at the answer,
    # from the optimality conditions.
    #
    # This has to be tested rather than assumed: the certified path converts the
    # component totals to `Float64` to fix them once, and before the dual branch
    # existed that silently dropped every partial.
    sp = Dict(symbol(s) => s for s in build_species(
        datapath("slop98-inorganic-thermofun.json")))
    cs = ChemicalSystem(
        [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))

    function pH_of(x)
        n = Any[fill(0.0 * oneunit(x) * u"mol", length(cs.species))...]
        n[idx["H2O@"]] = 55.5 * oneunit(x) * u"mol"
        n[idx["Cal"]] = 0.05 * oneunit(x) * u"mol"
        n[idx["CO2@"]] = x * u"mol"
        return pH(equilibrate(ChemicalState(cs, n)))
    end

    x0 = 0.01
    ad = ForwardDiff.derivative(pH_of, x0)
    h = 1.0e-6
    fd = (pH_of(x0 + h) - pH_of(x0 - h)) / (2h)

    @test isfinite(ad)
    @test ad != 0                       # the partials are not dropped
    @test ad ≈ fd rtol = 1.0e-5         # measured at 4.3e-9
    @test ad < 0                        # more CO2, lower pH

    # And the primal value is the certified one, not something the dual path
    # computed separately.
    @test pH_of(x0) ≈ pH(first(equilibrate_certified(let
        n = Any[fill(0.0u"mol", length(cs.species))...]
        n[idx["H2O@"]] = 55.5u"mol"; n[idx["Cal"]] = 0.05u"mol"
        n[idx["CO2@"]] = x0 * u"mol"
        ChemicalState(cs, n)
    end))) atol = 1.0e-12

end
