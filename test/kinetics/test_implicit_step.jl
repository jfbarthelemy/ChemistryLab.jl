# The implicit kinetic step of Leal et al. (2017).
#
# One Gibbs minimization per step, with the reaction extents as unknowns of the
# same system and the rate evaluated at the end-of-step composition. The
# assertions below are on properties the formulation must have, not on values
# fitted to it: an analytic answer for a constant rate, and the impossibility of
# crossing saturation at any step size.

using LinearAlgebra

@testsection "Implicit kinetic step" begin

    sp = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json")
            )
    )
    cs = ChemicalSystem(
        [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    A = Matrix{Float64}(cs.SM.A)
    i_Cal = findfirst(s -> symbol(s) == "Cal", cs.species)

    rxn = Reaction(
        OrderedDict(cs["Cal"] => 1.0),
        OrderedDict(cs["Ca+2"] => 1.0, cs["CO3-2"] => 1.0);
        symbol = "calcite dissolution",
    )

    function calcite_water()
        st = ChemicalState(cs)
        set_quantity!(st, "Cal", 1.0e-2u"mol")
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "H+", 1.0e-7u"mol")
        set_quantity!(st, "OH-", 1.0e-7u"mol")
        return st
    end

    @testsection "a constant rate gives the analytic answer" begin
        k = 1.0e-6
        kr = KineticReaction(
            cs, rxn, KineticFunc((T, P, t, n, lna, n0) -> k, NamedTuple(), u"mol/s"),
        )
        kss = KineticStepSolver(cs, DiluteSolutionModel(), [kr])

        # `K` is restricted to the non-aqueous participants, so on this reaction it
        # is the single entry −1 and `M = KᵀK = 1`. That is what makes `Δξ` the
        # reaction progress rather than `M` times it.
        @test size(kss.K) == (length(cs.species), 1)
        @test kss.K[i_Cal, 1] == -1.0
        @test kss.M ≈ [1.0]

        b0 = A * Float64[ustrip(us"mol", x) for x in calcite_water().n]
        for Δt in (10.0, 100.0, 1000.0)
            q = Ref(Float64[])
            s1 = kinetic_step(kss, calcite_water(), Δt; parameters = q)
            n1 = Float64[ustrip(us"mol", x) for x in s1.n]
            @test q[][1] ≈ Δt * k rtol = 1.0e-10          # Δξ = Δt·r
            @test n1[i_Cal] ≈ 1.0e-2 - k * Δt atol = 1.0e-13
            @test maximum(abs, A * n1 - b0) < 1.0e-10     # elements conserved
        end
    end

    @testsection "a rate reading the activities cannot cross saturation" begin
        # `r = k(1 − Ω)`: the feedback path from the aqueous equilibrium back into
        # the rate law. An EXPLICIT step of 1e6 s at k = 1e-5 would dissolve 10 mol,
        # a thousand times the calcite present. The implicit step cannot, because
        # the rate is evaluated at the composition it is producing.
        k = 1.0e-5
        stoich = Float64.(
            KineticReaction(
                cs, rxn,
                KineticFunc((T, P, t, n, lna, n0) -> 0.0, NamedTuple(), u"mol/s"),
            ).stoich,
        )
        R = 8.31446261815324
        GT = [
            ustrip(us"J/mol", s[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true)) /
                (R * 298.15) for s in cs.species
        ]
        kr = KineticReaction(
            cs, rxn,
            KineticFunc(
                (T, P, t, n, lna, n0) -> begin
                    Ω = saturation_ratio(
                        stoich, [lna[symbol(s)] for s in cs.species], GT,
                    )
                    k * (1 - Ω)
                end, NamedTuple(), u"mol/s",
            ),
        )
        kss = KineticStepSolver(cs, DiluteSolutionModel(), [kr])
        des = DualEquilibriumSolver(cs, DiluteSolutionModel())
        n_eq = Float64[
            ustrip(us"mol", x) for x in SciMLBase.solve(des, calcite_water()).n
        ]

        dissolved = Float64[]
        for Δt in (1.0, 100.0, 1.0e4, 1.0e6)
            n1 = Float64[
                ustrip(us"mol", x) for x in
                    kinetic_step(kss, calcite_water(), Δt).n
            ]
            # Never past equilibrium, however large the step.
            @test n1[i_Cal] >= n_eq[i_Cal] - 1.0e-12
            # And never more than the mineral present.
            @test n1[i_Cal] > 0
            push!(dissolved, 1.0e-2 - n1[i_Cal])
        end
        @test issorted(dissolved)                       # monotone in Δt
        # The long step approaches the equilibrium answer FROM BELOW and does not
        # reach it: `r = k(1 − Ω)` vanishes only at Ω = 1, so a finite step leaves
        # Ω = 0.999989 and the last 4e-6 relative undissolved. That is the rate
        # law's behavior, not a tolerance — hence the one-sided assertion.
        @test dissolved[end] < 1.0e-2 - n_eq[i_Cal]
        @test dissolved[end] ≈ 1.0e-2 - n_eq[i_Cal] rtol = 1.0e-4
        # An explicit step would have dissolved a thousand times too much.
        @test k * 1.0e6 > 1000 * dissolved[end]
    end

    @testsection "several reactions may share a mineral" begin
        CEM = Dict(
            symbol(s) => s for s in build_species(
                    datapath("cemdata18-thermofun.json")
                )
        )
        inp = split("C3A Gp H2O@ ettringite monosulphate12 Portlandite")
        spc = speciation(collect(values(CEM)), inp; aggregate_state = [AS_AQUEOUS])
        cs2 = ChemicalSystem(spc, CEMDATA_PRIMARIES)
        S(x) = cs2[x]
        r_aft = Reaction(
            OrderedDict(S("C3A") => 1.0, S("Gp") => 3.0, S("H2O@") => 26.0),
            OrderedDict(S("ettringite") => 1.0); symbol = "C3A -> AFt",
        )
        r_afm = Reaction(
            OrderedDict(S("C3A") => 1.0, S("Gp") => 1.0, S("H2O@") => 10.0),
            OrderedDict(S("monosulphate12") => 1.0); symbol = "C3A -> AFm",
        )
        f(k) = KineticFunc((T, P, t, n, lna, n0) -> k, NamedTuple(), u"mol/s")
        k1, k2 = 2.0e-6, 5.0e-7
        krs = [KineticReaction(cs2, r_aft, f(k1)), KineticReaction(cs2, r_afm, f(k2))]

        kss = KineticStepSolver(cs2, DiluteSolutionModel(), krs)
        # The two pathways must be distinguishable, or `Kᵀn` could not tell their
        # extents apart and the linear system would be singular.
        @test rank(kss.K) == 2
        @test size(kss.M) == (2, 2)

        st = ChemicalState(cs2)
        set_quantity!(st, "C3A", 1.0e-2u"mol")
        set_quantity!(st, "Gp", 3.0e-2u"mol")
        set_quantity!(st, "H2O@", 1.0u"kg")
        Δt = 1000.0
        q = Ref(Float64[])
        kinetic_step(kss, deepcopy(st), Δt; parameters = q)
        # `Δξ = Δt·M·r`, and with two pathways `M` is not the identity: the
        # measured extents are exactly that product, cross terms included.
        @test q[] ≈ Δt * (kss.M * [k1, k2]) rtol = 1.0e-8
    end

    @testsection "what is refused, and why" begin
        k = 1.0e-6
        kr = KineticReaction(
            cs, rxn, KineticFunc((T, P, t, n, lna, n0) -> k, NamedTuple(), u"mol/s"),
        )

        # No reaction at all: that is a plain equilibrium.
        @test_throws ArgumentError KineticStepSolver(
            cs, DiluteSolutionModel(), KineticReaction[],
        )

        @test_throws ArgumentError KineticStepSolver(
            cs, DiluteSolutionModel(), [kr]; coupling = :nonsense,
        )

        # Two reactions differing only in aqueous products give the same row twice.
        kr2 = KineticReaction(
            cs,
            Reaction(
                OrderedDict(cs["Cal"] => 1.0),
                OrderedDict(cs["Ca+2"] => 1.0, cs["HCO3-"] => 1.0, cs["OH-"] => -1.0);
                symbol = "calcite, bicarbonate route",
            ),
            KineticFunc((T, P, t, n, lna, n0) -> k, NamedTuple(), u"mol/s"),
        )
        @test_throws ArgumentError KineticStepSolver(
            cs, DiluteSolutionModel(), [kr, kr2],
        )
    end

end

@testsection "a solid solution in the kinetic step" begin

    # The coupling the package exists for: a mineral dissolving under a rate law
    # while the products — here a four-member C-S-H solid solution — are decided
    # by thermodynamics at the end of every step.
    substances = build_species(datapath("cemdata18-thermofun.json"))
    members = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD"]
    spc = speciation(
        substances, vcat("Portlandite", "C3S", members); aggregate_state = [AS_AQUEOUS],
    )
    ss = [
        x for x in build_solid_solutions(
                datapath("solid_solutions.toml"), Dict(symbol(s) => s for s in spc),
            ) if x.name == "CSHQ"
    ]
    cs2 = ChemicalSystem(spc, CEMDATA_PRIMARIES; solid_solutions = ss)
    S(x) = cs2[x]

    # C3S + 3 H2O -> 3 Ca+2 + SiO2@ + 6 OH-, balanced in matter and in charge.
    rxn2 = Reaction(
        OrderedDict(S("C3S") => 1.0, S("H2O@") => 3.0),
        OrderedDict(S("Ca+2") => 3.0, S("SiO2@") => 1.0, S("OH-") => 6.0);
        symbol = "C3S dissolution",
    )
    k = 2.0e-5
    kr = KineticReaction(
        cs2, rxn2, KineticFunc((T, P, t, n, lna, n0) -> k, NamedTuple(), u"mol/s"),
    )
    kss = KineticStepSolver(cs2, DiluteSolutionModel(), [kr])

    names2 = symbol.(cs2.species)
    i_C3S = findfirst(==("C3S"), names2)
    function fresh2()
        st = ChemicalState(cs2)
        set_quantity!(st, "C3S", 1.0u"mol")
        set_quantity!(st, "H2O@", 40.0u"mol")
        return st
    end

    totals = Float64[]
    for Δt in (1.0e3, 1.0e4)
        q = Ref(Float64[])
        cert = Ref{Any}(nothing)
        s1 = kinetic_step(kss, fresh2(), Δt; parameters = q, certificate = cert)
        n1 = Float64[ustrip(us"mol", x) for x in s1.n]
        amounts = [n1[findfirst(==(nm), names2)] for nm in members]

        # The extent and the mineral follow the rate law exactly.
        @test q[][1] ≈ k * Δt rtol = 1.0e-8
        @test n1[i_C3S] ≈ 1.0 - k * Δt atol = 1.0e-12

        # The solid solution FORMS, and over all four end-members — which is what
        # distinguishes a solution from four separate pure phases.
        @test all(>(1.0e-6), amounts)

        # And the step is certified — on the AUGMENTED problem, whose reactivity
        # rows carry the multipliers that make a kinetically held mineral's
        # stationarity satisfiable. Certified against the unconstrained problem it
        # would report a residual of 7 RT, because C3S is supersaturated by
        # construction: that is what being held back means.
        @test cert[].optimal
        @test cert[].stationarity < 1.0e-10
        @test cert[].balance < 1.0e-10

        push!(totals, sum(amounts))
    end

    # Ten times the step, ten times the C-S-H.
    @test totals[2] / totals[1] ≈ 10 rtol = 1.0e-2

    # Without the warm start the mixing phase is never admitted: the tangent-plane
    # test is made at a composition where the phase is absent, and it fails. This
    # is a property of the cold start, not of the kinetics — a plain equilibrium
    # from the same guess fails identically.
    cold = kinetic_step(kss, fresh2(), 1.0e3; warm_start = false)
    n_cold = Float64[ustrip(us"mol", x) for x in cold.n]
    @test count(>(1.0e-6), [n_cold[findfirst(==(nm), names2)] for nm in members]) < 4

end

@testsection "the step length can be chosen from a local error estimate" begin

    # What Reaktoro's kinetics does not have. Its options add a single initial
    # step to the equilibrium options; the step is the caller's and nothing
    # reports what taking it cost. The scheme is backward Euler, so a step ten
    # times too large is ten times less accurate, silently.
    sp = Dict(
        symbol(x) => x for x in build_species(
                datapath("slop98-inorganic-thermofun.json")
            )
    )
    cs = ChemicalSystem(
        [sp[x] for x in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    i_Cal = findfirst(x -> symbol(x) == "Cal", cs.species)
    function calcite_water()
        st = ChemicalState(cs)
        set_quantity!(st, "Cal", 1.0e-2u"mol")
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "H+", 1.0e-7u"mol")
        set_quantity!(st, "OH-", 1.0e-7u"mol")
        return st
    end

    rxn2 = Reaction(
        OrderedDict(cs["Cal"] => 1.0),
        OrderedDict(cs["Ca+2"] => 1.0, cs["CO3-2"] => 1.0); symbol = "calcite",
    )
    k = 1.0e-5
    stoich = Float64.(
        KineticReaction(
            cs, rxn2,
            KineticFunc((T, P, t, n, lna, n0) -> 0.0, NamedTuple(), u"mol/s"),
        ).stoich,
    )
    R = 8.31446261815324
    GT = [
        ustrip(us"J/mol", s[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true)) /
            (R * 298.15) for s in cs.species
    ]
    kr = KineticReaction(
        cs, rxn2,
        KineticFunc(
            (T, P, t, n, lna, n0) -> begin
                Ω = saturation_ratio(stoich, [lna[symbol(s)] for s in cs.species], GT)
                k * (1 - Ω)
            end, NamedTuple(), u"mol/s",
        ),
    )
    kss2 = KineticStepSolver(cs, DiluteSolutionModel(), [kr])
    cal(s) = ustrip(us"mol", s.n[i_Cal])

    # Reference: 256 equal steps to t = 400 s.
    TEND = 400.0
    st = calcite_water()
    for _ in 1:256
        st = kinetic_step(kss2, st, TEND / 256)
    end
    ref = cal(st)

    # One step of the whole interval, for comparison: first order, so it is out
    # by about 1.2e-6 mol.
    one_step = cal(kinetic_step(kss2, calcite_water(), TEND))
    @test abs(one_step - ref) > 1.0e-7

    function march(reltol)
        s = calcite_water()
        t, h, taken = 0.0, TEND, Float64[]
        while t < TEND - 1.0e-9 && length(taken) < 200
            s, used, h = kinetic_step_adaptive(
                kss2, s, min(h, TEND - t); t = t, reltol = reltol,
            )
            t += used
            push!(taken, used)
        end
        return (cal(s), taken)
    end

    v3, steps3 = march(1.0e-3)
    v4, steps4 = march(1.0e-4)

    # The march completes — it does not stall. That is the property the first
    # version of this controller lacked: with the tolerance taken relative to
    # `Δξ`, which vanishes with the step, the estimated error grew as the step
    # shrank and it halved to the floor without advancing.
    @test sum(steps3) ≈ TEND rtol = 1.0e-9
    @test sum(steps4) ≈ TEND rtol = 1.0e-9
    @test length(steps3) < 200 && length(steps4) < 200

    # A tighter tolerance takes more steps and lands closer.
    @test length(steps4) > length(steps3)
    @test abs(v4 - ref) < abs(one_step - ref)
    @test abs(v4 - ref) < 1.0e-10        # measured at 2.5e-13

    # And the step size is chosen, not fixed: it varies over two orders of
    # magnitude within one march.
    @test maximum(steps4) / minimum(steps4) > 10

    # An estimate is available to the caller.
    e = Ref(0.0)
    kinetic_step_adaptive(
        kss2, calcite_water(), TEND; reltol = 1.0e-4, error_estimate = e,
    )
    @test e[] > 0
    @test isfinite(e[])

end

@testsection "the adaptive step survives a case the stiff ODE route gets wrong" begin

    # Measured, and the reason the tutorial recommends the implicit route for a
    # rate law that reads the solution. On calcite dissolving under `r = k(1 − Ω)`
    # over 1e5 s with k = 1e-4:
    #
    #   Rodas5P and the default polyalgorithm : extent = −457 mol, retcode Success
    #   Tsit5, explicit                       : correct, 85 626 steps, 519 s
    #   one implicit step of 1e5 s            : wrong, and its certificate says so
    #   kinetic_step_adaptive                 : correct to eight digits, 7 steps
    #
    # A stiff method needs a Jacobian, the residual carries a re-speciation the
    # Jacobian does not see, and the error control built on it reports success on
    # nonsense. Only the last two rows are asserted here; the ODE route's failure
    # is guarded by a warning in the extension rather than pinned by a test, since
    # the answer it returns is not a number worth encoding.
    sp = Dict(
        symbol(x) => x for x in build_species(
                datapath("slop98-inorganic-thermofun.json")
            )
    )
    cs = ChemicalSystem(
        [sp[x] for x in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    i_Cal = findfirst(x -> symbol(x) == "Cal", cs.species)
    i_Ca = findfirst(x -> symbol(x) == "Ca+2", cs.species)
    rxn3 = Reaction(
        OrderedDict(cs["Cal"] => 1.0),
        OrderedDict(cs["Ca+2"] => 1.0, cs["CO3-2"] => 1.0); symbol = "calcite",
    )
    k = 1.0e-4
    stoich = Float64.(
        KineticReaction(
            cs, rxn3,
            KineticFunc((T, P, t, n, lna, n0) -> 0.0, NamedTuple(), u"mol/s"),
        ).stoich,
    )
    R = 8.31446261815324
    GT = [
        ustrip(us"J/mol", x[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true)) /
            (R * 298.15) for x in cs.species
    ]
    kr = KineticReaction(
        cs, rxn3,
        KineticFunc(
            (T, P, t, n, lna, n0) -> begin
                Ω = saturation_ratio(stoich, [lna[symbol(x)] for x in cs.species], GT)
                k * (1 - Ω)
            end, NamedTuple(), u"mol/s",
        ),
    )
    kss3 = KineticStepSolver(cs, DiluteSolutionModel(), [kr])
    des3 = DualEquilibriumSolver(cs, DiluteSolutionModel())

    function big()
        st = ChemicalState(cs)
        set_quantity!(st, "Cal", 5.0e-2u"mol")
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "H+", 1.0e-7u"mol")
        set_quantity!(st, "OH-", 1.0e-7u"mol")
        return st
    end

    eq = SciMLBase.solve(des3, big())
    Ca_eq = ustrip(us"mol", eq.n[i_Ca])
    Cal_eq = ustrip(us"mol", eq.n[i_Cal])

    # One step of the whole interval is wrong — and reports it.
    cert = Ref{Any}(nothing)
    one = kinetic_step(kss3, big(), 1.0e5; certificate = cert)
    @test !cert[].optimal
    @test ustrip(us"mol", one.n[i_Cal]) < 0.5 * Cal_eq      # it dissolved everything

    # The adaptive march gets there, choosing its own steps.
    TEND = 1.0e5
    st = big()
    t, h, taken = 0.0, TEND, Float64[]
    while t < TEND - 1.0e-6 && length(taken) < 300
        st, used, h = kinetic_step_adaptive(kss3, st, min(h, TEND - t); t = t, reltol = 1.0e-4)
        t += used
        push!(taken, used)
    end
    @test sum(taken) ≈ TEND rtol = 1.0e-9
    @test length(taken) < 30                                # measured at 7
    @test ustrip(us"mol", st.n[i_Ca]) ≈ Ca_eq rtol = 1.0e-6
    @test ustrip(us"mol", st.n[i_Cal]) ≈ Cal_eq rtol = 1.0e-8

end

@testsection "coupling = :species imposes the assemblage" begin

    # Two C3A pathways with distinct solid products, the form Lavergne et al.
    # (2018) uses. With one constraint per REACTION the minimization rearranges
    # the solids and the monosulfoaluminate ends at zero; with one per kinetic
    # SPECIES every product follows the coefficients written down.
    #
    # The bug that made this impossible is worth naming: a pinning row for a
    # product that starts absent has one positive entry and a zero budget, which
    # `degenerate_components` read as "this component is absent from the system".
    # It pinned that row's multiplier and declared the species dead — stationarity
    # residual 458, extents 4.3 times short. `conservation_rows` now restricts the
    # criterion to the rows that conserve something.
    CEM = Dict(
        symbol(x) => x for x in build_species(
                datapath("cemdata18-thermofun.json")
            )
    )
    spc = speciation(
        collect(values(CEM)),
        split("C3A Gp H2O@ ettringite monosulphate12 Portlandite");
        aggregate_state = [AS_AQUEOUS],
    )
    cs4 = ChemicalSystem(spc, CEMDATA_PRIMARIES)
    S(x) = cs4[x]
    r_aft = Reaction(
        OrderedDict(S("C3A") => 1.0, S("Gp") => 3.0, S("H2O@") => 26.0),
        OrderedDict(S("ettringite") => 1.0); symbol = "C3A -> AFt",
    )
    r_afm = Reaction(
        OrderedDict(S("C3A") => 1.0, S("Gp") => 1.0, S("H2O@") => 10.0),
        OrderedDict(S("monosulphate12") => 1.0); symbol = "C3A -> AFm",
    )
    f(k) = KineticFunc((T, P, t, n, lna, n0) -> k, NamedTuple(), u"mol/s")
    k1, k2 = 2.0e-6, 5.0e-7
    krs = [KineticReaction(cs4, r_aft, f(k1)), KineticReaction(cs4, r_afm, f(k2))]

    names4 = symbol.(cs4.species)
    idx4 = Dict(
        nm => findfirst(==(nm), names4) for nm in
            ("C3A", "Gp", "ettringite", "monosulphate12")
    )
    function fresh4()
        st = ChemicalState(cs4)
        set_quantity!(st, "C3A", 1.0e-2u"mol")
        set_quantity!(st, "Gp", 3.0e-2u"mol")
        set_quantity!(st, "H2O@", 1.0u"kg")
        return st
    end

    Δt = 1000.0
    kss_sp = KineticStepSolver(cs4, DiluteSolutionModel(), krs; coupling = :species)
    q = Ref(Float64[])
    cert_sp = Ref{Any}(nothing)
    s_sp = kinetic_step(kss_sp, fresh4(), Δt; parameters = q, certificate = cert_sp)
    n_sp = Float64[ustrip(us"mol", x) for x in s_sp.n]

    # The extents ARE the reaction progress: no `M` factor, because pinning the
    # species makes `Δξ` the progress itself.
    @test q[] ≈ Δt .* [k1, k2] rtol = 1.0e-6

    # And every species follows the stoichiometry, the products included.
    @test n_sp[idx4["C3A"]] ≈ 1.0e-2 - (k1 + k2) * Δt rtol = 1.0e-6
    @test n_sp[idx4["Gp"]] ≈ 3.0e-2 - (3k1 + k2) * Δt rtol = 1.0e-6
    @test n_sp[idx4["ettringite"]] ≈ k1 * Δt rtol = 1.0e-6
    @test n_sp[idx4["monosulphate12"]] ≈ k2 * Δt rtol = 1.0e-6

    # The element balance, and it is certified.
    #
    # Held by a linear row instead — a Lagrange multiplier that must reach the
    # mineral's own chemical potential, of order 10²-10³ in RT units — this same
    # step stalled at 6.1e-7 mol whatever `maxit`, `tol` or the number of
    # active-set updates. Eliminating the pinned species instead, so their amounts
    # are computed from the extents and the budget they leave is handed to a plain
    # equilibrium over what remains, brings it to 2e-14.
    A4 = Matrix{Float64}(cs4.SM.A)
    b4 = A4 * Float64[ustrip(us"mol", x) for x in fresh4().n]
    resid = maximum(abs, A4 * n_sp - b4)
    @test resid < 1.0e-11

    @test cert_sp[].optimal
    @test cert_sp[].balance < 1.0e-11

    # Both routes report a certificate of the SAME shape, whichever produced it.
    # `:reactions` fills it from `kkt_certificate` and `:species` from
    # `optimality_certificate`, and those name the balance differently.
    @test hasproperty(cert_sp[], :balance) && hasproperty(cert_sp[], :stationarity)

    # `:reactions` on the same problem: certified, and the assemblage is
    # thermodynamic instead — the monosulfoaluminate does not form.
    kss_rx = KineticStepSolver(cs4, DiluteSolutionModel(), krs)
    cert = Ref{Any}(nothing)
    s_rx = kinetic_step(kss_rx, fresh4(), Δt; certificate = cert)
    n_rx = Float64[ustrip(us"mol", x) for x in s_rx.n]
    @test cert[].optimal
    @test maximum(abs, A4 * n_rx - b4) < 1.0e-9
    @test n_rx[idx4["monosulphate12"]] < 1.0e-9
    @test n_rx[idx4["C3A"]] < n_sp[idx4["C3A"]]     # the two disagree, as they must

end
