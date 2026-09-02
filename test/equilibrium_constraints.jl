# Equilibrium under something other than fixed (T, P).
#
# The prescribed property is an unknown of the solver's own system, so these tests
# check three things: that the property really is satisfied at the answer, that
# the temperature or pressure found is the one an ordinary fixed-(T, P) solve
# reproduces, and that the answer agrees with a published number where one exists.

@testsection "Equilibrium constraints" begin

    sp = Dict(symbol(s) => s for s in build_species(
        datapath("slop98-inorganic-thermofun.json")))

    @testsection "adiabatic neutralization gives the published enthalpy" begin
        # H⁺ + OH⁻ → H₂O releases 55.8 kJ/mol. Nothing here is fitted to that: the
        # enthalpies of formation come from the database and the temperature is an
        # unknown of the system.
        cs = ChemicalSystem([sp[s] for s in split("H2O@ H+ OH-")], ["H2O@", "H+", "Zz"])
        des = DualEquilibriumSolver(cs, DiluteSolutionModel())

        function mix(n)
            st = ChemicalState(cs)
            set_quantity!(st, "H2O@", 1.0u"kg")
            set_quantity!(st, "H+", n * u"mol")
            set_quantity!(st, "OH-", n * u"mol")
            set_temperature!(st, 298.15u"K")
            return st
        end

        ΔT = Float64[]
        for n in (0.05, 0.10, 0.25, 0.50)
            st0 = mix(n)
            n0 = Float64[ustrip(us"mol", x) for x in st0.n]
            H0 = ChemistryLab._total_enthalpy(cs, n0, 298.15, 1.0e5)

            eq = SciMLBase.solve(des, mix(n); constraint = Adiabatic())
            Tf = ustrip(us"K", temperature(eq))
            nf = Float64[ustrip(us"mol", x) for x in eq.n]

            # The enthalpy is conserved — that is what "adiabatic" means, and it
            # holds to machine precision because it is an equation of the system.
            Hf = ChemistryLab._total_enthalpy(cs, nf, Tf, 1.0e5)
            @test abs((Hf - H0) / H0) < 1.0e-12

            # The heat released, measured at the initial temperature so it is the
            # enthalpy of reaction and not a heat-capacity effect.
            H_iso = ChemistryLab._total_enthalpy(cs, nf, 298.15, 1.0e5)
            @test (H_iso - H0) / n / 1000 ≈ -55.85 atol = 0.1

            push!(ΔT, Tf - 298.15)
        end

        # Temperature rises, and close to linearly in the amount neutralized:
        # 0.5 mol in 1 kg of water gives 6.6 K, against 27.9 kJ over 4.18 kJ/K.
        @test all(>(0), ΔT)
        @test ΔT[end] ≈ 6.62 atol = 0.1
        @test ΔT[end] / ΔT[1] ≈ 10 rtol = 0.05      # slightly sublinear
    end

    @testsection "the adiabatic answer is the fixed-T answer at the T it found" begin
        cs = ChemicalSystem(
            [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
            ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
        )
        des = DualEquilibriumSolver(cs, DiluteSolutionModel())
        function acidified()
            st = ChemicalState(cs)
            set_quantity!(st, "Cal", 1.0e-2u"mol")
            set_quantity!(st, "H2O@", 1.0u"kg")
            set_quantity!(st, "H+", 1.0e-2u"mol")
            set_quantity!(st, "OH-", 1.0e-10u"mol")
            set_temperature!(st, 298.15u"K")
            return st
        end

        eq_ad = SciMLBase.solve(des, acidified(); constraint = Adiabatic())
        Tf = temperature(eq_ad)

        st2 = acidified()
        set_temperature!(st2, Tf)
        eq_fixed = SciMLBase.solve(des, st2)

        n_ad = [ustrip(us"mol", x) for x in eq_ad.n]
        n_fx = [ustrip(us"mol", x) for x in eq_fixed.n]
        @test maximum(abs, n_ad .- n_fx) < 1.0e-12
        @test pH(eq_ad) ≈ pH(eq_fixed) atol = 1.0e-6
    end

    @testsection "FixedEnthalpy and Adiabatic agree" begin
        cs = ChemicalSystem([sp[s] for s in split("H2O@ H+ OH-")], ["H2O@", "H+", "Zz"])
        des = DualEquilibriumSolver(cs, DiluteSolutionModel())
        st = ChemicalState(cs)
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "H+", 0.1u"mol")
        set_quantity!(st, "OH-", 0.1u"mol")
        n0 = Float64[ustrip(us"mol", x) for x in st.n]
        H0 = ChemistryLab._total_enthalpy(cs, n0, 298.15, 1.0e5)

        a = SciMLBase.solve(des, deepcopy(st); constraint = Adiabatic())
        b = SciMLBase.solve(des, deepcopy(st); constraint = FixedEnthalpy(H0 * u"J"))
        @test temperature(a) ≈ temperature(b) rtol = 1.0e-8
    end

    @testsection "a volume constraint on a condensed system is refused, with the lever" begin
        # Not a solver limitation: the molar volumes of water and of the minerals
        # in these databases are exactly pressure-independent, so the volume of a
        # condensed system is fixed by its composition and the pressure has no
        # purchase on it. Diverging silently would be the wrong answer.
        cs = ChemicalSystem([sp[s] for s in split("H2O@ H+ OH-")], ["H2O@", "H+", "Zz"])
        des = DualEquilibriumSolver(cs, DiluteSolutionModel())
        st = ChemicalState(cs)
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "H+", 0.25u"mol")
        set_quantity!(st, "OH-", 0.25u"mol")

        @test_throws ArgumentError SciMLBase.solve(
            des, deepcopy(st); constraint = SealedVolume(),
        )

        n0 = Float64[ustrip(us"mol", x) for x in st.n]
        lever = ChemistryLab._pressure_lever(cs, n0, 298.15, 1.0e5)
        @test lever < 1.0e-4                     # measured at 1.0e-6
        @test ChemistryLab._total_volume(cs, n0, 298.15, 1.0e5) ≈
            ChemistryLab._total_volume(cs, n0, 298.15, 1.0e7) rtol = 1.0e-3
    end

    @testsection "a constraint needs the dual route, and says so" begin
        solid_only = ChemicalSystem(
            [Species("NaCl"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)],
            ["NaCl"],
        )
        st = ChemicalState(solid_only, [0.05u"mol"])
        @test_throws ArgumentError equilibrate_certified(st; constraint = Adiabatic())
    end

end

@testsection "Prescribed pH by implicit titrant" begin

    # The second vehicle: a prescribed chemical potential is a COLUMN added to the
    # conservation matrix — an unknown amount of titrant the system may draw on —
    # not a parameter. `Aq` carries it, so the linear rows read
    # `A x − A[:, titrant] q = b`, and the answer includes how much titrant it
    # took, which is what a titration measures.
    sp = Dict(symbol(s) => s for s in build_species(
        datapath("slop98-inorganic-thermofun.json")))
    cs = ChemicalSystem(
        [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    des = DualEquilibriumSolver(cs, DiluteSolutionModel())
    i_H = findfirst(s -> symbol(s) == "H+", cs.species)
    i_Cal = findfirst(s -> symbol(s) == "Cal", cs.species)

    function calcite_water()
        st = ChemicalState(cs)
        set_quantity!(st, "Cal", 1.0e-2u"mol")
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "H+", 1.0e-7u"mol")
        set_quantity!(st, "OH-", 1.0e-7u"mol")
        return st
    end

    free = SciMLBase.solve(des, calcite_water())
    @test pH(free) > 9                      # calcite in water is basic

    amounts = Float64[]
    dissolved = Float64[]
    for target in (6.0, 7.0, 8.0, 9.0)
        q = Ref(Float64[])
        eq = SciMLBase.solve(
            des, calcite_water(); constraint = FixedpH(target), parameters = q,
        )
        n = Float64[ustrip(us"mol", x) for x in eq.n]

        # The prescribed quantity is the ACTIVITY, and it is held exactly.
        p2 = ChemistryLab._build_params(eq)
        lnaH = des.lna(n, p2)[i_H]
        @test -lnaH / log(10) ≈ target atol = 1.0e-4

        # `pH(state)` is concentration-based, so it differs by the dilute model's
        # ρ ≈ 1 kg/L approximation — about 0.0013, and documented on `FixedpH`.
        @test abs(pH(eq) - target) < 0.01

        push!(amounts, q[][1])
        push!(dissolved, 1.0e-2 - n[i_Cal])
    end

    # Acid has to be ADDED to bring a basic solution down, and more of it the
    # lower the target. Both are physics, not tolerances.
    @test all(>(0), amounts)
    @test issorted(amounts; rev = true)
    @test issorted(dissolved; rev = true)
    # At pH 6 the calcite is gone entirely.
    @test dissolved[1] ≈ 1.0e-2 rtol = 1.0e-6

    @testsection "a titrant that is not in the system is refused" begin
        @test_throws ArgumentError SciMLBase.solve(
            des, calcite_water(); constraint = FixedpH(7.0; titrant = "Na+"),
        )
    end

    @testsection "FixedActivity is the general form" begin
        q = Ref(Float64[])
        eq = SciMLBase.solve(
            des, calcite_water();
            constraint = FixedActivity("H+", 10.0^-7), parameters = q,
        )
        n = Float64[ustrip(us"mol", x) for x in eq.n]
        lnaH = des.lna(n, ChemistryLab._build_params(eq))[i_H]
        @test -lnaH / log(10) ≈ 7.0 atol = 1.0e-4
        @test q[][1] > 0
    end

end
