@testsection "ChemicalState" begin
    # ── Shared fixtures ──────────────────────────────────────────────────────
    h2o = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    hplus = Species("H+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    oh = Species("OH-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    nacl = Species("NaCl"; aggregate_state = AS_CRYSTAL)

    cs_basic = ChemicalSystem([h2o, hplus, oh])

    @testsection "Construction — default T and P" begin
        state = ChemicalState(cs_basic)
        @test ustrip(temperature(state)) ≈ 298.15
        @test isapprox(ustrip(pressure(state)), 1.0e5; rtol = 1.0e-4)   # 1 bar in Pa
        @test length(state.n) == 3
        @test all(iszero.(ustrip.(state.n)))
    end

    @testsection "Construction with explicit n (moles)" begin
        state = ChemicalState(cs_basic, [55.5u"mol", 1.0e-7u"mol", 1.0e-7u"mol"])
        @test ustrip(moles(state, "H2O")) ≈ 55.5
        @test ustrip(moles(state, "H+")) ≈ 1.0e-7
    end

    @testsection "Construction with mass units" begin
        cs = ChemicalSystem([h2o, nacl])
        # 5.844 g NaCl ≈ 0.1 mol
        state = ChemicalState(cs, [55.5u"mol", 5.844u"g"])
        @test ustrip(moles(state, "H2O")) ≈ 55.5
        @test isapprox(ustrip(moles(state, "NaCl")), 0.1; rtol = 1.0e-3)
    end

    @testsection "Temperature and pressure accessors" begin
        state = ChemicalState(cs_basic; T = 350.0u"K", P = 2u"bar")
        @test ustrip(temperature(state)) ≈ 350.0
        @test isapprox(ustrip(pressure(state)), 2.0e5; rtol = 1.0e-4)
    end

    @testsection "set_temperature! and set_pressure!" begin
        state = ChemicalState(cs_basic)
        set_temperature!(state, 400.0u"K")
        @test ustrip(temperature(state)) ≈ 400.0

        set_pressure!(state, 5u"bar")
        @test isapprox(ustrip(pressure(state)), 5.0e5; rtol = 1.0e-4)
    end

    @testsection "set_quantity! by species" begin
        state = ChemicalState(cs_basic)
        set_quantity!(state, h2o, 10.0u"mol")
        @test ustrip(moles(state, "H2O")) ≈ 10.0
    end

    @testsection "set_quantity! by symbol" begin
        state = ChemicalState(cs_basic, [55.5u"mol", 0.0u"mol", 0.0u"mol"])
        set_quantity!(state, "H+", 1.0e-4u"mol")
        @test ustrip(moles(state, "H+")) ≈ 1.0e-4
    end

    @testsection "moles by phase (NamedTuple)" begin
        state = ChemicalState(cs_basic, [55.5u"mol", 1.0e-7u"mol", 1.0e-7u"mol"])
        n_phases = moles(state)
        @test hasproperty(n_phases, :liquid)
        @test hasproperty(n_phases, :total)
        @test ustrip(n_phases.liquid) ≈ 55.5 + 2.0e-7
        @test ustrip(n_phases.solid) ≈ 0.0
    end

    @testsection "mass accessors" begin
        state = ChemicalState(cs_basic, [1.0u"mol", 0.0u"mol", 0.0u"mol"])
        m = mass(state)
        @test hasproperty(m, :total)
        # 1 mol H2O ≈ 18.015 g
        @test isapprox(ustrip(uconvert(us"g", m.total)), 18.015; rtol = 1.0e-3)

        # Mass of a specific species
        m_h2o = mass(state, "H2O")
        @test isapprox(ustrip(uconvert(us"g", m_h2o)), 18.015; rtol = 1.0e-3)
    end

    @testsection "volume accessors (without V⁰)" begin
        # Without molar volumes, volume is a NamedTuple with zero values
        state = ChemicalState(cs_basic, [55.5u"mol", 1.0e-7u"mol", 1.0e-7u"mol"])
        v = volume(state)
        @test hasproperty(v, :total)
        # volume(state, species) returns nothing when V⁰ is absent
        @test isnothing(volume(state, h2o))
    end

    @testsection "pH and pOH — acidic system" begin
        # 1e-4 mol H+ in ~1L water → pH ≈ 4
        state = ChemicalState(cs_basic, [55.5u"mol", 1.0e-4u"mol", 0.0u"mol"])
        ph = pH(state)
        if !isnothing(ph)
            @test isapprox(ph, 4.0; atol = 0.1)
        end
    end

    @testsection "pH returns nothing without H+ or OH-" begin
        cs_no_ions = ChemicalSystem([h2o])
        state = ChemicalState(cs_no_ions, [55.5u"mol"])
        @test isnothing(pH(state))
        @test isnothing(pOH(state))
    end

    @testsection "porosity and saturation return NaN without V⁰" begin
        state = ChemicalState(cs_basic, [55.5u"mol", 0.0u"mol", 0.0u"mol"])
        # Without V⁰ all volumes are zero → porosity = NaN
        @test isnan(porosity(state))
        @test isnan(saturation(state))
    end

    @testsection "copy — independent mutable state" begin
        state = ChemicalState(cs_basic, [55.5u"mol", 1.0e-7u"mol", 1.0e-7u"mol"])
        clone = copy(state)

        set_quantity!(clone, "H2O", 10.0u"mol")
        @test ustrip(moles(state, "H2O")) ≈ 55.5   # original unchanged
        @test ustrip(moles(clone, "H2O")) ≈ 10.0   # clone updated

        # System is shared
        @test clone.system === state.system
    end

    @testsection "_update_derived! recalculates after mutation" begin
        state = ChemicalState(cs_basic, [55.5u"mol", 1.0e-7u"mol", 1.0e-7u"mol"])
        n_before = moles(state).liquid

        set_quantity!(state, "H2O", 100.0u"mol")
        n_after = moles(state).liquid

        @test ustrip(n_after) > ustrip(n_before)
    end

    @testsection "scaling: Base.* and Base./" begin
        state = ChemicalState(cs_basic, [2.0u"mol", 1.0u"mol", 0.0u"mol"])

        # state * α returns a new state
        s2 = state * 3.0
        @test ustrip(moles(s2, "H2O")) ≈ 6.0
        @test ustrip(moles(s2, "H+")) ≈ 3.0
        @test ustrip(moles(state, "H2O")) ≈ 2.0   # original unchanged

        # α * state
        s3 = 2.0 * state
        @test ustrip(moles(s3, "H2O")) ≈ 4.0

        # state / α
        s_half = state / 2.0
        @test ustrip(moles(s_half, "H2O")) ≈ 1.0

        # divide by zero raises
        @test_throws ErrorException state / 0.0

        # derived quantities are recomputed (total moles consistent)
        @test ustrip(moles(s2).total) ≈ ustrip(moles(state).total) * 3.0
    end

    @testsection "scaling: rescale! (mol)" begin
        state = ChemicalState(cs_basic, [2.0u"mol", 1.0u"mol", 0.0u"mol"])
        rescale!(state, 6.0u"mol")
        @test ustrip(us"mol", moles(state).total) ≈ 6.0
        # Ratios are preserved
        @test ustrip(moles(state, "H2O")) ≈ 4.0
        @test ustrip(moles(state, "H+")) ≈ 2.0
    end

    @testsection "scaling: rescale! (kg)" begin
        state = ChemicalState(cs_basic, [1.0u"mol", 0.0u"mol", 0.0u"mol"])
        m0 = ustrip(us"kg", mass(state).total)   # ≈ 0.018015 kg
        rescale!(state, 1.0u"kg")
        @test isapprox(ustrip(us"kg", mass(state).total), 1.0; rtol = 1.0e-6)
        # Factor is the inverse of m0
        @test isapprox(ustrip(moles(state, "H2O")), 1.0 / m0; rtol = 1.0e-6)
    end

    @testsection "scaling: rescale! (volume)" begin
        h2o_v = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
        # V⁰ must be a callable with signature (T, P; unit) — use NumericFunc
        h2o_v[:V⁰] = NumericFunc((T, P) -> 18.07e-6, 1u"m^3/mol")
        cs_v = ChemicalSystem([h2o_v])
        state_v = ChemicalState(cs_v, [1.0u"mol"])
        # volume of 1 mol H2O ≈ 18.07e-6 m³; rescale to 1 m³
        rescale!(state_v, 1.0u"m^3")
        @test isapprox(ustrip(us"m^3", volume(state_v).total), 1.0; rtol = 1.0e-5)
    end

    @testsection "scaling: rescale! error cases" begin
        state = ChemicalState(cs_basic, [2.0u"mol", 0.0u"mol", 0.0u"mol"])

        # Bad dimension
        @test_throws ErrorException rescale!(state, 1.0u"s")

        # Zero state
        s_zero = ChemicalState(cs_basic)   # all n = 0
        @test_throws ErrorException rescale!(s_zero, 1.0u"mol")
    end

    @testsection "ForwardDiff — _build_params / _build_n0" begin
        using ForwardDiff

        # Build minimal species with synthetic ΔₐG⁰ functions
        h2o_ad = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
        h2o_ad[:ΔₐG⁰] = SymbolicFunc(
            :(G0 + a * T); G0 = -237000.0, a = -70.0,
            units = [:T => "K", :G0 => "J/mol", :a => "J/(mol*K)"]
        )
        h2o_ad[:M] = 0.018u"kg/mol"
        h2o_ad[:V⁰] = SymbolicFunc(1.8e-5u"m^3/mol")

        hp_ad = Species("H+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
        hp_ad[:ΔₐG⁰] = SymbolicFunc(
            :(G0 + b * T); G0 = 0.0, b = 0.0,
            units = [:T => "K", :G0 => "J/mol", :b => "J/(mol*K)"]
        )
        hp_ad[:M] = 0.001u"kg/mol"

        cs_ad = ChemicalSystem([h2o_ad, hp_ad])
        state_ad = ChemicalState(cs_ad; T = 298.15u"K", P = 1.0e5u"Pa")
        set_quantity!(state_ad, "H2O", 55.5u"mol")
        set_quantity!(state_ad, "H+", 1.0e-7u"mol")

        # _build_params must return finite dimensionless values
        p = ChemistryLab._build_params(state_ad)
        @test all(isfinite, p.ΔₐG⁰overT)
        @test isfinite(p.T)
        @test isfinite(p.P)

        # _build_n0 must return the dimensionless mole vector
        n0 = ChemistryLab._build_n0(state_ad)
        @test length(n0) == 2
        @test n0[1] ≈ 55.5
        @test n0[2] ≈ 1.0e-7
    end
end

@testset "enthalpy, heat capacity, and the lazy thermodynamic functions" begin

    # `s[:Cp⁰]` used to return the not-found value 0 on any species whose
    # thermodynamic functions had not been forced yet, because `getindex` did not
    # do the lazy completion that `getproperty` does. It returned an `Int64`, so
    # the caller found out only when it tried to evaluate it.
    sp = build_species(datapath("cemdata18-thermofun.json"))
    qtz = first(s for s in sp if symbol(s) == "Qtz")
    @test !(qtz[:Cp⁰] isa Integer)          # i.e. not the 0 fallback
    @test qtz[:Cp⁰](T = 293.15, unit = true) isa AbstractQuantity
    # quartz near 25 °C, ≈ 44 J/(mol·K)
    @test 40 < ustrip(us"J/K/mol", qtz[:Cp⁰](T = 293.15, unit = true)) < 48

    dict = Dict(symbol(s) => s for s in sp)
    cs = ChemicalSystem(
        [dict[s] for s in ("H2O@", "H+", "OH-", "Ca+2", "Portlandite")],
        ["H2O@", "H+", "Ca+2", "Zz"],
    )
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 55.5u"mol")
    set_quantity!(st, "Portlandite", 1.0u"mol")

    @test isempty(missing_enthalpy(st))

    # Additive by construction, and the pieces are the tabulated values.
    H = ustrip(us"J", enthalpy(st))
    H_w = ustrip(us"J", enthalpy(st, "H2O@"))
    H_p = ustrip(us"J", enthalpy(st, "Portlandite"))
    # Additive over ALL species: the state also carries trace H⁺, OH⁻ and Ca²⁺,
    # so the sum is not just the two set by hand.
    @test H ≈ sum(ustrip(us"J", enthalpy(st, nm)) for nm in symbol.(cs.species)) rtol = 1.0e-12
    @test H ≈ H_w + H_p rtol = 1.0e-6
    @test H_w / 55.5 ≈ -285_830 rtol = 5.0e-3      # liquid water
    @test H_p ≈ -984_675 rtol = 5.0e-3             # portlandite

    C = ustrip(us"J/K", heat_capacity(st))
    @test C ≈ sum(ustrip(us"J/K", heat_capacity(st, nm)) for nm in symbol.(cs.species)) rtol = 1.0e-12
    @test 4100 < C < 4400                          # 55.5 × 75.3 + 87.5

    # Doubling the amounts doubles both: they are extensive.
    st2 = copy(st)
    st2.n .*= 2
    @test ustrip(us"J", enthalpy(st2)) ≈ 2H rtol = 1.0e-12
    @test ustrip(us"J/K", heat_capacity(st2)) ≈ 2C rtol = 1.0e-12

    # A species with no enthalpy of formation is reported rather than silently
    # dropped from the balance.
    bare = ChemicalSystem([Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)])
    @test "H2O" in missing_enthalpy(ChemicalState(bare, [55.5u"mol"]))

end
