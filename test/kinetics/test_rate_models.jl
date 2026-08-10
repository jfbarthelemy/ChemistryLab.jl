using ForwardDiff

# ── arrhenius_rate_constant ───────────────────────────────────────────────────

@testset "arrhenius_rate_constant" begin

    k = arrhenius_rate_constant(5.0e-4, 62000.0)

    @test k isa NumericFunc

    # At T_ref = 298.15, k(T_ref) = k₀
    @test isapprox(k(; T = 298.15), 5.0e-4; rtol = 1.0e-10)

    # Arrhenius: higher T → higher k
    @test k(; T = 350.0) > k(; T = 298.15)
    @test k(; T = 250.0) < k(; T = 298.15)

    # AD smoke-test: dk/dT should be finite and positive
    dkdT = ForwardDiff.derivative(T -> arrhenius_rate_constant(5.0e-4, 62000.0)(; T = T), 298.15)
    @test isfinite(dkdT)
    @test dkdT > 0

    # AD w.r.t. Ea: at T > T_ref, larger Ea increases k
    dk_dEa = ForwardDiff.derivative(Ea -> arrhenius_rate_constant(5.0e-4, Ea)(; T = 350.0), 62000.0)
    @test isfinite(dk_dEa)
    @test dk_dEa > 0

    # AD w.r.t. k₀
    dk_dk0 = ForwardDiff.derivative(k₀ -> arrhenius_rate_constant(k₀, 62000.0)(; T = 298.15), 5.0e-4)
    @test isapprox(dk_dk0, 1.0; rtol = 1.0e-10)

end

# ── KINETICS_RATE_FACTORIES ───────────────────────────────────────────────────

@testset "KINETICS_RATE_FACTORIES" begin

    @test haskey(KINETICS_RATE_FACTORIES, :arrhenius)

    factory = KINETICS_RATE_FACTORIES[:arrhenius]
    k = factory(; k₀ = 1.0e-6, Ea = 40000.0, T_ref = 298.15, R_gas = 8.31446)
    @test k isa AbstractFunc
    @test isapprox(k(; T = 298.15), 1.0e-6; rtol = 1.0e-10)

end

# ── saturation_ratio ──────────────────────────────────────────────────────────

@testset "saturation_ratio" begin

    stoich = [-1.0, 1.0]
    lna = [0.0, 0.0]
    G_over_T = [0.0, 0.0]
    Ω = saturation_ratio(stoich, lna, G_over_T)
    @test isapprox(Ω, 1.0; rtol = 1.0e-12)

    # Ω > 1: excess A → forward reaction is favored
    lna2 = [-1.0, 0.0]
    Ω2 = saturation_ratio(stoich, lna2, G_over_T)
    @test Ω2 > 1.0

    # AD smoke-test
    dΩ = ForwardDiff.derivative(x -> saturation_ratio(stoich, [x, 0.0], G_over_T), 0.0)
    @test isfinite(dΩ)

end

# ── StateView ─────────────────────────────────────────────────────────────────

@testset "StateView" begin

    data = [1.0, 2.0, 3.0]
    index = Dict("A" => 1, "B" => 2, "C" => 3)
    sv = StateView(data, index)

    @test sv isa StateView
    @test sv["A"] ≈ 1.0
    @test sv["B"] ≈ 2.0
    @test sv["C"] ≈ 3.0

    @test haskey(sv, "A")
    @test haskey(sv, "B")
    @test !haskey(sv, "D")

    # StateView is a thin wrapper — mutation in data is reflected
    data[2] = 99.0
    @test sv["B"] ≈ 99.0

    # AD smoke-test: ForwardDiff through StateView lookup
    function f_sv(x)
        d = [x, 2.0, 3.0]
        sv2 = StateView(d, index)
        return sv2["A"]^2
    end
    dfdx = ForwardDiff.derivative(f_sv, 5.0)
    @test isfinite(dfdx)
    @test isapprox(dfdx, 2 * 5.0; rtol = 1.0e-12)

end

# ── KineticFunc ───────────────────────────────────────────────────────────────

@testset "KineticFunc" begin

    # Simple closure: r = T * 1e-10
    f = (T, _P, _t, _n, _lna, _n0) -> T * 1.0e-10
    kf = KineticFunc(f, (T = 298.15u"K",), u"mol/s")

    @test kf isa KineticFunc

    # Calling convention: (T, P, t, n, lna, n0)
    index = Dict("X" => 1)
    sv = StateView([1.0], index)
    r = kf(298.15, 1.0e5, 0.0, sv, sv, sv)
    @test isapprox(r, 298.15 * 1.0e-10; rtol = 1.0e-12)

    # AD smoke-test w.r.t. T
    drdT = ForwardDiff.derivative(T -> kf(T, 1.0e5, 0.0, sv, sv, sv), 298.15)
    @test isfinite(drdT)
    @test isapprox(drdT, 1.0e-10; rtol = 1.0e-12)

end

# ── PK_PARAMS_* are NamedTuples ───────────────────────────────────────────────

@testset "PK_PARAMS are NamedTuples" begin

    for params in (PK_PARAMS_C3S, PK_PARAMS_C2S, PK_PARAMS_C3A, PK_PARAMS_C4AF)
        @test params isa NamedTuple
        @test haskey(params, :K₁)
        @test haskey(params, :N₁)
        @test haskey(params, :K₂)
        @test haskey(params, :N₂)
        @test haskey(params, :K₃)
        @test haskey(params, :N₃)
        @test haskey(params, :B)
        @test haskey(params, :Ea)
        @test haskey(params, :T_ref)
    end

end

# ── parrot_killoh ─────────────────────────────────────────────────────────────

@testset "parrot_killoh" begin

    # ── Construction ─────────────────────────────────────────────────────────
    pk = parrot_killoh(PK_PARAMS_C3S, "C3S")
    @test pk isa KineticFunc

    # ── Positivity and monotonicity ───────────────────────────────────────────
    # Build a minimal StateView with C3S moles
    n0 = 1.0      # mol
    T_K = 293.15  # K (T_ref → Aₜ = 1)

    index = Dict("C3S" => 1)
    n_sv = StateView([0.99 * n0], index)   # α ≈ 0.01
    n0_sv = StateView([n0], index)
    lna_sv = StateView([0.0], index)

    r_early = pk(T_K, 1.0e5, 0.0, n_sv, lna_sv, n0_sv)
    @test r_early > 0
    @test isfinite(r_early)

    # α = 0.50
    n_sv_mid = StateView([0.5 * n0], index)
    r_mid = pk(T_K, 1.0e5, 0.0, n_sv_mid, lna_sv, n0_sv)
    @test r_mid > 0
    @test isfinite(r_mid)

    # α ≈ 0.99
    n_sv_late = StateView([0.01 * n0], index)
    r_late = pk(T_K, 1.0e5, 0.0, n_sv_late, lna_sv, n0_sv)
    @test r_late > 0

    # ── Temperature dependence ────────────────────────────────────────────────
    r_low_T = pk(273.15, 1.0e5, 0.0, n_sv_mid, lna_sv, n0_sv)
    r_high_T = pk(323.15, 1.0e5, 0.0, n_sv_mid, lna_sv, n0_sv)
    @test r_high_T > r_low_T

    # ── α_max limit ───────────────────────────────────────────────────────────
    pk_max = parrot_killoh(PK_PARAMS_C3S, "C3S"; α_max = 0.85)
    # At n ≈ 0 (α ≈ 1.0 > α_max = 0.85), rate should be very small
    n_sv_zero = StateView([1.0e-15], index)
    r_at_zero = pk_max(T_K, 1.0e5, 0.0, n_sv_zero, lna_sv, n0_sv)
    @test r_at_zero < r_early

    # ── All 4 predefined clinker phases ──────────────────────────────────────
    for (params, name) in (
            (PK_PARAMS_C3S, "C3S"),
            (PK_PARAMS_C2S, "C2S"),
            (PK_PARAMS_C3A, "C3A"),
            (PK_PARAMS_C4AF, "C4AF"),
        )
        pk_k = parrot_killoh(params, name)
        idx = Dict(name => 1)
        n_k = StateView([0.5], idx)
        n0_k = StateView([1.0], idx)
        lna_k = StateView([0.0], idx)
        r_k = pk_k(293.15, 1.0e5, 0.0, n_k, lna_k, n0_k)
        @test r_k > 0
        @test isfinite(r_k)
    end

    # ── AD smoke-tests ────────────────────────────────────────────────────────
    # w.r.t. T
    dr_dT = ForwardDiff.derivative(
        T -> pk(T, 1.0e5, 0.0, n_sv_mid, lna_sv, n0_sv),
        293.15,
    )
    @test isfinite(dr_dT)
    @test dr_dT > 0   # higher T → faster rate

    # w.r.t. n_current (via StateView)
    function r_of_n(n_curr)
        sv = StateView([n_curr], index)
        return pk(T_K, 1.0e5, 0.0, sv, lna_sv, n0_sv)
    end
    dr_dn = ForwardDiff.derivative(r_of_n, 0.5)
    @test isfinite(dr_dn)
    # Decreasing n → decreasing α → changing rate (sign depends on regime)

end

# ── Surface area models ───────────────────────────────────────────────────────

@testset "SurfaceArea models" begin

    fixed = FixedSurfaceArea(0.5)
    @test surface_area(fixed, 0.01, 0.1) ≈ 0.5
    @test surface_area(fixed, 0.0, 0.1) ≈ 0.5

    bet = BETSurfaceArea(90.0)
    @test surface_area(bet, 0.01, 0.1) ≈ 90.0 * 0.01 * 0.1
    @test surface_area(bet, 0.0, 0.1) ≈ 0.0

    dA = ForwardDiff.derivative(n -> surface_area(bet, n, 0.1), 0.01)
    @test isfinite(dA)
    @test dA > 0

end

# ── KineticReaction constructors ──────────────────────────────────────────────

@testset "KineticReaction convenience constructor" begin

    calcite = Species("CaCO3"; symbol = "Calcite", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    h2o = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    ca2p = Species("Ca+2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    cs = ChemicalSystem([calcite, h2o, ca2p])
    n_sp = length(cs.species)

    pk = parrot_killoh(PK_PARAMS_C3S, "Calcite")

    # ── Constructor finds species by phreeqc formula ───────────────────────────
    kr = KineticReaction(cs, "Calcite", pk)
    @test kr isa KineticReaction
    @test kr.reaction isa AbstractReaction
    @test kr.idx_mineral == 1
    @test length(kr.stoich) == n_sp
    @test kr.stoich[1] ≈ -1.0
    @test all(iszero, kr.stoich[2:end])
    @test kr.rate_fn === pk
    @test kr.heat_per_mol === nothing

    # ── heat_per_mol is stored correctly ──────────────────────────────────────
    kr_heat = KineticReaction(cs, "Calcite", pk; heat_per_mol = 12_500.0)
    @test kr_heat.heat_per_mol isa Float64
    @test kr_heat.heat_per_mol ≈ 12_500.0

    # heat_per_mol as Quantity → converted to J/mol
    kr_hq = KineticReaction(cs, "Calcite", pk; heat_per_mol = 12.5u"kJ/mol")
    @test kr_hq.heat_per_mol isa Float64
    @test isapprox(kr_hq.heat_per_mol, 12_500.0; rtol = 1.0e-6)

    # ── ArgumentError for unknown species ─────────────────────────────────────
    @test_throws ArgumentError KineticReaction(cs, "Quartz", pk)

    # ── Low-level constructor ─────────────────────────────────────────────────
    rxn = Reaction([calcite, ca2p]; symbol = "calcite dissolution")
    kr_low = KineticReaction(rxn, pk, 1, [-1.0, 0.0, 1.0])
    @test kr_low.reaction isa AbstractReaction
    @test kr_low.idx_mineral == 1
    @test kr_low.stoich[3] ≈ 1.0

    # Low-level: validate idx_mineral > 0
    @test_throws ArgumentError KineticReaction(rxn, pk, 0, [-1.0, 0.0, 1.0])
    # Low-level: validate non-empty stoich
    @test_throws ArgumentError KineticReaction(rxn, pk, 1, Float64[])

    # ── Explicit-Reaction constructor from ChemicalSystem ──────────────────────
    kr_from_rxn = KineticReaction(cs, rxn, pk)
    @test kr_from_rxn.reaction === rxn
    @test kr_from_rxn.idx_mineral == 1

    # ── KineticsProblem via kinetic_species API ──────────────────────────────
    # Need 6 species for 2 nullspace reactions (4 atoms → 6 - 4 = 2 reactions)
    co3 = Species("CO3-2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    hplus = Species("H+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    oh = Species("OH-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    cs_kin = ChemicalSystem(
        [calcite, h2o, ca2p, co3, hplus, oh];
        kinetic_species = Dict("Calcite" => pk)
    )
    kp = KineticsProblem(cs_kin, ChemicalState(cs_kin), (0.0, 1.0))
    @test kp isa KineticsProblem
    @test length(kp.idx_kinetic) == 1
    u0 = build_u0(kp)
    @test length(u0) == 1

end

# ── KineticReaction from annotated Reaction ───────────────────────────────────

@testset "KineticReaction from annotated Reaction" begin

    calcite = Species("CaCO3"; symbol = "Calcite", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    h2o = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    ca2p = Species("Ca+2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    cs = ChemicalSystem([calcite, h2o, ca2p])

    rxn = Reaction([calcite, ca2p]; symbol = "calcite dissolution")

    # ── Missing :rate raises ArgumentError ────────────────────────────────────
    @test_throws ArgumentError KineticReaction(cs, rxn)

    # ── Attach a KineticFunc via :rate ────────────────────────────────────────
    pk = parrot_killoh(PK_PARAMS_C3S, "Calcite")
    rxn[:rate] = pk
    kr = KineticReaction(cs, rxn)
    @test kr isa KineticReaction
    @test kr.rate_fn isa KineticFunc
    @test kr.rate_fn === pk
    @test kr.idx_mineral == 1
    @test kr.heat_per_mol === nothing

    # Stoich: Calcite is reactant (-1), Ca+2 is product (+1), H2O = 0
    @test kr.stoich[1] < 0
    @test kr.stoich[3] > 0
    @test kr.stoich[2] == 0

    # ── Attach a plain callable → auto-wrapped in KineticFunc ─────────────────
    rxn2 = Reaction([calcite, ca2p]; symbol = "plain callable test")
    plain_fn = (T, P, t, n, lna, n0) -> 1.0e-9
    rxn2[:rate] = plain_fn
    kr2 = KineticReaction(cs, rxn2)
    @test kr2.rate_fn isa KineticFunc   # wrapped automatically

    # ── :heat_per_mol picked up correctly ─────────────────────────────────────
    rxn3 = Reaction([calcite, ca2p]; symbol = "heat test")
    rxn3[:rate] = pk
    rxn3[:heat_per_mol] = 12_500.0
    kr3 = KineticReaction(cs, rxn3)
    @test kr3.heat_per_mol isa Float64
    @test kr3.heat_per_mol ≈ 12_500.0

    # ── Build KineticsProblem from kinetic_species API ──────────────────────────
    co3 = Species("CO3-2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    hplus = Species("H+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    oh = Species("OH-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    cs_kin = ChemicalSystem(
        [calcite, h2o, ca2p, co3, hplus, oh];
        kinetic_species = Dict("Calcite" => pk)
    )
    state0 = ChemicalState(cs_kin)
    kp = KineticsProblem(cs_kin, state0, (0.0, 1.0))
    @test kp isa KineticsProblem
    @test length(kp.idx_kinetic) == 1

end
