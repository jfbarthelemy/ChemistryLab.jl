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
    # `u = [nₖ, ξ]`: one kinetic species and one reaction, with no equilibrium
    # solver and no calorimeter. The extents are carried alongside the moles.
    u0 = build_u0(kp)
    @test length(u0) == 1 + length(kp.kinetic_reactions)
    @test u0[1] == 0.0                       # no calcite in the empty state
    @test u0[end] == 0.0                     # ξ(t₀) = 0 by definition

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

# ── Correction factors ────────────────────────────────────────────────────────

@testset "blaine_factor / humidity_factor / powers_alpha_max" begin

    # ── Blaine ────────────────────────────────────────────────────────────────
    @test blaine_factor(385u"m^2/kg") == 1.0
    @test blaine_factor(385.0) == 1.0                    # plain Real = SI
    @test blaine_factor(770u"m^2/kg") ≈ 2.0
    @test blaine_factor(400u"m^2/kg"; blaine_ref = 400u"m^2/kg") == 1.0
    # Waller additions are referred to 400, not 385
    @test blaine_factor(2000u"m^2/kg"; blaine_ref = 400u"m^2/kg") ≈ 5.0

    # ── Relative humidity ─────────────────────────────────────────────────────
    # Hydration stops at and below 80 % RH; the limit from above is ≈ 0.0951.
    @test humidity_factor(0.5) == 0.0
    @test humidity_factor(0.8) == 0.0
    @test humidity_factor(0.8 + 1.0e-9) ≈ ((0.8 - 0.55) / 0.45)^4 rtol = 1.0e-6
    # Lavergne et al. (2018) quote "0.095" for this limit
    @test isapprox(((0.8 - 0.55) / 0.45)^4, 0.095; atol = 5.0e-4)
    @test humidity_factor(1.0) ≈ 1.0
    # Strictly increasing above the cut
    hs = 0.81:0.01:1.0
    @test issorted(humidity_factor.(hs))

    # ── Powers water-availability cap ─────────────────────────────────────────
    @test powers_alpha_max(0.5) == 1.0
    @test powers_alpha_max(0.42) ≈ 1.0
    @test powers_alpha_max(0.32) ≈ 0.32 / 0.42
    @test powers_alpha_max(0.21) ≈ 0.5

    # AD through the cap (calibration on w/c must not need finite differences)
    @test ForwardDiff.derivative(powers_alpha_max, 0.3) ≈ 1 / 0.42

end

# ── parrot_killoh_avrami ──────────────────────────────────────────────────────

@testset "parrot_killoh_avrami" begin

    index = Dict("C3S" => 1)
    lna_sv = StateView([0.0], index)
    n0_sv = StateView([1.0], index)
    T_ref = 293.15
    at(pk, α; T = T_ref, t = 3600.0) =
        pk(T, 1.0e5, t, StateView([1.0 - α], index), lna_sv, n0_sv)

    pk = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S")
    @test pk isa KineticFunc

    # ── The rate is strictly positive at α = 0 ────────────────────────────────
    # Without the PK_AVRAMI_SEED floor the Avrami branch vanishes there, α ≡ 0
    # solves the ODE and hydration never starts.
    @test at(pk, 0.0) > 0
    @test all(isfinite(at(pk, α)) && at(pk, α) > 0 for α in (0.0, 1.0e-8, 0.01, 0.5, 0.99))

    # ── Decreasing once past the nucleation peak, and vanishing at α → 1 ──────
    @test at(pk, 0.9) < at(pk, 0.5) < at(pk, 0.2)
    @test at(pk, 0.999) < 1.0e-6

    # ── Arrhenius ─────────────────────────────────────────────────────────────
    @test at(pk, 0.3; T = 313.15) > at(pk, 0.3) > at(pk, 0.3; T = 278.15)

    # ── Which branch limits which phase ───────────────────────────────────────
    # Lavergne et al. (2018), §1.1.1: with the 1984 parameters "there is no
    # nucleation-growth step for C2S, and no diffusion-controlled step for C3S".
    # Re-derive the three branches independently of the implementation.
    function branches(p, ξ)
        # `ustrip` on a DynamicQuantities quantity already returns SI (s⁻¹ here)
        k₁ = ustrip(p.k₁); n₁ = p.n₁
        k₂ = ustrip(p.k₂)
        k₃ = ustrip(p.k₃); n₃ = p.n₃
        om = max(1 - ξ, 1.0e-12)
        ξa = max(ξ, ChemistryLab.PK_AVRAMI_SEED)
        return (
            (k₁ / n₁) * om * (-log(1 - ξa))^(1 - n₁),          # Avrami
            k₂ * om^(2 / 3) / max(1 - om^(1 / 3), 1.0e-12),    # Jander
            k₃ * om^n₃,                                        # power law
        )
    end
    active(p) = unique(argmin(collect(branches(p, ξ))) for ξ in 0.0:0.01:0.99)

    @test active(PK84_PARAMS_C2S) == [3]           # power law only
    @test sort(active(PK84_PARAMS_C3S)) == [1, 3]  # Avrami then power law, no Jander
    @test 2 in active(PK84_PARAMS_C3A)             # aluminates do see diffusion
    @test 2 in active(PK84_PARAMS_C4AF)

    # ── The implementation agrees with that independent derivation ────────────
    for α in (0.05, 0.2, 0.5, 0.8)
        @test at(pk, α) ≈ minimum(branches(PK84_PARAMS_C3S, α)) rtol = 1.0e-10
    end

    # ── Blaine and humidity corrections compose multiplicatively ──────────────
    pk_fine = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S"; blaine = 770u"m^2/kg")
    @test at(pk_fine, 0.3) ≈ 2 * at(pk, 0.3) rtol = 1.0e-10

    # ── The same oracle, applied to the DEPRECATED smoothed variant ──────────
    #
    # This is the test that would have caught it. `parrot_killoh` with
    # `PK_PARAMS_*` is limited by its diffusion branch for EVERY phase from
    # α ≈ 0.028 on, because its prefactor `3K₃/N₃ = 0.0018 d⁻¹` is 28 times
    # smaller than the canonical `k₂ = 0.05 d⁻¹`. The rate then integrates in
    # closed form, `α = α_max[1 − (1 − √(2K₃t/N₃))³]`, which does not depend on
    # the phase at all: all four reach 0.2386 at seven days, against the 0.61
    # the cement literature reports for a CEM I at w/c = 0.40.
    #
    # The behavior is pinned rather than fixed: the variant is deprecated and
    # its attribution withdrawn, so what matters is that the discrepancy stays
    # visible instead of being rediscovered from a demo that looks plausible.
    function branches_smoothed(p, ξ)
        K₁ = ustrip(p.K₁); N₁ = p.N₁
        K₂ = ustrip(p.K₂); N₂ = p.N₂
        K₃ = ustrip(p.K₃); N₃ = p.N₃
        om = 1 - ξ
        return (
            (K₁ / N₁) * om^N₁ / (1 + p.B * ξ^N₃),                       # nucleation-growth
            K₂ * om^N₂,                                                 # interaction
            3 * K₃ * om^(2 / 3) / (N₃ * max(1 - om^(1 / 3), 1.0e-10)),  # diffusion (Jander)
        )
    end

    # `min(max(r_NG, r_I), r_D)`: diffusion caps the other two.
    smoothed_is_diffusion_limited(p, ξ) = let b = branches_smoothed(p, ξ)
        b[3] <= max(b[1], b[2])
    end

    # Measured, over the interval a seven-day run traverses: THREE of the four
    # phases are diffusion-limited there and therefore land on the same
    # phase-blind Jander value. C₄AF is the exception — its own
    # nucleation-growth branch is slower than diffusion, so it limits instead,
    # and C₄AF comes out LOWER still (0.193 against 0.239 in the run that
    # exposed this). Both failures have the same cause and neither is physical.
    for p in (PK_PARAMS_C3S, PK_PARAMS_C2S, PK_PARAMS_C3A)
        @test all(smoothed_is_diffusion_limited(p, ξ) for ξ in 0.06:0.01:0.30)
    end
    @test !smoothed_is_diffusion_limited(PK_PARAMS_C3S, 0.001)
    @test !any(smoothed_is_diffusion_limited(PK_PARAMS_C4AF, ξ) for ξ in 0.06:0.01:0.30)

    # The closed-form Jander solution, and the fact that it is phase-blind:
    # K₃ = 0.0024 d⁻¹ and N₃ = 4 are identical across all four parameter sets,
    # even though K₁ spans a factor of 18.
    jander_alpha(p, t_days) = let
        K₃ = ustrip(us"1/d", p.K₃)
        1 - (1 - sqrt(2 * K₃ * t_days / p.N₃))^3
    end
    for p in (PK_PARAMS_C3S, PK_PARAMS_C2S, PK_PARAMS_C3A)
        @test jander_alpha(p, 7.0) ≈ 0.250524 atol = 1.0e-5
    end
    @test maximum(ustrip(us"1/d", p.K₁) for p in
        (PK_PARAMS_C3S, PK_PARAMS_C2S, PK_PARAMS_C3A, PK_PARAMS_C4AF)) /
        minimum(ustrip(us"1/d", p.K₁) for p in
        (PK_PARAMS_C3S, PK_PARAMS_C2S, PK_PARAMS_C3A, PK_PARAMS_C4AF)) > 18.0

    # Phase-blind: one distinct value across the four parameter sets, because
    # K₃ and N₃ are identical in all of them.
    @test length(unique(round(jander_alpha(p, 7.0); digits = 6) for p in
        (PK_PARAMS_C3S, PK_PARAMS_C2S, PK_PARAMS_C3A, PK_PARAMS_C4AF))) == 1

    pk_dry = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S"; humidity = 0.7)
    @test at(pk_dry, 0.3) == 0.0                                  # hydration stopped
    pk_h = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S"; humidity = 0.9)
    @test at(pk_h, 0.3) ≈ humidity_factor(0.9) * at(pk, 0.3) rtol = 1.0e-10

    # A time-dependent humidity is accepted as a callable
    pk_ht = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S"; humidity = t -> t < 100 ? 0.99 : 0.7)
    @test at(pk_ht, 0.3; t = 10.0) > 0
    @test at(pk_ht, 0.3; t = 1000.0) == 0.0

    # ── α_max cap ─────────────────────────────────────────────────────────────
    pk_cap = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S"; α_max = powers_alpha_max(0.32))
    @test at(pk_cap, 0.76) < at(pk_cap, 0.3)      # already near the cap
    @test at(pk_cap, 0.8) >= 0

    # ── AD through the rate ───────────────────────────────────────────────────
    d = ForwardDiff.derivative(
        α -> pk(T_ref, 1.0e5, 3600.0, StateView([1.0 - α], index), lna_sv, n0_sv), 0.4
    )
    @test isfinite(d)
    dT = ForwardDiff.derivative(
        T -> pk(T, 1.0e5, 3600.0, StateView([0.6], index), lna_sv, n0_sv), T_ref
    )
    @test isfinite(dT) && dT > 0

    # ── The two variants are distinct objects with non-transferable parameters ─
    @test parrot_killoh(PK_PARAMS_C3S, "C3S")(T_ref, 1.0e5, 3600.0, StateView([0.6], index), lna_sv, n0_sv) !=
        at(pk, 0.4)

end

# ── waller ────────────────────────────────────────────────────────────────────

@testset "waller" begin

    index = Dict("FA" => 1)
    lna_sv = StateView([0.0], index)
    n0_sv = StateView([1.0], index)
    T_ref = 293.15
    day = 86400.0
    at(w, α, t; T = T_ref) = w(T, 1.0e5, t, StateView([1.0 - α], index), lna_sv, n0_sv)

    w = waller(WALLER_PARAMS_FLY_ASH, "FA")
    @test w isa KineticFunc

    # ── Positive and finite everywhere, including the singular α = 0 point ────
    @test all(
        isfinite(at(w, α, t)) && at(w, α, t) > 0
            for (α, t) in ((0.0, 1.0), (0.0, day), (0.01, day), (0.5, 80day), (0.95, 500day))
    )

    # ── The closed form integrates back to the sigmoid it came from ───────────
    # α(t) = 1/(1+(τ/t)^n); integrate α̇(α) by RK4 from t = 1 h and compare.
    τ, n_w = 80 * day, 0.7
    sigmoid(t) = 1 / (1 + (τ / t)^n_w)
    t0 = 3600.0
    α = sigmoid(t0)
    t = t0
    nsteps = 40_000
    h = (200day - t0) / nsteps
    f(a, tt) = at(w, a, tt)
    for _ in 1:nsteps
        k1 = f(α, t); k2 = f(α + h / 2 * k1, t + h / 2)
        k3 = f(α + h / 2 * k2, t + h / 2); k4 = f(α + h * k3, t + h)
        α += h / 6 * (k1 + 2k2 + 2k3 + k4)
        t += h
    end
    @test isapprox(α, sigmoid(200day); rtol = 1.0e-3)

    # ── α(τ) = 1/2 for the shipped fly-ash parameters ─────────────────────────
    @test isapprox(sigmoid(τ), 0.5; atol = 1.0e-12)

    # ── Monotone decay of the rate with the degree of reaction ────────────────
    @test at(w, 0.9, 100day) < at(w, 0.5, 100day) < at(w, 0.1, 100day)

    # ── Pozzolanic reactions are more thermo-activated than the clinker ───────
    ratio_fa = at(w, 0.3, 30day; T = 313.15) / at(w, 0.3, 30day)
    pk = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S")
    idx3 = Dict("C3S" => 1)
    r(T) = pk(T, 1.0e5, 30day, StateView([0.7], idx3), StateView([0.0], idx3), StateView([1.0], idx3))
    @test ratio_fa > r(313.15) / r(293.15)

    # ── Blaine is referred to the addition's own reference (400, not 385) ─────
    w_fine = waller(WALLER_PARAMS_FLY_ASH, "FA"; blaine = 800u"m^2/kg")
    @test at(w_fine, 0.3, 30day) ≈ 2 * at(w, 0.3, 30day) rtol = 1.0e-10

    # Silica fume shares the kinetics; its reactivity is carried by the fineness
    w_sf = waller(WALLER_PARAMS_SILICA_FUME, "FA"; blaine = 2000u"m^2/kg")
    @test at(w_sf, 0.3, 30day) ≈ 5 * at(w, 0.3, 30day) rtol = 1.0e-10

    # ── Slag reacts more slowly and incompletely ──────────────────────────────
    w_slag = waller(WALLER_PARAMS_SLAG, "FA"; α_max = 0.9)
    @test at(w_slag, 0.3, 30day) < at(w, 0.3, 30day)

    # ── AD ────────────────────────────────────────────────────────────────────
    @test isfinite(
        ForwardDiff.derivative(
            a -> w(T_ref, 1.0e5, 30day, StateView([1.0 - a], index), lna_sv, n0_sv), 0.3
        )
    )

end
