using OrderedCollections

# A minimal cement-like system built from CEMDATA18, reused by every testset here.
const _VF_DATA = datapath("cemdata18-thermofun.json")

function _vf_system()
    substances = build_species(_VF_DATA)
    syms = split("C3S C2S C3A C4AF Portlandite Jennite ettringite C3AH6 C3FH6 H2O@")
    sp = speciation(substances, syms; aggregate_state = [AS_AQUEOUS])
    return ChemicalSystem(sp, CEMDATA_PRIMARIES)
end

@testset "volume_fractions" begin

    cs = _vf_system()

    state = ChemicalState(cs)
    set_quantity!(state, "C3S", 1.0u"mol")
    set_quantity!(state, "Portlandite", 2.0u"mol")
    set_quantity!(state, "H2O@", 10.0u"mol")

    # ── Default normalization: fractions of the CURRENT volume, summing to 1 ──
    f = volume_fractions(state)
    @test f isa OrderedDict{String, Float64}
    @test isapprox(sum(values(f)), 1.0; atol = 1.0e-12)
    # Species with a zero amount are absent, not zero-valued
    @test !haskey(f, "C2S")
    @test haskey(f, "C3S") && haskey(f, "Portlandite") && haskey(f, "H2O@")

    V_tot = volume(state).total

    # Solid and solvent fractions are positive and dominate
    @test f["C3S"] > 0 && f["Portlandite"] > 0 && f["H2O@"] > 0
    @test f["C3S"] + f["Portlandite"] + f["H2O@"] > 0.999

    # Aqueous solutes may contribute a NEGATIVE fraction (electrostriction).
    # They are kept, which is what makes the balance close against `volume`.
    @test haskey(f, "OH-")
    @test f["OH-"] < 0
    @test isapprox(f["OH-"], ustrip(volume(state, "OH-") / V_tot); rtol = 1.0e-12)

    # Consistent with the per-species `volume` accessor already in the package
    @test isapprox(f["Portlandite"], ustrip(volume(state, "Portlandite") / V_tot); rtol = 1.0e-12)

    # `reference = state` itself must reproduce the default up to the void entry
    f_self = volume_fractions(state; reference = state)
    @test isapprox(f_self["void"], 0.0; atol = 1.0e-10)
    for k in keys(f)
        @test isapprox(f_self[k], f[k]; rtol = 1.0e-12)
    end

    # ── Zero reference volume is an error, not a NaN ──────────────────────────
    @test_throws ArgumentError volume_fractions(ChemicalState(cs))

end

@testset "chemical_shrinkage and the sealed-volume void" begin

    cs = _vf_system()

    # C3S + 5.3 H → C1.7SH4 + 1.3 CH is contracting. CEMDATA18 has no C1.7SH4,
    # so use the balanced Jennite reaction shipped with the package scripts:
    #   C3S + 103/30 H2O → Jennite + 4/3 Portlandite
    state0 = ChemicalState(cs)
    set_quantity!(state0, "C3S", 1.0u"mol")
    set_quantity!(state0, "H2O@", 5.0u"mol")

    ξ = 0.4
    state = ChemicalState(cs)
    set_quantity!(state, "C3S", (1.0 - ξ)u"mol")
    set_quantity!(state, "H2O@", (5.0 - ξ * 103 / 30)u"mol")
    set_quantity!(state, "Jennite", (ξ * 1.0)u"mol")
    set_quantity!(state, "Portlandite", (ξ * 4 / 3)u"mol")

    δV = chemical_shrinkage(state, state0)
    # Hydration of alite contracts: products occupy less than reactants + water.
    @test ustrip(δV) > 0                       # SI: m³

    # The void closes the balance exactly, by construction
    f = volume_fractions(state; reference = state0)
    @test isapprox(sum(values(f)), 1.0; atol = 1.0e-12)
    @test isapprox(f["void"], ustrip(δV / volume(state0).total); rtol = 1.0e-10)

    # Without a reference the void is absent and the fractions renormalize
    f_open = volume_fractions(state)
    @test !haskey(f_open, "void")
    @test isapprox(sum(values(f_open)), 1.0; atol = 1.0e-12)
    @test f_open["Portlandite"] > f["Portlandite"]   # same matter, smaller denominator

    # A custom void key is honored
    f_named = volume_fractions(state; reference = state0, void_key = "chemical shrinkage")
    @test haskey(f_named, "chemical shrinkage")
    @test !haskey(f_named, "void")

    # Zero extent → zero shrinkage
    @test isapprox(ustrip(chemical_shrinkage(state0, state0)), 0.0; atol = 1.0e-20)

end

@testset "volume_fractions with groups" begin

    cs = _vf_system()
    state0 = ChemicalState(cs)
    set_quantity!(state0, "C3S", 1.0u"mol")
    set_quantity!(state0, "C2S", 0.5u"mol")
    set_quantity!(state0, "H2O@", 6.0u"mol")

    state = ChemicalState(cs)
    set_quantity!(state, "C3S", 0.4u"mol")
    set_quantity!(state, "C2S", 0.3u"mol")
    set_quantity!(state, "H2O@", 3.0u"mol")
    set_quantity!(state, "Jennite", 0.8u"mol")
    set_quantity!(state, "Portlandite", 0.9u"mol")

    groups = [
        "anhydrous" => ["C3S", "C2S", "C3A", "C4AF"],
        "C-S-H" => "Jennite",
        "CH" => "Portlandite",
        "water" => "H2O@",
    ]

    g = volume_fractions(state, groups; reference = state0)
    @test isapprox(sum(values(g)), 1.0; atol = 1.0e-12)

    # A group is exactly the sum of its members
    f = volume_fractions(state; reference = state0)
    @test isapprox(g["anhydrous"], f["C3S"] + f["C2S"]; rtol = 1.0e-12)
    @test isapprox(g["C-S-H"], f["Jennite"]; rtol = 1.0e-12)

    # Groups declared but absent from the state are kept at zero — the caller
    # asked for that phase and a missing key would break downstream indexing
    g2 = volume_fractions(state, ["AFt" => "ettringite", "CH" => "Portlandite"]; reference = state0)
    @test g2["AFt"] == 0.0

    # Unassigned species are collected, never silently dropped
    partial = volume_fractions(state, ["CH" => "Portlandite"]; reference = state0)
    @test haskey(partial, "other")
    @test isapprox(sum(values(partial)), 1.0; atol = 1.0e-12)

    # ── Double assignment is an error, not a double count ─────────────────────
    @test_throws ArgumentError volume_fractions(
        state, ["a" => ["C3S", "Jennite"], "b" => "Jennite"]; reference = state0
    )

    # Without a reference, no void group appears
    g_open = volume_fractions(state, groups)
    @test !haskey(g_open, "void")
    @test isapprox(sum(values(g_open)), 1.0; atol = 1.0e-12)

end

@testset "missing_molar_volumes" begin

    cs = _vf_system()
    state = ChemicalState(cs)
    set_quantity!(state, "C3S", 1.0u"mol")
    # Every CEMDATA18 species used here carries V⁰
    @test isempty(missing_molar_volumes(state))

    # A hand-built species without V⁰ is invisible to the volume balance, and
    # `missing_molar_volumes` is what makes that visible.
    sp_x = Species(
        "CaAl2Si2O8"; symbol = "GGBS", aggregate_state = AS_CRYSTAL,
        properties = Dict{Symbol, Any}(:M => 0.095u"kg/mol")
    )
    cs2 = ChemicalSystem(vcat(cs.species, sp_x), CEMDATA_PRIMARIES)
    state2 = ChemicalState(cs2)
    set_quantity!(state2, "C3S", 1.0u"mol")
    set_quantity!(state2, "GGBS", 1.0u"mol")

    @test "GGBS" in missing_molar_volumes(state2)
    @test !("C3S" in missing_molar_volumes(state2))
    # It contributes nothing to the fractions — which is exactly why the warning matters
    @test !haskey(volume_fractions(state2), "GGBS")

    # A species with no V⁰ but also no matter present is not reported
    state3 = ChemicalState(cs2)
    set_quantity!(state3, "C3S", 1.0u"mol")
    @test isempty(missing_molar_volumes(state3))

end

@testset "porosity of a setting binder" begin

    # The one-argument `porosity` is `(V_liquid + V_gas)/V_total`, which is wrong
    # for a hydrating cement on both ends: the denominator shrinks with the
    # reactions, and the empty porosity left by the chemical shrinkage is not a
    # species and so is invisible. The two-argument method fixes both.
    cs = _vf_system()

    state0 = ChemicalState(cs)
    set_quantity!(state0, "C3S", 1.0u"mol")
    set_quantity!(state0, "H2O@", 5.0u"mol")

    # The Jennite reaction, contracting (see the testset above).
    ξ = 0.4
    st = ChemicalState(cs)
    set_quantity!(st, "C3S", (1.0 - ξ)u"mol")
    set_quantity!(st, "H2O@", (5.0 - ξ * 103 / 30)u"mol")
    set_quantity!(st, "Jennite", (ξ * 1.0)u"mol")
    set_quantity!(st, "Portlandite", (ξ * 4 / 3)u"mol")

    p = porosity(st, state0)
    @test p.liquid > 0
    @test p.void > 0                                   # sealed hydration contracts
    @test p.total ≈ p.liquid + p.void

    # ── Consistent with `volume_fractions` on the same reference ─────────────
    # This is the contract: the porosity is the water plus the void of the very
    # same balance, so a micromechanical model fed by either agrees.
    f = volume_fractions(st, ["water" => "H2O@"]; reference = state0)
    @test isapprox(p.void, f["void"]; rtol = 1.0e-10)
    # `f["water"]` is the solvent only; `p.liquid` is the whole aqueous phase,
    # solutes included, so it is the larger of the two by the (negative) solute
    # contribution — tiny here, but the ordering is the point.
    @test p.liquid <= f["water"] + 1.0e-6

    # ── It is strictly larger than the one-argument version ──────────────────
    # Both errors push the same way: the naive ratio understates the porosity.
    @test p.total > porosity(st)

    # ── At the reference state the void vanishes ─────────────────────────────
    p0 = porosity(state0, state0)
    @test isapprox(p0.void, 0.0; atol = 1.0e-10)
    @test isapprox(p0.total, porosity(state0); rtol = 1.0e-10)

    # ── Saturation counts the empty porosity ────────────────────────────────
    s = saturation(st, state0)
    @test 0 < s < 1
    @test isapprox(s, p.liquid / p.total; rtol = 1.0e-12)
    # A sealed paste desaturates as it hydrates, though no water ever leaves it.
    @test s < saturation(state0, state0)

    # ── Guard rail ───────────────────────────────────────────────────────────
    @test_throws ArgumentError porosity(st, ChemicalState(cs))

end
