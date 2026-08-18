using LinearAlgebra

# The dual solver exists because the interior-point route does not reach
# stationarity on a cement, and — more to the point — cannot tell you whether it
# has. These tests check the two things that matter: that the certificate is a
# real proof, and that the certified answer is the right one.

const DUAL_DATA = joinpath(pkgdir(ChemistryLab), "data", "cemdata18-thermofun.json")

# Reaktoro 2.13.0, Cemdata18, ideal activities; 55.5 mol H2O, 0.05 mol calcite,
# 0.01 mol CO2. The same reference as `equilibrium_reference.jl`.
const DUAL_REAKTORO = Dict(
    "H2O@" => 55.4961389, "H+" => 3.69286738e-7, "OH-" => 2.70633096e-8,
    "CO2@" => 0.00613892838, "HCO3-" => 0.00738915525, "CO3-2" => 9.37864261e-7,
    "Ca+2" => 0.00352901991, "CaOH+" => 1.58551676e-9,
    "Ca(CO3)@" => 5.54793457e-6, "Ca(HCO3)+" => 0.000332647348,
    "Cal" => 0.0461327832,
)

function _dual_calcite_system()
    dict = Dict(symbol(s) => s for s in build_species(DUAL_DATA))
    syms = collect(keys(DUAL_REAKTORO))
    return ChemicalSystem([dict[s] for s in syms], ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
end

@testset "the certificate is a proof, not an opinion" begin

    cs = _dual_calcite_system()
    names = symbol.(cs.species)
    n = Any[fill(0.0u"mol", length(cs.species))...]
    n[findfirst(==("H2O@"), names)] = 55.5u"mol"
    n[findfirst(==("Cal"), names)] = 0.05u"mol"
    n[findfirst(==("CO2@"), names)] = 0.01u"mol"
    st = ChemicalState(cs, n)

    A = Float64.(cs.SM.A)
    b = A * [ustrip(us"mol", x) for x in st.n]

    des = DualEquilibriumSolver(cs, DiluteSolutionModel())
    ipm = equilibrate(st, OptimaOptimizer())
    dual = SciMLBase.solve(des, ipm; b = b)

    cert = optimality_certificate(des, dual; b = b)
    @test cert.optimal
    @test cert.stationarity < 1.0e-9        # interior species are stationary
    @test cert.balance < 1.0e-9             # matter is conserved
    @test cert.worst_supersaturation <= 1.0e-8   # nothing absent is supersaturated

    # The certificate must also be able to say NO. The interior-point answer to
    # the same problem does not satisfy the KKT conditions, and the certificate
    # detects it — which is the whole reason for having one.
    cert_ipm = optimality_certificate(des, ipm; b = b)
    @test !cert_ipm.optimal
    @test cert_ipm.stationarity > 1.0e-3

end

@testset "the certified answer matches Reaktoro" begin

    cs = _dual_calcite_system()
    names = symbol.(cs.species)
    n = Any[fill(0.0u"mol", length(cs.species))...]
    n[findfirst(==("H2O@"), names)] = 55.5u"mol"
    n[findfirst(==("Cal"), names)] = 0.05u"mol"
    n[findfirst(==("CO2@"), names)] = 0.01u"mol"
    st = ChemicalState(cs, n)
    b = Float64.(cs.SM.A) * [ustrip(us"mol", x) for x in st.n]

    des = DualEquilibriumSolver(cs, DiluteSolutionModel())
    dual = SciMLBase.solve(des, equilibrate(st, OptimaOptimizer()); b = b)

    # EVERY species to 1 %, where the interior-point method needs 5 % and misses
    # one species by 147 %.
    for (k, nm) in enumerate(names)
        @test ustrip(us"mol", dual.n[k]) ≈ DUAL_REAKTORO[nm] rtol = 1.0e-2
    end

    # `CaOH+` in particular. `equilibrium_reference.jl` records this species as
    # `@test_broken`: the interior-point answer is 1.6 times too large, because a
    # trace species contributes almost nothing to the Gibbs energy and its
    # stationarity is the first thing an interior-point method gives up on. The
    # certified answer has it right, and that is not luck — stationarity is
    # imposed on every interior species alike, whatever it weighs.
    k = findfirst(==("CaOH+"), names)
    @test ustrip(us"mol", dual.n[k]) ≈ DUAL_REAKTORO["CaOH+"] rtol = 1.0e-2

end

@testset "phase stability is decided, not assumed" begin

    # Calcite in pure water. Portlandite is available and must NOT appear; the
    # active set has to reject it. Getting this wrong is not benign: admitting it
    # makes two phases simultaneously saturated, which over-determines the element
    # potentials and leaves the Newton grinding on infeasible rows.
    substances = build_species(DUAL_DATA)
    sp = speciation(substances, ["Cal", "Portlandite"]; aggregate_state = [AS_AQUEOUS])
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)

    st = ChemicalState(cs)
    set_quantity!(st, "Cal", 0.01u"mol")
    set_quantity!(st, "H2O@", 55.5u"mol")
    b = Float64.(cs.SM.A) * [ustrip(us"mol", x) for x in st.n]

    des = DualEquilibriumSolver(cs, DiluteSolutionModel())
    dual = SciMLBase.solve(des, equilibrate(st, OptimaOptimizer()); b = b)
    cert = optimality_certificate(des, dual; b = b)

    @test cert.optimal
    names = symbol.(cs.species)
    @test ustrip(us"mol", dual.n[findfirst(==("Portlandite"), names)]) < 1.0e-20
    @test ustrip(us"mol", dual.n[findfirst(==("Cal"), names)]) > 1.0e-3

    # Calcite in pure water sits near pH 9.9. The interior-point answer to the
    # same problem is pH 7.0 — not an imprecision but a wrong answer, and one
    # nothing in that solver's output reveals.
    @test 9.5 < pH(dual) < 10.3

end

@testset "a system without a solvent is refused" begin

    # The element potentials are fixed by the mass-action laws of the aqueous
    # species. Without them there is nothing to determine `y`, and the right
    # response is to say so rather than to return a number.
    dict = Dict(symbol(s) => s for s in build_species(DUAL_DATA))
    cs = ChemicalSystem([dict["Cal"], dict["Portlandite"]])
    @test_throws ArgumentError DualEquilibriumSolver(cs, DiluteSolutionModel())

end

@testset "a vanished component is recognised by sign, not by zero" begin

    # The criterion itself is linear algebra and is tested in OptimaSolver. What
    # belongs here is its chemical consequence: with iron, sulfur and aluminium
    # absent from a calcite–water system, their species must come back at exactly
    # zero while the carbonate system is untouched.
    #
    # Getting the criterion wrong is not benign. `b_k = 0` forces the species of
    # component `k` to vanish only when the non-zero entries of that row share a
    # sign. The `H+` row carries `+1` for `H+` and `−1` for `OH-`, so its zero
    # total is the ordinary state of pure water; treating it as degenerate kills
    # the whole acid–base system and returns pH 7.000 with the calcite
    # undissolved, instead of pH 9.90.
    substances = build_species(DUAL_DATA)
    sp = speciation(substances, ["Cal", "Portlandite"]; aggregate_state = [AS_AQUEOUS])
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    st = ChemicalState(cs)
    set_quantity!(st, "Cal", 0.01u"mol")
    set_quantity!(st, "H2O@", 55.5u"mol")
    bv = Float64.(cs.SM.A) * [ustrip(us"mol", x) for x in st.n]

    des = DualEquilibriumSolver(cs, DiluteSolutionModel())
    dual = SciMLBase.solve(des, equilibrate(st, OptimaOptimizer()); b = bv)
    names = symbol.(cs.species)

    for nm in ("Ca+2", "CO3-2", "OH-")
        k = findfirst(==(nm), names)
        k === nothing && continue
        @test ustrip(us"mol", dual.n[k]) > 0     # the carbonate system survives
    end
    @test pH(dual) > 9.0                          # and it is not pure water

    @test optimality_certificate(des, dual; b = bv).optimal

    # And the reduction itself, on a system that HAS an absent component: bring
    # ettringite in, which puts sulfur and aluminium among the components, and
    # start with none of either. Every species carrying them must be excluded —
    # they are zero by the CONSTRAINT, and the potential of a component nobody
    # supplies is determined by nothing, so their saturation indices are
    # meaningless and testing them reported violations of several hundred RT.
    sp2 = speciation(substances, ["Cal", "ettringite"]; aggregate_state = [AS_AQUEOUS])
    cs2 = ChemicalSystem(sp2, CEMDATA_PRIMARIES)
    st2 = ChemicalState(cs2)
    set_quantity!(st2, "Cal", 0.01u"mol")
    set_quantity!(st2, "H2O@", 55.5u"mol")
    bv2 = Float64.(cs2.SM.A) * [ustrip(us"mol", x) for x in st2.n]

    des2 = DualEquilibriumSolver(cs2, DiluteSolutionModel())
    dual2 = SciMLBase.solve(des2, equilibrate(st2, OptimaOptimizer()); b = bv2)
    cert2 = optimality_certificate(des2, dual2; b = bv2)

    @test cert2.n_absent_component > 0
    @test cert2.optimal
    names2 = symbol.(cs2.species)
    k = findfirst(==("ettringite"), names2)
    @test ustrip(us"mol", dual2.n[k]) < 1.0e-20   # no sulfur, hence no ettringite

end

@testset "solid solutions are refused, not silently mishandled" begin

    # A pure phase has unit activity, hence `gᵢ = uᵢ` and a bound-constrained
    # variable that can be EXACTLY zero. An end-member of a solid solution
    # cannot: its activity goes to −∞ as its mole fraction goes to zero, so it is
    # never exactly absent while the phase exists, and the active set belongs at
    # the level of the PHASE, decided by a tangent-plane test. Treating an
    # end-member as a pure phase returns a composition that looks reasonable and
    # is not the minimum — which a caller cannot see, hence the explicit refusal.
    substances = build_species(DUAL_DATA)
    toml = joinpath(pkgdir(ChemistryLab), "data", "solid_solutions.toml")
    members = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD"]
    sp = speciation(substances, vcat("Portlandite", members); aggregate_state = [AS_AQUEOUS])
    ss = [x for x in build_solid_solutions(toml, Dict(symbol(s) => s for s in sp)) if x.name == "CSHQ"]
    @test length(ss) == 1

    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)
    @test_throws ArgumentError DualEquilibriumSolver(cs, DiluteSolutionModel())

    # Without the declaration the same species are ordinary pure phases and the
    # solver accepts them: the refusal is about the solution, not the species.
    cs_plain = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    @test DualEquilibriumSolver(cs_plain, DiluteSolutionModel()) isa DualEquilibriumSolver

end
