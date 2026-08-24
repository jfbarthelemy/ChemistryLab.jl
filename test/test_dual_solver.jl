using LinearAlgebra

# The dual solver exists because the interior-point route does not reach
# stationarity on a cement, and — more to the point — cannot tell you whether it
# has. These tests check the two things that matter: that the certificate is a
# real proof, and that the certified answer is the right one.

const DUAL_DATA = datapath("cemdata18-thermofun.json")

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
    # belongs here is its chemical consequence: with iron, sulfur and aluminum
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
    # ettringite in, which puts sulfur and aluminum among the components, and
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

@testset "a solid solution is solved as a mixing phase" begin

    # A pure phase has unit activity: `gᵢ = uᵢ` when present, and it can be
    # exactly zero — a bound-constrained variable with an active set. An
    # end-member of a solid solution cannot: its activity goes to −∞ as its mole
    # fraction goes to zero, so it is never exactly absent while the phase
    # exists. The active set belongs at the level of the PHASE, admitted by
    # Michelsen's tangent-plane measure rather than by the sign of a saturation
    # index.
    #
    # The difference is not presentational. Declared as separate pure phases, the
    # four CSHQ end-members give a composition that is not certifiable and that
    # puts everything in ONE of them; declared as a solution, the same system is
    # certified and distributes over all four, which is what a solid solution is.
    substances = build_species(DUAL_DATA)
    toml = datapath("solid_solutions.toml")
    members = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD"]
    sp = speciation(substances, vcat("Portlandite", members); aggregate_state = [AS_AQUEOUS])
    ss = [x for x in build_solid_solutions(toml, Dict(symbol(s) => s for s in sp)) if x.name == "CSHQ"]
    @test length(ss) == 1

    # Portlandite plus silica at C/S = 1.7, electrically neutral by construction.
    function solve_cs(cs)
        st = ChemicalState(cs)
        set_quantity!(st, "Portlandite", 1.7u"mol")
        set_quantity!(st, "SiO2@", 1.0u"mol")
        set_quantity!(st, "H2O@", 40.0u"mol")
        b = Float64.(cs.SM.A) * [ustrip(us"mol", x) for x in st.n]
        des = DualEquilibriumSolver(cs, DiluteSolutionModel())
        d = SciMLBase.solve(des, equilibrate(st, OptimaOptimizer()); b = b)
        return d, optimality_certificate(des, d; b = b), des
    end

    cs_ss = ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)
    d_ss, cert_ss, des_ss = solve_cs(cs_ss)

    @test cert_ss.optimal
    @test cert_ss.stationarity < 1.0e-9
    @test cert_ss.balance < 1.0e-9

    names = symbol.(cs_ss.species)
    amounts = [ustrip(us"mol", d_ss.n[findfirst(==(nm), names)]) for nm in members]
    @test all(>(0), amounts)                       # a solution mixes
    @test count(>(1.0e-3), amounts) >= 3           # and not into a single corner
    @test 12.0 < pH(d_ss) < 12.8                   # a C-S-H pore solution

    # The end-members are carried as a mixing phase, so none is left among the
    # bound-constrained variables where it would be a pure phase.
    ss_idx = Set(findfirst(==(nm), names) for nm in members)
    @test length(des_ss.ss_groups) == 1
    @test Set(des_ss.ss_groups[1]) == ss_idx
    @test isempty(intersect(Set(des_ss.idx_pure), ss_idx))
    @test findfirst(==("Portlandite"), names) in des_ss.idx_pure

    # The same species WITHOUT the declaration are ordinary pure phases, and that
    # is a different — and worse — problem.
    #
    # Four pure phases of nearly the same composition compete, and the free energy
    # awards everything to the winner: 1.0 mol in a single end-member and zero in
    # the other three, because nothing in that model pays for mixing. The optimum
    # is a corner of the feasible set, and the solver does not reach it exactly —
    # the answer comes back with a stationarity residual of order 1 and is
    # correctly refused by the certificate. The solution model is not a
    # presentational choice: it is what makes the problem solvable.
    cs_pure = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    d_pure, cert_pure, _ = solve_cs(cs_pure)
    npure = [ustrip(us"mol", d_pure.n[findfirst(==(nm), symbol.(cs_pure.species))]) for nm in members]
    @test !cert_pure.optimal
    # Strictly fewer end-members than the solution model occupies, and the count
    # is not asserted exactly: nothing in that model pays for mixing, so which
    # corner of the feasible set it settles in is a property of the solver's path,
    # not of the thermodynamics. What IS a property of the thermodynamics is that
    # it settles in a corner at all, while the solution spreads over all four.
    @test count(>(1.0e-3), npure) < count(>(1.0e-3), amounts)
    # No mass-balance assertion on this answer: it is NOT an equilibrium — the
    # certificate says so — and its element balance does not close. Asserting a
    # conservation law on a composition that fails the constraints would be
    # asserting a property of a non-solution.
    @test pH(d_pure) ≈ pH(d_ss) rtol = 1.0e-2      # the pore solution is still right

end

@testset "solve_certified returns a proved answer, whichever start produced it" begin

    # The certificate is an oracle, so offering several starting compositions and
    # keeping a PROVED answer is rigorous rather than opportunistic: the proof does
    # not depend on having predicted which start would win. It is also necessary —
    # measured on an LC³ equilibrium, the interior-point back ends each certify a
    # degree of reaction the other misses.
    cs = _dual_calcite_system()
    names = symbol.(cs.species)
    n = Any[fill(0.0u"mol", length(cs.species))...]
    n[findfirst(==("H2O@"), names)] = 55.5u"mol"
    n[findfirst(==("Cal"), names)] = 0.05u"mol"
    n[findfirst(==("CO2@"), names)] = 0.01u"mol"
    good = ChemicalState(cs, n)
    b = Float64.(cs.SM.A) * [ustrip(us"mol", x) for x in good.n]

    des = DualEquilibriumSolver(cs, DiluteSolutionModel())

    # A start with the carbonate budget in the wrong species entirely: everything
    # dissolved, nothing precipitated.
    bad_n = Any[fill(1.0e-12u"mol", length(cs.species))...]
    bad_n[findfirst(==("H2O@"), names)] = 55.5u"mol"
    bad = ChemicalState(cs, bad_n)

    # Bad start FIRST, so the function has to move on rather than report it.
    eq, cert = solve_certified(des, [bad, good]; b = b)
    @test cert.optimal
    @test cert.stationarity < 1.0e-9
    @test cert.balance < 1.0e-9
    # pH 6.43, not the 9.9 of calcite in PURE water: this system carries 0.01 mol
    # of CO₂, which acidifies it, and 6.43 is Reaktoro's own answer for it
    # (`H⁺ = 3.69e-7` in `DUAL_REAKTORO` above).
    @test 6.3 < pH(eq) < 6.6

    # Both starts land on the same point, which is what convexity says must happen
    # — and a solver returning different answers from different starts would be
    # stopping short of stationarity, not finding different minima.
    eq_from_good, _ = solve_certified(des, [good]; b = b)
    for k in eachindex(names)
        @test ustrip(us"mol", eq.n[k]) ≈ ustrip(us"mol", eq_from_good.n[k]) rtol = 1.0e-5
    end

    # With only the bad start it must still return something, with the certificate
    # that says what it is — never a silent claim.
    eq_b, cert_b = solve_certified(des, [bad]; b = b)
    @test eq_b isa ChemicalState
    @test cert_b isa NamedTuple
    @test haskey(cert_b, :optimal)

    # A single good start behaves exactly like `solve` followed by the certificate.
    eq_g, cert_g = solve_certified(des, [good]; b = b)
    @test cert_g.optimal
    @test ustrip(us"mol", eq_g.n[findfirst(==("Cal"), names)]) ≈
        ustrip(us"mol", eq.n[findfirst(==("Cal"), names)]) rtol = 1.0e-6

end
