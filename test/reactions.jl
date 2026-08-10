using ChemistryLab
using Test

@testsection "Reactions" begin
    @testsection "Constructor from equation" begin
        r = Reaction("H2O = H+ + OH-")
        @test reactants(r)[Species("H2O")] == 1
        @test products(r)[Species("H+")] == 1
        @test products(r)[Species("OH-")] == 1
        @test r.equal_sign == '='
    end

    @testsection "Arrow variants" begin
        # Forward arrow
        r_fwd = Reaction("H2O → H+ + OH-")
        @test r_fwd.equal_sign == '→'
        @test reactants(r_fwd)[Species("H2O")] == 1

        # Equilibrium arrow
        r_eq = Reaction("H2O ↔ H+ + OH-")
        @test r_eq.equal_sign == '↔'

        # Double arrow
        r_dbl = Reaction("H2O ⇌ H+ + OH-")
        @test r_dbl.equal_sign == '⇌'
    end

    @testsection "Getindex and properties" begin
        r = Reaction("CaCO3 = Ca+2 + CO3-2")
        @test r[Species("CaCO3")] == -1
        @test r[Species("Ca+2")] == 1
        r[:testprop] = "abc"
        @test r[:testprop] == "abc"
        @test r.testprop == "abc"
    end

    @testsection "Charge balance" begin
        # Balanced reaction: charge on left == charge on right
        r_balanced = Reaction("H2O = H+ + OH-")
        @test charge(r_balanced) == 0

        # Unbalanced: net charge ≠ 0
        r_unbalanced = Reaction("Fe+2 = Fe+3")
        @test charge(r_unbalanced) != 0
    end

    @testsection "Reaction iteration" begin
        r = Reaction("CaCO3 = Ca+2 + CO3-2")
        all_reac = collect(reactants(r))
        all_prod = collect(products(r))
        @test length(all_reac) == 1
        @test length(all_prod) == 2
        # Each item is a Pair{Species, Number}
        @test all(p -> p isa Pair, all_reac)
        @test all(p -> p isa Pair, all_prod)
    end

    @testsection "Simplify reaction" begin
        reac = Dict(Species("OH⁻") => -1, Species("H2O") => -1)
        prod = Dict(Species("H2O") => 1, Species("H⁺") => 1)
        r = Reaction(reac, prod)
        rs = simplify_reaction(r)
        @test length(reactants(rs)) == 1
        @test length(products(rs)) == 2
        @test haskey(reactants(rs), Species("OH⁻"))
        @test haskey(products(rs), Species("H⁺"))
        # H2O appears only once per side after simplification (not canceled here
        # because it has coeff -1 in reac and +1 in prod, creating net zero which
        # simplify_reaction keeps on the products side with coeff 1)
        @test !haskey(reactants(rs), Species("H2O"))
    end

    @testsection "Addition and subtraction" begin
        r1 = Reaction("H2O = H+ + OH-")
        r2 = Reaction("2H2O = H2 + 2OH⁻")
        rsum = r1 + r2
        @test reactants(rsum)[Species("H2O")] == 3

        rsub = r2 - r1
        # r2 - r1 keeps both H2O terms unsimplified: H2O appears in reactants with
        # coefficient 2 (from r2). Use simplify_reaction to obtain the net coefficient.
        rs = simplify_reaction(rsub)
        @test haskey(reactants(rs), Species("H2O")) || haskey(products(rs), Species("H2O"))
    end

    @testsection "split/merge species by stoich" begin
        s = Dict(Species("CaCl2") => -1, Species("Ca+2") => 1, Species("Cl-") => 2)
        reac, prod = ChemistryLab.split_species_by_stoich(s)
        @test haskey(reac, Species("CaCl2"))
        @test haskey(prod, Species("Ca+2"))
        merged = ChemistryLab.merge_species_by_stoich(reac, prod)
        @test haskey(merged, Species("CaCl2")) && haskey(merged, Species("Ca+2"))
    end

    @testsection "scale_stoich!" begin
        s = Dict(Species("O") => 2, Species("H") => 4)
        ChemistryLab.scale_stoich!(s)
        @test s[Species("O")] == 4 && s[Species("H")] == 8
    end
end
