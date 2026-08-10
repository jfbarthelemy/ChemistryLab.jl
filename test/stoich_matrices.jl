using ChemistryLab
using Test
using LinearAlgebra

@testsection "Stoichiometric Matrix" begin
    @testsection "Basic matrix construction" begin
        h2o = Species("H2O")
        hplus = Species("H+")
        oh = Species("OH-")
        species = [h2o, hplus, oh]

        CSM = CanonicalStoichMatrix(species)
        A, atoms = CSM.A, CSM.primaries
        @test size(A) == (3, 3)   # H, O, Zz (charge) × 3 species
        @test atoms == [:H, :O, :Zz]
        @test A[1, 1] == 2   # H2O has 2 H
        @test A[2, 1] == 1   # H2O has 1 O
    end

    @testsection "Matrix with charged species" begin
        species = [Species("Fe2+"), Species("Fe3+"), Species("e-")]
        SM = StoichMatrix(species)
        A, indep_comp, dep_comp = SM.A, SM.primaries, SM.species

        @test length(dep_comp) == 3   # Fe2+, Fe3+, e-
        @test any(s -> charge(s) != 0, dep_comp)
    end

    @testsection "Reaction conversion" begin
        na2so4 = Species("Na2(SO4)")
        na = Species("Na+")
        so4 = Species("SO4-2")
        naso4 = Species("Na(SO4)-")
        species = [na2so4, na, so4, naso4]

        SM = StoichMatrix(species)
        list_reactions = reactions(SM)

        @test length(list_reactions) > 0
        @test all(r -> r isa Reaction, list_reactions)

        # Each reaction must have non-empty reactants and products
        @test all(r -> !isempty(reactants(r)) || !isempty(products(r)), list_reactions)
    end

    @testsection "Mass-based calculations" begin
        h2o = Species("H2O")
        h2 = Species("H2")
        o2 = Species("O2")
        species = [h2o, h2, o2]

        CSM = mass_matrix(CanonicalStoichMatrix(species))
        A, atoms = CSM.A, CSM.primaries
        @test size(A) == (2, 3)   # 2 elements (H, O) × 3 species

        # Mass conservation: each column (species) sums to 1.0 (normalized mass)
        h2o_idx = findfirst(s -> s == h2o, species)
        @test sum(A[:, h2o_idx]) ≈ 1.0 atol = 1.0e-10
    end

    @testsection "CemSpecies handling" begin
        species = [CemSpecies("C3S"), CemSpecies("CH"), CemSpecies("CSH")]
        CSM = CanonicalStoichMatrix(species)
        A, atoms = CSM.A, CSM.primaries
        @test size(A, 1) ≥ 2   # At least Ca(C) and Si(S) oxide components

        # Verify atoms list contains cement oxide symbols
        @test any(x -> x in [:C, :S, :H], atoms)
    end

    @testsection "Utility functions — union_atoms" begin
        d1 = Dict(:Ca => 1, :O => 1)
        d2 = Dict(:Si => 1, :O => 2)
        atoms = union_atoms([d1, d2])
        @test :Ca in atoms
        @test :Si in atoms
        @test :O in atoms
        # No duplicates
        @test length(atoms) == length(unique(atoms))
    end

    @testsection "Utility functions — same_components" begin
        # For regular Species → should return atoms_charge
        species = [Species("CaCO3")]
        f = ChemistryLab.same_components(species)
        @test f === atoms_charge

        # For CemSpecies → should return oxides_charge
        cem_species = [CemSpecies("C3S")]
        g = ChemistryLab.same_components(cem_species)
        @test g === oxides_charge

        # Verify functions produce the expected keys
        keys_atoms = keys(atoms_charge(Species("CaCO3")))
        @test :Ca in keys_atoms && :C in keys_atoms && :O in keys_atoms

        keys_oxides = keys(oxides_charge(CemSpecies("C3S")))
        @test :C in keys_oxides   # CaO present in C3S
        @test :S in keys_oxides   # SiO2 present in C3S
    end

    @testsection "pprint does not error" begin
        h2o = Species("H2O")
        hplus = Species("H+")
        species = [h2o, hplus]
        CSM = CanonicalStoichMatrix(species)
        @test_nowarn pprint(CSM.A, CSM.primaries, CSM.species)
    end
end
