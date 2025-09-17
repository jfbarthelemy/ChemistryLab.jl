same_components(::Vector) = atoms_charge
same_components(::Vector{<:CemSpecies}) = oxides_charge

item_order(::Vector) = ATOMIC_ORDER
item_order(::Vector{<:CemSpecies}) = OXIDE_ORDER

union_atoms(atom_dicts::Vector{<:Dict}, order_vec = ATOMIC_ORDER) = sort!(collect(union(keys.(atom_dicts)...)), by=k -> findfirst(==(k), order_vec))

function stoich_matrix(species::Vector; display = true)
    involved_atoms_dicts = same_components(species).(species)
    involved_atoms = union_atoms(involved_atoms_dicts, item_order(species))
    T = promote_type(valtype.(involved_atoms_dicts)...)
    stoich_matrix = zeros(T, length(involved_atoms), length(species))
    for (j, atoms) in enumerate(involved_atoms_dicts)
        for (i, atom) in enumerate(involved_atoms)
            stoich_matrix[i, j] = get(atoms, atom, zero(T))
        end
    end
    species_names = name.(species)
    if display
        pretty_table(
           stoich_matrix,
           column_labels = species_names,
           row_labels = involved_atoms,
           style= TextTableStyle(; table_border = crayon"green")
           )
    end
    return stoich_matrix, involved_atoms, species_names
end

function stoich_matrix(species::Vector, candidate_primaries::Vector; display = true)

    vec_components = same_components(union(species,candidate_primaries))

    species = deepcopy(species)
    num_initial_species = length(species)
    initial_involved_atoms = union_atoms(vec_components.(species), item_order(species))
    candidate_primaries = deepcopy(candidate_primaries)

    for x in candidate_primaries
        if x ∉ species && all(k -> first(k) ∈ initial_involved_atoms || first(k) == :Zz, vec_components(x))
            push!(species, x)
        end
    end

    initial_involved_atoms = union_atoms(vec_components.(species), item_order(species))

    SpType(::Vector) = Species
    SpType(::Vector{<:CemSpecies}) = CemSpecies
    Zz = SpType(species)("Zz")
    charged = :Zz ∈ initial_involved_atoms
    if charged
        if Zz ∉ species push!(species, Zz) end
        if Zz ∉ candidate_primaries push!(candidate_primaries, Zz) end
    end

    M, involved_atoms, species_names = stoich_matrix(species; display = false)
    redox = charged && rank(M[:, 1:end-1]; rtol=1.e-6) != rank(M[1:end-1, 1:end-1]; rtol=1.e-6)

    if !redox && charged
        pop!(species)
        pop!(species_names)
        M = M[1:end-1, 1:end-1]
    end

    cols_candidates = [findfirst(==(name(x)), species_names) for x in candidate_primaries]
    filter!(x-> x !== nothing, cols_candidates)
    M_subset = M[:, cols_candidates]
    F = qr(M_subset, Val(true))
    r = Int(rank(M_subset; rtol=1.e-6))
    pivot_idx = F.p[1:r]
    independent_cols_indices = sort(cols_candidates[pivot_idx])
    sort!(independent_cols_indices, by = x->species_names[x] !== "H2O@" && species_names[x] !== "H2O" && species_names[x] !== "H")

    M_indep = M[:, independent_cols_indices]
    A = stoich_coef_round.(pinv(M_indep)*M)

    indep_comp = species_names[independent_cols_indices]
    dep_comp = species_names[1:num_initial_species]
    A = A[:, 1:num_initial_species]

    if redox && "Zz" ∈ dep_comp
        A = A[:, 1:end-1]
        dep_comp = dep_comp[1:end-1]
    end

    if display
        pretty_table(
           A,
           column_labels = dep_comp,
           row_labels = indep_comp,
           style= TextTableStyle(; table_border = crayon"green")
           )
    end

    return A, indep_comp, dep_comp 
end
