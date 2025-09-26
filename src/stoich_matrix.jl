same_components(::Vector{<:AbstractSpecies}) = atoms_charge
same_components(::Vector{<:CemSpecies}) = oxides_charge

item_order(::Vector{<:AbstractSpecies}) = ATOMIC_ORDER
item_order(::Vector{<:CemSpecies}) = OXIDE_ORDER

union_atoms(atom_dicts::Vector{<:Dict}, order_vec = ATOMIC_ORDER) = sort!(collect(union(keys.(atom_dicts)...)), by=k -> findfirst(==(k), order_vec))

function print_stoich_matrix(A::AbstractMatrix, indep_comp_names::Vector, dep_comp_names::Vector)
    pretty_table(
        A,
        column_labels=dep_comp_names,
        row_labels=indep_comp_names,
        style=TextTableStyle(; table_border=crayon"green")
    )
end

function stoich_matrix_to_equations(A::AbstractMatrix, indep_comp_names::Vector, dep_comp_names::Vector; scaling=1, display=true, equal_sign='=')
    eqns = String[]
    for (j, sp) in enumerate(dep_comp_names)
        if sp in indep_comp_names
            if display
                println(rpad("$(sp)", 11), "| $sp $(equal_sign) $sp")
            end
        else
            coeffs = Dict(zip(indep_comp_names, -A[:, j]))
            coeffs[sp] = 1
            eqn = format_equation(coeffs; scaling=scaling, equal_sign=equal_sign)
            push!(eqns, eqn)
            if display
                println(rpad("$(sp)", 11), "| ", eqn)
            end
        end
    end
    return eqns
end

function stoich_matrix(species::Vector{<:AbstractSpecies}; display = true)
    involved_atoms_dicts = same_components(species).(species)
    involved_atoms = union_atoms(involved_atoms_dicts, item_order(species))
    T = promote_type(valtype.(involved_atoms_dicts)...)
    A = zeros(T, length(involved_atoms), length(species))
    for (j, atoms) in enumerate(involved_atoms_dicts)
        for (i, atom) in enumerate(involved_atoms)
            A[i, j] = get(atoms, atom, zero(T))
        end
    end
    if display print_stoich_matrix(A, involved_atoms, symbol.(species)) end
    return A, involved_atoms
end

function stoich_matrix(s::Vector{<:AbstractSpecies}, candidate_primaries::Vector{<:AbstractSpecies}; display = true)

    vec_components = same_components(union(s,candidate_primaries))

    S = promote_type(typeof.(s)..., typeof.(candidate_primaries)...)

    species = S[]
    append!(species, s)
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

    M, involved_atoms = stoich_matrix(species; display=false)
    redox = charged && rank(M[:, 1:end-1]; rtol=1.e-6) != rank(M[1:end-1, 1:end-1]; rtol=1.e-6)

    if !redox && charged
        pop!(species)
        M = M[1:end-1, 1:end-1]
    end

    cols_candidates = [findfirst(y -> y == x, species) for x in candidate_primaries]
    filter!(x-> x !== nothing, cols_candidates)
    M_subset = M[:, cols_candidates]
    F = qr(M_subset, Val(true))
    r = Int(rank(M_subset; rtol=1.e-6))
    pivot_idx = F.p[1:r]
    independent_cols_indices = sort(cols_candidates[pivot_idx])
    sort!(independent_cols_indices, by = x->symbol(species[x]) !== "H2O@" && symbol(species[x]) !== "H2O" && symbol(species[x]) !== "H₂O" && symbol(species[x]) !== "H")

    M_indep = M[:, independent_cols_indices]
    A = stoich_coef_round.(pinv(M_indep)*M)

    indep_comp = species[independent_cols_indices]
    dep_comp = species[1:num_initial_species]
    A = A[:, 1:num_initial_species]

    if redox && Zz ∈ dep_comp
        A = A[:, 1:end-1]
        dep_comp = dep_comp[1:end-1]
    end

    if display print_stoich_matrix(A, symbol.(indep_comp), symbol.(dep_comp)) end

    return A, indep_comp, dep_comp 
end

