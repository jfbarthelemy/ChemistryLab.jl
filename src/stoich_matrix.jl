same_components(::Vector{<:AbstractSpecies}) = atoms_charge
same_components(::Vector{<:CemSpecies}) = oxides_charge

item_order(::Vector{<:AbstractSpecies}) = ATOMIC_ORDER
item_order(::Vector{<:CemSpecies}) = OXIDE_ORDER

union_atoms(atom_dicts::Vector{<:Dict}, order_vec = ATOMIC_ORDER) = sort!(collect(union(keys.(atom_dicts)...)), by=k -> findfirst(==(k), order_vec))

function print_stoich_matrix(A::AbstractMatrix, indep_comp_names::Vector, dep_comp_names::Vector)
    hl_p = TextHighlighter(
        (data, i, j) -> (data[i, j] > 0),
        crayon"bold light_red"
    )
    hl_n = TextHighlighter(
        (data, i, j) -> (data[i, j] < 0),
        crayon"bold light_blue"
    )
    hl_z = TextHighlighter(
        (data, i, j) -> (data[i, j] == 0),
        crayon"conceal"
    )
    pretty_table(
        A,
        column_labels = dep_comp_names,
        row_labels    = indep_comp_names,
        highlighters  = [hl_p, hl_n, hl_z],
        style         = TextTableStyle(;
                            row_label = crayon"magenta bold",
                            first_line_column_label = crayon"cyan bold",
                            table_border = crayon"green bold")
    )
end

function stoich_matrix_to_equations(A::AbstractMatrix, indep_comp_names::AbstractVector, dep_comp_names::AbstractVector; scaling=1, display=true, equal_sign='=')
    eqns = String[]
    pad = 11
    for (j, sp) in enumerate(dep_comp_names)
        if sp in indep_comp_names
            if display
                println(rpad("$(sp)", pad), "| $(colored_formula(sp)) $(string(COL_PAR(string(equal_sign)))) $(colored_formula(sp))")
            end
        else
            coeffs = Dict(zip(indep_comp_names, -A[:, j]))
            coeffs[sp] = 1
            eqn = format_equation(coeffs; scaling=scaling, equal_sign=equal_sign)
            push!(eqns, eqn)
            if display
                println(rpad("$(sp)", pad), "| ", colored_equation(eqn))
            end
        end
    end
    return eqns
end

function stoich_matrix_to_reactions(A::AbstractMatrix, indep_comp_names::AbstractVector{<:AbstractSpecies}, dep_comp_names::AbstractVector{<:AbstractSpecies}; scaling=1, display=true, equal_sign='=')
    eqns = Reaction[]
    pad = 11
    for (j, sp) in enumerate(dep_comp_names)
        if sp in indep_comp_names
            if display
                println(rpad("$(symbol(sp))", pad), "| $(colored(sp)) $(string(COL_PAR(string(equal_sign)))) $(colored(sp))")
            end
        else
            coeffs = Dict(zip(indep_comp_names, -A[:, j]))
            coeffs[sp] = 1
            eqn = scaling*Reaction(coeffs)
            push!(eqns, eqn)
            if display
                println(rpad("$(symbol(sp))", pad), "| ", colored(eqn))
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

function stoich_matrix(vs::Vector{<:AbstractSpecies}, candidate_primaries::Vector{<:AbstractSpecies}; display = true, involve_all_atoms = false)

    safe_rank(A; rtol=1e-6) = try rank(A, rtol=rtol) catch; rank(A) end

    all_species = union(vs,candidate_primaries)
    vec_components = same_components(all_species)

    S = promote_type(typeof.(vs)..., typeof.(candidate_primaries)...)

    species = S[]
    append!(species, vs)
    num_initial_species = length(species)
    initial_involved_atoms = involve_all_atoms ? union_atoms(vec_components.(all_species), item_order(species)) : union_atoms(vec_components.(species), item_order(species))
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
    redox = charged && safe_rank(M[:, 1:end-1]) != safe_rank(M[1:end-1, 1:end-1])

    if !redox && charged
        pop!(species)
        M = M[1:end-1, 1:end-1]
    end

    cols_candidates = [findfirst(y -> y == x, species) for x in candidate_primaries]
    filter!(x-> !isnothing(x), cols_candidates)
    M_subset = M[:, cols_candidates]
    if size(M_subset,1) >= size(M_subset, 2)
        independent_cols_indices = cols_candidates
    else
        F = qr(M_subset, Val(true))
        r = Int(safe_rank(M_subset))
        pivot_idx = F.p[1:r]
        independent_cols_indices = sort(cols_candidates[pivot_idx])
    end
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

