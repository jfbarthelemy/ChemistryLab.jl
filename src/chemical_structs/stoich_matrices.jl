# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using LinearAlgebra
using OrderedCollections
using PrettyTables

"""
        StoichMatrix{T,P}

Container holding a stoichiometric matrix `A` together with the
`primaries` (independent components) and the full `species` vector.

# Fields

    - `A`: matrix (components × species) of stoichiometric coefficients.
    - `primaries`: vector of independent components (Symbols or `Species`).
    - `species`: vector of `Species` corresponding to the columns of `A`.
    - `N`: nullspace matrix — columns are vectors of the kernel of `A`.
    Built analytically when primaries ⊆ species (trivial case); otherwise computed via
    RREF over ℚ with integer coefficients (cleared denominators, reduced by GCD).

# Examples

```julia
julia> species = [Species("H₂O"), Species("OH⁻"), Species("H⁺")]
3-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 H⁺ {H⁺} [H⁺ ◆ H+]

julia> SM = StoichMatrix(species)
┌─────┬─────┬─────┬────┐
│     │ H₂O │ OH⁻ │ H⁺ │
├─────┼─────┼─────┼────┤
│ H₂O │   1 │     │  1 │
│ OH⁻ │     │   1 │ -1 │
└─────┴─────┴─────┴────┘
```
"""
struct StoichMatrix{
        T <: Number, P, V <: AbstractVector{P}, M <: AbstractMatrix{T}, S <: AbstractSpecies,
    } <: AbstractMatrix{T}
    A::M
    primaries::V
    species::Vector{S}
    N::M
end

Base.eltype(::StoichMatrix{T}) where {T} = T

primtype(::StoichMatrix{T, P}) where {T, P} = P

# Lossless conversion helpers — avoid the overflow in rationalize(BigInt, ::Rational{Int64}).
_to_qbig(x::Integer) = Rational{BigInt}(BigInt(x))
_to_qbig(x::Rational) = Rational{BigInt}(BigInt(numerator(x)), BigInt(denominator(x)))
# Use tol = 1e-3 to recover the "intended" rational (e.g. 5//3, 21//10) from
# Float64 approximations produced by stoich_coef_round, instead of tol = 1e-9
# which would give huge denominators and break A*N = 0 for the original A.
function _to_qbig(x::AbstractFloat)
    r = rationalize(x; tol = 1.0e-3)
    return Rational{BigInt}(BigInt(numerator(r)), BigInt(denominator(r)))
end
# Fallback for other concrete Number types (e.g. ForwardDiff.Dual): extract Float64 value.
# Symbolic types (Symbolics.Num) are excluded upstream by _is_rationalizable.
_to_qbig(x::Number) = _to_qbig(Float64(x))

# Convert A to a BigInt integer matrix.
# For pure integer input (common case), this is just BigInt.(A).
# For mixed/rational input, clear fractional denominators by multiplying by their LCM.
_to_bigint_matrix(A::AbstractMatrix{<:Integer}) = BigInt.(A)
function _to_bigint_matrix(A::AbstractMatrix)
    B_rat = _to_qbig.(A)
    D = mapreduce(denominator, lcm, B_rat; init = one(BigInt))
    return [numerator(x) * (D ÷ denominator(x)) for x in B_rat]
end

# Exact rational null space via Bareiss algorithm (forward pass) + back-substitution.
#
# Bareiss's theorem guarantees that `÷ prev` is exact ONLY for forward (downward)
# elimination. Clearing rows above the pivot (full RREF) does NOT satisfy this property
# and silently corrupts B via truncation. We therefore do a forward-only pass (upper
# triangular form) and extract null vectors by rational back-substitution.
#
# Returns Matrix{Rational{BigInt}} of shape (n, dim_ker).
function _rational_nullspace(A::AbstractMatrix)
    m, n = size(A)
    B = _to_bigint_matrix(A)
    pivot_cols = Int[]
    prev = one(BigInt)
    row = 1
    # Forward-only Bareiss: produces upper-triangular form with exact BigInt entries.
    for col in 1:n
        i = findfirst(!iszero, B[row:end, col])
        isnothing(i) && continue
        i += row - 1
        i != row && (B[[row, i], :] = B[[i, row], :])
        piv = B[row, col]
        for k in (row + 1):m           # forward only — no backward clearing
            bkc = B[k, col]
            iszero(bkc) && continue
            for j in 1:n
                B[k, j] = (piv * B[k, j] - bkc * B[row, j]) ÷ prev
            end
        end
        push!(pivot_cols, col)
        prev = piv
        row += 1
        row > m && break
    end
    free_cols = setdiff(1:n, pivot_cols)
    isempty(free_cols) && return zeros(Rational{BigInt}, n, 0)
    # Rational back-substitution on the upper-triangular Bareiss matrix.
    # For each free column fc: set v[fc]=1, then solve for v[pivot_cols] back-to-front.
    rk = length(pivot_cols)
    N_rat = zeros(Rational{BigInt}, n, length(free_cols))
    for (j, fc) in enumerate(free_cols)
        v = zeros(Rational{BigInt}, n)
        v[fc] = one(Rational{BigInt})
        for r in rk:-1:1
            pc = pivot_cols[r]
            s = -Rational{BigInt}(B[r, fc])
            for r2 in (r + 1):rk
                s -= Rational{BigInt}(B[r, pivot_cols[r2]]) * v[pivot_cols[r2]]
            end
            v[pc] = s / Rational{BigInt}(B[r, pc])
        end
        N_rat[:, j] = v
    end
    return N_rat
end

# ── Stage 2: conversion from rational nullspace to concrete types ─────────────

# Convert Rational{BigInt} → concrete Real without precision loss where possible.
# Int for integers, Rational{Int} for small denominators, Float64 as last resort.
function _to_stoich_real(r::Rational{BigInt})
    n, d = numerator(r), denominator(r)
    if isone(d)
        return typemin(Int) <= n <= typemax(Int) ? Int(n) : stoich_coef_round(Float64(r))
    elseif d < 10 && typemin(Int) <= n <= typemax(Int) && typemin(Int) <= d
        return Rational{Int}(Int(n), Int(d))
    else
        return stoich_coef_round(Float64(r))
    end
end

# Detect whether all entries in a column are (near-)integer.
_is_int_entry(x::Integer) = true
_is_int_entry(x::Rational) = isone(denominator(x))
_is_int_entry(x::AbstractFloat) = isapprox(x, round(x); atol = 1.0e-9)

# Check whether a matrix element type supports exact rational conversion.
# Symbolic types (Symbolics.Num <: Real) are excluded — their nullspace
# cannot be computed numerically and must remain empty.
_is_rationalizable(::AbstractMatrix{<:Union{Integer, AbstractFloat, Rational}}) = true
_is_rationalizable(::AbstractMatrix) = false

# Clear all denominators + GCD + sign normalize on a pre-computed rational nullspace.
# Returns Matrix{BigInt}.
function _integer_from_rational(N_rat::Matrix{Rational{BigInt}})
    n, k = size(N_rat)
    k == 0 && return zeros(BigInt, n, 0)
    N_out = zeros(BigInt, n, k)
    for j in 1:k
        v = N_rat[:, j]
        d = mapreduce(denominator, lcm, v; init = one(BigInt))
        vi = numerator.(v .* d)
        g = mapreduce(abs, gcd, vi; init = zero(BigInt))
        iszero(g) && continue
        vi = vi .÷ g
        lead = findfirst(!iszero, vi)
        !isnothing(lead) && vi[lead] < 0 && (vi = -vi)
        N_out[:, j] = vi
    end
    return N_out
end

# Smart conversion: integer where A is integer, fractional otherwise.
# Returns Matrix{Real} with concrete subtypes (Int, Rational{Int}, Float64) — never BigInt.
function _optimal_from_rational(N_rat::Matrix{Rational{BigInt}}, A::AbstractMatrix)
    n, k = size(N_rat)
    k == 0 && return zeros(Int, n, 0)
    is_int_col = [all(_is_int_entry, view(A, :, j)) for j in axes(A, 2)]
    N_out = Matrix{Real}(undef, n, k)
    for j in 1:k
        v = N_rat[:, j]
        int_idx = [i for i in 1:n if is_int_col[i] && !iszero(v[i])]
        d = if !isempty(int_idx)
            mapreduce(i -> denominator(v[i]), lcm, int_idx; init = one(BigInt))
        else
            mapreduce(denominator, lcm, v; init = one(BigInt))
        end
        v = v .* d
        g = if !isempty(int_idx)
            mapreduce(i -> abs(numerator(v[i])), gcd, int_idx; init = zero(BigInt))
        else
            mapreduce(x -> abs(numerator(x)), gcd, v; init = zero(BigInt))
        end
        !iszero(g) && g > 1 && (v = v .// g)
        lead = findfirst(!iszero, v)
        !isnothing(lead) && v[lead] < 0 && (v = -v)
        for i in 1:n
            N_out[i, j] = _to_stoich_real(v[i])
        end
    end
    return N_out
end

# Thin wrappers: compute rational nullspace then convert.
_integer_nullspace(A::AbstractMatrix) = _integer_from_rational(_rational_nullspace(A))
_optimal_nullspace(A::AbstractMatrix) = _optimal_from_rational(_rational_nullspace(A), A)

# ── Stage 3: kinetic species diagonalization ─────────────────────────────────

# Diagonalize rows of the rational nullspace corresponding to kinetic species.
# After this, each kinetic species row has exactly one non-zero entry.
# In-place column operations preserve A*N = 0 exactly (Rational{BigInt} arithmetic).
function _diagonalize_kinetic_rows!(
        N_rat::Matrix{Rational{BigInt}}, kinetic_indices::AbstractVector{Int}
    )
    claimed = Set{Int}()
    for k in kinetic_indices
        # Find an unclaimed pivot column where N_rat[k, j] ≠ 0
        j = nothing
        for c in axes(N_rat, 2)
            c ∈ claimed && continue
            iszero(N_rat[k, c]) && continue
            j = c
            break
        end
        isnothing(j) && throw(
            ArgumentError(
                "Kinetic species at index $k does not participate in any " *
                    "remaining independent reaction (row is zero in unclaimed columns)."
            )
        )
        push!(claimed, j)
        # Eliminate row k from all other columns
        for c in axes(N_rat, 2)
            c == j && continue
            iszero(N_rat[k, c]) && continue
            factor = N_rat[k, c] // N_rat[k, j]
            N_rat[:, c] .-= factor .* N_rat[:, j]
        end
    end
    return nothing
end

# Resolve kinetic species names or objects to row indices in the species vector.
function _resolve_kinetic_indices(
        kinetic_species::AbstractVector{<:AbstractSpecies},
        species::AbstractVector{<:AbstractSpecies},
    )
    return [
        begin
            idx = findfirst(s -> s == ks, species)
            isnothing(idx) &&
                throw(ArgumentError("Kinetic species \"$(symbol(ks))\" not found in species list."))
            idx
        end for ks in kinetic_species
    ]
end

function _resolve_kinetic_indices(
        kinetic_species::AbstractVector{<:AbstractString},
        species::AbstractVector{<:AbstractSpecies},
    )
    return [
        begin
            idx = findfirst(s -> symbol(s) == name, species)
            isnothing(idx) &&
                throw(ArgumentError("Kinetic species \"$name\" not found in species list."))
            idx
        end for name in kinetic_species
    ]
end

function StoichMatrix(
        A::M, primaries::Union{Vector{Symbol}, Vector{S}}, species::AbstractVector{S};
        kinetic_species = nothing,
    ) where {M <: AbstractMatrix, S <: AbstractSpecies}
    indices_in = [findfirst(x -> x == p, species) for p in primaries]
    if any(x -> isnothing(x), indices_in)
        # General case: primaries are not contained in species (e.g., primaries are atoms).
        if _is_rationalizable(A)
            N_rat = _rational_nullspace(A)
            if !isnothing(kinetic_species)
                kin_idx = _resolve_kinetic_indices(kinetic_species, species)
                _diagonalize_kinetic_rows!(N_rat, kin_idx)
            end
            N_opt = _optimal_from_rational(N_rat, A)
            N = similar(A, size(N_opt)...)
            N .= N_opt
        else
            # Symbolic matrix (e.g. Symbolics.Num entries): nullspace is undefined
            # numerically — return empty nullspace as in the analytical branch.
            N = similar(A, 0, 0)
        end
        return StoichMatrix(A, primaries, Vector(species), N)
    else
        # Trivial case: primaries ⊆ species and A[:,indices_in] = I — analytical formula.
        # Each non-primary species already has exactly one column → diagonalization
        # is automatic. Just validate kinetic species are not primaries.
        if !isnothing(kinetic_species)
            kin_idx = _resolve_kinetic_indices(kinetic_species, species)
            for k in kin_idx
                k ∈ indices_in && throw(
                    ArgumentError(
                        "Kinetic species \"$(symbol(species[k]))\" cannot be a primary. " *
                            "Primaries appear in every reaction (identity block of A)."
                    )
                )
            end
        end
        p, n = size(A)
        N = similar(A, n, n - p)
        indices_out = setdiff(eachindex(species), indices_in)
        N[indices_out, :] .= Diagonal(ones(eltype(A), length(indices_out)))
        N[indices_in, :] .= -A[:, indices_out]
        return StoichMatrix(A, primaries, Vector(species), N)
    end
end

for f in (:size, :length, :getindex, :ndims, :iterate)
    @eval Base.$f(SM::StoichMatrix, args...; kwargs...) = $f(SM.A, args...; kwargs...)
end

Base.show(io::IO, SM::StoichMatrix) = show(io, SM.A)

function Base.show(::IO, ::MIME"text/plain", SM::StoichMatrix)
    (; A, primaries, species) = SM
    column_labels = try
        symbol.(species)
    catch
        species
    end
    row_labels = try
        symbol.(primaries)
    catch
        primaries
    end
    formatters = [(v, i, j) -> iszero(v) ? "" : AnsiTextCell(string(v))]
    pretty_table(
        A; column_labels = column_labels, row_labels = row_labels, formatters = formatters
    )
    return nothing
end

function apply(f::Function, SM::StoichMatrix, args...; kwargs...)
    return StoichMatrix(f(SM.A, args...; kwargs...), SM.primaries, SM.species)
end

"""
    pprint(A::AbstractMatrix, indep_comp_names::AbstractVector, dep_comp_names::AbstractVector)

Print a stoichiometric matrix with colored formatting.

# Arguments

  - `A`: stoichiometric matrix.
  - `indep_comp_names`: row labels (independent components).
  - `dep_comp_names`: column labels (dependent components).

Uses text highlighters to color positive (red), negative (blue), and zero (concealed) values.
"""
function pprint(
        A::AbstractMatrix,
        indep_comp_names::AbstractVector,
        dep_comp_names::AbstractVector;
        row_label = :symbol,
        col_label = :symbol,
        label = :identity,
        kwargs...,
    )
    if label != :identity
        row_label = label
        col_label = label
    end
    column_labels = try
        eval(col_label).(dep_comp_names)
    catch
        dep_comp_names
    end
    row_labels = try
        eval(col_label).(indep_comp_names)
    catch
        indep_comp_names
    end
    hl_p = TextHighlighter((data, i, j) -> (data[i, j] > 0), crayon"bold light_red")
    hl_n = TextHighlighter((data, i, j) -> (data[i, j] < 0), crayon"bold light_blue")
    formatters = [(v, i, j) -> iszero(v) ? "" : AnsiTextCell(string(v))]
    return try
        pretty_table(
            A;
            column_labels = column_labels,
            row_labels = row_labels,
            formatters = formatters,
            highlighters = [hl_p, hl_n],
            style = TextTableStyle(;
                row_label = crayon"magenta bold",
                first_line_column_label = crayon"cyan bold",
                table_border = crayon"green bold",
            ),
            kwargs...,
        )
    catch
        pretty_table(
            A;
            column_labels = column_labels,
            row_labels = row_labels,
            formatters = formatters,
            style = TextTableStyle(;
                row_label = crayon"magenta bold",
                first_line_column_label = crayon"cyan bold",
                table_border = crayon"green bold",
            ),
            kwargs...,
        )
    end
end

function pprint(
        SM::StoichMatrix; row_label = :symbol, col_label = :symbol, label = :identity, kwargs...
    )
    (; A, primaries, species) = SM
    return pprint(
        A,
        primaries,
        species;
        row_label = row_label,
        col_label = col_label,
        label = label,
        kwargs...,
    )
end

"""
        CanonicalStoichMatrix(species)

Construct a StoichMatrix from a list of species i.e. the matrix giving the stoichiometric
coefficients of the species in columns with respect to the atoms in rows.

# Arguments

  - `species`: list of species.

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> CSM = CanonicalStoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌────┬─────┬────┬─────┬─────┬───────┬───────┐
│    │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├────┼─────┼────┼─────┼─────┼───────┼───────┤
│  C │     │    │     │   1 │     1 │     1 │
│  H │   2 │  1 │   1 │     │     1 │       │
│  O │   1 │    │   1 │   2 │     3 │     3 │
│ Zz │     │  1 │  -1 │     │    -1 │    -2 │
└────┴─────┴────┴─────┴─────┴───────┴───────┘
```
"""
function CanonicalStoichMatrix(species::AbstractVector{<:AbstractSpecies})
    involved_atoms_dicts = same_components(species).(species)
    involved_atoms = union_atoms(involved_atoms_dicts, item_order(species))
    T = promote_type(valtype.(involved_atoms_dicts)...)
    A = zeros(T, length(involved_atoms), length(species))
    for (j, atoms) in enumerate(involved_atoms_dicts)
        for (i, atom) in enumerate(involved_atoms)
            A[i, j] = get(atoms, atom, zero(T))
        end
    end
    if _is_rationalizable(A)
        N_int = _integer_nullspace(A)
        N = similar(A, size(N_int)...)
        N .= N_int
    else
        N = similar(A, 0, 0)
    end
    return StoichMatrix(A, involved_atoms, Vector(species), N)
end

"""
        StoichMatrix(species, candidate_primaries=species; involve_all_atoms=true)

Construct a StoichMatrix from a list of species and a list of candidate primary species
(by default the list of species itself).

# Arguments

  - `species`: list of species.
  - `candidate_primaries`: list of candidate primary species (default: `species`).
  - `involve_all_atoms`: if true the algorithm is allowed to use species
    of `candidate_primaries` containing atoms which are not in `species` (default: true).

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │   1 │     │     1 │     1 │
│  H⁺ │     │  1 │  -1 │     │    -1 │    -2 │
│ CO₂ │     │    │     │   1 │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> SM.N # nullspace (columns are vectors of the kernel of the stoichiometric matrix)
6×3 Matrix{Int64}:
 -1  -1  -1
  1   1   2
  1   0   0
  0  -1  -1
  0   1   0
  0   0   1

julia> SM.A * SM.N
3×3 Matrix{Int64}:
 0  0  0
 0  0  0
 0  0  0
```
"""
function StoichMatrix(
        species::AbstractVector{<:AbstractSpecies},
        candidate_primaries::AbstractVector{<:AbstractSpecies} = species;
        involve_all_atoms = true,
        optimize_primaries = false,
        kinetic_species = nothing,
    )
    safe_rank(A; rtol = 1.0e-6) =
    try
        rank(A; rtol = rtol)
    catch
        try
            rank(A)
        catch
            min(size(A)...)
        end
    end
    safe_pinv(A) =
    try
        pinv(A)
    catch
        inv(A)
    end

    all_species = union(species, candidate_primaries)
    vec_components = same_components(all_species)

    S = promote_type(typeof.(species)..., typeof.(candidate_primaries)...)

    newspecies = S[]
    append!(newspecies, species)
    num_initial_species = length(newspecies)
    initial_involved_atoms = if involve_all_atoms
        union_atoms(vec_components.(all_species), item_order(newspecies))
    else
        union_atoms(vec_components.(newspecies), item_order(newspecies))
    end
    candidate_primaries = copy(candidate_primaries)

    SpType(::AbstractVector) = Species
    SpType(::AbstractVector{<:CemSpecies}) = CemSpecies
    Zz = SpType(newspecies)("Zz")
    charged = :Zz ∈ initial_involved_atoms
    if charged
        if Zz ∉ newspecies
            push!(newspecies, Zz)
            num_initial_species += 1
        end
        if Zz ∉ candidate_primaries
            push!(candidate_primaries, Zz)
        end
    end

    for x in candidate_primaries
        idx = findfirst(y -> x == y, newspecies)
        if isnothing(idx) &&
                all(k -> first(k) ∈ initial_involved_atoms || first(k) == :Zz, vec_components(x))
            push!(newspecies, x)
        end
    end

    CSM = CanonicalStoichMatrix(newspecies)
    M, involved_atoms = CSM.A, CSM.primaries
    redox =
        charged &&
        safe_rank(M[:, begin:end .!= num_initial_species]) !=
        safe_rank(M[1:(end - 1), begin:end .!= num_initial_species])

    if !redox && charged
        deleteat!(newspecies, num_initial_species)
        M = M[1:(end - 1), 1:end .!= num_initial_species]
        num_initial_species -= 1
    end

    cols_candidates = [findfirst(y -> y == x, newspecies) for x in candidate_primaries]
    cols_candidates = cols_candidates[.!isnothing.(cols_candidates)]
    # Exclude kinetic species from primary selection — they must stay dependent
    # so they can each be isolated to a single nullspace column.
    if !isnothing(kinetic_species)
        kin_idx_in_new = _resolve_kinetic_indices(kinetic_species, newspecies)
        filter!(c -> c ∉ kin_idx_in_new, cols_candidates)
    end
    M_subset = M[:, cols_candidates]

    r = Int(safe_rank(M_subset))
    if optimize_primaries
        F = qr(M_subset, Val(true))
        pivot_idx = F.p[1:r]
        independent_cols_indices = sort(cols_candidates[pivot_idx])
    else
        pivot_idx = Int[]
        current_matrix = zeros(eltype(M_subset), size(M_subset, 1), 0)
        for j in axes(M_subset, 2)
            candidate = hcat(current_matrix, M_subset[:, j])
            if Int(safe_rank(candidate)) > size(current_matrix, 2)
                push!(pivot_idx, j)
                current_matrix = candidate
                if length(pivot_idx) == r
                    break
                end
            end
        end
        independent_cols_indices = cols_candidates[pivot_idx]
    end

    sort!(
        independent_cols_indices;
        by = x ->
        symbol(newspecies[x]) !== "H2O@" &&
            symbol(newspecies[x]) !== "H2O" &&
            symbol(newspecies[x]) !== "H₂O" &&
            symbol(newspecies[x]) !== "H",
    )
    M_indep = M[:, independent_cols_indices]
    M_indep = promote_type(typeof.(M_indep)...).(M_indep)
    A = stoich_coef_round.(safe_pinv(M_indep) * M)

    indep_comp = newspecies[independent_cols_indices]
    dep_comp = newspecies[1:num_initial_species]
    A = A[:, 1:num_initial_species]

    if redox && Zz ∈ dep_comp # && Zz ∉ indep_comp
        A = A[:, 1:(end - 1)]
        dep_comp = dep_comp[1:(end - 1)]
        num_initial_species -= 1
    end

    zero_rows = all(iszero.(A); dims = 2)[:, 1]
    A = A[.!zero_rows, :]
    indep_comp = indep_comp[.!zero_rows]

    return StoichMatrix(A, indep_comp, dep_comp; kinetic_species = kinetic_species)
end

function StoichMatrix(
        species::AbstractDict,
        candidate_primaries = species;
        involve_all_atoms = true,
        optimize_primaries = false,
        kinetic_species = nothing,
    )
    gather_species(d::AbstractDict{T, S} where {T, S <: AbstractSpecies}) = collect(values(d))
    gather_species(d::AbstractDict{S, T} where {S <: AbstractSpecies, T}) = collect(keys(d))
    gather_species(d) = d
    return StoichMatrix(
        gather_species(species),
        gather_species(candidate_primaries);
        involve_all_atoms = involve_all_atoms,
        optimize_primaries = optimize_primaries,
        kinetic_species = kinetic_species,
    )
end

function StoichMatrix(
        species,
        candidate_primaries = species;
        involve_all_atoms = true,
        optimize_primaries = false,
        kinetic_species = nothing,
    )
    return StoichMatrix(
        collect(species),
        collect(candidate_primaries);
        involve_all_atoms = involve_all_atoms,
        optimize_primaries = optimize_primaries,
        kinetic_species = kinetic_species,
    )
end

"""
        pull_primaries(SM::StoichMatrix)

Construct a StoichMatrix by reordering the species such that the primaries appear first
(works only if primaries are contained within species).

# Arguments

  - `SM`: StoichMatrix.

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │   1 │     │     1 │     1 │
│  H⁺ │     │  1 │  -1 │     │    -1 │    -2 │
│ CO₂ │     │    │     │   1 │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> pull_primaries(SM)
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ CO₂ │ OH⁻ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │     │   1 │     1 │     1 │
│  H⁺ │     │  1 │     │  -1 │    -1 │    -2 │
│ CO₂ │     │    │   1 │     │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> pull_primaries(SM).N
6×3 Matrix{Int64}:
 -1  -1  -1
  1   1   2
  0  -1  -1
  1   0   0
  0   1   0
  0   0   1
```
"""
function pull_primaries(SM::StoichMatrix)
    (; A, primaries, species) = SM
    indices_in = [findfirst(x -> x == p, species) for p in primaries]
    if any(x -> isnothing(x), indices_in)
        return SM
    else
        indices_out = setdiff(eachindex(species), indices_in)
        indices = [indices_in; indices_out]
        return StoichMatrix(
            A[:, indices], primaries, species[indices], [-A[:, indices_out]; I]
        )
    end
end

"""
        push_primaries(SM::StoichMatrix)

Construct a StoichMatrix by reordering the species such that the primaries appear at the end
(works only if primaries are contained within species).

# Arguments

  - `SM`: StoichMatrix.

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │   1 │     │     1 │     1 │
│  H⁺ │     │  1 │  -1 │     │    -1 │    -2 │
│ CO₂ │     │    │     │   1 │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> push_primaries(SM)
┌─────┬─────┬───────┬───────┬─────┬────┬─────┐
│     │ OH⁻ │ HCO₃⁻ │ CO₃²⁻ │ H₂O │ H⁺ │ CO₂ │
├─────┼─────┼───────┼───────┼─────┼────┼─────┤
│ H₂O │   1 │     1 │     1 │   1 │    │     │
│  H⁺ │  -1 │    -1 │    -2 │     │  1 │     │
│ CO₂ │     │     1 │     1 │     │    │   1 │
└─────┴─────┴───────┴───────┴─────┴────┴─────┘

julia> push_primaries(SM).N
6×3 Matrix{Int64}:
  1   0   0
  0   1   0
  0   0   1
 -1  -1  -1
  1   1   2
  0  -1  -1
```
"""
function push_primaries(SM::StoichMatrix)
    (; A, primaries, species) = SM
    indices_in = [findfirst(x -> x == p, species) for p in primaries]
    if any(x -> isnothing(x), indices_in)
        return SM
    else
        indices_out = setdiff(eachindex(species), indices_in)
        indices = [indices_out; indices_in]
        return StoichMatrix(
            A[:, indices], primaries, species[indices], typeof(A)([I; -A[:, indices_out]])
        )
    end
end

"""
        mass_matrix(SM::StoichMatrix)

Construct a StoichMatrix with mass correspondence instead of molar stoichiometry.

# Arguments

  - `SM`: StoichMatrix.

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │   1 │     │     1 │     1 │
│  H⁺ │     │  1 │  -1 │     │    -1 │    -2 │
│ CO₂ │     │    │     │   1 │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> mass_matrix(SM)
┌─────┬─────┬─────┬────────────┬─────┬────────────┬────────────┐
│     │ H₂O │  H⁺ │        OH⁻ │ CO₂ │      HCO₃⁻ │      CO₃²⁻ │
├─────┼─────┼─────┼────────────┼─────┼────────────┼────────────┤
│ H₂O │ 1.0 │     │    1.05927 │     │    0.29525 │    0.30021 │
│  H⁺ │     │ 1.0 │ -0.0592697 │     │ -0.0165203 │ -0.0335955 │
│ CO₂ │     │     │            │ 1.0 │    0.72127 │   0.733386 │
└─────┴─────┴─────┴────────────┴─────┴────────────┴────────────┘

julia> CSM = CanonicalStoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌────┬─────┬────┬─────┬─────┬───────┬───────┐
│    │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├────┼─────┼────┼─────┼─────┼───────┼───────┤
│  C │     │    │     │   1 │     1 │     1 │
│  H │   2 │  1 │   1 │     │     1 │       │
│  O │   1 │    │   1 │   2 │     3 │     3 │
│ Zz │     │  1 │  -1 │     │    -1 │    -2 │
└────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> mass_matrix(CSM)
┌────┬──────────┬─────┬───────────┬──────────┬───────────┬──────────┐
│    │      H₂O │  H⁺ │       OH⁻ │      CO₂ │     HCO₃⁻ │    CO₃²⁻ │
├────┼──────────┼─────┼───────────┼──────────┼───────────┼──────────┤
│  C │          │     │           │ 0.272921 │   0.19685 │ 0.200157 │
│  H │ 0.111907 │ 1.0 │ 0.0592697 │          │ 0.0165203 │          │
│  O │ 0.888093 │     │   0.94073 │ 0.727079 │   0.78663 │ 0.799843 │
│ Zz │          │     │           │          │           │          │
└────┴──────────┴─────┴───────────┴──────────┴───────────┴──────────┘
```
"""
function mass_matrix(SM::StoichMatrix{T, S}) where {T, S <: AbstractSpecies}
    Mspecies = ustrip.(getproperty.(SM.species, :M))
    Mprimaries = ustrip.(getproperty.(SM.primaries, :M))
    return StoichMatrix(Mprimaries .* SM.A .* inv.(Mspecies)', SM.primaries, SM.species)
end

function mass_matrix(SM::StoichMatrix{T, Symbol}) where {T}
    Mspecies = ustrip.(getproperty.(SM.species, :M))
    S = promote_type(root_type.(typeof.(SM.species))...)
    Mprimaries = ustrip.(getproperty.(S.(string.(SM.primaries)), :M))
    return StoichMatrix(Mprimaries .* SM.A .* inv.(Mspecies)', SM.primaries, SM.species)
end

# Internal: build reactions from nullspace columns (general case).
# Each column of N gives stoichiometric coefficients for the species;
# negative = reactant, positive = product. Symbol is auto-detected.
function _reactions_from_nullspace(
        species::AbstractVector{<:AbstractSpecies}, N::AbstractMatrix
    )
    rxns = Reaction[]
    for V in eachcol(N)
        all(iszero, V) && continue
        dict = OrderedDict(sp => v for (sp, v) in zip(species, V) if !iszero(v))
        isempty(dict) && continue
        push!(rxns, Reaction(dict; symbol = nothing))
    end
    return unique!(rxns)
end

"""
        reactions(SM::StoichMatrix)

Construct the list of non-trivial reactions from a StoichMatrix whose primaries
are species. Uses `push_primaries` when primaries are contained in species;
falls back to direct nullspace extraction otherwise.

# Arguments

  - `SM`: StoichMatrix.

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │   1 │     │     1 │     1 │
│  H⁺ │     │  1 │  -1 │     │    -1 │    -2 │
│ CO₂ │     │    │     │   1 │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> reactions(SM)
3-element Vector{Reaction{Species{Int64}, Int64, Species{Int64}, Int64, Int64}}:
 H₂O = OH⁻ + H⁺
 H₂O + CO₂ = HCO₃⁻ + H⁺
 H₂O + CO₂ = CO₃²⁻ + 2H⁺
```
"""
# Species-primary matrices: use push_primaries when possible, fall back to general extraction.
function reactions(SM::StoichMatrix)
    if !isempty(SM.N)
        pSM = push_primaries(SM)
        if pSM === SM
            # push_primaries could not reorder — fall back to general extraction
            return _reactions_from_nullspace(SM.species, SM.N)
        end
        return [
            Reaction(OrderedDict(zip(pSM.species, V)); symbol = symbol(pSM.species[j])) for
                (j, V) in enumerate(eachcol(pSM.N))
        ]
    else
        lr = unique!(
            [
                Reaction(
                    merge(
                        +, OrderedDict(SM.species[j] => 1), OrderedDict(zip(SM.primaries, -V))
                    );
                    symbol = symbol(SM.species[j]),
                ) for (j, V) in enumerate(eachcol(SM.A))
            ]
        )
        return lr[.!isempty.(lr)]
    end
end

"""
        pprint(reactions::AbstractVector{<:Reaction})

Pretty print a list of reactions.

# Arguments

  - `SM`: StoichMatrix.

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │   1 │     │     1 │     1 │
│  H⁺ │     │  1 │  -1 │     │    -1 │    -2 │
│ CO₂ │     │    │     │   1 │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> pprint(reactions(SM))
  OH⁻ │ H₂O = OH⁻ + H⁺
HCO₃⁻ │ H₂O + CO₂ = HCO₃⁻ + H⁺
CO₃²⁻ │ H₂O + CO₂ = CO₃²⁻ + 2H⁺
```
"""
function pprint(reactions::AbstractVector{<:Reaction}; kwargs...)
    pad = maximum(length.(symbol.(reactions)))
    for r in reactions
        println(lpad("$(symbol(r))", pad), " │ ", colored(r))
    end
    return
end

const oxides_as_species = [Species(d; symbol = string(k)) for (k, d) in CEMENT_TO_MENDELEEV]

const Aoxides, atoms_in_oxides =
    getfield.(Ref(CanonicalStoichMatrix(oxides_as_species)), [:A, :primaries])

const order_atom_in_oxides = Dict(atom => i for (i, atom) in enumerate(atoms_in_oxides))
