# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
#
# Volume fractions of a ChemicalState, and the chemical-shrinkage void that
# closes the balance in a sealed system.
#
# The bookkeeping is the one of Lavergne et al. (2018), Eqs. (8)-(9): all
# amounts refer to a FIXED reference volume of fresh material. Hydration
# reactions consume more volume than they create, so the sum of the species
# volumes falls below that reference; the deficit is a gas-filled void, and
# ignoring it makes the fractions sum to less than one.

using OrderedCollections

# ── chemical_shrinkage ───────────────────────────────────────────────────────

"""
    chemical_shrinkage(state::ChemicalState, state₀::ChemicalState) -> Quantity

Volume lost to the chemical reactions between the reference state `state₀` and
`state`, ``δV = V_0 - V``, positive for the usual contraction (Le Chatelier).

Both states must share the same [`ChemicalSystem`](@ref). Species without a
standard molar volume `V⁰` contribute nothing, exactly as in
[`volume`](@ref) — see [`missing_molar_volumes`](@ref) to find out which those
are before trusting the number.

# Examples

```julia
δV = chemical_shrinkage(state_28d, state_0)
ustrip(δV / volume(state_0).total)          # relative contraction, dimensionless
```

See also: [`volume_fractions`](@ref), [`missing_molar_volumes`](@ref).
"""
function chemical_shrinkage(state::ChemicalState, state₀::ChemicalState)
    state.system === state₀.system || state.system == state₀.system ||
        throw(ArgumentError("both states must refer to the same ChemicalSystem"))
    return volume(state₀).total - volume(state).total
end

# ── missing_molar_volumes ────────────────────────────────────────────────────

"""
    missing_molar_volumes(state::ChemicalState; atol = 0.0) -> Vector{String}

Symbols of the species that are present in `state` (amount in mol strictly above
`atol`) but carry no standard molar volume `V⁰`, and are therefore invisible to
[`volume`](@ref), [`porosity`](@ref) and [`volume_fractions`](@ref).

A non-empty result means the volume balance is incomplete. Databases differ in
coverage — CEMDATA18 supplies `V⁰` for the cement phases but a species added by
hand, or one taken from an aqueous-only dataset, may not have it.

# Examples

```julia
isempty(missing_molar_volumes(state)) || @warn "incomplete volume balance" missing_molar_volumes(state)
```
"""
function missing_molar_volumes(state::ChemicalState; atol::Real = 0.0)
    out = String[]
    for (i, sp) in enumerate(state.system.species)
        _has_molar_volume(sp) && continue
        _primal(ustrip(us"mol", state.n[i])) > atol || continue
        push!(out, symbol(sp))
    end
    return out
end

# ── volume_fractions ─────────────────────────────────────────────────────────

"""
    volume_fractions(state::ChemicalState; reference = nothing, void_key = "void")
        -> OrderedDict{String, Float64}

Volume fraction of every species carrying a standard molar volume, keyed by
species symbol. Species with a zero amount are omitted.

An individual fraction may be **negative**: aqueous solutes have negative
standard partial molar volumes, electrostriction contracting the solvent around
the ion (`V⁰(OH⁻) ≈ -4.7 cm³/mol`). Those contributions are kept, so that this
function and [`volume`](@ref) always agree; they are negligible at the
concentrations of a pore solution, and grouping folds them back into the liquid.

Two normalizations, selected by `reference`:

  - `reference = nothing` (default) — fractions are relative to the *current*
    total volume of the state, and sum to 1 by construction. This is the right
    choice for a closed, fully-saturated system.
  - `reference::ChemicalState` — fractions are relative to
    `volume(reference).total`, held fixed. They then sum to less than 1, and the
    deficit is returned under `void_key` as the chemical-shrinkage void. This is
    the sealed-curing convention of Lavergne et al. (2018): the specimen keeps
    its volume while the reactions consume some, and the resulting empty
    porosity is a phase of the microstructure, not a rounding error.

Passing `reference = state` itself is equivalent to the default.

# Examples

```julia
f = volume_fractions(state_28d; reference = state_0)
f["Portlandite"], f["void"]
sum(values(f))                       # 1.0
```

See also: [`volume_fractions(state, groups)`](@ref), [`chemical_shrinkage`](@ref),
[`missing_molar_volumes`](@ref), [`porosity`](@ref).
"""
function volume_fractions(
        state::ChemicalState; reference = nothing, void_key::AbstractString = "void"
    )
    T = temperature(state)
    P = pressure(state)
    V_ref = reference === nothing ? volume(state).total : volume(reference).total
    _primal(ustrip(us"m^3", V_ref)) > 0 ||
        throw(ArgumentError("reference volume is zero: no fractions can be defined"))

    out = OrderedDict{String, Float64}()
    V_sum = 0.0
    for (i, sp) in enumerate(state.system.species)
        _has_molar_volume(sp) || continue
        # Filter on the AMOUNT, not on the sign of the contribution: aqueous
        # solutes have negative standard partial molar volumes (electrostriction
        # contracts the solvent around an ion — V⁰(OH⁻) ≈ -4.7 cm³/mol). Dropping
        # them would leave `volume` and `volume_fractions` disagreeing, and the
        # fractions summing to slightly more than one.
        iszero(_primal(ustrip(us"mol", state.n[i]))) && continue
        Vᵢ = state.n[i] * _molar_volume(sp)(T = T, P = P; unit = true)
        fᵢ = _primal(ustrip(Vᵢ / V_ref))
        out[symbol(sp)] = fᵢ
        V_sum += fᵢ
    end

    if reference !== nothing
        void = 1.0 - V_sum
        # A negative deficit means the reactions EXPANDED the solid+liquid beyond
        # the reference volume — physically possible (e.g. delayed ettringite) but
        # not representable as a void. Report it rather than clamp it silently.
        void < -1.0e-10 && @warn "volume expanded beyond the reference state; " *
            "the sealed-volume convention does not apply" excess = -void
        out[void_key] = max(void, 0.0)
    end
    return out
end

"""
    volume_fractions(state::ChemicalState, groups; kwargs...)
        -> OrderedDict{String, Float64}

Volume fractions aggregated into named families.

`groups` is any iterable of `name => symbols` pairs, where `symbols` is a
species symbol or a collection of them — the granularity a mean-field
homogenization scheme consumes, where "C-S-H" is one phase rather than four
solid-solution end members.

Keyword arguments are those of [`volume_fractions(state)`](@ref); when
`reference` is given, the chemical-shrinkage void is appended under `void_key`
as a group of its own.

Species that appear in no group are collected under `other_key` when their total
fraction exceeds `atol` **in magnitude**, and dropped otherwise — silently
discarding matter would defeat the purpose of a volume balance, and the leftover
may be negative when only aqueous solutes are unassigned. A species listed in
two groups raises an error rather than being counted twice.

# Examples

```julia
groups = [
    "anhydrous" => ["C3S", "C2S", "C3A", "C4AF"],
    "C-S-H"     => ["CSHQ-TobH", "CSHQ-TobD", "CSHQ-JenH", "CSHQ-JenD"],
    "CH"        => "Portlandite",
    "AFt"       => "ettringite",
    "water"     => "H2O@",
]
f = volume_fractions(state, groups; reference = state_0)
```

See also: [`volume_fractions(state)`](@ref).
"""
function volume_fractions(
        state::ChemicalState, groups;
        reference = nothing, void_key::AbstractString = "void",
        other_key::AbstractString = "other", atol::Real = 1.0e-12
    )
    per_species = volume_fractions(state; reference = reference, void_key = void_key)

    # Map species symbol → group name, rejecting any double assignment.
    owner = Dict{String, String}()
    for (gname, syms) in groups
        for s in (syms isa AbstractString ? (syms,) : syms)
            haskey(owner, s) && throw(
                ArgumentError(
                    "species \"$s\" is listed in both groups " *
                        "\"$(owner[s])\" and \"$(gname)\""
                )
            )
            owner[s] = String(gname)
        end
    end

    out = OrderedDict{String, Float64}()
    for (gname, _) in groups
        out[String(gname)] = 0.0
    end
    leftover = 0.0
    for (sym, f) in per_species
        if reference !== nothing && sym == void_key
            continue                       # appended last, keeps its own key
        elseif haskey(owner, sym)
            out[owner[sym]] += f
        else
            leftover += f
        end
    end

    # Test the MAGNITUDE: a leftover can be negative (aqueous solutes), and
    # dropping a negative one would push the total above 1.
    abs(leftover) > atol && (out[other_key] = leftover)
    reference !== nothing && (out[void_key] = per_species[void_key])
    return out
end
