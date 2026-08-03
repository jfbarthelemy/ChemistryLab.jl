# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using Crayons
using DynamicQuantities
using ForwardDiff

"""Terminal color for charge display (cyan bold)."""
const COL_CHARGE = crayon"cyan bold"
"""Terminal color for parentheses (magenta bold)."""
const COL_PAR = crayon"magenta bold"
"""Terminal color for internal stoichiometric coefficients (red bold)."""
const COL_STOICH_INT = crayon"red bold"
"""Terminal color for external stoichiometric coefficients (yellow bold)."""
const COL_STOICH_EXT = crayon"yellow bold"

"""
    _primal(x::Real) -> Real

Extract the primal (non-derivative) value of `x` for use in control-flow decisions.

For plain `Real` numbers this is a no-op.  For forward-mode AD types such as
`ForwardDiff.Dual` (which are subtypes of `Real`) this returns `x` itself, but
since those types overload comparison operators to compare only the primal value,
writing `_primal(x) > threshold` makes the intent explicit and ensures that
branching decisions are always taken on physical values rather than derivatives.

Adding a method `_primal(x::YourDualType) = x.value` enables support for
non-standard AD frameworks that do not overload comparison operators.

# Examples
```julia
julia> _primal(3.14)
3.14
```
"""
@inline _primal(x::Real) = x
@inline _primal(x::ForwardDiff.Dual) = ForwardDiff.value(x)

"""
    print_title(title; crayon=:none, indent="", style=:none)

Print a formatted title with optional styling and indentation.

# Arguments

  - `title`: title string to display.
  - `crayon`: Crayon for colored output (default `:none`).
  - `indent`: left indentation string (default `""`).
  - `style`: formatting style - `:none`, `:underline`, or `:box` (default `:none`).

# Examples

```julia
julia> print_title("Test Title"; style=:none)
Test Title

julia> print_title("Section"; style=:underline, indent="  ")
  Section
  ───────
```
"""
function print_title(title; crayon = :none, indent = "", style = :none)
    draw = crayon != :none ? x -> println(crayon(x)) : println
    if style == :underline
        width = length(title)
        draw(indent * "$title")
        draw(indent * "─"^width)
    elseif style == :box
        width = length(title) + 4
        draw(indent * "┌" * "─"^width * "┐")
        draw(indent * "│  $title  │")
        draw(indent * "└" * "─"^width * "┘")
    else
        draw(indent * "$title")
    end
    return print(crayon"reset")
end

"""
    root_type(T) -> Type

Get the root type wrapper from a potentially parameterized type.
"""
root_type(T) = T isa UnionAll ? root_type(T.body) : T.name.wrapper

"""
    safe_ustrip(unit::UnionAbstractQuantity, q::UnionAbstractQuantity) -> Number
    safe_ustrip(unit::UnionAbstractQuantity, q) -> Number

Safely remove units from `q`, converting to `unit` when their dimensions match.
If `q` is a plain numeric value (no units) it is returned unchanged.

When the dimensions of `unit` and `q` differ (e.g. `unit` is dimensionless but
`q` carries a physical dimension) the conversion is skipped and `q` is stripped
in its native unit.  This prevents a `DimensionError` in contexts such as
`ThermoFactory` where a dimensionless placeholder unit `u"1"` may be passed
alongside a quantity with physical dimensions.

# Examples

```julia
julia> safe_ustrip(u"m", 5u"m")
5.0

julia> safe_ustrip(u"cm", 1u"m")
100.0

julia> safe_ustrip(u"m", 3.0)
3.0
```
"""
function safe_ustrip(unit::UnionAbstractQuantity, q::UnionAbstractQuantity)
    # Avoid comparing dimension types directly (SymbolicDimensions vs Dimensions
    # comparison raises an error in DynamicQuantities ≤ 1.11).  Instead, attempt
    # the direct conversion, then try uconvert + ustrip as fallback.
    try
        return ustrip(unit, q)
    catch
        try
            return ustrip(uconvert(unit, q))
        catch
            return ustrip(q)
        end
    end
end

safe_ustrip(::UnionAbstractQuantity, q) = ustrip(q)

"""
    safe_uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, q::UnionAbstractQuantity{<:Any, <:Dimensions}) -> UnionAbstractQuantity
    safe_uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, q) -> Any

Convert `q` to the unit type represented by `qout` using `uconvert` when `q`
is a quantity with compatible dimensions. If `q` is a plain numeric value,
return `q` unchanged.

# Arguments

  - `qout` : unitful quantity target (e.g., u"m", u"kg").
  - `q` : a quantity or a plain numeric value.

# Returns

  - Converted quantity (of type `qout`) or the original `q` when `q` has no units.

# Examples

```julia
julia> safe_uconvert(us"m", 5u"cm")
0.05 m

julia> safe_uconvert(us"m", 3.0)
3.0
```
"""
safe_uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, q::UnionAbstractQuantity{<:Any, <:Dimensions}) = uconvert(qout, q)

safe_uconvert(::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, q) = q

"""
    _ensure_unit(unit, x::Real) -> Quantity
    _ensure_unit(unit, x) -> Quantity

Ensure a value is a `Quantity` with the given `unit`.

  - Plain `Real` → interpreted as SI, wrapped: `x * unit`.
  - `Quantity` → converted to `unit` via [`safe_uconvert`](@ref).
"""
_ensure_unit(unit, x::Real) = x * unit
_ensure_unit(unit, x) = safe_uconvert(unit, x)

"""
    safe_uparse(x::AbstractString) -> AbstractQuantity
    safe_uparse(x::AbstractQuantity) -> AbstractQuantity

Parse a unit string into a quantity via `uparse`, or return an existing
quantity unchanged. Acts as an idempotent wrapper.
"""
safe_uparse(x::AbstractString) = uparse(x)

safe_uparse(x::AbstractQuantity) = x

"""
    force_uconvert(qout::UnionAbstractQuantity, q) -> AbstractQuantity

Strip units from `q` (converting to `qout` dimensions when possible) and
re-attach the unit of `qout`. Unlike `safe_uconvert`, this always returns
a quantity with the units of `qout`, even when the input is dimensionless.
"""
force_uconvert(qout::UnionAbstractQuantity, q) = safe_ustrip(qout, q) * qout
