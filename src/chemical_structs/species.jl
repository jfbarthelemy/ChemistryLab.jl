# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using OrderedCollections
using PrettyTables

"""
    @enum AggregateState

Enumeration for species aggregate states.

# Values

  - `AS_UNDEF`: undefined state.
  - `AS_AQUEOUS`: aqueous solution.
  - `AS_CRYSTAL`: crystalline solid.
  - `AS_GAS`: gas phase.
"""
@enum AggregateState AS_UNDEF AS_AQUEOUS AS_CRYSTAL AS_GAS

"""
    @enum Class

Enumeration for species chemical classes.

# Values

  - `SC_UNDEF`: undefined class.
  - `SC_AQSOLVENT`: aqueous solvent.
  - `SC_AQSOLUTE`: aqueous solute.
  - `SC_COMPONENT`: component.
  - `SC_GASFLUID`: gas or fluid.
  - `SC_SSENDMEMBER`: end-member of a solid solution phase.
"""
@enum Class SC_UNDEF SC_AQSOLVENT SC_AQSOLUTE SC_COMPONENT SC_GASFLUID SC_SSENDMEMBER

"""
    abstract type AbstractSpecies end

Abstract base type for all chemical species representations.

All concrete species types (`Species`, `CemSpecies`) inherit from this type.
"""
abstract type AbstractSpecies end

const PropertyType = Union{
    Number,
    AbstractVector{<:Number},
    AbstractVector{<:Pair{Symbol}},
    Function,
    AbstractFunc,
    AbstractString,
    Missing,
}

"""
    Base.isequal(s1::AbstractSpecies, s2::AbstractSpecies) -> Bool

Compare two species for equality based on formula, aggregate state, and class.

# Examples

```jldoctest
julia> s1 = Species("H2O"; aggregate_state=AS_AQUEOUS);

julia> s2 = Species("H₂O"; aggregate_state=AS_AQUEOUS);

julia> s1 == s2
true
```
"""
function Base.isequal(s1::AbstractSpecies, s2::AbstractSpecies)
    return isequal(formula(s1), formula(s2)) &&
        isequal(aggregate_state(s1), aggregate_state(s2)) &&
        isequal(class(s1), class(s2))
end
==(s1::AbstractSpecies, s2::AbstractSpecies) = isequal(s1, s2)

"""
    Base.hash(s::AbstractSpecies, h::UInt) -> UInt

Compute hash for a species based on symbol, formula, aggregate state, and class.
"""
function Base.hash(s::AbstractSpecies, h::UInt)
    return hash(symbol(s), hash(formula(s), hash(aggregate_state(s), hash(class(s), h))))
end

"""
    name(s::AbstractSpecies) -> String

Return the name of the species.

# Examples

```jldoctest
julia> s1 = Species("H2O"; aggregate_state=AS_AQUEOUS);

julia> s1.name == "H2O"
true
```
"""
name(s::AbstractSpecies) = s.name

"""
    symbol(s::AbstractSpecies) -> String

Return the symbol of the species.
"""
symbol(s::AbstractSpecies) = s.symbol

"""
    formula(s::AbstractSpecies) -> Formula

Return the Formula object associated with the species.

# Examples

```jldoctest
julia> s1 = Species("H2O"; aggregate_state=AS_AQUEOUS);

julia> formula(s1) == Formula("H2O")
true
```
"""
formula(s::AbstractSpecies) = s.formula

"""
    atoms(s::AbstractSpecies) -> OrderedDict{Symbol,Number}

Return the atomic composition (element => coefficient) of the species.
"""
atoms(s::AbstractSpecies) = composition(formula(s))

"""
    charge(s::AbstractSpecies) -> Int8

Return the formal charge of the species.

# Examples

```jldoctest
julia> s1 = Species("Ca(HSiO3)+");

julia> charge(s1) == 1
true
```
"""
charge(s::AbstractSpecies) = charge(formula(s))

"""
    aggregate_state(s::AbstractSpecies) -> AggregateState

Return the aggregate state of the species.

# Examples

```jldoctest
julia> s1 = Species("H2O"; aggregate_state=AS_AQUEOUS);

julia> aggregate_state(s1) == AS_AQUEOUS
true
```
"""
aggregate_state(s::AbstractSpecies) = s.aggregate_state

"""
    class(s::AbstractSpecies) -> Class

Return the chemical class of the species.
"""
class(s::AbstractSpecies) = s.class

"""
    properties(s::AbstractSpecies) -> OrderedDict{Symbol,PropertyType}

Return the properties dictionary of the species.
"""
properties(s::AbstractSpecies) = s.properties

"""
    check_mendeleev(s::AbstractSpecies) -> Bool

Validate that all element symbols in the species exist in the periodic table.
"""
check_mendeleev(s::AbstractSpecies) = check_mendeleev(formula(s))

"""
    mendeleev_filter(s::AbstractSpecies) -> Union{AbstractSpecies,Nothing}

Return the species if valid according to Mendeleev check, otherwise `nothing`.
"""
mendeleev_filter(s::AbstractSpecies) = check_mendeleev(s) ? s : nothing

"""
    atoms_charge(s::AbstractSpecies) -> OrderedDict{Symbol,Number}

Return atomic composition including the charge as a :Zz key if non-zero.

# Examples

```jldoctest
julia> s = Species("Ca+2");

julia> atoms_charge(s)
OrderedDict{Symbol, Int64} with 2 entries:
  :Ca => 1
  :Zz => 2
```
"""
function atoms_charge(s::AbstractSpecies)
    z = charge(s)
    if iszero(z)
        return atoms(s)
    else
        ac = copy(atoms(s))
        ac[:Zz] = z
        return ac
    end
end

"""
    Base.getindex(s::AbstractSpecies, i::Symbol) -> Any

Access species components, atoms, or properties by symbol key.

Returns 0 if the key is not found.

# Examples

```jldoctest
julia> s = Species("H2O");

julia> s[:H]
2

julia> s[:N]
0
```
"""
function Base.getindex(s::AbstractSpecies, i::Symbol)
    coef = get(components(s), i, get(atoms(s), i, get(properties(s), i, nothing)))
    if isnothing(coef)
        # The thermodynamic functions are built on demand, and `getproperty`
        # already knows that. `getindex` did not, so `s[:Cp⁰]` returned the
        # not-found value 0 on any species whose functions had not been forced
        # yet — silently, and as an `Int64` that blows up only when the caller
        # tries to evaluate it. Returning 0 is right for a missing ATOM, which is
        # what that fallback is for; it is never right for a property that the
        # species can produce.
        if i in [:Cp⁰, :ΔₐH⁰, :S⁰, :ΔₐG⁰, :V⁰]
            complete_thermo_functions!(s)
            haskey(properties(s), i) && return properties(s)[i]
        end
        # println("$(i) not found in $(root_type(typeof(s))) $(colored(s))")
        return 0
    end
    return coef
end

"""
    Base.setindex!(s::AbstractSpecies, value, i::Symbol)

Set a property value for the species.
"""
Base.setindex!(s::AbstractSpecies, value, i::Symbol) = setindex!(properties(s), value, i)

"""
    Base.getproperty(s::AbstractSpecies, sym::Symbol) -> Any

Access species fields or registered properties.

Throws an error if the symbol is neither a field nor a property.
"""
function Base.getproperty(s::AbstractSpecies, sym::Symbol)
    if sym in fieldnames(typeof(s))
        return getfield(s, sym)
    else
        if !haskey(properties(s), sym) && sym in [:Cp⁰, :ΔₐH⁰, :S⁰, :ΔₐG⁰, :V⁰]
            complete_thermo_functions!(s)
        end
        return properties(s)[sym]
    end
end

"""
    Base.haskey(s::AbstractSpecies, sym::Symbol) -> Bool

Check if a property key exists in the species properties dictionary.
"""
function Base.haskey(s::AbstractSpecies, sym::Symbol)
    if haskey(properties(s), sym)
        return true
    else
        if sym in [:Cp⁰, :ΔₐH⁰, :S⁰, :ΔₐG⁰, :V⁰]
            complete_thermo_functions!(s)
            return haskey(properties(s), sym)
        else
            return false
        end
    end
end

"""
    Base.setproperty!(s::AbstractSpecies, sym::Symbol, value)

Set a property value, preventing direct modification of structural fields.

Throws an error if attempting to modify a structural field directly.
"""
function Base.setproperty!(s::AbstractSpecies, sym::Symbol, value)
    if !ismissing(value)
        if sym in fieldnames(typeof(s))
            error(
                "Cannot modify field '$sym' directly. Use constructor or dedicated methods."
            )
        else
            properties(s)[sym] = value
        end
    end
    return s
end

"""
    ordered_dict_with_default(gen, key_type, val_type) -> OrderedDict

Create an OrderedDict from a generator, ensuring proper typing even when empty.
"""
function ordered_dict_with_default(gen, key_type, val_type)
    d = OrderedDict(gen)
    isempty(d) && (d = OrderedDict{key_type, val_type}())
    return d
end

"""
    struct Species{T<:Number} <: AbstractSpecies

Standard chemical species representation using atomic composition.

# Fields

  - `name::String`: human-readable name.
  - `symbol::String`: species symbol.
  - `formula::Formula{T}`: chemical formula with stoichiometric coefficients.
  - `aggregate_state::AggregateState`: physical state.
  - `class::Class`: chemical class.
  - `properties::OrderedDict{Symbol,PropertyType}`: thermodynamic and other properties.

# Examples

```jldoctest
julia> s = Species("H2O"; name="Water", aggregate_state=AS_AQUEOUS);

julia> atoms(s)
OrderedDict{Symbol, Int64} with 2 entries:
  :H => 2
  :O => 1
```
"""
struct Species{T <: Number} <: AbstractSpecies
    name::String
    symbol::String
    formula::Formula{T}
    aggregate_state::AggregateState
    class::Class
    properties::OrderedDict{Symbol, PropertyType}
end

"""
    expr(s::Species) -> String

Return the original expression string of the species formula.

# Examples

```jldoctest
julia> expr(Species("H2O"; name="Water", aggregate_state=AS_AQUEOUS))
"H2O"

julia> expr(Species("H2O"; name="Water", aggregate_state=AS_AQUEOUS)) == expr(Formula("H2O"))
true
```
"""
expr(s::Species) = expr(formula(s))

"""
    phreeqc(s::Species) -> String

Return the PHREEQC-compatible representation of the species formula.
"""
phreeqc(s::Species) = phreeqc(formula(s))

"""
    unicode(s::Species) -> String

Return the Unicode pretty representation of the species formula.
"""
unicode(s::Species) = unicode(formula(s))

"""
    colored(s::Species) -> String

Return the colored terminal representation of the species formula.
"""
colored(s::Species) = colored(formula(s))

"""
    components(s::Species) -> OrderedDict{Symbol,Number}

Return the components of a Species (atomic composition with charge).
"""
components(s::Species) = atoms_charge(s)

"""
    mainformula(s::Species) -> Formula

Return the main formula representation for the species.
"""
mainformula(s::Species) = s.formula

"""
    Species(formula::Formula; name, symbol, aggregate_state, class, properties) -> Species

Construct a Species from a Formula object.

# Arguments

  - `formula`: Formula object with atomic composition.
  - `name`: species name (default: formula expression).
  - `symbol`: species symbol (default: formula expression).
  - `aggregate_state`: physical state (default: AS_UNDEF).
  - `class`: chemical class (default: SC_UNDEF).
  - `properties`: property dictionary (default: empty OrderedDict).

# Examples

```julia
julia> f = Formula("NaCl");

julia> s = Species(f; name="Sodium chloride", aggregate_state=AS_CRYSTAL);

julia> name(s)
"Sodium chloride"
```
"""
function Species(
        formula::Formula;
        name = expr(formula),
        symbol = expr(formula),
        aggregate_state = AS_UNDEF,
        class = SC_UNDEF,
        properties::AbstractDict = OrderedDict{Symbol, PropertyType}(),
    )
    atoms = composition(formula)
    if !haskey(properties, :M) && check_mendeleev(formula)
        properties[:M] = calculate_molar_mass(atoms)
    end
    return Species{valtype(atoms)}(
        name,
        symbol,
        formula,
        aggregate_state,
        class,
        OrderedDict{Symbol, PropertyType}(k => v for (k, v) in properties),
    )
end

"""
    Species(; expr, name, symbol, aggregate_state, class, properties) -> Species

Construct a Species from keyword arguments.

# Arguments

  - `expr`: formula string to parse (default: "").
  - `name`: species name (default: expr).
  - `symbol`: species symbol (default: expr).
  - `aggregate_state`: physical state (default: AS_UNDEF).
  - `class`: chemical class (default: SC_UNDEF).
  - `properties`: property dictionary (default: empty OrderedDict).
"""
function Species(;
        expr::AbstractString = "",
        name = expr,
        symbol = expr,
        aggregate_state = AS_UNDEF,
        class = SC_UNDEF,
        properties::AbstractDict = OrderedDict{Symbol, PropertyType}(),
    )
    return Species(
        Formula(expr);
        name = name,
        symbol = symbol,
        aggregate_state = aggregate_state,
        class = class,
        properties = properties,
    )
end

"""
    Species(f::AbstractString; name, symbol, aggregate_state, class, properties) -> Species

Construct a Species from a formula string.

# Arguments

  - `f`: formula string to parse.
  - `name`: species name (default: f).
  - `symbol`: species symbol (default: f).
  - `aggregate_state`: physical state (default: AS_UNDEF).
  - `class`: chemical class (default: SC_UNDEF).
  - `properties`: property dictionary (default: empty OrderedDict).

# Examples

```julia
julia> s = Species("Ca+2"; aggregate_state=AS_AQUEOUS);

julia> charge(s)
2
```
"""
function Species(
        f::AbstractString;
        name = f,
        symbol = f,
        aggregate_state = AS_UNDEF,
        class = SC_UNDEF,
        properties::AbstractDict = OrderedDict{Symbol, PropertyType}(),
    )
    return Species(;
        expr = f,
        name = name,
        symbol = symbol,
        aggregate_state = aggregate_state,
        class = class,
        properties = properties,
    )
end

"""
    Species(atoms::AbstractDict{Symbol,T}, charge=0; name, symbol, aggregate_state, class, properties) where {T} -> Species

Construct a Species from an atomic composition dictionary.

# Arguments

  - `atoms`: dictionary mapping element symbols to stoichiometric coefficients.
  - `charge`: formal charge (default 0).
  - `name`: species name (default: computed from formula).
  - `symbol`: species symbol (default: name).
  - `aggregate_state`: physical state (default: AS_UNDEF).
  - `class`: chemical class (default: SC_UNDEF).
  - `properties`: property dictionary (default: empty OrderedDict).
"""
function Species(
        atoms::AbstractDict{Symbol, T},
        charge = 0;
        name = "",
        symbol = "",
        aggregate_state = AS_UNDEF,
        class = SC_UNDEF,
        properties::AbstractDict = OrderedDict{Symbol, PropertyType}(),
    ) where {T}
    formula = Formula(atoms, charge)
    if length(name) == 0
        name = unicode(formula)
    end
    if length(symbol) == 0
        symbol = name
    end
    return Species(
        formula;
        name = name,
        symbol = symbol,
        aggregate_state = aggregate_state,
        class = class,
        properties = properties,
    )
end

"""
    Species(atoms::Pair{Symbol,T}...; name, symbol, aggregate_state, class, properties) where {T} -> Species

Construct a Species from element => coefficient pairs.

# Examples

```julia
julia> s = Species(:H => 2, :O => 1; name="Water");

julia> atoms(s)
OrderedDict{Symbol, Int64} with 2 entries:
  :H => 2
  :O => 1
```
"""
function Species(
        atoms::Pair{Symbol, T}...;
        name = "",
        symbol = "",
        aggregate_state = AS_UNDEF,
        class = SC_UNDEF,
        properties::AbstractDict = OrderedDict{Symbol, PropertyType}(),
    ) where {T}
    return Species(
        OrderedDict(atoms...);
        name = name,
        symbol = symbol,
        aggregate_state = aggregate_state,
        class = class,
        properties = properties,
    )
end

"""
    Species{T}(s::Species; kwargs...) where {T} -> Species{T}

Construct a Species with a specific coefficient type from another Species.
"""
function Species{T}(
        s::Species;
        name = name(s),
        symbol = symbol(s),
        aggregate_state = aggregate_state(s),
        class = class(s),
        properties = properties(s),
    ) where {T}
    return Species(
        convert(T, formula(s));
        name = name,
        symbol = symbol,
        aggregate_state = aggregate_state,
        class = class,
        properties = properties,
    )
end

function Species{T}(s::Species{T}; kwargs...) where {T}
    if isempty(kwargs)
        return s
    end
    return Species(
        convert(T, formula(s));
        name = name(s),
        symbol = symbol(s),
        aggregate_state = aggregate_state(s),
        class = class(s),
        properties = properties(s),
        kwargs...,
    )
end

"""
    Species(s::Species; kwargs...) -> Species

Copy constructor for Species with optional field overrides.
"""
function Species(s::Species; kwargs...)
    if isempty(kwargs)
        return s
    end
    return Species(
        formula(s);
        name = name(s),
        symbol = symbol(s),
        aggregate_state = aggregate_state(s),
        class = class(s),
        properties = properties(s),
        kwargs...,
    )
end

"""
    Base.show(io::IO, s::Species)

Compact single-line representation of a Species.
"""
function Base.show(io::IO, s::Species)
    return print(io, symbol(s), " {", name(s), "} [", formula(s), "]")
end

"""
    Base.show(io::IO, ::MIME"text/plain", s::Species)

Detailed multi-line REPL display for Species.
"""
function Base.show(io::IO, ::MIME"text/plain", s::Species)
    complete_thermo_functions!(s)
    pad = 15
    println(io, typeof(s))
    if name(s) != formula(s) && length(name(s)) > 0
        println(io, lpad("name", pad), ": ", name(s))
    end
    if symbol(s) != formula(s) && length(symbol(s)) > 0
        println(io, lpad("symbol", pad), ": ", symbol(s))
    end
    # println(io, lpad("formula", pad), ": ", colored_formula(expr(s)), " | ", colored_formula(phreeqc(s)), " | ", colored_formula(unicode(s)))
    print_formula(io, formula(s), "formula", pad)
    println(io, lpad("atoms", pad), ": ", join(["$k => $v" for (k, v) in atoms(s)], ", "))
    println(io, lpad("charge", pad), ": ", charge(s))
    println(io, lpad("aggregate_state", pad), ": ", aggregate_state(s))
    pr = length(properties(s)) > 0 ? println : print
    pr(io, lpad("class", pad), ": ", class(s))
    return if length(properties(s)) > 0
        print(
            io,
            lpad("properties", pad),
            ": ",
            join(["$k = $v" for (k, v) in properties(s)], "\n" * repeat(" ", pad + 2)),
        )
    end
end

"""
    pprint(s::Species)

Pretty-print a Species to standard output using the same multi-line layout
as the MIME "text/plain" show method.

# Arguments

  - `s` : Species instance to print.

# Returns

  - `nothing` (side-effect: formatted output to stdout).
"""
function pprint(s::Species)
    pad = 15
    println(typeof(s))
    if name(s) != formula(s) && length(name(s)) > 0
        println(lpad("name", pad), ": ", name(s))
    end
    if symbol(s) != formula(s) && length(symbol(s)) > 0
        println(lpad("symbol", pad), ": ", symbol(s))
    end
    pprint_formula(formula(s), "formula", pad)
    println(lpad("atoms", pad), ": ", join(["$k => $v" for (k, v) in atoms(s)], ", "))
    println(lpad("charge", pad), ": ", charge(s))
    println(lpad("aggregate_state", pad), ": ", aggregate_state(s))
    pr = length(properties(s)) > 0 ? println : print
    pr(lpad("class", pad), ": ", class(s))
    if length(properties(s)) > 0
        print(
            lpad("properties", pad),
            ": ",
            join(["$k = $v" for (k, v) in properties(s)], "\n" * repeat(" ", pad + 2)),
        )
    end
    return println()
end

"""
    struct CemSpecies{T<:Number,S<:Number} <: AbstractSpecies

Cement chemistry species representation using oxide notation.

# Fields

  - `name::String`: human-readable name.
  - `symbol::String`: species symbol.
  - `formula::Formula{T}`: atomic composition formula.
  - `cemformula::Formula{S}`: oxide notation formula.
  - `aggregate_state::AggregateState`: physical state.
  - `class::Class`: chemical class.
  - `properties::OrderedDict{Symbol,PropertyType}`: thermodynamic and other properties.

# Examples

```julia
julia> s = CemSpecies("C3A"; name="Tricalcium aluminate");

julia> oxides(s)
OrderedDict{Symbol, Int64} with 2 entries:
  :C => 3
  :A => 1
```
"""
struct CemSpecies{T <: Number, S <: Number} <: AbstractSpecies
    name::String
    symbol::String
    formula::Formula{T}
    cemformula::Formula{S}
    aggregate_state::AggregateState
    class::Class
    properties::OrderedDict{Symbol, PropertyType}
end

"""
    cemformula(s::CemSpecies) -> Formula

Return the oxide notation formula of the cement species.
"""
cemformula(s::CemSpecies) = s.cemformula

"""
    mainformula(s::CemSpecies) -> Formula

Return the main formula representation (oxide notation) for the cement species.
"""
mainformula(s::CemSpecies) = s.cemformula

"""
    expr(s::CemSpecies) -> String

Return the expression string of the cement formula.
"""
expr(s::CemSpecies) = expr(cemformula(s))

"""
    phreeqc(s::CemSpecies) -> String

Return the PHREEQC-compatible representation of the cement formula.
"""
phreeqc(s::CemSpecies) = phreeqc(cemformula(s))

"""
    unicode(s::CemSpecies) -> String

Return the Unicode representation of the cement formula.
"""
unicode(s::CemSpecies) = unicode(cemformula(s))

"""
    colored(s::CemSpecies) -> String

Return the colored terminal representation of the cement formula.
"""
colored(s::CemSpecies) = colored(cemformula(s))

"""
    oxides(s::CemSpecies) -> OrderedDict{Symbol,Number}

Return the oxide composition of the cement species.
"""
oxides(s::CemSpecies) = composition(cemformula(s))

oxides(s::Species) = oxides(CemSpecies(s))

"""
    oxides_charge(s::CemSpecies) -> OrderedDict{Symbol,Number}

Return oxide composition including the charge as a :Zz key if non-zero.
"""
function oxides_charge(s::CemSpecies)
    z = charge(s)
    if iszero(z)
        return oxides(s)
    else
        ac = copy(oxides(s))
        ac[:Zz] = z
        return ac
    end
end

"""
    components(s::CemSpecies) -> OrderedDict{Symbol,Number}

Return the components of a CemSpecies (oxide composition with charge).
"""
components(s::CemSpecies) = oxides_charge(s)

"""
    CemSpecies(cemformula::Formula; name, symbol, aggregate_state, class, properties) -> CemSpecies

Construct a CemSpecies from an oxide formula.

# Arguments

  - `cemformula`: Formula object in oxide notation.
  - `name`: species name (default: formula expression).
  - `symbol`: species symbol (default: formula expression).
  - `aggregate_state`: physical state (default: AS_UNDEF).
  - `class`: chemical class (default: SC_UNDEF).
  - `properties`: property dictionary (default: empty OrderedDict).
"""
function CemSpecies(
        cemformula::Formula;
        name = expr(cemformula),
        symbol = expr(cemformula),
        aggregate_state = AS_UNDEF,
        class = SC_UNDEF,
        properties::AbstractDict = OrderedDict{Symbol, PropertyType}(),
    )
    formula = Formula(to_mendeleev(composition(cemformula)), charge(cemformula))
    atoms = composition(formula)
    if !haskey(properties, :M) && check_mendeleev(formula)
        properties[:M] = calculate_molar_mass(atoms)
    end
    return CemSpecies{valtype(atoms), valtype(composition(cemformula))}(
        name, symbol, formula, cemformula, aggregate_state, class, properties
    )
end

"""
    CemSpecies(; expr, name, symbol, aggregate_state, class, properties) -> CemSpecies

Construct a CemSpecies from keyword arguments with an oxide formula string.
"""
function CemSpecies(;
        expr::AbstractString = "",
        name = expr,
        symbol = expr,
        aggregate_state = AS_UNDEF,
        class = SC_UNDEF,
        properties::AbstractDict = OrderedDict{Symbol, PropertyType}(),
    )
    return CemSpecies(
        Formula(expr);
        name = name,
        symbol = symbol,
        aggregate_state = aggregate_state,
        class = class,
        properties = properties,
    )
end

"""
    CemSpecies(f::AbstractString; name, symbol, aggregate_state, class, properties) -> CemSpecies

Construct a CemSpecies from an oxide formula string.

# Examples

```julia
julia> s = CemSpecies("C3S"; name="Alite");

julia> oxides(s)
OrderedDict{Symbol, Int64} with 2 entries:
  :C => 3
  :S => 1
```
"""
function CemSpecies(
        f::AbstractString;
        name = f,
        symbol = f,
        aggregate_state = AS_UNDEF,
        class = SC_UNDEF,
        properties::AbstractDict = OrderedDict{Symbol, PropertyType}(),
    )
    return CemSpecies(;
        expr = f,
        name = name,
        symbol = symbol,
        aggregate_state = aggregate_state,
        class = class,
        properties = properties,
    )
end

"""
    CemSpecies(oxides::AbstractDict{Symbol,T}, charge=0; name, symbol, aggregate_state, class, properties) where {T} -> CemSpecies

Construct a CemSpecies from an oxide composition dictionary.

# Arguments

  - `oxides`: dictionary mapping oxide symbols to stoichiometric coefficients.
  - `charge`: formal charge (default 0).
  - `name`: species name (default: computed from formula).
  - `symbol`: species symbol (default: name).
  - `aggregate_state`: physical state (default: AS_UNDEF).
  - `class`: chemical class (default: SC_UNDEF).
  - `properties`: property dictionary (default: empty OrderedDict).
"""
function CemSpecies(
        oxides::AbstractDict{Symbol, T},
        charge = 0;
        name = "",
        symbol = "",
        aggregate_state = AS_UNDEF,
        class = SC_UNDEF,
        properties::AbstractDict = OrderedDict{Symbol, PropertyType}(),
    ) where {T}
    cemformula = Formula(oxides, charge; order = OXIDE_ORDER)
    if length(name) == 0
        name = unicode(cemformula)
    end
    if length(symbol) == 0
        symbol = name
    end
    return CemSpecies(
        cemformula;
        name = name,
        symbol = symbol,
        aggregate_state = aggregate_state,
        class = class,
        properties = properties,
    )
end

"""
    CemSpecies(oxides::Pair{Symbol,T}...; name, symbol, aggregate_state, class, properties) where {T} -> CemSpecies

Construct a CemSpecies from oxide => coefficient pairs.

# Examples

```julia
julia> s = CemSpecies(:C => 3, :S => 2; name="C3S2");

julia> oxides(s)
OrderedDict{Symbol, Int64} with 2 entries:
  :C => 3
  :S => 2
```
"""
function CemSpecies(
        oxides::Pair{Symbol, T}...;
        name = "",
        symbol = "",
        aggregate_state = AS_UNDEF,
        class = SC_UNDEF,
        properties::AbstractDict = OrderedDict{Symbol, PropertyType}(),
    ) where {T}
    return CemSpecies(
        OrderedDict(oxides...);
        name = name,
        symbol = symbol,
        aggregate_state = aggregate_state,
        class = class,
        properties = properties,
    )
end

"""
    CemSpecies(s::Species; name, symbol, aggregate_state, class, properties) -> CemSpecies

Convert a Species to CemSpecies by decomposing into oxide notation.

Throws an error if the species cannot be decomposed into cement oxides.

# Arguments

  - `s`: source Species.
  - `name`: override name (default: keep original).
  - `symbol`: override symbol (default: keep original).
  - `aggregate_state`: override state (default: keep original).
  - `class`: override class (default: keep original).
  - `properties`: override properties (default: keep original).
"""
function CemSpecies(
        s::Species;
        name = name(s),
        symbol = symbol(s),
        aggregate_state = aggregate_state(s),
        class = class(s),
        properties = properties(s),
    )
    satoms = atoms(s)
    b = zeros(valtype(satoms), size(Aoxides, 1))
    for (atom, coef) in satoms
        if atom ∉ atoms_in_oxides
            error(
                "$(name) cannot be decomposed in cement oxides since $(atom) does not belong to cement atoms",
            )
        else
            b[order_atom_in_oxides[atom]] = coef
        end
    end
    x = stoich_coef_round.(Aoxides \ b)
    bcalc = stoich_coef_round.(Aoxides * x)
    if try
            isequal(bcalc, b) || isapprox(bcalc, b; rtol = 1.0e-4)
        catch
            false
        end
        oxides = OrderedDict(
            OXIDE_ORDER[i] => vx for (i, vx) in enumerate(x) if !iszero(vx)
        )
        return CemSpecies(
            oxides,
            charge(s);
            name = name,
            symbol = symbol,
            aggregate_state = aggregate_state,
            class = class,
            properties = properties,
        )
    else
        SM = StoichMatrix([s], oxides_as_species; pprint = false)
        A, indep_comp = SM.A, SM.primaries
        oxides = OrderedDict(Symbol(indep_comp[i].symbol) => A[i, 1] for i in 1:size(A, 1))
        if !isempty(oxides)
            cemspecies = CemSpecies(
                oxides,
                charge(s);
                name = name,
                symbol = symbol,
                aggregate_state = aggregate_state,
                class = class,
                properties = OrderedDict{Symbol, PropertyType}(
                    k => v for (k, v) in properties
                ),
            )
            if cemspecies == s
                return cemspecies
            end
        end
    end
    return error("$(name) cannot be decomposed in cement oxides")
end

"""
    Species(s::CemSpecies; name, symbol, aggregate_state, class, properties) -> Species

Convert a CemSpecies to Species using atomic composition.

# Arguments

  - `s`: source CemSpecies.
  - `name`: override name (default: keep original).
  - `symbol`: override symbol (default: keep original).
  - `aggregate_state`: override state (default: keep original).
  - `class`: override class (default: keep original).
  - `properties`: override properties (default: keep original).
"""
function Species(
        s::CemSpecies;
        name = name(s),
        symbol = symbol(s),
        aggregate_state = aggregate_state(s),
        class = class(s),
        properties = properties(s),
    )
    return Species{valtype(atoms(s))}(
        name,
        symbol,
        formula(s),
        aggregate_state,
        class,
        OrderedDict{Symbol, PropertyType}(k => v for (k, v) in properties),
    )
end

"""
    CemSpecies{S}(s::CemSpecies; kwargs...) where {S} -> CemSpecies{S}

Construct a CemSpecies with a specific coefficient type from another CemSpecies.
"""
function CemSpecies{S}(
        s::CemSpecies;
        name = name(s),
        symbol = symbol(s),
        aggregate_state = aggregate_state(s),
        class = class(s),
        properties = properties(s),
    ) where {S}
    return CemSpecies(
        convert(S, cemformula(s));
        name = name,
        symbol = symbol,
        aggregate_state = aggregate_state,
        class = class,
        properties = properties,
    )
end

function CemSpecies{S}(s::CemSpecies{S}; kwargs...) where {S}
    if isempty(kwargs)
        return s
    end
    return CemSpecies(
        convert(S, cemformula(s));
        name = name(s),
        symbol = symbol(s),
        aggregate_state = aggregate_state(s),
        class = class(s),
        properties = properties(s),
        kwargs...,
    )
end

"""
    CemSpecies{S,T}(s::CemSpecies; kwargs...) where {S,T} -> CemSpecies{S,T}

Construct a CemSpecies with specific coefficient types from another CemSpecies.
"""
function CemSpecies{S, T}(
        s::CemSpecies;
        name = name(s),
        symbol = symbol(s),
        aggregate_state = aggregate_state(s),
        class = class(s),
        properties = properties(s),
    ) where {S, T}
    return CemSpecies(
        convert(S, cemformula(s));
        name = name,
        symbol = symbol,
        aggregate_state = aggregate_state,
        class = class,
        properties = properties,
    )
end

function CemSpecies{S, T}(s::CemSpecies{S, T}; kwargs...) where {S, T}
    if isempty(kwargs)
        return s
    end
    return CemSpecies(
        convert(S, cemformula(s));
        name = name(s),
        symbol = symbol(s),
        aggregate_state = aggregate_state(s),
        class = class(s),
        properties = properties(s),
        kwargs...,
    )
end

"""
    CemSpecies(s::CemSpecies; kwargs...) -> CemSpecies

Copy constructor for CemSpecies with optional field overrides.
"""
function CemSpecies(s::CemSpecies; kwargs...)
    if isempty(kwargs)
        return s
    end
    return CemSpecies(
        cemformula(s);
        name = name(s),
        symbol = symbol(s),
        aggregate_state = aggregate_state(s),
        class = class(s),
        properties = properties(s),
        kwargs...,
    )
end

"""
    Base.show(io::IO, s::CemSpecies)

Compact single-line representation of a CemSpecies.
"""
function Base.show(io::IO, s::CemSpecies)
    return print(io, symbol(s), " {", name(s), "} [", cemformula(s), "]")
end

"""
    Base.show(io::IO, ::MIME"text/plain", s::CemSpecies)

Detailed multi-line REPL display for CemSpecies.
"""
function Base.show(io::IO, ::MIME"text/plain", s::CemSpecies)
    complete_thermo_functions!(s)
    pad = 15
    println(io, typeof(s))
    if name(s) != expr(s) && length(name(s)) > 0
        println(io, lpad("name", pad), ": ", name(s))
    end
    if symbol(s) != expr(s) && length(symbol(s)) > 0
        println(io, lpad("symbol", pad), ": ", symbol(s))
    end
    cf = cemformula(s)
    f = formula(s)
    # println(io, lpad("cemformula", pad), ": ", colored_formula(expr(cf)), " | ", colored_formula(phreeqc(cf)), " | ", colored_formula(unicode(cf)))
    print_formula(io, cf, "cemformula", pad)
    println(io, lpad("oxides", pad), ": ", join(["$k => $v" for (k, v) in oxides(s)], ", "))
    # println(io, lpad("formula", pad), ": ", colored_formula(expr(f)), " | ", colored_formula(phreeqc(f)), " | ", colored_formula(unicode(f)))
    print_formula(io, f, "formula", pad)
    println(io, lpad("atoms", pad), ": ", join(["$k => $v" for (k, v) in atoms(s)], ", "))
    println(io, lpad("charge", pad), ": ", charge(s))
    println(io, lpad("aggregate_state", pad), ": ", aggregate_state(s))
    pr = length(properties(s)) > 0 ? println : print
    pr(io, lpad("class", pad), ": ", class(s))
    return if length(properties(s)) > 0
        print(
            io,
            lpad("properties", pad),
            ": ",
            join(["$k = $v" for (k, v) in properties(s)], "\n" * repeat(" ", pad + 2)),
        )
    end
end

"""
    pprint(s::CemSpecies)

Pretty-print a CemSpecies to standard output using the same multi-line layout
as the MIME "text/plain" show method.

# Arguments

  - `s` : CemSpecies instance to print.

# Returns

  - `nothing` (side-effect: formatted output to stdout).
"""
function pprint(s::CemSpecies)
    pad = 15
    println(typeof(s))
    if name(s) != expr(s) && length(name(s)) > 0
        println(lpad("name", pad), ": ", name(s))
    end
    if symbol(s) != expr(s) && length(symbol(s)) > 0
        println(lpad("symbol", pad), ": ", symbol(s))
    end
    cf = cemformula(s)
    f = formula(s)
    # println(lpad("cemformula", pad), ": ", colored_formula(expr(cf)), " | ", colored_formula(phreeqc(cf)), " | ", colored_formula(unicode(cf)))
    pprint_formula(cf, "cemformula", pad)
    println(lpad("oxides", pad), ": ", join(["$k => $v" for (k, v) in oxides(s)], ", "))
    # println(lpad("formula", pad), ": ", colored_formula(expr(f)), " | ", colored_formula(phreeqc(f)), " | ", colored_formula(unicode(f)))
    pprint_formula(f, "formula", pad)
    println(lpad("atoms", pad), ": ", join(["$k => $v" for (k, v) in atoms(s)], ", "))
    println(lpad("charge", pad), ": ", charge(s))
    println(lpad("aggregate_state", pad), ": ", aggregate_state(s))
    pr = length(properties(s)) > 0 ? println : print
    pr(lpad("class", pad), ": ", class(s))
    if length(properties(s)) > 0
        print(
            lpad("properties", pad),
            ": ",
            join(["$k = $v" for (k, v) in properties(s)], "\n" * repeat(" ", pad + 2)),
        )
    end
    return println()
end

"""
    Base.promote_rule(::Type{Species}, ::Type{<:AbstractSpecies}) -> Type{Species}

Define promotion rule to convert AbstractSpecies to Species.
"""
Base.promote_rule(::Type{Species}, ::Type{<:AbstractSpecies}) = Species
Base.promote_rule(::Type{<:AbstractSpecies}, ::Type{Species}) = Species

Base.promote_rule(::Type{Species}, ::Type{Species{T}}) where {T} = Species
Base.promote_rule(::Type{Species{T}}, ::Type{Species}) where {T} = Species

Base.promote_rule(::Type{<:CemSpecies}, ::Type{Species{T}}) where {T} = Species
Base.promote_rule(::Type{Species{T}}, ::Type{<:CemSpecies}) where {T} = Species

"""
    apply(func::Function, s::S, args...; kwargs...) where {S<:AbstractSpecies} -> S

Apply a function element-wise to all numeric components and properties of a species.

# Arguments

  - `func`: function to apply.
  - `s`: source species.
  - `args...`: additional arguments for func.
  - `kwargs...`: keyword arguments (including potential overrides for name, symbol, etc.).

# Returns

  - New species with transformed values.

Handles Quantity types, attempting to preserve dimensions when possible.
"""
function apply(func::Function, s::S, args...; kwargs...) where {S <: AbstractSpecies}
    function tryfunc(v)
        return if v isa Quantity
            (
                try
                    func(ustrip(v), args...; kwargs...) *
                        func(dimension(v), args...; kwargs...)
                catch
                    try
                        func(ustrip(v), args...; kwargs...) * dimension(v)
                    catch
                        v
                    end
                end
            )
        else
            (
                try
                    func(v, args...; kwargs...)
                catch
                    v
                end
            )
        end
    end
    newcomponents = OrderedDict(k => tryfunc(v) for (k, v) in components(s))
    newSpecies = root_type(typeof(s))(
        newcomponents,
        tryfunc(charge(s));
        name = get(kwargs, :name, name(s)),
        symbol = get(kwargs, :symbol, symbol(s)),
        aggregate_state = get(kwargs, :aggregate_state, aggregate_state(s)),
        class = get(kwargs, :class, class(s)),
    )
    for (k, v) in properties(s)
        newSpecies[k] = tryfunc(v)
    end
    return newSpecies
end

"""
    find_species(s::AbstractString, species_list=nothing, S::Type{<:AbstractSpecies}=Species; aggregate_state=AS_UNDEF, class=SC_UNDEF) -> AbstractSpecies

Find or construct a species from a string identifier.

# Arguments

  - `s`: species identifier string (formula, symbol, or name).
  - `species_list`: optional list of species to search (default: nothing constructs new species).
  - `S`: species type to construct if not found (default: Species).
  - `aggregate_state`: filter by aggregate state (default: AS_UNDEF, no filter).
  - `class`: filter by chemical class (default: SC_UNDEF, no filter).

# Returns

  - Matching species from list, or newly constructed species if not found.

The function searches by symbol, PHREEQC format, Unicode format, formula expression,
and name. If multiple matches are found, a warning is displayed and the first match is returned.

# Examples

```julia
julia> species_list = [
           Species("H2O"; aggregate_state=AS_AQUEOUS), Species("H2O"; aggregate_state=AS_GAS)
       ];

julia> s = find_species("H2O", species_list; aggregate_state=AS_AQUEOUS);

julia> aggregate_state(s)
AS_AQUEOUS::AggregateState = 1
```
"""
function find_species(
        s::AbstractString,
        species_list = nothing,
        S::Type{<:AbstractSpecies} = Species;
        aggregate_state = AS_UNDEF,
        class = SC_UNDEF,
    )
    if isnothing(species_list)
        return S(s)
    else
        for crit in (symbol, phreeqc, unicode, expr ∘ mainformula, name)
            crit_vals = crit.(species_list)
            fil = species_list[
                .!isnothing.(species_list) .&& .!ismissing.(species_list) .&& ((s .== crit_vals) .|| (phreeqc_to_unicode(s) .== crit_vals) .|| (unicode_to_phreeqc(s) .== crit_vals)) .&& (aggregate_state .== AS_UNDEF) .|| (aggregate_state .== (x -> x.aggregate_state).(species_list)) .&& (class .== SC_UNDEF) .|| (
                    class .== (x -> x.class).(
                        species_list
                    )
                ),
            ]
            if length(fil) > 1
                println(crayon"red bold"("Several species correspond to $s:"))
                for x in fil
                    println("∙ ", x)
                end
                println(
                    crayon"red bold"(
                        "!!! In absence of more precision $(fil[1]) will be chosen !!!"
                    ),
                )
            end
            if length(fil) > 0
                return fil[1]
            end
        end
        comp_vals = composition.(mainformula.(species_list))
        fil = species_list[
            .!isnothing.(species_list) .&& .!ismissing.(species_list) .&& (comp_vals .== Ref(parse_formula(s))) .&& (aggregate_state .== AS_UNDEF) .|| (aggregate_state .== (x -> x.aggregate_state).(species_list)) .&& (class .== SC_UNDEF) .|| (
                class .== (x -> x.class).(
                    species_list
                )
            ),
        ]
        if length(fil) > 1
            println(crayon"red bold"("Several species correspond to $s:"))
            for x in fil
                println("∙ ", x)
            end
            println(
                crayon"red bold"(
                    "!!! In absence of more precision $(fil[1]) will be chosen !!!"
                ),
            )
        end
        if length(fil) > 0
            return fil[1]
        end
        return S(s)
    end
end

"""
    complete_thermo_functions!(s::AbstractSpecies)

Populate thermodynamic properties (`Cp⁰`, `ΔₐH⁰`, `S⁰`, `ΔₐG⁰`, `V⁰`) from parameters in `s.properties`.

If `thermo_params` dictionary is present in `properties`, it initializes thermodynamic functions
using `:thermo_method` (e.g., `"cp_ft_equation"`, `"solute_hkf88_reaktoro"`) or scalar defaults.
New thermodynamic models are registered by dispatching `build_thermo_functions(Val(:model_name), params)`.
"""
function complete_thermo_functions!(s::AbstractSpecies)
    if haskey(properties(s), :thermo_params)
        params = s[:thermo_params]
        dict_params = Dict(params)
        s.Tref = dict_params[:T]
        s.Pref = dict_params[:P]
        if haskey(properties(s), :thermo_method)
            dtf = build_thermo_functions(Symbol(s[:thermo_method]), params)
            for (k, v) in dtf
                s[k] = v
            end
            delete!(s.properties, :thermo_method)
        else
            if !haskey(properties(s), :Cp⁰) &&
                    haskey(dict_params, :Cp⁰) &&
                    !ismissing(dict_params[:Cp⁰])
                dtf = build_thermo_functions(
                    :cp_ft_equation, [:a₀ => dict_params[:Cp⁰]; params]
                )
                for (k, v) in dtf
                    s[k] = v
                end
            end
        end
        if haskey(properties(s), :V_method)
            s[:V⁰] = SymbolicFunc(dict_params[:V⁰])
            delete!(s.properties, :V_method)
        else
            for k in [:V⁰]
                if !haskey(properties(s), k) &&
                        haskey(dict_params, k) &&
                        !ismissing(dict_params[k])
                    s[k] = SymbolicFunc(dict_params[k])
                end
            end
        end
        for k in [:Cp⁰, :ΔₐH⁰, :S⁰, :ΔₐG⁰, :V⁰]
            if haskey(dict_params, k) && !ismissing(dict_params[k])
                s[Symbol(k, "_Tref")] = dict_params[k]
                if !haskey(properties(s), k)
                    s[k] = SymbolicFunc(dict_params[k])
                end
            end
        end
        delete!(s.properties, :thermo_params)
    end
    return s
end

"""
    with_class(s::Species, c::Class) -> Species

Return a copy of `s` with its class set to `c`. All other fields (name, symbol,
formula, aggregate state, properties) are preserved unchanged.

Useful to requalify database species as `SC_SSENDMEMBER` before grouping them
into a [`SolidSolutionPhase`](@ref), since `Species` is immutable.

# Examples

```jldoctest
julia> s = Species("CaCO3"; aggregate_state=AS_CRYSTAL, class=SC_COMPONENT);

julia> class(s)
SC_COMPONENT::Class = 3

julia> s2 = with_class(s, SC_SSENDMEMBER);

julia> class(s2)
SC_SSENDMEMBER::Class = 5
```
"""
function with_class(s::Species{T}, c::Class) where {T}
    return Species{T}(s.name, s.symbol, s.formula, s.aggregate_state, c, s.properties)
end
