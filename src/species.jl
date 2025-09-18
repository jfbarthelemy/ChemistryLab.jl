abstract type AbstractSpecies end

==(s1::AbstractSpecies, s2::AbstractSpecies) = formula(s1) == formula(s2)

name(s::AbstractSpecies) = s.name
symbol(s::AbstractSpecies) = s.symbol
formula(s::AbstractSpecies) = s.formula
atoms(s::AbstractSpecies) = composition(formula(s))
charge(s::AbstractSpecies) = charge(formula(s))
properties(s::AbstractSpecies) = s.properties

function atoms_charge(s::AbstractSpecies)
    z = charge(s)
    if iszero(z)
        return atoms(s)
    else
        ac = deepcopy(atoms(s))
        ac[:Zz] = z
        return ac
    end
end

Base.getindex(s::AbstractSpecies, i::Symbol) = get(atoms(s), i, get(properties(s), i, 0))

Base.setindex!(s::AbstractSpecies, value, i::Symbol) = setindex!(properties(s), value, i)

function Base.getproperty(s::AbstractSpecies, sym::Symbol)
    if sym in fieldnames(typeof(s))
        return getfield(s, sym)
    elseif sym in keys(properties(s))
        return properties(s)[sym]
    else
        error("Symbol '$sym' is neither a field nor a registered property.")
    end
end

function Base.setproperty!(s::AbstractSpecies, sym::Symbol, value)
    if sym in fieldnames(typeof(s))
        error("Cannot modify field '$sym' directly. Use constructor or dedicated methods.")
    else
        properties(s)[sym] = value
    end
end

struct Species{T<:Number} <: AbstractSpecies
    name::String
    symbol::String
    formula::Formula{T}
    properties::Dict{Symbol,Number}
end

expr(s::Species) = expr(formula(s))
phreeqc(s::Species) = phreeqc(formula(s))
unicode(s::Species) = unicode(formula(s))
components(s::Species) = atoms_charge(s)

function Species(formula::Formula; name=expr(formula), symbol=expr(formula))
    atoms = composition(formula)
    properties = Dict(:molar_mass => calculate_molar_mass(atoms))
    return Species{valtype(atoms)}(name, symbol, formula, properties)
end

function Species(;expr::AbstractString="", name=expr, symbol=expr)
    formula = Formula(expr)
    return Species(formula; name=name, symbol=symbol)
end
Species(f::AbstractString; name=f, symbol=f) = Species(;expr=f, name=name, symbol=symbol)

function Species(atoms::AbstractDict{Symbol,T}, charge=0; name="", symbol="") where {T}
    formula = Formula(atoms, charge)
    properties = Dict(:molar_mass => calculate_molar_mass(atoms))
    return Species{valtype(atoms)}(name, symbol, formula, properties)
end

function Base.show(io::IO, s::Species)
    pad = 11
    println(io, typeof(s))
    if name(s) != formula(s) && length(name(s))>0
        println(io, lpad("name", pad), ": ", name(s))
    end
    if symbol(s) != formula(s) && length(symbol(s))>0
        println(io, lpad("symbol", pad), ": ", symbol(s))
    end
    println(io, lpad("formula", pad), ": ", expr(s), " | ", phreeqc(s), " | ", unicode(s))
    println(io, lpad("atoms", pad), ": ", join(["$k:$v" for (k, v) in atoms(s)], ", "))
    println(io, lpad("charge", pad), ": ", charge(s))
    if length(properties(s))>0 println(io, lpad("properties", pad), ": ", join(["$k = $v" for (k, v) in properties(s)], "\n"*repeat(" ", pad+2))) end
end

struct CemSpecies{T<:Number, S<:Number} <: AbstractSpecies
    name::String
    symbol::String
    formula::Formula{T}
    cemformula::Formula{S}
    properties::Dict{Symbol,Number}
end

cemformula(s::CemSpecies) = s.cemformula
expr(s::CemSpecies) = expr(cemformula(s))
phreeqc(s::CemSpecies) = phreeqc(cemformula(s))
unicode(s::CemSpecies) = unicode(cemformula(s))
oxides(s::CemSpecies) = composition(cemformula(s))
function oxides_charge(s::CemSpecies)
    z = charge(s)
    if iszero(z)
        return oxides(s)
    else
        ac = deepcopy(oxides(s))
        ac[:Zz] = z
        return ac
    end
end
components(s::CemSpecies) = oxides_charge(s)

function CemSpecies(cemformula::Formula; name=expr(cemformula), symbol=expr(cemformula))
    formula = Formula(to_mendeleev(composition(cemformula)), charge(cemformula))
    atoms = composition(formula)
    properties = Dict(:molar_mass => calculate_molar_mass(atoms))
    return CemSpecies{valtype(atoms),valtype(composition(cemformula))}(name, symbol, formula, cemformula, properties)
end

function CemSpecies(;expr::AbstractString="", name=expr, symbol=expr)
    cemformula = Formula(expr)
    return CemSpecies(cemformula; name=name, symbol=symbol)
end

CemSpecies(f::AbstractString; name=f, symbol=f) = CemSpecies(;expr=f, name=name, symbol=symbol)

function CemSpecies(oxides::AbstractDict{Symbol,T}, charge=0; name="", symbol="") where {T}
    cemformula = Formula(oxides, charge; order=OXIDE_ORDER)
    return CemSpecies(cemformula; name=name, symbol=symbol)
end

Species(s::CemSpecies) = Species{valtype(atoms(s))}(name(s), symbol(s), formula(s), properties(s))

Base.convert(::Type{<:Species}, s::CemSpecies) = Species(s)

function Base.show(io::IO, s::CemSpecies)
    pad = 11
    println(io, typeof(s))
    if name(s) != expr(s) && length(name(s))>0
        println(io, lpad("name", pad), ": ", name(s))
    end
    if symbol(s) != expr(s) && length(symbol(s))>0
        println(io, lpad("symbol", pad), ": ", symbol(s))
    end
    cf = cemformula(s)
    f = formula(s)
    println(io, lpad("cemformula", pad), ": ", expr(cf), " | ", phreeqc(cf), " | ", unicode(cf))
    println(io, lpad("oxides", pad), ": ", join(["$k:$v" for (k, v) in oxides(s)], ", "))
    println(io, lpad("formula", pad), ": ", expr(f), " | ", phreeqc(f), " | ", unicode(f))
    println(io, lpad("atoms", pad), ": ", join(["$k:$v" for (k, v) in atoms(s)], ", "))
    println(io, lpad("charge", pad), ": ", charge(s))
    if length(properties(s))>0 println(io, lpad("properties", pad), ": ", join(["$k = $v" for (k, v) in properties(s)], "\n"*repeat(" ", pad+2))) end
end

Base.promote_rule(::Type{Species}, ::Type{CemSpecies}) = Species
Base.promote_rule(::Type{CemSpecies}, ::Type{Species}) = Species
