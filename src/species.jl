struct AtomGroup{T<:Number} coef::T ; sym::Symbol end

*(n::Number, sym::Symbol) = AtomGroup(n, sym)

Base.convert(AtomGroup, sym::Symbol) = AtomGroup(1, sym)

AtomGroup(sym::Symbol) = AtomGroup(1, sym)

abstract type AbstractSpecies{T<:Number} end

struct Species{T} <: AbstractSpecies{T}
    name::String
    formula::String
    symbol::String
    phreeqc::String
    unicode::String
    atoms::Dict{Symbol,T}
    charge::Int8
    molar_mass::typeof(1.0g/mol)
    properties::Dict{Symbol,Number}
end

function Species(;formula::AbstractString="", name=formula, symbol=formula)
    atoms = parse_formula(formula)
    charge = extract_charge(formula)
    molar_mass = calculate_molar_mass(atoms)
    return Species{valtype(atoms)}(name, formula, symbol, unicode_to_phreeqc(formula), phreeqc_to_unicode(formula), atoms, charge, molar_mass, Dict{Symbol,Number}())
end

Species(f::AbstractString; name=f, symbol=f) = Species(;formula=f, name=name, symbol=symbol)

function Species(atoms::AbstractDict{Symbol,T}, charge=0; name="", symbol="") where {T}
    scharge = (:z, :z⁺, :e, :e⁻)
    formula = join(k in scharge ? "" : string(k)*string(isone(v) ? "" : v) for (k,v) in atoms)
    charge += get(atoms, :z, 0)
    charge += get(atoms, :z⁺, 0)
    charge -= get(atoms, :e, 0)
    charge -= get(atoms, :e⁻, 0) 
    if !iszero(charge)
        sign = charge<0 ? "-" : "+"
        abscharge = abs(charge)
        strch = isone(abscharge) ? "" : string(abscharge)
        formula *= sign * strch
    end
    if iszero(length(name))
        name = formula
    end
    if iszero(length(symbol))
        symbol = formula
    end
    newatoms = copy(atoms)
    for s in scharge pop!(newatoms, s, 0) end
    molar_mass = calculate_molar_mass(newatoms)
    return Species{T}(name, formula, symbol, unicode_to_phreeqc(formula), phreeqc_to_unicode(formula), newatoms, charge, molar_mass, Dict{Symbol,Number}())
end

name(s::Species) = s.name
formula(s::Species) = s.formula
symbol(s::Species) = s.symbol
phreeqc(s::Species) = s.phreeqc
unicode(s::Species) = s.unicode
atoms(s::Species) = s.atoms
charge(s::Species) = s.charge
molar_mass(s::Species) = s.molar_mass
properties(s::Species) = s.properties

==(s1::Species, s2::Species) = s1.atoms == s2.atoms && s1.charge == s2.charge

Base.getindex(s::AbstractSpecies{T}, i::Symbol) where {T} = get(atoms(s), i, get(properties(s), i, zero(T)))

Base.setindex!(s::AbstractSpecies, value, i::Symbol) = setindex!(properties(s), value, i)

function Base.getproperty(s::Species, sym::Symbol)
    if sym in fieldnames(Species)
        return getfield(s, sym)
    elseif sym in keys(s.properties)
        return s.properties[sym]
    else
        error("Symbol '$sym' is neither a field nor a registered property.")
    end
end

function Base.setproperty!(s::Species, sym::Symbol, value)
    if sym in fieldnames(typeof(s))
        error("Cannot modify field '$sym' directly. Use constructor or dedicated methods.")
    else
        s.properties[sym] = value
    end
end

function Base.show(io::IO, s::Species)
    pad = 15
    if name(s) != formula(s) && length(name(s))>0
        println(io, lpad("name", pad), ": ", name(s))
    end
    println(io, lpad("formula", pad), ": ", formula(s))
    if symbol(s) != formula(s) && length(symbol(s))>0
        println(io, lpad("symbol", pad), ": ", symbol(s))
    end
    println(io, lpad("Phreeqc formula", pad), ": ", phreeqc(s))
    println(io, lpad("Unicode formula", pad), ": ", unicode(s))
    println(io, lpad("atoms", pad), ": ", join(["$k:$v" for (k, v) in atoms(s)], ", "))
    println(io, lpad("charge", pad), ": ", charge(s))
    println(io, lpad("molar mass", pad), ": ", molar_mass(s))
    if length(properties(s))>0 println(io, lpad("properties", pad), ": ", join(["$k = $v" for (k, v) in properties(s)], "\n"*repeat(" ", pad+2))) end
end

function *(s::Species, x::T) where {T<:Number}
    atoms = Dict(k => x*v for (k,v) in s.atoms)
    return Species(atoms, s.charge; name=s.name, symbol=s.symbol)
end

function /(s::Species, x::T) where {T<:Number}
    atoms = Dict(k => v/x for (k,v) in s.atoms)
    return Species(atoms, s.charge; name=s.name, symbol=s.symbol)
end

function //(s::Species, x::T) where {T<:Number}
    atoms = Dict(k => v//x for (k,v) in s.atoms)
    return Species(atoms, s.charge; name=s.name, symbol=s.symbol)
end

function +(s::Species, atom::Symbol)
    atoms = copy(f.atoms)
    atoms[atom] = get(f.atoms, atom, 0) + 1
    return Species(atoms, f.charge; name="", symbol="")
end

function +(s::Species, atom::AtomGroup)
    atoms = copy(f.atoms)
    atoms[atom.sym] = get(f.atoms, atom.sym, 0) + atom.coef
    return Species(atoms, f.charge; name="", symbol="")
end

function +(a::AtomGroup{T}, b::AtomGroup{S}) where {T,S}
    atoms = Dict{Symbol, promote_type(T, S)}(a.sym => a.coef, b.sym => b.coef)
    return Species(atoms, 0)
end

+(a::AtomGroup, b::Symbol) = a + AtomGroup(b)

+(a::Symbol, b::AtomGroup) = AtomGroup(a) + b

+(a::Symbol, b::Symbol) = AtomGroup(a) + AtomGroup(b)
