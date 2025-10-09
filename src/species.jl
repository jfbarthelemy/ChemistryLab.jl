abstract type AbstractSpecies end

Base.isequal(s1::AbstractSpecies, s2::AbstractSpecies) = isequal(formula(s1), formula(s2)) && isequal(symbol(s1), symbol(s2))
==(s1::AbstractSpecies, s2::AbstractSpecies) = isequal(s1, s2)
Base.hash(s::AbstractSpecies, h::UInt) = hash(symbol(s), hash(formula(s), h))

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
    properties::OrderedDict{Symbol,Number}
end

expr(s::Species) = expr(formula(s))
phreeqc(s::Species) = phreeqc(formula(s))
unicode(s::Species) = unicode(formula(s))
colored(s::Species) = colored(formula(s))
components(s::Species) = atoms_charge(s)
mainformula(s::Species) = s.formula

function Species(formula::Formula; name=expr(formula), symbol=expr(formula))
    atoms = composition(formula)
    properties = OrderedDict(:molar_mass => calculate_molar_mass(atoms))
    return Species{valtype(atoms)}(name, symbol, formula, properties)
end

function Species(;expr::AbstractString="", name=expr, symbol=expr)
    formula = Formula(expr)
    return Species(formula; name=name, symbol=symbol)
end
Species(f::AbstractString; name=f, symbol=f) = Species(;expr=f, name=name, symbol=symbol)

function Species(atoms::AbstractDict{Symbol,T}, charge=0; name="", symbol="") where {T}
    formula = Formula(atoms, charge)
    properties = OrderedDict(:molar_mass => calculate_molar_mass(atoms))
    if length(name) == 0 name = unicode(formula) end
    if length(symbol) == 0 symbol = name end
    return Species{valtype(atoms)}(name, symbol, formula, properties)
end

function Species(atoms::Pair{Symbol,T}...; name="", symbol="") where {T}
    return Species(OrderedDict(atoms...); name=name, symbol=symbol)
end

function Base.convert(::Type{Species{T}}, s::Species; name=name(s), symbol=symbol(s)) where {T}
    return Species(convert(T, formula(s)); name=name, symbol=symbol)
end

function Species{T}(s::Species; name=name(s), symbol=symbol(s)) where {T}
    return convert(Species{T}, s; name=name, symbol=symbol)
end

function Species(s::Species; name=name(s), symbol=symbol(s))
    return Species(formula(s); name=name, symbol=symbol)
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
    # println(io, lpad("formula", pad), ": ", colored_formula(expr(s)), " | ", colored_formula(phreeqc(s)), " | ", colored_formula(unicode(s)))
    print_formula(io, formula(s), "formula", pad)
    println(io, lpad("atoms", pad), ": ", join(["$k:$v" for (k, v) in atoms(s)], ", "))
    println(io, lpad("charge", pad), ": ", charge(s))
    if length(properties(s))>0 print(io, lpad("properties", pad), ": ", join(["$k = $v" for (k, v) in properties(s)], "\n"*repeat(" ", pad+2))) end
end

struct CemSpecies{T<:Number, S<:Number} <: AbstractSpecies
    name::String
    symbol::String
    formula::Formula{T}
    cemformula::Formula{S}
    properties::OrderedDict{Symbol,Number}
end

cemformula(s::CemSpecies) = s.cemformula
mainformula(s::CemSpecies) = s.cemformula
expr(s::CemSpecies) = expr(cemformula(s))
phreeqc(s::CemSpecies) = phreeqc(cemformula(s))
unicode(s::CemSpecies) = unicode(cemformula(s))
colored(s::CemSpecies) = colored(cemformula(s))
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
    properties = OrderedDict(:molar_mass => calculate_molar_mass(atoms))
    return CemSpecies{valtype(atoms),valtype(composition(cemformula))}(name, symbol, formula, cemformula, properties)
end

function CemSpecies(; expr::AbstractString="", name=expr, symbol=expr)
    cemformula = Formula(expr)
    return CemSpecies(cemformula; name=name, symbol=symbol)
end

CemSpecies(f::AbstractString; name=f, symbol=f) = CemSpecies(; expr=f, name=name, symbol=symbol)

function CemSpecies(oxides::AbstractDict{Symbol,T}, charge=0; name="", symbol="") where {T}
    cemformula = Formula(oxides, charge; order=OXIDE_ORDER)
    if length(name) == 0 name = unicode(cemformula) end
    if length(symbol) == 0 symbol = name end
    return CemSpecies(cemformula; name=name, symbol=symbol)
end

function CemSpecies(oxides::Pair{Symbol,T}...; name="", symbol="") where {T}
    return CemSpecies(OrderedDict(oxides...); name=name, symbol=symbol)
end

function CemSpecies(s::Species)
    candidate_primaries = [Species(d; symbol=string(k)) for (k,d) in cement_to_mendeleev]
    A, indep_comp, dep_comp = stoich_matrix([s], candidate_primaries; display=false)
    oxides = OrderedDict(Symbol(symbol(indep_comp[i])) => A[i,1] for i in 1:size(A, 1))
    return CemSpecies(oxides, charge(s); name=name(s), symbol=symbol(s))
end

Species(s::CemSpecies) = Species{valtype(atoms(s))}(name(s), symbol(s), formula(s), properties(s))

Base.convert(::Type{<:Species}, s::CemSpecies) = Species(s)

function Base.convert(::Type{CemSpecies{S}}, s::CemSpecies; name=name(s), symbol=symbol(s)) where {S}
    return CemSpecies(convert(S, cemformula(s)); name=name, symbol=symbol)
end

function CemSpecies{S}(s::CemSpecies; name=name(s), symbol=symbol(s)) where {S}
    return convert(CemSpecies{S}, s; name=name, symbol=symbol)
end

function CemSpecies{S,T}(s::CemSpecies; name=name(s), symbol=symbol(s)) where {S,T}
    return convert(CemSpecies{S}, s; name=name, symbol=symbol)
end

function CemSpecies(s::CemSpecies; name=name(s), symbol=symbol(s))
    return CemSpecies(cemformula(s); name=name, symbol=symbol)
end

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
    # println(io, lpad("cemformula", pad), ": ", colored_formula(expr(cf)), " | ", colored_formula(phreeqc(cf)), " | ", colored_formula(unicode(cf)))
    print_formula(io, cf, "cemformula", pad)
    println(io, lpad("oxides", pad), ": ", join(["$k:$v" for (k, v) in oxides(s)], ", "))
    # println(io, lpad("formula", pad), ": ", colored_formula(expr(f)), " | ", colored_formula(phreeqc(f)), " | ", colored_formula(unicode(f)))
    print_formula(io, f, "formula", pad)
    println(io, lpad("atoms", pad), ": ", join(["$k:$v" for (k, v) in atoms(s)], ", "))
    println(io, lpad("charge", pad), ": ", charge(s))
    if length(properties(s))>0 print(io, lpad("properties", pad), ": ", join(["$k = $v" for (k, v) in properties(s)], "\n"*repeat(" ", pad+2))) end
end

Base.promote_rule(::Type{Species}, ::Type{<:AbstractSpecies}) = Species
Base.promote_rule(::Type{<:AbstractSpecies}, ::Type{Species}) = Species

Base.promote_rule(::Type{Species}, ::Type{Species{T}}) where {T} = Species
Base.promote_rule(::Type{Species{T}}, ::Type{Species}) where {T} = Species

Base.promote_rule(::Type{<:CemSpecies}, ::Type{Species{T}}) where {T} = Species
Base.promote_rule(::Type{Species{T}}, ::Type{<:CemSpecies}) where {T} = Species

function apply(func::Function, s::S, args... ; kwargs...) where {S<:AbstractSpecies}
    tryfunc(v) = try func(ustrip(v), args... ; kwargs...) * unit(v) catch; v end
    newcomponents = OrderedDict(k => tryfunc(v) for (k,v) ∈ components(s))
    newSpecies = root_type(typeof(s))(newcomponents, tryfunc(charge(s)); name=get(kwargs, :name, name(s)), symbol=get(kwargs, :symbol, symbol(s)))
    for (k,v) in properties(s)
        newSpecies[k] = tryfunc(v)
    end
    return newSpecies
end
