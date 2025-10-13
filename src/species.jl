@enum AggregateState AS_UNDEF AS_AQUEOUS AS_CRYSTAL AS_GAS

@enum Class SC_UNDEF SC_AQSOLVENT SC_AQSOLUTE SC_COMPONENT SC_GAS_FLUID

abstract type AbstractSpecies end

Base.isequal(s1::AbstractSpecies, s2::AbstractSpecies) = isequal(formula(s1), formula(s2)) && isequal(aggregate_state(s1), aggregate_state(s2)) && isequal(class(s1), class(s2))
==(s1::AbstractSpecies, s2::AbstractSpecies) = isequal(s1, s2)
Base.hash(s::AbstractSpecies, h::UInt) = hash(symbol(s), hash(formula(s), hash(aggregate_state(s), hash(class(s), h))))

name(s::AbstractSpecies) = s.name
symbol(s::AbstractSpecies) = s.symbol
formula(s::AbstractSpecies) = s.formula
atoms(s::AbstractSpecies) = composition(formula(s))
charge(s::AbstractSpecies) = charge(formula(s))
aggregate_state(s::AbstractSpecies) = s.aggregate_state
class(s::AbstractSpecies) = s.class
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

Base.getindex(s::AbstractSpecies, i::Symbol) = get(atoms(s), i, get(properties(s), i, nothing))

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
    return s
end

const PropertyType = Union{Number,AbstractVector{<:Number},Function,AbstractString}

struct Species{T<:Number} <: AbstractSpecies
    name::String
    symbol::String
    formula::Formula{T}
    aggregate_state::AggregateState
    class::Class
    properties::OrderedDict{Symbol,PropertyType}
end

expr(s::Species) = expr(formula(s))
phreeqc(s::Species) = phreeqc(formula(s))
unicode(s::Species) = unicode(formula(s))
colored(s::Species) = colored(formula(s))
components(s::Species) = atoms_charge(s)
mainformula(s::Species) = s.formula

function Species(formula::Formula; name=expr(formula), symbol=expr(formula), aggregate_state=AS_UNDEF, class=SC_UNDEF, properties=OrderedDict{Symbol,PropertyType}())
    atoms = composition(formula)
    if !haskey(properties, :molar_mass) properties[:molar_mass] = calculate_molar_mass(atoms) end
    return Species{valtype(atoms)}(name, symbol, formula, aggregate_state, class, properties)
end

Species(; expr::AbstractString="", name=expr, symbol=expr, aggregate_state=AS_UNDEF, class=SC_UNDEF, properties=OrderedDict{Symbol,PropertyType}()) = Species(Formula(expr); name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)

Species(f::AbstractString; name=f, symbol=f, aggregate_state=AS_UNDEF, class=SC_UNDEF, properties=OrderedDict{Symbol,PropertyType}()) = Species(; expr=f, name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)

function Species(atoms::AbstractDict{Symbol,T}, charge=0; name="", symbol="", aggregate_state=AS_UNDEF, class=SC_UNDEF, properties=OrderedDict{Symbol,PropertyType}()) where {T}
    formula = Formula(atoms, charge)
    if length(name) == 0
        name = unicode(formula)
    end
    if length(symbol) == 0
        symbol = name
    end
    return Species(formula; name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)
end

function Species(atoms::Pair{Symbol,T}...; name="", symbol="", aggregate_state=AS_UNDEF, class=SC_UNDEF, properties=OrderedDict{Symbol,PropertyType}()) where {T}
    return Species(OrderedDict(atoms...); name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)
end

function Base.convert(::Type{Species{T}}, s::Species; name=name(s), symbol=symbol(s), aggregate_state=aggregate_state(s), class=class(s), properties=properties(s)) where {T}
    return Species(convert(T, formula(s)); name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)
end

function Species{T}(s::Species; name=name(s), symbol=symbol(s), aggregate_state=aggregate_state(s), class=class(s), properties=properties(s)) where {T}
    return convert(Species{T}, s; name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)
end

function Species(s::Species; name=name(s), symbol=symbol(s), aggregate_state=aggregate_state(s), class=class(s), properties=properties(s))
    return Species(formula(s); name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)
end

function Base.show(io::IO, s::Species)
    print(io, name(s), " [", formula(s), "]")
end

function Base.show(io::IO, ::MIME"text/plain", s::Species)
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
    println(io, lpad("class", pad), ": ", class(s))
    if length(properties(s)) > 0
        print(io, lpad("properties", pad), ": ", join(["$k = $v" for (k, v) in properties(s)], "\n" * repeat(" ", pad + 2)))
    end
end

struct CemSpecies{T<:Number,S<:Number} <: AbstractSpecies
    name::String
    symbol::String
    formula::Formula{T}
    cemformula::Formula{S}
    aggregate_state::AggregateState
    class::Class
    properties::OrderedDict{Symbol,PropertyType}
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

function CemSpecies(cemformula::Formula; name=expr(cemformula), symbol=expr(cemformula), aggregate_state=AS_UNDEF, class=SC_UNDEF, properties=OrderedDict{Symbol,PropertyType}())
    formula = Formula(to_mendeleev(composition(cemformula)), charge(cemformula))
    atoms = composition(formula)
    if !haskey(properties, :molar_mass) properties[:molar_mass] = calculate_molar_mass(atoms) end
    return CemSpecies{valtype(atoms),valtype(composition(cemformula))}(name, symbol, formula, cemformula, aggregate_state, class, properties)
end

CemSpecies(; expr::AbstractString="", name=expr, symbol=expr, aggregate_state=AS_UNDEF, class=SC_UNDEF, properties=OrderedDict{Symbol,PropertyType}()) = CemSpecies(Formula(expr); name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)

CemSpecies(f::AbstractString; name=f, symbol=f, aggregate_state=AS_UNDEF, class=SC_UNDEF, properties=OrderedDict{Symbol,PropertyType}()) = CemSpecies(; expr=f, name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)

function CemSpecies(oxides::AbstractDict{Symbol,T}, charge=0; name="", symbol="", aggregate_state=AS_UNDEF, class=SC_UNDEF, properties=OrderedDict{Symbol,PropertyType}()) where {T}
    cemformula = Formula(oxides, charge; order=OXIDE_ORDER)
    if length(name) == 0
        name = unicode(cemformula)
    end
    if length(symbol) == 0
        symbol = name
    end
    return CemSpecies(cemformula; name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)
end

function CemSpecies(oxides::Pair{Symbol,T}...; name="", symbol="", aggregate_state=AS_UNDEF, class=SC_UNDEF, properties=OrderedDict{Symbol,PropertyType}()) where {T}
    return CemSpecies(OrderedDict(oxides...); name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)
end

function CemSpecies(s::Species; name=name(s), symbol=symbol(s), aggregate_state=aggregate_state(s), class=class(s), properties=properties(s))
    satoms = atoms(s)
    b = zeros(valtype(satoms), size(Aoxides, 1))
    for (atom, coef) in satoms
        if atom ∉ atoms_in_oxides
            error("$(name) cannot be decomposed in cement oxides since $(atom) does not belong to cement atoms")
        else
            b[order_atom_in_oxides[atom]] = coef
        end
    end
    x = stoich_coef_round.(Aoxides\b)
    bcalc = stoich_coef_round.(Aoxides*x)
    if try isequal(bcalc, b) || isapprox(bcalc, b; rtol=1.e-4) catch; false end
        oxides = OrderedDict(OXIDE_ORDER[i] => vx for (i,vx) in enumerate(x) if !iszero(vx))
        return CemSpecies(oxides, charge(s); name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)
    else
        A, indep_comp, dep_comp = stoich_matrix([s], oxides_as_species; display=false)
        oxides = OrderedDict(Symbol(indep_comp[i].symbol) => A[i, 1] for i in 1:size(A, 1))
        cemspecies = CemSpecies(oxides, charge(s); name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)
        if cemspecies == s
            return cemspecies
        else
            error("$(name) cannot be decomposed in cement oxides")
        end
    end
end

Species(s::CemSpecies; name=name(s), symbol=symbol(s), aggregate_state=aggregate_state(s), class=class(s), properties=properties(s)) = Species{valtype(atoms(s))}(name, symbol, formula(s), aggregate_state, class, properties)

Base.convert(::Type{<:Species}, s::CemSpecies; name=name(s), symbol=symbol(s), aggregate_state=aggregate_state(s), class=class(s), properties=properties(s)) = Species(s; name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)

function Base.convert(::Type{CemSpecies{S}}, s::CemSpecies; name=name(s), symbol=symbol(s), aggregate_state=aggregate_state(s), class=class(s), properties=properties(s)) where {S}
    return CemSpecies(convert(S, cemformula(s)); name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)
end

function CemSpecies{S}(s::CemSpecies; name=name(s), symbol=symbol(s), aggregate_state=aggregate_state(s), class=class(s), properties=properties(s)) where {S}
    return convert(CemSpecies{S}, s; name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)
end

function CemSpecies{S,T}(s::CemSpecies; name=name(s), symbol=symbol(s), aggregate_state=aggregate_state(s), class=class(s), properties=properties(s)) where {S,T}
    return convert(CemSpecies{S}, s; name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)
end

function CemSpecies(s::CemSpecies; name=name(s), symbol=symbol(s), aggregate_state=aggregate_state(s), class=class(s), properties=properties(s))
    return CemSpecies(cemformula(s); name=name, symbol=symbol, aggregate_state=aggregate_state, class=class, properties=properties)
end

function Base.show(io::IO, s::CemSpecies)
    print(io, name(s), " [", cemformula(s), "]")
end

function Base.show(io::IO, ::MIME"text/plain", s::CemSpecies)
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
    println(io, lpad("class", pad), ": ", class(s))
    if length(properties(s)) > 0
        print(io, lpad("properties", pad), ": ", join(["$k = $v" for (k, v) in properties(s)], "\n" * repeat(" ", pad + 2)))
    end
end

Base.promote_rule(::Type{Species}, ::Type{<:AbstractSpecies}) = Species
Base.promote_rule(::Type{<:AbstractSpecies}, ::Type{Species}) = Species

Base.promote_rule(::Type{Species}, ::Type{Species{T}}) where {T} = Species
Base.promote_rule(::Type{Species{T}}, ::Type{Species}) where {T} = Species

Base.promote_rule(::Type{<:CemSpecies}, ::Type{Species{T}}) where {T} = Species
Base.promote_rule(::Type{Species{T}}, ::Type{<:CemSpecies}) where {T} = Species

function apply(func::Function, s::S, args...; kwargs...) where {S<:AbstractSpecies}
    tryfunc(v) = v isa Quantity ? (
        try func(ustrip(v), args...; kwargs...) * func(unit(v), args...; kwargs...); catch; try func(ustrip(v), args...; kwargs...) * unit(v); catch; v; end; end
        ) : (
        try func(v, args...; kwargs...); catch; v; end
        )
    newcomponents = OrderedDict(k => tryfunc(v) for (k, v) ∈ components(s))
    newSpecies = root_type(typeof(s))(newcomponents, tryfunc(charge(s)); name=get(kwargs, :name, name(s)), symbol=get(kwargs, :symbol, symbol(s)), aggregate_state=get(kwargs, :aggregate_state, aggregate_state(s)), class=get(kwargs, :class, class(s)))
    for (k, v) in properties(s)
        newSpecies[k] = tryfunc(v)
    end
    return newSpecies
end
