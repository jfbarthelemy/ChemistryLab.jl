struct AtomGroup{T<:Number} coef::T ; sym::Symbol end

*(n::Number, sym::Symbol) = AtomGroup(stoich_coef_round(n), sym)
*(sym::Symbol, n::Number) = AtomGroup(stoich_coef_round(n), sym)
/(sym::Symbol, n::Number) = AtomGroup(stoich_coef_round(inv(n)), sym)
//(sym::Symbol, n::Number) = AtomGroup(stoich_coef_round(1//n), sym)

Base.convert(::Type{AtomGroup}, sym::Symbol) = AtomGroup(1, sym)

AtomGroup(sym::Symbol) = AtomGroup(1, sym)

struct Formula{T<:Number}
    expr::String
    phreeqc::String
    unicode::String
    composition::Dict{Symbol,T}
    charge::Int8
end

stoichtype(::Formula{T}) where {T} = T
expr(f::Formula) = f.expr
phreeqc(f::Formula) = f.phreeqc
unicode(f::Formula) = f.unicode
composition(f::Formula) = f.composition
charge(f::Formula) = f.charge

function Formula(expr::AbstractString="")
    if expr == "Zz" || expr == "Zz+" || expr == "Zz⁺"
        composition = Dict{Symbol,Int}()
        charge = 1
    elseif expr == "e" || expr == "e-" || expr == "e⁻"
        composition = Dict{Symbol,Int}()
        charge = -1
    else
        composition = parse_formula(expr)
        charge = extract_charge(expr)
    end
    return Formula{valtype(composition)}(expr, unicode_to_phreeqc(expr), phreeqc_to_unicode(expr), composition, charge)
end

function Formula(composition::AbstractDict{Symbol,T}, charge=0; order=ATOMIC_ORDER) where {T<:Number}
    # 1. Filter out charge symbols and convert to Vector for sorting
    filtered_keys = setdiff(keys(composition), [:Zz, :Zz⁺, :e, :e⁻])
    sorted_keys = sort(collect(filtered_keys), by=k -> findfirst(==(k), order))

    # 2. Build the formula string
    expr_parts = String[]
    for k in sorted_keys
        v = composition[k]
        push!(expr_parts, string(k) * (isone(v) ? "" : string(v)))
    end
    expr = join(expr_parts, "")

    # 3. Handle charge (e.g., Ca²⁺, SO₄²⁻)
    charge += get(composition, :Zz, 0) + get(composition, :Zz⁺, 0)
    charge -= get(composition, :e, 0) + get(composition, :e⁻, 0)
    if !iszero(charge)
        sign = charge < 0 ? "-" : "+"
        abscharge = abs(charge)
        strch = isone(abscharge) ? "" : string(abscharge)
        expr *= sign * strch
    end

    # 4. Clean composition (remove charge keys)
    newcomposition = copy(composition)
    for s in [:Zz, :Zz⁺, :e, :e⁻]
        pop!(newcomposition, s, nothing)
    end

    return Formula{T}(expr, unicode_to_phreeqc(expr), phreeqc_to_unicode(expr), newcomposition, charge)
end


function Formula(f::Formula)
    return Formula(composition(f))
end

Base.getindex(f::Formula{T}, i::Symbol) where {T} = get(composition(f), i, zero(T))
Base.length(f::Formula) = length(composition(f))

==(f1::Formula, f2::Formula) = composition(f1) == composition(f2) && charge(f1) == charge(f2)
Base.isequal(f1::Formula, f2::Formula) = f1 == f2
Base.hash(f::Formula, h::UInt) = hash(composition(f), hash(charge(f), h))

function Base.show(io::IO, s::Formula)
    pad = 11
    println(io, typeof(s))
    println(io, lpad("formula", pad), ": ", colored_formula(expr(s)), " | ", colored_formula(phreeqc(s)), " | ", colored_formula((unicode(s))))
    println(io, lpad("composition", pad), ": ", join(["$k:$v" for (k, v) in composition(s)], ", "))
    println(io, lpad("charge", pad), ": ", charge(s))
end

function *(f::Formula, x::T) where {T<:Number}
    composition = Dict(k => x*v for (k,v) in f.composition)
    return Formula(composition, f.charge)
end

function /(f::Formula, x::T) where {T<:Number}
    composition = Dict(k => v/x for (k,v) in f.composition)
    return Formula(composition, f.charge)
end

function //(f::Formula, x::T) where {T<:Number}
    composition = Dict(k => v//x for (k,v) in f.composition)
    return Formula(composition, f.charge)
end

function +(f::Formula, atom::Symbol)
    composition = copy(f.composition)
    composition[atom] = get(f.composition, atom, 0) + 1
    return Formula(composition, f.charge)
end

function +(f::Formula, atom::AtomGroup)
    composition = copy(f.composition)
    composition[atom.sym] = get(f.composition, atom.sym, 0) + atom.coef
    return Formula(composition, f.charge)
end

function +(a::AtomGroup{T}, b::AtomGroup{S}) where {T,S}
    composition = a.sym == b.sym ? Dict{Symbol, promote_type(T, S)}(a.sym => a.coef+b.coef) : Dict{Symbol, promote_type(T, S)}(a.sym => a.coef, b.sym => b.coef)
    return Formula(composition, 0)
end

+(a::AtomGroup, b::Symbol) = a + AtomGroup(b)

+(a::Symbol, b::AtomGroup) = AtomGroup(a) + b

+(a::Symbol, b::Symbol) = AtomGroup(a) + AtomGroup(b)

function Base.convert(T::Type{<:Number}, f::Formula)
    newcomposition = Dict(k => convert(T, v) for (k,v) ∈ composition(f))
    return Formula{T}(expr(f), phreeqc(f), unicode(f), newcomposition, charge(f))
end

function Base.map(func::Function, f::Formula, args... ; kwargs...)
        newcomposition = Dict(k => func(v, args... ; kwargs...) for (k,v) ∈ composition(f))
    return Formula{valtype(newcomposition)}(expr(f), phreeqc(f), unicode(f), newcomposition, charge(f))
end
