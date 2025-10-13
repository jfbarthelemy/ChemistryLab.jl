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
    colored::String
    composition::OrderedDict{Symbol,T}
    charge::Int8
end

stoichtype(::Formula{T}) where {T} = T
expr(f::Formula) = f.expr
phreeqc(f::Formula) = f.phreeqc
unicode(f::Formula) = f.unicode
colored(f::Formula) = f.colored
composition(f::Formula) = f.composition
charge(f::Formula) = f.charge

function Formula(expr::AbstractString="")
    if expr == "Zz" || expr == "Zz+" || expr == "Zz⁺"
        composition = OrderedDict{Symbol,Int}()
        charge = 1
    elseif expr == "e" || expr == "e-" || expr == "e⁻"
        composition = OrderedDict{Symbol,Int}()
        charge = -1
    else
        composition = parse_formula(expr)
        charge = extract_charge(expr)
    end
    phreeqc_expr = unicode_to_phreeqc(expr)
    unicode_expr = phreeqc_to_unicode(replace(expr, r"\|\-?\d+\|" => ""))
    return Formula{valtype(composition)}(expr, phreeqc_expr, unicode_expr, colored_formula(unicode_expr), composition, charge)
end

function Formula(composition::AbstractDict{Symbol,T}, charge=0; order=ATOMIC_ORDER) where {T<:Number}
    charge_symbols = [:Zz, :Zz⁺, :e, :e⁻]
    # 1. Filter out charge symbols and convert to Vector for sorting
    filtered_keys = setdiff(keys(composition), charge_symbols)
    sorted_keys = sort(collect(filtered_keys), by=k -> findfirst(==(k), order))

    # 2. Build the formula string
    expr_parts = String[]
    uni_parts = String[]
    col_expr_parts = String[]
    for k in sorted_keys
        v = stoich_coef_round(composition[k])
        if !iszero(v)
            strv0 = get(dict_frac_unicode, v, string(v))
            strv = replace(strv0, " "=>"", "*"=>"")
            strvuni = strv
            colstrv = strv
            if occursin("+", strv0) || occursin("-", strv0) || occursin("*", strv0)
                strv = "(" * strv *")"
            end
            if any(x->x ∉ keys(dict_all_normal_to_sub), strv)
                strvuni = strv
                colstrv = string(COL_STOICH_INT(strv))
            else
                strvuni = all_normal_to_sub(strvuni)
                colstrv = string(COL_STOICH_INT(strvuni))
            end
            push!(expr_parts, string(k) * (isone(v) ? "" : strv))
            push!(uni_parts, string(k) * (isone(v) ? "" : strvuni))
            push!(col_expr_parts, string(k) * (isone(v) ? "" : colstrv))
        end
    end
    expr = join(expr_parts, "")
    uni = join(uni_parts, "")
    col_expr = join(col_expr_parts, "")

    # 3. Handle charge (e.g., Ca²⁺, SO₄²⁻)
    if iszero(charge)
        charge = get(composition, :Zz, 0) + get(composition, :Zz⁺, 0)
        charge -= get(composition, :e, 0) + get(composition, :e⁻, 0)
    end
    if !iszero(charge)
        sign = charge < 0 ? "-" : "+"
        abscharge = abs(charge)
        strch = isone(abscharge) ? "" : string(abscharge)
        expr *= sign * strch
        uni *= sign * strch
        col_expr *= string(COL_CHARGE(normal_to_super(sign * strch)))
    end

    # 4. Clean composition (remove charge keys)
    newcomposition = OrderedDict(k => v for (k,v) in composition if k ∉ charge_symbols && !iszero(v))

    return Formula{T}(expr, unicode_to_phreeqc(expr), uni, col_expr, newcomposition, charge)
end

function Formula(f::Formula)
    return Formula(composition(f))
end

Base.getindex(f::Formula{T}, i::Symbol) where {T} = get(composition(f), i, zero(T))
Base.length(f::Formula) = length(composition(f))

Base.isequal(f1::Formula, f2::Formula) = isequal(composition(f1), composition(f2)) && isequal(charge(f1), charge(f2))
==(f1::Formula, f2::Formula) = isequal(f1, f2)
Base.hash(f::Formula, h::UInt) = hash(composition(f), hash(charge(f), h))

function Base.show(io::IO, f::Formula)
    print(io, join(unique!([expr(f),  phreeqc(f), unicode(f), colored(f)]), " ∙ "))
end

print_formula(io::IO, f::Formula, title::String, pad::Int) = println(io, lpad(title, pad), ": ", join(unique!([expr(f),  phreeqc(f), unicode(f), colored(f)]), " ∙ "))

function Base.show(io::IO, ::MIME"text/plain", f::Formula)
    pad = 11
    println(io, typeof(f))
    print_formula(io, f, "formula", pad)
    println(io, lpad("composition", pad), ": ", join(["$k => $v" for (k, v) in composition(f)], ", "))
    println(io, lpad("charge", pad), ": ", charge(f))
end

function *(f::Formula, x::T) where {T<:Number}
    composition = OrderedDict(k => x*v for (k,v) in f.composition)
    return Formula(composition, f.charge)
end

function /(f::Formula, x::T) where {T<:Number}
    composition = OrderedDict(k => v/x for (k,v) in f.composition)
    return Formula(composition, f.charge)
end

function //(f::Formula, x::T) where {T<:Number}
    composition = OrderedDict(k => v//x for (k,v) in f.composition)
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
    composition = a.sym == b.sym ? OrderedDict{Symbol, promote_type(T, S)}(a.sym => a.coef+b.coef) : Dict{Symbol, promote_type(T, S)}(a.sym => a.coef, b.sym => b.coef)
    return Formula(composition, 0)
end

+(a::AtomGroup, b::Symbol) = a + AtomGroup(b)

+(a::Symbol, b::AtomGroup) = AtomGroup(a) + b

+(a::Symbol, b::Symbol) = AtomGroup(a) + AtomGroup(b)

function Base.convert(T::Type{<:Number}, f::Formula)
    newcomposition = OrderedDict(k => convert(T, v) for (k,v) ∈ composition(f))
    return Formula{T}(expr(f), phreeqc(f), unicode(f), colored(f), newcomposition, charge(f))
end

function apply(func::Function, f::Formula, args... ; kwargs...)
    tryfunc(v) = v isa Quantity ? (
        try func(ustrip(v), args...; kwargs...) * func(unit(v), args...; kwargs...); catch; try func(ustrip(v), args...; kwargs...) * unit(v); catch; v; end; end
        ) : (
        try func(v, args...; kwargs...); catch; v; end
        )
    newcomposition = OrderedDict(k => tryfunc(v) for (k,v) ∈ composition(f))
    return Formula{valtype(newcomposition)}(expr(f), phreeqc(f), unicode(f), colored(f), newcomposition, tryfunc(charge(f)))
end
