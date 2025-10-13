struct Reaction{SR<:AbstractSpecies, TR<:Number, SP<:AbstractSpecies, TP<:Number}
    equation::String
    colored::String
    reactants::OrderedDict{SR, TR}
    products::OrderedDict{SP, TP}
    equal_sign::Char
    properties::OrderedDict{Symbol,PropertyType}
end

equation(r::Reaction) = r.equation
colored(r::Reaction) = r.colored
reactants(r::Reaction) = r.reactants
products(r::Reaction) = r.products
equal_sign(r::Reaction) = r.equal_sign
properties(r::Reaction) = r.properties

Base.getindex(r::Reaction, i::Symbol) = get(properties(r), i, nothing)

Base.setindex!(r::Reaction, value, i::Symbol) = setindex!(properties(r), value, i)

function Base.getproperty(r::Reaction, sym::Symbol)
    if sym in fieldnames(typeof(r))
        return getfield(r, sym)
    elseif sym in keys(properties(r))
        return properties(r)[sym]
    else
        error("Symbol '$sym' is neither a field nor a registered property.")
    end
end

function Base.setproperty!(r::Reaction, sym::Symbol, value)
    if sym in fieldnames(typeof(r))
        error("Cannot modify field '$sym' directly. Use constructor or dedicated methods.")
    else
        properties(r)[sym] = value
    end
    return r
end

function ordered_dict_with_default(gen, key_type, val_type)
    d = OrderedDict(gen)
    isempty(d) && (d = OrderedDict{key_type, val_type}())
    return d
end

function remove_zeros(d::AbstractDict)
    for (k,v) in d
        if iszero(v)
            pop!(d, k)
        end
    end
    return d
end

function Reaction(equation::AbstractString, S::Type{<:AbstractSpecies} = Species)
    reactants, products, equal_sign = parse_equation(equation)
    return Reaction(equation,
                    colored_equation(equation),
                    ordered_dict_with_default((S(k) => stoich_coef_round(v) for (k, v) in reactants if !iszero(v)), S, Number),
                    ordered_dict_with_default((S(k) => stoich_coef_round(v) for (k, v) in products if !iszero(v)), S, Number),
                    equal_sign,
                    OrderedDict{Symbol,PropertyType}()
                    )
end

CemReaction(equation::AbstractString) = Reaction(equation, CemSpecies)

function split_species_by_stoich(species_stoich::AbstractDict{S, T}; side::Symbol=:sign) where {S<:AbstractSpecies, T<:Number}
    reactants = OrderedDict{S,T}()
    products = OrderedDict{S,T}()
    for (species, coef) in species_stoich
        if !iszero(coef)
            if try (coef < 0 && side == :sign) || side == :reactants catch; false end
                reactants[species] = -stoich_coef_round(coef)
            else
                products[species] = stoich_coef_round(coef)
            end
        end
    end
    return reactants, products
end

function merge_species_by_stoich(reactants::AbstractDict{SR,TR}, products::AbstractDict{SP,TP}) where {SR<:AbstractSpecies, TR<:Number, SP<:AbstractSpecies, TP<:Number}
    return merge(+, ordered_dict_with_default((species => -stoich_coef_round(coef) for (species, coef) in reactants), SR, TR),
                    ordered_dict_with_default((species => stoich_coef_round(coef) for (species, coef) in products), SP, TP))
end

function format_side(side::AbstractDict{S, T}) where {S<:AbstractSpecies, T<:Number}
    equation = String[]
    coleq = String[]
    ch = 0
    Zz = root_type(S)("Zz")
    for (species, coef) in side
        if !iszero(coef) && species != Zz
            coeff_str = isone(coef) ? "" : string(stoich_coef_round(coef))
            coeff_str = add_parentheses_if_needed(coeff_str)
            coeff_str = replace(coeff_str, " "=>"", "*"=>"")
            # if occursin("+", coeff_str) || occursin("-", coeff_str) || occursin("*", coeff_str)
            #     coeff_str = "(" * coeff_str *")"
            # end
            push!(equation, coeff_str*unicode(species))
            push!(coleq, string(COL_STOICH_EXT(coeff_str))*colored(species))
            ch += coef*charge(species)
        end
    end
    if isempty(equation) equation = "∅" ; coleq = "∅" end
    return join(equation, " + "), join(coleq, " + "), ch
end

function Reaction(reactants::AbstractDict{SR, TR}, products::AbstractDict{SP, TP} ; equal_sign='=', properties::AbstractDict{Symbol,PropertyType}=OrderedDict{Symbol,PropertyType}(), side::Symbol=:none) where {SR<:AbstractSpecies, TR<:Number, SP<:AbstractSpecies, TP<:Number}

    if side != :none
        reactants, products = split_species_by_stoich(merge_species_by_stoich(reactants, products); side=side)
    end

    sreac, creac, charge_left = format_side(reactants)
    sprod, cprod, charge_right = format_side(products)
    
    charge_diff = charge_right - charge_left

    if !isapprox(charge_diff, 0, atol=1e-4)
        needed_e = charge_diff < 0 ? -stoich_coef_round(charge_diff) : stoich_coef_round(charge_diff)
        e_term = needed_e == 1 ? "e⁻" : "$needed_e" * "e⁻"
        ce_term = needed_e == 1 ? "e⁻" : string(COL_STOICH_EXT(add_parentheses_if_needed("$needed_e"))) * "e⁻"

        if charge_diff < 0
            # Add e- to the left (reactants)
            sreac = isempty(sreac) ? e_term : "$sreac + $e_term"
            creac = isempty(creac) ? e_term : "$creac + $ce_term"
        else
            # Add e- to the right (products)
            sprod = isempty(sprod) ? e_term : "$sprod + $e_term"
            cprod = isempty(cprod) ? e_term : "$cprod + $ce_term"
        end
    end

    equation = sreac * " " * string(equal_sign) * " " * sprod
    colored = creac * " " * string(COL_PAR(string(equal_sign))) * " " * cprod

    return Reaction(equation,
                    colored,
                    reactants,
                    products,
                    equal_sign,
                    properties
                    )
end

function Reaction(species_stoich::AbstractDict{S, T}; equal_sign::Char='=', properties::AbstractDict{Symbol,PropertyType}=OrderedDict{Symbol,PropertyType}(), side::Symbol=:sign) where {S<:AbstractSpecies, T<:Number}
    reactants, products = split_species_by_stoich(species_stoich; side=side)
    return Reaction(reactants, products; equal_sign=equal_sign, properties=properties)
end

Base.convert(::Type{Reaction}, s::S) where {S<:AbstractSpecies} = Reaction(OrderedDict(s => 1))
Base.convert(::Type{Reaction{U,T}}, s::S) where {U<:AbstractSpecies, T<:Number, S<:AbstractSpecies} = Reaction(OrderedDict(s => 1))
Reaction(s::S) where {S<:AbstractSpecies} = Reaction(OrderedDict(s => 1))
Reaction{U,T}(s::S) where {U<:AbstractSpecies, T<:Number, S<:AbstractSpecies} = Reaction(OrderedDict(s => 1))

Reaction(r::R; equal_sign=r.equal_sign, properties=r.properties, side::Symbol=:none) where {R<:Reaction} = Reaction(reactants(r), products(r); equal_sign=equal_sign, side=side, properties=properties)

function simplify_reaction(r::Reaction)
    reac = remove_zeros(deepcopy(reactants(r)))
    prod = remove_zeros(deepcopy(products(r)))
    common_species = intersect(keys(reac), keys(prod))
    for species in common_species
        coef = prod[species] - reac[species]
        if iszero(coef)
            pop!(reac, species)
            pop!(prod, species)
        elseif try coef>0 catch; true end
            prod[species] = coef
            pop!(reac, species)
        else
            reac[species] = -coef
            pop!(prod, species)
        end
    end
    return Reaction(reac, prod; equal_sign=equal_sign(r), properties=properties(r))
end

function scale_stoich!(species_stoich::AbstractDict{<:AbstractSpecies, <:Number})
    v = values(species_stoich)
    if all(x -> x isa Integer || x isa Rational, v)
        mult = gcd(v...)
        for k in keys(species_stoich) species_stoich[k] *= mult end
    end
end

function build_species_stoich(species::AbstractVector{<:AbstractSpecies}; scaling=1, auto_scale=false)
    A, indep_comp, dep_comp = stoich_matrix(species[1:1], species[2:end]; display=false, involve_all_atoms=true)
    S, T = promote_type(typeof.(indep_comp)..., typeof.(dep_comp)...), eltype(A)
    species_stoich = OrderedDict{S,T}()
    species_stoich[dep_comp[1]] = -scaling
    for (i, s) in enumerate(indep_comp)
        species_stoich[s] = A[i,1]*scaling
    end
    if auto_scale scale_stoich!(species_stoich) end
    return species_stoich
end

function Reaction(species::AbstractVector{<:AbstractSpecies}; equal_sign='=', properties::AbstractDict{Symbol,PropertyType}=OrderedDict{Symbol,PropertyType}(), scaling=1, auto_scale=false, side::Symbol=:sign)
    species_stoich = build_species_stoich(species; scaling=scaling, auto_scale=auto_scale)
    return Reaction(species_stoich; equal_sign=equal_sign, properties=properties, side=side)
end

function Reaction(reac::AbstractVector{<:AbstractSpecies}, prod::AbstractVector{<:AbstractSpecies}; equal_sign='=', properties::AbstractDict{Symbol,PropertyType}=OrderedDict{Symbol,PropertyType}(), scaling=1, auto_scale=false, side::Symbol=:none)
    species = [reac ; prod]
    species_stoich = build_species_stoich(species; scaling=scaling, auto_scale=auto_scale)
    S, T = keytype(species_stoich), valtype(species_stoich)
    if side != :none
        return Reaction(species_stoich; equal_sign=equal_sign, properties=properties, side=side)
    else
        return Reaction(ordered_dict_with_default((k=>-v for (k,v) in species_stoich if k in reac), S, T), 
                        ordered_dict_with_default((k=>v for (k,v) in species_stoich if k in prod), S, T)
                        ; 
                        equal_sign=equal_sign, properties=properties, side=:none)
    end
end

*(ν::Number, s::AbstractSpecies) = Reaction(OrderedDict(s => ν))

*(ν::Number, r::Reaction{SR,TR,SP,TP}) where {SR<:AbstractSpecies, TR<:Number, SP<:AbstractSpecies, TP<:Number} = Reaction(ordered_dict_with_default((k => ν*v for (k,v) in reactants(r)), SR, TR), ordered_dict_with_default((k => ν*v for (k,v) in products(r)), SP, TP); equal_sign=r.equal_sign, properties=r.properties)

-(s::AbstractSpecies) = Reaction(OrderedDict(s => -1))

-(r::Reaction) = Reaction(products(r), reactants(r); equal_sign=r.equal_sign, properties=r.properties)

function +(s::S1, t::S2) where {S1<:AbstractSpecies, S2<:AbstractSpecies}
    S = promote_type(S1, S2)
    s == t ? Reaction(OrderedDict(S(s) => 2)) : Reaction(OrderedDict(S(s) => 1, S(t) => 1))
end

function -(s::S1, t::S2) where {S1<:AbstractSpecies, S2<:AbstractSpecies}
    S = promote_type(S1, S2)
    s == t ? Reaction(OrderedDict{S,Number}()) : Reaction(OrderedDict(S(s) => 1, S(t) => -1))
end

function add_stoich(d1::AbstractDict{S1,T1}, d2::AbstractDict{S2,T2}) where {S1<:AbstractSpecies, T1<:Number, S2<:AbstractSpecies, T2<:Number}
    S = promote_type(S1, S2)
    T = promote_type(T1, T2)
    d = OrderedDict{S,T}()
    for (k, v) in d1
        d[k] = get(d, k, 0) + v
    end
    for (k, v) in d2
        d[k] = get(d, k, 0) + v
    end
    return d
end

function +(r::R, s::S) where {R<:Reaction, S<:AbstractSpecies}
    return Reaction(reactants(r), add_stoich(products(r), OrderedDict(s => 1)); equal_sign=r.equal_sign, properties=r.properties)
end

function -(r::R, s::S) where {R<:Reaction, S<:AbstractSpecies}
    return Reaction(reactants(r), add_stoich(products(r), OrderedDict(s => -1)); equal_sign=r.equal_sign, properties=r.properties)
end

+(s::S, r::R) where {S<:AbstractSpecies, R<:Reaction} = +(r,s)

-(s::S, r::R) where {S<:AbstractSpecies, R<:Reaction} = +(s,-r)

function +(r::R, u::U) where {R<:Reaction, U<:Reaction}
    return Reaction(add_stoich(reactants(r), reactants(u)), add_stoich(products(r), products(u)); equal_sign=r.equal_sign, properties=merge(properties(r), properties(u)))
end

function -(r::R, u::U) where {R<:Reaction, U<:Reaction}
    return Reaction(add_stoich(reactants(r), products(u)), add_stoich(products(r), reactants(u)); equal_sign=r.equal_sign, properties=merge(properties(r), properties(u)))
end

const EQUAL_OPS = union(fwd_arrows[2:end], bwd_arrows[2:end], double_arrows, pure_rate_arrows, equal_signs[2:end])

for OP in Symbol.(EQUAL_OPS)
    @eval $OP(r,s) = Reaction(-Reaction(r)+Reaction(s); equal_sign=first(string($OP)), side=:sign)
end

function Base.show(io::IO, r::Reaction)
    print(io, colored(r))
    # if length(properties(r)) > 0
    #     pad = 11
    #     print(io, "  [ ", join(["$k = $v" for (k, v) in properties(r)], " ; "), " ]")
    # end
end

function Base.show(io::IO, ::MIME"text/plain", r::Reaction)
    println(io, colored(r))
    pad = 10
    if length(reactants(r)) > 0
        println(io, lpad("reactants", pad), ": ", join(["$(name(k))→$v" for (k, v) in reactants(r)], ", "))
        # println(io, lpad("reactants", pad), ": ", join(["$(name(k))→$v" for (k, v) in reactants(r)], "\n" * repeat(" ", pad + 2)))
    end
    if length(products(r)) > 0
        println(io, lpad("products", pad), ": ", join(["$(name(k))→$v" for (k, v) in products(r)], ", "))
        # println(io, lpad("products", pad), ": ", join(["$(name(k))→$v" for (k, v) in products(r)], "\n" * repeat(" ", pad + 2)))
    end
    if length(properties(r)) > 0
        print(io, lpad("properties", pad), ": ", join(["$k = $v" for (k, v) in properties(r)], "\n" * repeat(" ", pad + 2)))
    end
end

function apply(func::Function, r::Reaction{SR, TR, SP, TP}, args... ; kwargs...) where {SR<:AbstractSpecies, TR<:Number, SP<:AbstractSpecies, TP<:Number}
    tryfunc(v) = v isa Quantity ? (
        try func(ustrip(v), args...; kwargs...) * func(unit(v), args...; kwargs...); catch; try func(ustrip(v), args...; kwargs...) * unit(v); catch; v; end; end
        ) : (
        try func(v, args...; kwargs...); catch; v; end
        )
    reac = OrderedDict{SR, TR}(apply(func, k, args... ; kwargs..., name="", symbol ="") => tryfunc(v) for (k,v) ∈ reactants(r))
    prod = OrderedDict{SP, TP}(apply(func, k, args... ; kwargs..., name="", symbol ="") => tryfunc(v) for (k,v) ∈ products(r))
    newReaction = Reaction(reac, prod; equal_sign=get(kwargs, :equal_sign, equal_sign(r)), 
                                       properties=get(kwargs, :properties, properties(r)),
                                       side=get(kwargs, :side, :none))
        for (k, v) in properties(r)
        newReaction[k] = tryfunc(v)
    end
    return newReaction

end
