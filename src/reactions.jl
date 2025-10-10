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

function Reaction(equation::AbstractString)
    reactants, products, equal_sign = parse_equation(equation)
    return Reaction(equation,
                    colored_equation(equation),
                    OrderedDict(Species(k) => stoich_coef_round(v) for (k, v) in reactants),
                    OrderedDict(Species(k) => stoich_coef_round(v) for (k, v) in products),
                    equal_sign,
                    OrderedDict{Symbol,PropertyType}()
                    )
end

function CemReaction(equation::AbstractString)
    reactants, products, equal_sign = parse_equation(equation)
    return Reaction(equation,
                    colored_equation(equation),
                    OrderedDict(CemSpecies(k) => stoich_coef_round(v) for (k, v) in reactants),
                    OrderedDict(CemSpecies(k) => stoich_coef_round(v) for (k, v) in products),
                    equal_sign,
                    OrderedDict{Symbol,PropertyType}()
                    )
end

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

function merge_species_by_stoich(reactants::AbstractDict{<:AbstractSpecies, <:Number}, products::AbstractDict{<:AbstractSpecies, <:Number})
    return merge(+, OrderedDict(species => -stoich_coef_round(coef) for (species, coef) in reactants),
                    OrderedDict(species => stoich_coef_round(coef) for (species, coef) in products))
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

function Reaction(reactants::AbstractDict{SR, TR}, products::AbstractDict{SP, TP}, properties::AbstractDict{Symbol,PropertyType}=OrderedDict{Symbol,PropertyType}() ; equal_sign='=', side::Symbol=:sign) where {SR<:AbstractSpecies, TR<:Number, SP<:AbstractSpecies, TP<:Number}

    if side != :sign
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

function Reaction(species_stoich::AbstractDict{S, T}, properties::AbstractDict{Symbol,PropertyType}=OrderedDict{Symbol,PropertyType}(); equal_sign::Char='=', side::Symbol=:sign) where {S<:AbstractSpecies, T<:Number}
    reactants, products = split_species_by_stoich(species_stoich; side=side)
    return Reaction(reactants, products, properties; equal_sign=equal_sign)
end

Base.convert(::Type{Reaction}, s::S) where {S<:AbstractSpecies} = Reaction(OrderedDict(s => 1))
Base.convert(::Type{Reaction{U,T}}, s::S) where {U<:AbstractSpecies, T<:Number, S<:AbstractSpecies} = Reaction(OrderedDict(s => 1))
Reaction(s::S) where {S<:AbstractSpecies} = Reaction(OrderedDict(s => 1))
Reaction{U,T}(s::S) where {U<:AbstractSpecies, T<:Number, S<:AbstractSpecies} = Reaction(OrderedDict(s => 1))

Reaction(r::R; equal_sign=r.equal_sign, side::Symbol=:none) where {R<:Reaction} = equal_sign == r.equal_sign && side == :none ? r : Reaction(reactants(r), products(r), properties(r); equal_sign=equal_sign, side=side)

function simplify_reaction(r::Reaction)
    reac = deepcopy(reactants(r))
    prod = deepcopy(products(r))
    common_species = intersect(keys(reac), keys(prod))
    for species in common_species
        coef = prod[species] - reac[species]
        if try coef>0 catch; true end
            prod[species] = coef
            pop!(reac, species)
        else
            reac[species] = -coef
            pop!(prod, species)
        end
    end
    return Reaction(reac, prod, properties(r); equal_sign=equal_sign(r))
end

function scale_stoich!(species_stoich::AbstractDict{<:AbstractSpecies, <:Number})
    v = values(species_stoich)
    if all(x -> x isa Integer || x isa Rational, v)
        mult = gcd(v...)
        for k in keys(species_stoich) species_stoich[k] *= mult end
    end
end

function Reaction(species::AbstractVector{<:AbstractSpecies}, properties::AbstractDict{Symbol,PropertyType}=OrderedDict{Symbol,PropertyType}(); scaling=1, equal_sign='=', side::Symbol=:sign, auto_scale=false)
    A, indep_comp, dep_comp = stoich_matrix(species[1:1], species[2:end]; display=false, involve_all_atoms=true)
    species_stoich = OrderedDict{promote_type(typeof.(indep_comp)..., typeof.(dep_comp)...),eltype(A)}()
    species_stoich[dep_comp[1]] = -scaling
    for (i, s) in enumerate(indep_comp)
        species_stoich[s] = A[i,1]*scaling
    end
    if auto_scale scale_stoich!(species_stoich) end
    return Reaction(species_stoich, properties; equal_sign=equal_sign, side=side)
end

function Reaction(reac::AbstractVector{<:AbstractSpecies}, prod::AbstractVector{<:AbstractSpecies}, properties::AbstractDict{Symbol,PropertyType}=OrderedDict{Symbol,PropertyType}(); scaling=1, equal_sign='=', side::Symbol=:sign, auto_scale=false)
    species = [reac ; prod]
    A, indep_comp, dep_comp = stoich_matrix(species[1:1], species[2:end]; display=false, involve_all_atoms=true)
    species_stoich = OrderedDict{promote_type(typeof.(indep_comp)..., typeof.(dep_comp)...),eltype(A)}()
    species_stoich[dep_comp[1]] = -scaling
    for (i, s) in enumerate(indep_comp)
        species_stoich[s] = A[i,1]*scaling
    end
    if auto_scale scale_stoich!(species_stoich) end
    if side != :sign
        return Reaction(species_stoich; equal_sign=equal_sign, side=side)
    else
        return Reaction(OrderedDict(k=>-v for (k,v) in species_stoich if k in reac), OrderedDict(k=>v for (k,v) in species_stoich if k in prod), properties; equal_sign=equal_sign)
    end
end

*(ν::Number, s::AbstractSpecies) = Reaction(OrderedDict(s => ν))

*(ν::Number, r::Reaction) = Reaction(OrderedDict(k => ν*v for (k,v) in reactants(r)), OrderedDict(k => ν*v for (k,v) in products(r)), properties(r); equal_sign=r.equal_sign)

-(s::AbstractSpecies) = Reaction(OrderedDict(s => -1))

-(r::Reaction) = Reaction(products(r), reactants(r), properties(r); equal_sign=r.equal_sign)

function +(s::S1, t::S2) where {S1<:AbstractSpecies, S2<:AbstractSpecies}
    S = promote_type(S1, S2)
    s == t ? Reaction(OrderedDict(S(s) => 2)) : Reaction(OrderedDict(S(s) => 1, S(t) => 1))
end

function -(s::S1, t::S2) where {S1<:AbstractSpecies, S2<:AbstractSpecies}
    S = promote_type(S1, S2)
    s == t ? Reaction(OrderedDict()) : Reaction(OrderedDict(S(s) => 1, S(t) => -1))
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
    return Reaction(reactants(r), add_stoich(products(r), OrderedDict(s => 1)), properties(r); equal_sign=r.equal_sign)
end

function -(r::R, s::S) where {R<:Reaction, S<:AbstractSpecies}
    return Reaction(reactants(r), add_stoich(products(r), OrderedDict(s => -1)), properties(r); equal_sign=r.equal_sign)
end

+(s::S, r::R) where {S<:AbstractSpecies, R<:Reaction} = +(r,s)

-(s::S, r::R) where {S<:AbstractSpecies, R<:Reaction} = +(s,-r)

function +(r::R, u::U) where {R<:Reaction, U<:Reaction}
    return Reaction(add_stoich(reactants(r), reactants(u)), add_stoich(products(r), products(u)), merge(properties(r), properties(u)); equal_sign=r.equal_sign)
end

function -(r::R, u::U) where {R<:Reaction, U<:Reaction}
    return Reaction(add_stoich(reactants(r), products(u)), add_stoich(products(r), reactants(u)), merge(properties(r), properties(u)); equal_sign=r.equal_sign)
end

const EQUAL_OPS = union(fwd_arrows[2:end], bwd_arrows[2:end], double_arrows, pure_rate_arrows, equal_signs[2:end])

for OP in Symbol.(EQUAL_OPS)
    @eval $OP(r,s) = Reaction(-Reaction(r)+Reaction(s); equal_sign=first(string($OP)), side=:sign)
end

function Base.show(io::IO, r::Reaction)
    print(io, colored(r))
    if length(properties(r)) > 0
        pad = 11
        print(io, "  [ ", join(["$k = $v" for (k, v) in properties(r)], " ;\t"), " ]")
    end
end

function apply(func::Function, r::Reaction{SR, TR, SP, TP}, args... ; kwargs...) where {SR<:AbstractSpecies, TR<:Number, SP<:AbstractSpecies, TP<:Number}
    tryfunc(v) =
        try
            func(ustrip(v), args...; kwargs...) * func(unit(v), args...; kwargs...)
        catch
            try
                func(v, args...; kwargs...)
            catch
                v
            end
        end
    reac = OrderedDict{SR, TR}(apply(func, k, args... ; kwargs..., name="", symbol ="") => tryfunc(v) for (k,v) ∈ reactants(r))
    prod = OrderedDict{SP, TP}(apply(func, k, args... ; kwargs..., name="", symbol ="") => tryfunc(v) for (k,v) ∈ products(r))
    newReaction = Reaction(reac, prod; equal_sign=equal_sign(r))
        for (k, v) in properties(r)
        newReaction[k] = tryfunc(v)
    end
    return newReaction

end

# function multiply_by_least_integer(v)
#     denominators = [denominator(x) for x in v if x isa Rational]
#     lcm_denominator = isempty(denominators) ? 1 : foldl(lcm, denominators)
#     result = lcm_denominator .* v
#     return result, lcm_denominator
# end
