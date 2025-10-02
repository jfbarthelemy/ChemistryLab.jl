struct Reaction{SR<:AbstractSpecies, TR<:Number, SP<:AbstractSpecies, TP<:Number}
    equation::String
    colored::String
    reactants::Dict{SR, TR}
    products::Dict{SP, TP}
    equal_sign::Union{Nothing,Char}
end

equation(r::Reaction) = r.equation
colored(r::Reaction) = r.colored
reactants(r::Reaction) = r.reactants
products(r::Reaction) = r.products
equal_sign(r::Reaction) = r.equal_sign

function Reaction(equation::AbstractString)
    reactants, products, equal_sign = parse_equation(equation)
    return Reaction(equation,
                    colored_equation(equation),
                    Dict(Species(k) => stoich_coef_round(v) for (k, v) in reactants),
                    Dict(Species(k) => stoich_coef_round(v) for (k, v) in products),
                    equal_sign,
                    )
end

function CemReaction(equation::AbstractString)
    reactants, products, equal_sign = parse_equation(equation)
    return Reaction(equation,
                    colored_equation(equation),
                    Dict(CemSpecies(k) => stoich_coef_round(v) for (k, v) in reactants),
                    Dict(CemSpecies(k) => stoich_coef_round(v) for (k, v) in products),
                    equal_sign,
                    )
end

function split_species_by_stoich(species_stoich::AbstractDict{S, T}) where {S<:AbstractSpecies, T<:Number}
    reactants = Dict{S,T}()
    products = Dict{S,T}()
    for (species, coef) in species_stoich
        if !iszero(coef)
            if coef < 0
                reactants[species] = -stoich_coef_round(coef)
            else
                products[species] = stoich_coef_round(coef)
            end
        end
    end
    return reactants, products
end

function merge_species_by_stoich(reactants::AbstractDict{<:AbstractSpecies, <:Number}, products::AbstractDict{<:AbstractSpecies, <:Number})
    return merge(Dict(species => -stoich_coef_round(coef) for (species, coef) in reactants),
                    Dict(species => stoich_coef_round(coef) for (species, coef) in products))
end

function format_side(side::AbstractDict{S, T}) where {S<:AbstractSpecies, T<:Number}
    equation = String[]
    coleq = String[]
    ch = 0
    for (species, coef) in side
        if !iszero(coef)
            coeff_str = isone(coef) ? "" : string(stoich_coef_round(coef))
            if occursin("+", coeff_str) || occursin("-", coeff_str) || occursin("*", coeff_str)
                coeff_str = "(" * coeff_str *")"
            end
            push!(equation, coeff_str*unicode(species))
            push!(coleq, string(COL_STOICH_EXT(coeff_str))*colored(species))
            ch += coef*charge(species)
        end
    end
    if isempty(equation) equation = "∅" ; coleq = "∅" end
    return join(equation, " + "), join(coleq, " + "), ch
end

function Reaction(reactants::AbstractDict{SR, TR}, products::AbstractDict{SP, TP}, equal_sign=nothing) where {SR<:AbstractSpecies, TR<:Number, SP<:AbstractSpecies, TP<:Number}
    sreac, creac, charge_left = format_side(reactants)
    sprod, cprod, charge_right = format_side(products)
    if isnothing(equal_sign) equal_sign = '=' end
    
    charge_diff = charge_right - charge_left

    if !isapprox(charge_diff, 0, atol=1e-6)
        needed_e = charge_diff < 0 ? -stoich_coef_round(charge_diff) : stoich_coef_round(charge_diff)
        e_term = needed_e == 1 ? "e⁻" : "$needed_e" * "e⁻"

        if charge_diff < 0
            # Add e- to the left (reactants)
            sreac = isempty(sreac) ? e_term : "$sreac + $e_term"
            creac = isempty(creac) ? e_term : "$creac + $e_term"
        else
            # Add e- to the right (products)
            sprod = isempty(sprod) ? e_term : "$sprod + $e_term"
            cprod = isempty(cprod) ? e_term : "$cprod + $e_term"
        end
    end

    equation = sreac * " " * string(equal_sign) * " " * sprod
    colored = creac * " " * string(COL_PAR(string(equal_sign))) * " " * cprod

    return Reaction{SR,TR,SP,TP}(equation,
                             colored,
                             reactants,
                             products,
                             equal_sign,
                             )
end

function Reaction(species_stoich::Dict{S, T}, equal_sign::Union{Nothing,Char}=nothing) where {S<:AbstractSpecies, T<:Number}
    reactants, products = split_species_by_stoich(species_stoich)
    return Reaction(reactants, products, equal_sign)
end

Base.convert(::Type{Reaction}, s::S) where {S<:AbstractSpecies} = Reaction(Dict(s => 1))
Base.convert(::Type{Reaction{U,T}}, s::S) where {U<:AbstractSpecies, T<:Number, S<:AbstractSpecies} = Reaction(Dict(s => 1))
Reaction(s::S) where {S<:AbstractSpecies} = Reaction(Dict(s => 1))
Reaction{U,T}(s::S) where {U<:AbstractSpecies, T<:Number, S<:AbstractSpecies} = Reaction(Dict(s => 1))

Reaction(r::R) where {R<:Reaction} = r
Reaction(r::R, equal_sign) where {R<:Reaction} = equal_sign == r.equal_sign ? r : Reaction(reactants(r), products(r), equal_sign)

function Reaction(species::Vector{<:AbstractSpecies}; scaling=1, equal_sign='=')
    A, indep_comp, dep_comp = stoich_matrix(species[1:1], species[2:end]; display=false, involve_all_atoms=true)
    species_stoich = Dict{promote_type(typeof.(indep_comp)..., typeof.(dep_comp)...),eltype(A)}()
    for (i, s) in enumerate(indep_comp)
        species_stoich[s] = A[i,1]*scaling
    end
    species_stoich[dep_comp[1]] = -scaling
    return Reaction(species_stoich, equal_sign)
end

function Reaction(reac::Vector{<:AbstractSpecies}, prod::Vector{<:AbstractSpecies}; scaling=1, equal_sign='=')
    species = [reac ; prod]
    A, indep_comp, dep_comp = stoich_matrix(species[1:1], species[2:end]; display=false, involve_all_atoms=true)
    species_stoich = Dict{promote_type(typeof.(indep_comp)..., typeof.(dep_comp)...),eltype(A)}()
    for (i, s) in enumerate(indep_comp)
        species_stoich[s] = A[i,1]*scaling
    end
    species_stoich[dep_comp[1]] = -scaling
    return Reaction(Dict(k=>-v for (k,v) in species_stoich if k in reac), Dict(k=>v for (k,v) in species_stoich if k in prod), equal_sign)
end

*(ν::Number, s::AbstractSpecies) = Reaction(Dict(s => ν))

*(ν::Number, r::Reaction) = Reaction(Dict(k => ν*v for (k,v) in reactants(r)), Dict(k => ν*v for (k,v) in products(r)), r.equal_sign)

-(s::AbstractSpecies) = Reaction(Dict(s => -1))

-(r::Reaction) = Reaction(products(r), reactants(r), r.equal_sign)

function +(s::S1, t::S2) where {S1<:AbstractSpecies, S2<:AbstractSpecies}
    S = promote_type(S1, S2)
    s == t ? Reaction(Dict(S(s) => 2)) : Reaction(Dict(S(s) => 1, S(t) => 1))
end

function -(s::S1, t::S2) where {S1<:AbstractSpecies, S2<:AbstractSpecies}
    S = promote_type(S1, S2)
    s == t ? Reaction(Dict()) : Reaction(Dict(S(s) => 1, S(t) => -1))
end

function add_stoich(d1::AbstractDict{S1,T1}, d2::AbstractDict{S2,T2}) where {S1<:AbstractSpecies, T1<:Number, S2<:AbstractSpecies, T2<:Number}
    S = promote_type(S1, S2)
    T = promote_type(T1, T2)
    d = Dict{S,T}()
    for (k, v) in d1
        d[k] = get(d, k, 0) + v
    end
    for (k, v) in d2
        d[k] = get(d, k, 0) + v
    end
    return d
end

function +(r::R, s::S) where {R<:Reaction, S<:AbstractSpecies}
    return Reaction(reactants(r), add_stoich(products(r), Dict(s => 1)), r.equal_sign)
end

function -(r::R, s::S) where {R<:Reaction, S<:AbstractSpecies}
    return Reaction(reactants(r), add_stoich(products(r), Dict(s => -1)), r.equal_sign)
end

+(s::S, r::R) where {S<:AbstractSpecies, R<:Reaction} = +(r,s)

-(s::S, r::R) where {S<:AbstractSpecies, R<:Reaction} = +(s,-r)

function +(r::R, u::U) where {R<:Reaction, U<:Reaction}
    return Reaction(add_stoich(reactants(r), reactants(u)), add_stoich(products(r), products(u)), r.equal_sign)
end

function -(r::R, u::U) where {R<:Reaction, U<:Reaction}
    return Reaction(add_stoich(reactants(r), products(u)), add_stoich(products(r), reactants(u)), r.equal_sign)
end

const EQUAL_OPS = union(fwd_arrows[2:end], bwd_arrows[2:end], double_arrows, pure_rate_arrows, equal_signs[2:end])

for OP in Symbol.(EQUAL_OPS)
    @eval $OP(r,s) = Reaction(-Reaction(r)+Reaction(s), first(string($OP)))
end

function Base.show(io::IO, r::Reaction)
    print(io, colored(r))
end

function Base.map(func::Function, r::Reaction, args... ; kwargs...)
        reac = Dict(k => func(v, args... ; kwargs...) for (k,v) ∈ reactants(r))
        prod = Dict(k => func(v, args... ; kwargs...) for (k,v) ∈ products(r))
    return Reaction(reac, prod, equal_sign(r))
end
