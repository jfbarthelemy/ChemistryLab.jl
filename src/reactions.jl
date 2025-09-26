struct Reaction{S<:AbstractSpecies, T<:Number}
    equation::String
    species_stoich::Dict{S, T}
    equal_sign::Union{Nothing,Char}
end

function Reaction(equation::AbstractString)
    reactants, equal_sign, products = parse_equation(equation)
    return Reaction(equation,
                    merge(Dict(Species(k) => -stoich_coef_round(v) for (k, v) in reactants),
                    Dict(Species(k) => stoich_coef_round(v) for (k, v) in products)),
                    equal_sign,
                    )
end

function Reaction(species_stoich::Dict{S, T}, equal_sign=nothing) where {S<:AbstractSpecies, T<:Number}
    equation = format_equation(Dict(symbol(k) => v for (k,v) in species_stoich); equal_sign=equal_sign)
    return Reaction{S,T}(equation,
                    species_stoich,
                    equal_sign,
                    )
end

function CemReaction(equation::AbstractString)
    reactants, equal_sign, products = parse_equation(equation)
    return Reaction(equation,
                    merge(Dict(CemSpecies(k) => -stoich_coef_round(v) for (k, v) in reactants),
                    Dict(CemSpecies(k) => stoich_coef_round(v) for (k, v) in products)),
                    equal_sign,
                    )
end

Base.convert(::Type{Reaction}, s::S) where {S<:AbstractSpecies} = Reaction(Dict(s => 1))
Base.convert(::Type{Reaction{U,T}}, s::S) where {U<:AbstractSpecies, T<:Number, S<:AbstractSpecies} = Reaction(Dict(s => 1))
Reaction(s::S) where {S<:AbstractSpecies} = Reaction(Dict(s => 1))
Reaction{U,T}(s::S) where {U<:AbstractSpecies, T<:Number, S<:AbstractSpecies} = Reaction(Dict(s => 1))

Reaction(r::R) where {R<:Reaction} = r
Reaction(r::R, equal_sign) where {R<:Reaction} = equal_sign == r.equal_sign ? r : Reaction(r.species_stoich, equal_sign)

function Reaction(species::S, candidate_primaries::Vector; scaling=1, equal_sign='=') where {S<:AbstractSpecies}
    if !isa(species, Vector) species = [species] end
    A, indep_comp, dep_comp = stoich_matrix(species, candidate_primaries; display=false)
    species_stoich = Dict{promote_type(typeof.(indep_comp)..., typeof.(dep_comp)...),eltype(A)}()
    for (i, s) in enumerate(indep_comp)
        species_stoich[s] = -A[i,1]*scaling
    end
    species_stoich[dep_comp[1]] = scaling
    return Reaction(species_stoich, equal_sign)
end

species_list(r::Reaction) = collect(keys(r.species_stoich))
stoich_list(r::Reaction) = collect(values(r.species_stoich))

*(ν::Number, s::AbstractSpecies) = Reaction(Dict(s => ν))

*(ν::Number, r::Reaction) = Reaction(Dict(k => ν*v for (k,v) in r.species_stoich), r.equal_sign)

-(s::AbstractSpecies) = Reaction(Dict(s => -1))

-(r::Reaction) = Reaction(Dict(k => -v for (k,v) in r.species_stoich), r.equal_sign)

function +(s::S1, t::S2) where {S1<:AbstractSpecies, S2<:AbstractSpecies}
    S = promote_type(S1, S2)
    s == t ? Reaction(Dict(S(s) => 2)) : Reaction(Dict(S(s) => 1, t => 1))
end

function +(r::R, s::S) where {R<:Reaction, S<:AbstractSpecies}
    dict = r.species_stoich
    species_stoich = Dict{promote_type(keytype(dict), S), valtype(dict)}()
    for (k, v) in dict
        species_stoich[k] = get(species_stoich, k, 0) + v
    end
    species_stoich[s] = get(species_stoich, s, 0) + 1
    return Reaction(species_stoich, r.equal_sign)
end

+(s::S, r::R) where {S<:AbstractSpecies, R<:Reaction} = +(r,s)

function +(r::R, u::U) where {R<:Reaction, U<:Reaction}
    dict1 = r.species_stoich
    dict2 = u.species_stoich
    S = promote_type(keytype(dict1), keytype(dict2))
    T = promote_type(valtype(dict1), valtype(dict2))
    species_stoich = Dict{S, T}()
    for (k, v) in dict1
        species_stoich[k] = get(species_stoich, k, 0) + v
    end
    for (k, v) in dict2
        species_stoich[k] = get(species_stoich, k, 0) + v
    end
    return Reaction(species_stoich, r.equal_sign)
end

const EQUAL_OPS = union(fwd_arrows[2:end], bwd_arrows[2:end], double_arrows, pure_rate_arrows, equal_signs[2:end])

for OP in Symbol.(EQUAL_OPS)
    @eval $OP(r,s) = Reaction(-Reaction(r)+Reaction(s), first(string($OP)))
end

function Base.show(io::IO, r::Reaction)
    print(io, r.equation)
end
