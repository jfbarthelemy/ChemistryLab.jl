using ChemistryLab
using DynamicQuantities
using Symbolics
using Plots

# Formula
fgen = Formula("A1//2B3C0.4")
convert(Float64, fgen)
ATOMIC_ORDER # provides the order of atoms in formulas

# Species
H₂O = Species("H₂O"; name = "Water", symbol = "H₂O@", aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
HSO₄⁻ = Species("HSO₄⁻"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
CO₂ = Species(Dict(:C => 1, :O => 2); name = "Carbon dioxide", symbol = "CO₂⤴", aggregate_state = AS_GAS, class = SC_GASFLUID)
species = [H₂O, HSO₄⁻, CO₂] ;
CSM = CanonicalStoichMatrix(species)
pprint(CSM; label = :name)
pprint(CSM; label = :symbol)
pprint(CSM; label = :formula)

water_without_name_symbol = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
water_without_name_symbol == H₂O # true since atoms, aggregate_state and class are equal despite instances are different
vapour = Species("H₂O"; name = "Vapour", symbol = "H₂O⤴", aggregate_state = AS_GAS, class = SC_GASFLUID)
vapour == H₂O # false since aggregate_state or class are different despite atoms are identical

# CemSpecies
OXIDE_ORDER # provides the order of oxides in cement formulas
C3S = CemSpecies("C3S"; name = "Alite", symbol = "C₃S", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
C2S = CemSpecies("C₂S"; name = "Belite", symbol = "C₂S", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
C3A = CemSpecies("C3A"; name = "Aluminate", symbol = "C₃A", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
C4AF = CemSpecies(Dict(:C => 4, :A => 1, :F => 1); name = "Ferrite", symbol = "C₄AF", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
cemspecies = [C3S, C2S, C3A, C4AF]
CSM = CanonicalStoichMatrix(cemspecies)
pprint(CSM; label = :name)
pprint(CSM; label = :symbol)
pprint(CSM; label = :formula)
# conversion CemSpecies → Species always possible
spC3S = Species(C3S)
spC3S == C3S # true since atoms, aggregate_state and class are identical
# conversion Species → CemSpecies possible only if the species can be decomposed in cement oxides
spCH = Species("Ca(OH)2"; name = "Portlandite", symbol = "CH", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
CH = CemSpecies(spCH)
spCH == CH # true
spCH == CemSpecies("CH") # false since aggregate_state and class are undef
spCH == CemSpecies("CH"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT) # true even though names and/or symbols do not coincide
try
    CemSpecies(Species("Ca(OH)"))
catch
    "ERROR: Ca(OH) cannot be decomposed in cement oxides"
end
CemSpecies(Species("CaCO3"; name = "Calcite", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)) # ok here

# Thermofun input
println("LOADING DATABASES...")
df_elements, df_substances, df_reactions = read_thermofun_database(datapath("cemdata18-merged.json"))
species_vec = build_species(df_substances)
dict_species = Dict(symbol(s) => s for s in species_vec)
dict_reactions = Dict(symbol(r) => r for r in build_reactions(df_reactions, species_vec))

df_elements_psi, df_substances_psi, df_reactions_psi = read_thermofun_database(datapath("psinagra-12-07-thermofun.json"))
species_psi = build_species(df_substances_psi)
dict_species_psi = Dict(symbol(s) => s for s in species_psi)
dict_reactions_psi = Dict(symbol(r) => r for r in build_reactions(df_reactions_psi, species_psi))

df_elements_aq17, df_substances_aq17, df_reactions_aq17 = read_thermofun_database(datapath("aq17-thermofun.json"))
species_aq17 = build_species(df_substances_aq17)
dict_species_aq17 = Dict(symbol(s) => s for s in species_aq17)
dict_reactions_aq17 = Dict(symbol(r) => r for r in build_reactions(df_reactions_aq17, species_aq17))

df_elements_orga, df_substances_orga, df_reactions_orga = read_thermofun_database(datapath("slop98-organic-thermofun.json"))
species_orga = build_species(df_substances_orga)
dict_species_orga = Dict(symbol(s) => s for s in species_orga)
dict_reactions_orga = Dict(symbol(r) => r for r in build_reactions(df_reactions_orga, species_orga))

df_elements_inorga, df_substances_inorga, df_reactions_inorga = read_thermofun_database(datapath("slop98-inorganic-thermofun.json"))
species_inorga = build_species(df_substances_inorga)
dict_species_inorga = Dict(symbol(s) => s for s in species_inorga)
dict_reactions_inorga = Dict(symbol(r) => r for r in build_reactions(df_reactions_inorga, species_inorga))

# Extraction of primaries from .dat
df_primaries = extract_primary_species(datapath("CEMDATA18-31-03-2022-phaseVol.dat"))

# Construction of stoich matrix with species from database
function get_secondaries(all_species, atomlist::AbstractVector{<:Symbol}, aggregate_states = [AS_AQUEOUS], excluded_species = [])
    return filter(
        x -> x.aggregate_state ∈ aggregate_states
            && all(in.(keys(atoms(x)), Ref(atomlist)))
            && x ∉ excluded_species,
        collect(values(all_species))
    )
end
function get_secondaries(all_species, species_list, aggregate_states = [AS_AQUEOUS], excluded_species = [])
    return get_secondaries(all_species, union_atoms(atoms.(collect(values(species_list)))), aggregate_states, excluded_species)
end
given_species = collect(values(build_species(df_substances, split("C2S Portlandite Jennite H2O@"))))
secondaries = get_secondaries(dict_species, given_species, [AS_AQUEOUS], [dict_species["H2@"], dict_species["O2@"]])
species = unique([given_species; secondaries])
candidate_primaries = [s == "Zz" ? Species("Zz") : dict_species[s] for s in df_primaries.symbol]
SM = StoichMatrix(species, candidate_primaries)
pprint(SM)
list_reactions = reactions(SM)
pprint(list_reactions)

# Same calculation with filtered databases
df_given_species = filter(row -> row.symbol in split("C2S Portlandite Jennite H2O@"), df_substances)
involved_atoms = union_atoms(parse_formula.(df_given_species.formula))
df_secondaries = filter(
    row -> last(only(row.aggregate_state)) == "AS_AQUEOUS"
        && all(in.(union_atoms([parse_formula(row.formula)]), Ref(involved_atoms)))
        && row.symbol ∉ split("H2@ O2@"), df_substances
)
df_union = unique(vcat(df_given_species, df_secondaries))
# Same as above in one command
df_union = get_compatible_species(
    df_substances, split("C2S Portlandite Jennite H2O@");
    aggregate_states = [AS_AQUEOUS], exclude_species = split("H2@ O2@"), union = true
)
all_species = build_species(df_union)
dict_all_species = Dict(symbol(s) => s for s in all_species)
candidate_primaries = collect(skipmissing(get.(Ref(dict_all_species), CEMDATA_PRIMARIES, missing)))
candidate_primaries = get.(Ref(dict_all_species), intersect(keys(dict_all_species), CEMDATA_PRIMARIES), missing)
candidate_primaries = [dict_all_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_all_species, s)]

SM = StoichMatrix(dict_all_species, candidate_primaries)
pprint(SM)
list_reactions = reactions(SM)
pprint(list_reactions)
# Filtering reactions with involved species from database
df_involved_reactions = filter(row -> all(in.([x.symbol for x in row.reactants], Ref(expr.(values(dict_all_species))))), df_reactions)
dict_involved_reactions = build_reactions(df_involved_reactions, all_species)


df_hydrates_clinker = get_compatible_species(df_substances, split("C3S C2S C3A C4AF Gp H2O@"); aggregate_states = [AS_CRYSTAL])


# Construction of stoich matrix with aqueous species from database
aqueous_species = filter(x -> last(x).aggregate_state == AS_AQUEOUS, dict_species)
species = collect(values(aqueous_species))
candidate_primaries = [s == "Zz" ? Species("Zz") : dict_species[s] for s in df_primaries.symbol]
SM = StoichMatrix(species, candidate_primaries)
pprint(SM)
list_reactions = reactions(SM)
pprint(list_reactions)

# Conversion to cement notation with species from database
CemSpecies(H₂O)
CemSpecies(vapour)
Jennite = dict_species["Jennite"]
cemJennite = CemSpecies(Jennite)
Jennite == cemJennite

# Equation parsing
equation = "13H⁺ + NO₃⁻ + CO₃²⁻ + 10e⁻ = 6H₂O@ + HCN@"
## simple parsing giving Dicts
reac, prod, equal_sign = parse_equation(equation)
## construction of a Reaction struct from a string
r = Reaction(equation)
## simplification of a Reaction
r = Reaction("2H₂O → H⁺ + OH⁻ + H₂O")
println(simplify_reaction(r))
## construction of a Reaction of only CemSpecies by CemReaction
eqC3S = "C₃S + 5.3H = 1.3CH + C₁.₇SH₄"
rC3S = CemReaction(eqC3S)
## construction of a Reaction by operations on Species
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
CSH = CemSpecies("C1.7SH4")
r = C3S + 5.3H ↔ 1.3CH + CSH
## construction of a Reaction by a balance calculation
r = Reaction([C3S, H, CH, CSH]; equal_sign = '→')
Reaction(CemSpecies.(["C3S", "H", "CH", "C1.8SH4"]))
for c_over_s in 1.5:0.1:2.0
    println(Reaction(CemSpecies.(["C3S", "H", "CH", "C$(c_over_s)SH4"])))
end

# Chen & Brouwers
@variables a b g
CSH = CemSpecies(Dict(:C => a, :S => one(a), :A => b, :H => g))
HT = CemSpecies("M₅AH₁₃")
HG = CemSpecies("C₆AFS₂H₈")
AFt = CemSpecies("C₆AS̄₃H₃₂")
ST = CemSpecies("C₂ASH₈")
AH = CemSpecies("C₄AH₁₃")
CSM = CanonicalStoichMatrix([CSH, HT, HG, AFt, ST, AH]); pprint(CSM)
A, ox = CSM.A, CSM.primaries
A = typeof(a).(A[1:(end - 1), 1:end]) # end-1 to remove the line corresponding to water H
oxides = (CemSpecies ∘ string).(ox[1:(end - 1)])
hydrates = [CSH, HT, HG, AFt, ST, AH]
pprint(A, symbol.(oxides), symbol.(hydrates))
pprint(inv(A), symbol.(hydrates), symbol.(oxides))
Mhyd = ustrip.(us"g/mol", getproperty.(hydrates, :M))
Mox = ustrip.(us"g/mol", getproperty.(oxides, :M))
B = Mox .* A .* inv.(Mhyd)'
# or directly
mCSM = mass_matrix(CSM)
B, ox = mCSM.A, mCSM.primaries
B = B[1:(end - 1), 1:end] # to remove the H line
pprint(B, "m_" .* symbol.(oxides), "m_" .* symbol.(hydrates))
pprint(substitute.(inv(B), Ref(Dict(a => 1.8, b => 1, g => 4))), "m_" .* symbol.(hydrates), "m_" .* symbol.(oxides))

# Alkane combustion with Symbolics
@variables n
CₙH₂ₙ₊₂ = Species(:C => n, :H => 2n + 2)
O₂, H₂O, CO₂ = Species.(split("O₂ H₂O CO₂")) ;
r = Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂])
for vn in 1:9
    print("n=$vn ⇒ "); println(colored(apply(substitute, r, n => vn)))
end
for vn in 1:9
    print("n=$vn ⇒ "); println(colored(apply(Symbolics.symbolic_to_float ∘ stoich_coef_round ∘ (x -> x isa Num ? x.val : x) ∘ substitute, r, n => vn)))
end

# Example from https://github.com/thermohub/chemicalfun
formulas = ["Ca+2", "Fe+2", "Fe|3|+3", "H+", "OH-", "SO4-2", "CaSO4@", "CaOH+", "FeO@", "HFe|3|O2@", "FeOH+", "Fe|3|OH+2", "H2O@", "FeS|-2|", "FeS|0|S|-2|", "S|4|O2"] ;
species = Species.(formulas) ;
candidate_primaries = species[1:6] ;
SM = StoichMatrix(species)
pprint(SM)
list_reactions = reactions(SM)
pprint(list_reactions)

lT = ((0:1:100) .+ 273.15)
@time plot(lT, (T -> dict_species["Jennite"].ΔₐH⁰(T = T)).(lT))

# Check which species involved in reactions have not been previously constructed in the list of substances (in this case they are built on-the-fly and don't have thermo properties)
for re in values(dict_reactions)
    for s in keys(re.reactants)
        if !haskey(s, :Cp⁰)
            println(colored(s), "  ∙  ", re.symbol)
        end
    end
end

# Check consistency of logKr at Tref and logKr_Tref of the database
for r in values(dict_reactions)
    println(colored(r), " → ", r.logKr(), " == ", r.logKr_Tref)
end

r = dict_reactions["Cal"]
Tref = 298.15u"K"
ΔᵣG⁰(T) = sum(ν * s.ΔₐG⁰(T = T) for (s, ν) in r)
logK = -ΔᵣG⁰(Tref) / (Constants.R * Tref) / log(10)
K(T) = 10^(-ΔᵣG⁰(T) / (Constants.R * T) / log(10))
pK(T) = ΔᵣG⁰(T) / (Constants.R * T) / log(10)
r.logKr()
plot(ustrip.(lT), ustrip.(ΔᵣG⁰.(lT)))
for (re, ν) in r.reactants
    println(re, " ΔₐG⁰=", re.ΔₐG⁰(T = Tref))
end
for (pr, ν) in r.products
    println(pr, " ΔₐG⁰=", pr.ΔₐG⁰(T = Tref))
end
for (s, ν) in r
    println(s, " ΔₐG⁰=", s.ΔₐG⁰(T = Tref))
end

Tref = 298.15
for r in values(dict_reactions)
    if true # abs((r.logKr_Tref -r.logKr())/r.logKr_Tref)>0.01
        println(collect(keys(r))[1].symbol)
        println(colored(r))
        println("     → logKr given at $(Tref) = ", r.logKr_Tref)
        println("     → logKr calculated at $(Tref) = ", r.logKr())
        try
            println("     → logKr=-ΔᵣG⁰(Tref)/(R Tref ln(10)) = ", -r.ΔᵣG⁰(T = Tref) / ustrip(Constants.R * Tref) / log(10))
        catch
            println("     → logKr=-ΔᵣG⁰(Tref)/(R Tref ln(10)) = XXXXXXXXXXXXXXXXXXXXX")
        end
    end
end

# Construction of a Reaction from a string and a list of species (possibly with thermodynamic properties from a database)
r = Reaction("NaOH → Na+ + OH-"; species_list = dict_species)
for re in keys(r.reactants)
    display(re)
end
for pr in keys(r.products)
    display(pr)
end

# Example from Leal et al. 2017
H₂O = dict_species["H2O@"]
H⁺ = dict_species["H+"]
OH⁻ = dict_species["OH-"]
CO₂ = dict_species["CO2@"]
HCO₃⁻ = dict_species["HCO3-"]
CO₃²⁻ = dict_species["CO3-2"]
for sp in (:H₂O, :H⁺, :OH⁻, :CO₂, :HCO₃⁻, :CO₃²⁻)
    symsp = String(sp)
    @eval $sp = Species($(String(sp)))
end
CSM = CanonicalStoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻]); pprint(CSM)
SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻]); pprint(SM)

# Example from Miron et al. 2023, https://doi.org/10.21105/joss.04624
p1 = plot(xlabel = "Temperature [°C]", ylabel = "Cp⁰ [J K⁻¹]", xticks = 0:50:250, yticks = -2000:200:2000)
for sp in getindex.(Ref(dict_species_aq17), split("Na+ Ca+2 SiO2@ CO3-2 OH-"))
    plot!(p1, θ -> sp.Cp⁰(T = 273.15 + θ), 0:0.1:250, label = unicode(sp))
end
p2 = plot(xlabel = "Temperature [°C]", ylabel = "log₁₀K⁰", xticks = 0:50:250, yticks = -20:2:10)
for lsp in [
        "Calcite Ca+2 CO3-2",
        "H2O@ H+ OH-",
        "NaCl@ Na+ Cl-",
        "Al+3 H2O@ AlOH+2 H+",
    ]
    rr = Reaction(getindex.(Ref(dict_species_aq17), split(lsp)))
    plot!(p2, θ -> rr.logK⁰(T = 273.15 + θ), 0:0.1:250, label = rr.equation)
end
plot(p1, p2, layout = (1, 2), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm, size = (900, 400))

# Vapor pressure of water
l = dict_species_aq17["H2O@"]
v = dict_species_aq17["H2O"]
T0 = 373.15
ΔHₗᵥ = v.ΔₐH⁰(T = T0) - l.ΔₐH⁰(T = T0)
Rankine(T) = 13.7 - 5120 / T
lnP(T) = ΔHₗᵥ / ustrip(Constants.R) * (1 / T0 - 1 / T)
lθ = 0:1:200
plot(xlabel = "T [°C]", ylabel = "ln(P/P0)")
plot!(θ -> Rankine(273.15 + θ), lθ, label = "Rankine")
plot!(θ -> lnP(273.15 + θ), lθ, label = "ΔHₗᵥ/Constants.R*(1/T0-1/T)")

# using JSON
# json_str = JSON.json(df_substances)           # DataFrame → chaîne JSON
# open("data/test.json", "w") do f        # écriture dans un fichier
#     JSON.print(f, json_str, 4)
# end
