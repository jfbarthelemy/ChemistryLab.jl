using Revise, ChemistryLab, Unicode
import Unitful: u, g, cm, K, J, mol, Quantity, uconvert, ustrip, unit, uparse
using SymPy
import Symbolics: @variables, substitute, value

# Formula
fgen = Formula("A1//2B3C0.4")
convert(Float64, fgen)
ATOMIC_ORDER # provides the order of atoms in formulas
fCO2 = :C+2*:O # works only with a part of real atoms registered in ATOMIC_ORDER
fHSO₄⁻ = :H+:S+4*:O+:e
fNa⁺ = :Na+:Zz

# Species
fH₂O = 2*:H + :O
H₂O = Species(fH₂O; name="Water", symbol="H₂O&", aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)
HSO₄⁻ = Species("HSO₄⁻"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)
CO₂ = Species(Dict(:C=>1, :O=>2); name="Carbon dioxide", symbol="CO₂⤴", aggregate_state=AS_GAS, class=SC_GAS_FLUID)
species = [H₂O, HSO₄⁻, CO₂] ;
A, atomlist = canonical_stoich_matrix(species; label=:name) ; # label only for display
A, atomlist = canonical_stoich_matrix(species; label=:symbol) ;
A, atomlist = canonical_stoich_matrix(species; label=:formula) ;

water_without_name_symbol = Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)
water_without_name_symbol == H₂O # true since atoms, aggregate_state and class are equal despite instances are different
vapour = Species(fH₂O; name="Vapour", symbol="H₂O⤴", aggregate_state=AS_GAS, class=SC_GAS_FLUID)
vapour == H₂O # false since aggregate_state or class are different despite atoms are identical

# CemSpecies
OXIDE_ORDER # provides the order of oxides in cement formulas
C3S = CemSpecies("C3S"; name="Alite", symbol="C₃S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C2S = CemSpecies("C₂S"; name="Belite", symbol="C₂S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C3A = CemSpecies("C3A"; name="Aluminate", symbol="C₃A", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C4AF = CemSpecies(Dict(:C=>4, :A=>1, :F=>1); name="Ferrite", symbol="C₄AF", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
cemspecies = [C3S, C2S, C3A, C4AF]
A, indep_comp = canonical_stoich_matrix(cemspecies; label=:name) ;
A, indep_comp = canonical_stoich_matrix(cemspecies; label=:symbol) ;
A, indep_comp = canonical_stoich_matrix(cemspecies; label=:formula) ;
 # conversion CemSpecies → Species always possible
spC3S = Species(C3S)
spC3S == C3S # true since atoms, aggregate_state and class are identical
 # conversion Species → CemSpecies possible only if the species can be decomposed in cement oxides
spCH = Species("Ca(OH)2"; name="Portlandite", symbol="CH", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
CH = CemSpecies(spCH)
spCH == CH # true
spCH == CemSpecies("CH") # false since aggregate_state and class are undef
spCH == CemSpecies("CH"; aggregate_state=AS_CRYSTAL, class=SC_COMPONENT) # true even though names and/or symbols do not coincide 
try CemSpecies(Species("Ca(OH)")) catch; "ERROR: Ca(OH) cannot be decomposed in cement oxides" end
CemSpecies(Species("CaCO3"; name="Calcite", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)) # ok here

# Thermofun cemdata18
df_elements, df_substances, df_reactions = read_thermofun("data/cemdata18-merged.json"; debug=false, with_units=true, all_properties=true) # debug only for conception phase (not to be put in the doc)
dict_species = Dict(zip(df_substances.symbol, df_substances.species))
# filter(p->!haskey(p.second, :Cp), dict_species)
# filter(p->haskey(p.second, :Cp) && !iszero(p.second.Cp.a3), dict_species)
df_primaries = extract_primary_species("data/CEMDATA18-31-03-2022-phaseVol.dat")

# Construction of stoich matrix with species from database
given_species = filter(row->row.symbol ∈ split("C3S Portlandite Jennite H2O@"), df_substances)
secondaries = filter(row->row.aggregate_state == "AS_AQUEOUS" 
                          && all(k->first(k) ∈ union_atoms(atoms.(given_species.species)), atoms(row.species))
                          && row.symbol ∉ split("H2@ O2@"),
                          df_substances)
all_species = unique(vcat(given_species, secondaries), :symbol)
# species = Species.(all_species.formula)
# candidate_primaries = Species.(df_primaries.formula)
species = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(all_species.formula, all_species.symbol)]
candidate_primaries = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(df_primaries.formula, df_primaries.symbol)]
A, indep_comp, dep_comp = stoich_matrix(species, candidate_primaries) ;

# Construction of stoich matrix with aqueous species from database
aqueous_species = filter(row->row.aggregate_state == "AS_AQUEOUS", df_substances)
species = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(aqueous_species.formula, aqueous_species.symbol)]
candidate_primaries = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(df_primaries.formula, df_primaries.symbol)]
A, indep_comp, dep_comp = stoich_matrix(species, candidate_primaries) ;
stoich_matrix_to_reactions(A, indep_comp, dep_comp) ;

# CemSpecies with Sym coef
â, b̂, ĝ = symbols("â b̂ ĝ", real=true)
ox = Dict(:C => â, :S => one(Sym), :A => b̂, :H => ĝ)
CSH = CemSpecies(ox; aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
numCSH = apply(N, apply(subs, CSH, â=>1.8, b̂=>1, ĝ=>5))
floatCSH = apply(x->convert(Float64, x), numCSH) # only coefficients of oxides are converted to Float64 here not those of atoms

# Conversion to cement notation with species from database
CemSpecies(H₂O)
CemSpecies(vapour)
df_Jennite = filter(row->row.symbol == "Jennite", df_substances)
Jennite = Species(df_Jennite.formula[1]; name=df_Jennite.name[1], symbol=df_Jennite.symbol[1], aggregate_state=eval(Meta.parse(df_Jennite.aggregate_state[1])), class=eval(Meta.parse(df_Jennite.class[1])))
cemJennite = CemSpecies(Jennite)
Jennite == cemJennite

# Extraction of properties
cemJennite.ΔfG = df_Jennite.ΔfG[1].values*uparse(df_Jennite.ΔfG[1].units)
cemJennite

# Equation parsing
equation = "13H⁺ + NO₃⁻ + CO₃²⁻ + 10e⁻ = 6H₂O@ + HCN@"
 ## simple parsing giving Dicts
reac, prod, equal_sign = parse_equation(equation)
 ## construction of a Reaction struct from a string
r = Reaction(equation)
 ## simplification of a Reaction
r = Reaction("2H₂O → H⁺ + OH⁻ +H₂O")
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
r = Reaction([C3S, H, CH, CSH]; equal_sign='→')
Reaction(CemSpecies.(["C3S", "H", "CH", "C1.8SH4"]))
for c_over_s in 1.5:0.1:2.
    println(Reaction(CemSpecies.(["C3S", "H", "CH", "C$(c_over_s)SH4"])))
end
 ## construction of a Reaction by a balance calculation with symbolic numbers
â, b̂, ĝ = symbols("â b̂ ĝ", real=true)
CSH = CemSpecies(Dict(:C => â, :S => one(â), :H => ĝ))
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
r = Reaction([CSH, C3S, H, CH]; equal_sign='→')
 ## application of a function to stoichimetric coefficients (here simplify)
r = apply(simplify, Reaction([C3S, H], [CH, CSH]; equal_sign='→'))
A, _, _ = stoich_matrix([C3S], [CSH, H, CH]; involve_all_atoms=true) ;
simplify.(A)

# Chen & Brouwers
CSH = CemSpecies(Dict(:C => â, :S => one(â), :A => b̂, :H => ĝ))
HT = CemSpecies("M₅AH₁₃")
HG = CemSpecies("C₆AFS₂H₈")
AFt = CemSpecies("C₆S̄₃H₃₂")
ST = CemSpecies("C₂ASH₈")
AH = CemSpecies("C₄AH₁₃")
A, ox = canonical_stoich_matrix([CSH, HT, HG, AFt, ST, AH]);
A = typeof(â).(A[1:end-1, 1:end]) # end-1 to remove the line corresponding to water H
oxides = (CemSpecies∘string).(ox[1:end-1])
hydrates = [CSH, HT, HG, AFt, ST, AH]
print_stoich_matrix(A, symbol.(oxides), symbol.(hydrates))
print_stoich_matrix(inv(A), symbol.(hydrates), symbol.(oxides))
Mhyd = ustrip.(getproperty.(hydrates, :molar_mass))
Mox = ustrip.(getproperty.(oxides, :molar_mass))
B = Mox .* A .* inv.(Mhyd)'
# or directly
B, ox = canonical_stoich_matrix([CSH, HT, HG, AFt, ST, AH]; mass=true) ;
B = B[1:end-1, 1:end] # to remove the H line
print_stoich_matrix(B, "m_" .* symbol.(oxides), "m_" .* symbol.(hydrates))
print_stoich_matrix(subs.(inv(B), Ref(Dict(â=>1.8, b̂=>1, ĝ=>4))), "m_" .* symbol.(hydrates), "m_" .* symbol.(oxides))

# Alkane combustion with SymPy
n = symbols("n", real=true)
CₙH₂ₙ₊₂ = Species(:C => n, :H => 2n+2)
O₂ = Species("O₂")
H₂O = Species("H₂O")
CO₂ = Species("CO₂")
r = Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂])
apply(factor, r)
println(2r)
for vn in 1:9 println("n=$vn ⇒ ", apply(subs, r, n=>vn)) end
println(Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂]; side=:products))
println(Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂]; side=:reactants))
@show r[O₂]

# Alkane combustion with Symbolics
@variables n
CₙH₂ₙ₊₂ = Species(:C => n, :H => 2n+2)
O₂ = Species("O₂")
H₂O = Species("H₂O")
CO₂ = Species("CO₂")
r = Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂])
for vn in 1:9 println("n=$vn ⇒ ", apply(substitute, r, n=>vn)) end
for vn in 1:9 println("n=$vn ⇒ ", apply(stoich_coef_round∘value∘substitute, r, n=>vn)) end

# Example from https://github.com/thermohub/chemicalfun
formulas = ["Ca+2", "Fe+2", "Fe|3|+3", "H+", "OH-", "SO4-2", "CaSO4@", "CaOH+", "FeO@", "HFe|3|O2@", "FeOH+", "Fe|3|OH+2", "H2O@",  "FeS|-2|", "FeS|0|S|-2|", "S|4|O2"] ;
species = Species.(formulas) ;
candidate_primaries = species[1:6] ;
A, indep_comp, dep_comp = stoich_matrix(species) ;
B, indep_comp, dep_comp = stoich_matrix(species; mass=true) ;
lr = stoich_matrix_to_reactions(A, indep_comp, dep_comp) ;

# Callable
 # with units (coefficient units should be consistent with the basis of functions provided in thermofun database)
cemJennite.Cp = ThermoFunction(:Cp, [210.0J/K/mol, 0.120J/mol/K^2, -3.07e6J*K/mol, 0.0J/mol/√K])
@show cemJennite.Cp ;
@show fieldnames(typeof(cemJennite.Cp)) ;
@show fieldnames(typeof(cemJennite.Cp)) ;
@show cemJennite.Cp(298.15K) ;
@show cemJennite.Cp() ; # application by default on Tref
 # same without units
cemJennite.Cp = ThermoFunction(:Cp, [210.0, 0.120, -3.07e6, 0.0])
@show cemJennite.Cp ;
@show fieldnames(typeof(cemJennite.Cp)) ;
@show fieldnames(typeof(cemJennite.Cp)) ;
@show cemJennite.Cp(298.15) ;
@show cemJennite.Cp() ; # application by default on Tref

using Plots
lT = ((0:1:100) .+ 273.15).*K
@time plot(lT, dict_species["Jennite"].Cp.(lT))

# Check which species involved in reactions have not been previously constructed in the list of substances (in this case they are built on-the-fly and don't have thermo properties)
for row in eachrow(df_reactions)
    re = row.reaction
    for s in keys(re.products)
        if !haskey(s, :Cp) println(s) end
    end
    for s in keys(re.reactants)
        if !haskey(s, :Cp) println(s) end
    end
end

# Check consistency of logKr at Tref and logKr0 of the database
for r in df_reactions.reaction println(r.logKr(), " == ", r.logKr0) end