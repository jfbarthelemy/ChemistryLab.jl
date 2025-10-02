using Revise, CementChemistry, Unicode
import Unitful: u, g, cm, K, J, mol, Quantity, uconvert
using PeriodicTable
using PhysicalConstants.CODATA2022

# Formula
fgen = Formula("A1//2B3C0.4")
convert(Float64, fgen)
ATOMIC_ORDER # provides the order of atoms in formulas
fCO2 = :C+2*:O # works only with a part of real atoms registered in ATOMIC_ORDER
fHSO₄⁻ = :H+:S+4*:O+:e
fNa⁺ = :Na+:Zz

# Species
fH₂O = 2*:H + :O
H₂O = Species(fH₂O)
HSO₄⁻ = Species("HSO₄⁻")
CO₂ = Species(Dict(:C=>1, :O=>2); name="CO₂")
species = [H₂O, HSO₄⁻, CO₂] ;
stoich_matrix(species) ;

# CemSpecies
OXIDE_ORDER # provides the order of oxides in cement formulas
C3S = CemSpecies("C3S")
C2S = CemSpecies("C2S")
C3A = CemSpecies("C3A")
C4AF = CemSpecies(Dict(:C=>4, :A=>1, :F=>1); name="C4AF")
cemspecies = [C3S, C2S, C3A, C4AF]
A, indep_comp = stoich_matrix(cemspecies) ;

# Thermofun cemdata18
df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("data/cemdata18-merged.json")
df_primaries = extract_primary_species("data/CEMDATA18-31-03-2022-phaseVol.dat")

# Construction of stoich matrix with species from database
given_species = filter(row->row.symbol ∈ split("C3S Portlandite Jennite H2O@"), df_substances)
secondaries = filter(row->row.aggregate_state == "AS_AQUEOUS" 
                          && all(k->first(k) ∈ union_atoms(given_species.atoms), row.atoms)
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
using SymPy
â, b̂, ĝ = symbols("â b̂ ĝ", real=true)
ox = Dict(:C => â, :S => one(Sym), :A => b̂, :H => ĝ)
CSH = CemSpecies(ox)
numCSH = CemSpecies(map(N, map(subs, cemformula(CSH), â=>1.8, b̂=>1, ĝ=>5)))
floatCSH = Species(convert(Float64, formula(numCSH)))

# Conversion to cement notation
H₂O = Species("H₂O")
CemSpecies(H₂O)
df_Jennite = filter(row->row.symbol == "Jennite", df_substances)
Jennite = Species(df_Jennite.formula[1]; name=df_Jennite.name[1], symbol=df_Jennite.symbol[1])
cemJennite = CemSpecies(Jennite)
println(unicode(Jennite), " ≡ ", unicode(cemJennite))

# Equation parsing
equation = "13H⁺ + NO₃⁻ + CO₃²⁻ + 10e⁻ = 6H₂O@ + HCN@"
 ## simple parsing giving Dicts
reac, prod, equal_sign = parse_equation(equation)
 ## construction of a Reaction struct from a string
r = Reaction(equation)
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
 ## application of a function to stoichimetric coefficients (here simplify)
r = map(simplify, Reaction([C3S, H], [CH, CSH]; equal_sign='→'))
Reaction(CemSpecies.(["C3S", "H", "CH", "C1.8SH4"]))
for c_over_s in 1.5:0.1:2.
    println(Reaction(CemSpecies.(["C3S", "H", "CH", "C$(c_over_s)SH4"])))
end
 ## construction of a Reaction by a balance calculation with symbolic numbers
â, b̂, ĝ = symbols("â b̂ ĝ", real=true)
CSH = CemSpecies(Dict(:C => â, :S => one(Sym), :H => ĝ))
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
r = Reaction([CSH, C3S, H, CH]; equal_sign='→')
stoich_matrix([C3S], [CSH, H, CH]; involve_all_atoms=true) ;

# Chen & Brouwers
CSH = CemSpecies(Dict(:C => â, :S => 1, :A => b̂, :H => ĝ))
HT = CemSpecies("M₅AH₁₃")
HG = CemSpecies("C₆AFS₂H₈")
AFt = CemSpecies("C₆S̄₃H₃₂")
ST = CemSpecies("C₂ASH₈")
AH = CemSpecies("C₄AH₁₃")
A, ox = stoich_matrix([CSH, HT, HG, AFt, ST, AH]);
A = Sym.(A[1:end-1, 1:end])
oxides = (CemSpecies∘string).(ox[1:end-1])
hydrates = [CSH, HT, HG, AFt, ST, AH]
print_stoich_matrix(A, symbol.(oxides), symbol.(hydrates))
print_stoich_matrix(inv(A), symbol.(hydrates), symbol.(oxides))
Mhyd = getproperty.(hydrates, :molar_mass)
Mox = getproperty.(oxides, :molar_mass)
B = Mox .* A .* inv.(Mhyd)'
print_stoich_matrix(B, "m_" .* symbol.(oxides), "m_" .* symbol.(hydrates))
print_stoich_matrix(subs.(inv(B), â=>1.8, b̂=>1, ĝ=>4), "m_" .* symbol.(hydrates), "m_" .* symbol.(oxides))
