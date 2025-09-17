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
CO₂ = Species(Dict(:C=>1, :O=>2))
species = [H₂O, HSO₄⁻, CO₂] ;
stoich_matrix(species) ;

# CemSpecies
OXIDE_ORDER # provides the order of oxides in cement formulas
C3S = CemSpecies("C3S")
C2S = CemSpecies("C2S")
C3A = CemSpecies("C3A")
C4AF = CemSpecies(Dict(:C=>4, :A=>1, :F=>1))
cemspecies = [C3S, C2S, C3A, C4AF]
stoich_matrix(cemspecies) ;

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
species = [Species(f; name=phreeqc_to_unicode(n)) for (f,n) in zip(all_species.formula, all_species.symbol)]
candidate_primaries = [Species(f; name=phreeqc_to_unicode(n)) for (f,n) in zip(df_primaries.formula, df_primaries.symbol)]
stoich_matrix(species, candidate_primaries) ;

# Construction of stoich matrix with aqueous species from database
aqueous_species = filter(row->row.aggregate_state == "AS_AQUEOUS", df_substances)
species = [Species(f; name=phreeqc_to_unicode(n)) for (f,n) in zip(aqueous_species.formula, aqueous_species.symbol)]
candidate_primaries = [Species(f; name=phreeqc_to_unicode(n)) for (f,n) in zip(df_primaries.formula, df_primaries.symbol)]
A, indep_comp, dep_comp = stoich_matrix(species, candidate_primaries) ;

# CemSpecies with Sym coef
using SymPy
â, b̂, ĝ = symbols("â b̂ ĝ", real=true)
ox = Dict(:C => â, :S => one(Sym), :A => b̂, :H => ĝ)
CSH = CemSpecies(ox)
numCSH = CemSpecies(map(N, map(subs, cemformula(CSH), â=>1.8, b̂=>1, ĝ=>5)))
floatCSH = Species(convert(Float64, formula(numCSH)))
