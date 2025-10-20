using CementChemistry, Unicode, DataFrames

# 1. Test formula parsing with charge and Unicode
fSO4 = Formula("SO₄²⁻")
println("Parsed SO₄²⁻: ", fSO4)

# 1.bis Test convert Formula
println("As Float64: ", convert(Float64, fSO4))

# 2. Test Species creation from formula string and dictionary
species_SO4 = Species("SO₄²⁻")
species_SO4 = Species(Dict(:S=>1, :O=>4), -2)
species_dict = Species(Dict(:Na=>1, :O=>1))
println("Species from string: ", species_SO4)
println("Species from dict: ", species_dict)

# 3. Test CemSpecies with custom oxide composition
custom_oxide = CemSpecies(Dict(:C=>2, :S=>1, :A=>1, :F=>0.5); name="CustomPhase")
println("Custom CemSpecies: ", custom_oxide)
println("Stoichiometry: ", cemformula(custom_oxide))

# 4. Test conversion between Phreeqc and Unicode notation
phreeqc = "CaSO4"
unicode = phreeqc_to_unicode(phreeqc)
println("Phreeqc to Unicode: ", unicode)
println("Unicode to Phreeqc: ", unicode_to_phreeqc(unicode))

# 5. Test filtering database for solid phases only
df_elements, df_substances, df_reactions = read_thermofun("../data/cemdata18-merged.json")
solid_species = filter(row->row.symbol == "M15SH", df_substances)
solid_species = filter(row->row.aggregate_state == "AS_CRYSTAL", df_substances)
println("Number of solid species: ", nrow(solid_species))

# 6. Test symbolic CemSpecies with more variables
using SymPy
x̂, ŷ = symbols("x̂ ŷ", real=true)
sym_oxide = CemSpecies(Dict(:C=>x̂, :S=>ŷ, :A=>1, :H=>2))
println("Symbolic CemSpecies: ", sym_oxide)
num_oxide = CemSpecies(map(N, map(subs, cemformula(sym_oxide), x̂=>3, ŷ=>2)))
println("Numeric CemSpecies: ", num_oxide)