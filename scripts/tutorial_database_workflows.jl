using ChemistryLab
using OrderedCollections

df_elements, df_substances, df_reactions = read_thermofun_database(datapath("cemdata18-merged.json"))
df_union = get_compatible_species(
    df_substances, split("C3S C2S C3A C4AF Portlandite Jennite H2O@");
    aggregate_states = [AS_AQUEOUS], exclude_species = split("H2@ O2@ FeOH+ Fe+2"), union = true
)
dict_all_species = Dict(symbol(s) => s for s in build_species(df_union))
candidate_primaries = [dict_all_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_all_species, s)]
SM = StoichMatrix(dict_all_species, candidate_primaries) ; pprint(SM)
list_reactions = reactions(SM) ; pprint(list_reactions)
for r in list_reactions
    display(r); println()
end

df_hydrates_clinker = get_compatible_species(df_substances, split("C3S C2S C3A C4AF Gp H2O@"); aggregate_states = [AS_CRYSTAL, AS_AQUEOUS])
dict_hydrates = Dict(symbol(s) => s for s in build_species(df_hydrates_clinker))
candidate_primaries = [dict_hydrates[s] for s in CEMDATA_PRIMARIES if haskey(dict_hydrates, s)]
@time SM = StoichMatrix(dict_hydrates, candidate_primaries) ; pprint(SM)
@time list_reactions = reactions(SM) ; pprint(list_reactions)

species = build_species(df_substances)
dict_species = Dict(symbol(s) => s for s in species)
dict_reactions = Dict(symbol(r) => r for r in build_reactions(df_reactions, species))

df_calcite = get_compatible_species(
    df_substances, split("Cal H2O@ CO2");
    aggregate_states = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"), union = true
)
dict_species_calcite = Dict(symbol(s) => s for s in build_species(df_calcite))
primaries = [dict_species_calcite[s] for s in split("H2O@ H+ CO2@ Ca+2")]
SM = StoichMatrix(values(dict_species_calcite), primaries); pprint(SM)
list_reactions = reactions(SM) ; pprint(list_reactions)
for r in list_reactions
    display(r); println()
end
dict_reactions_calcite = Dict(r.symbol => r for r in list_reactions)
