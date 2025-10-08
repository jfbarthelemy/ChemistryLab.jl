# Examples

## Get Stoichiometric Matrix from a list of species

Let's imagine that we want to form the stochiometric matrix of a list of solid and water species. For that, we need to read the database from which these species originate and retrieve the list of primary species from that database.

```julia
using CementChemistry
df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("../../../data/cemdata18-merged.json")
df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat")
```
See [`CementChemistry.parse_cemdata18_thermofun`](@ref) and [`CementChemistry.extract_primary_species`](@ref)

It is then necessary to identify the list of secondary species likely to appear during the reactions.

```julia
given_species = filter(row -> row.symbol ∈ split("C3S Portlandite Jennite H2O@"), df_substances)
secondaries = filter(row -> row.aggregate_state == "AS_AQUEOUS" &&
                          all(k -> first(k) ∈ union_atoms(given_species.atoms), row.atoms) &&
                          row.symbol ∉ split("H2@ O2@"),
                          df_substances)
```

We can then deduce the primary species concerned by the reaction.

```julia
all_species = unique(vcat(given_species, secondaries), :symbol)
species = [Species(f; symbol = phreeqc_to_unicode(n)) for (f, n) in zip(all_species.formula, all_species.symbol)]
candidate_primaries = [Species(f; symbol = phreeqc_to_unicode(n)) for (f, n) in zip(df_primaries.formula, df_primaries.symbol)]
```

And construct the stoichiometric matrix

```@setup example1
    using CementChemistry #hide
    df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("../../../data/cemdata18-merged.json") #hide
    df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat") #hide

    given_species = filter(row -> row.symbol ∈ split("C3S Portlandite Jennite H2O@"), df_substances) #hide
    secondaries = filter(row -> row.aggregate_state == "AS_AQUEOUS" &&
                            all(k -> first(k) ∈ union_atoms(given_species.atoms), row.atoms) &&
                            row.symbol ∉ split("H2@ O2@"),
                            df_substances) #hide


    all_species = unique(vcat(given_species, secondaries), :symbol) #hide
    species = [Species(f; symbol = phreeqc_to_unicode(n)) for (f, n) in zip(all_species.formula, all_species.symbol)] #hide
    candidate_primaries = [Species(f; symbol = phreeqc_to_unicode(n)) for (f, n) in zip(df_primaries.formula, df_primaries.symbol)] #hide
```

```@example example1
A, indep_comp, dep_comp = stoich_matrix(species, candidate_primaries)

using PrettyTables #hide
```

## Get Stoichiometric Matrix from a database file

```@example example2
using CementChemistry
df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("../../../data/cemdata18-merged.json")
df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat")
aqueous_species = filter(row->row.aggregate_state == "AS_AQUEOUS", df_substances)
species = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(aqueous_species.formula, aqueous_species.symbol)]
candidate_primaries = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(df_primaries.formula, df_primaries.symbol)]
A, indep_comp, dep_comp = stoich_matrix(species, candidate_primaries) ;

using PrettyTables #hide
```

All the reactions of the species contained in the database can thus be reconstructed. Here, only ionic species are listed given the choice to only read ionic species in the database ("AS_AQUEOUS").

```@example example2
stoich_matrix_to_reactions(A, indep_comp, dep_comp) ;
```

---

The exercise can also be done on solid species. In this case, the data filter is carried out using the keyword "AS_CRYSTAL", in accordance with the terminology adopted in Thermofun.

```@setup example3
using CementChemistry
df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("../../../data/cemdata18-merged.json")
df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat")
aqueous_species = filter(row->row.aggregate_state == "AS_CRYSTAL", df_substances)
species = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(aqueous_species.formula, aqueous_species.symbol)]
candidate_primaries = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(df_primaries.formula, df_primaries.symbol)]
```

```@example example3
A, indep_comp, dep_comp = stoich_matrix(species, candidate_primaries) ; #hide

using PrettyTables #hide
```

```@example example3
stoich_matrix_to_reactions(A, indep_comp, dep_comp) ; #hide
```

---

Or with gases ("AS_GAS")

```@setup example4
using CementChemistry
df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("../../../data/cemdata18-merged.json")
df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat")
aqueous_species = filter(row->row.aggregate_state == "AS_GAS", df_substances)
species = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(aqueous_species.formula, aqueous_species.symbol)]
candidate_primaries = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(df_primaries.formula, df_primaries.symbol)]
```

```@example example4
A, indep_comp, dep_comp = stoich_matrix(species, candidate_primaries) ; #hide

using PrettyTables #hide
```

```@example example4
stoich_matrix_to_reactions(A, indep_comp, dep_comp) ; #hide
```