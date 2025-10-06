# Stoichiometric Matrix

From the definition of species, it is possible to construct a stoichiometric matrix that establishes the relationship between species and chemical elements for species or oxides for cement species. This is called canonical decomposition.

```@setup database_stoichiometry
    using CementChemistry
    import Pkg; Pkg.add("PrettyTables")
```

## Stochiometric matrix for species

Any species can be described as a linear combination of chemical elements. A species vector can be expressed as a function of the chemical elements on which they depend. This dependence leads to the creation of a stochiometric matrix.


```@example
using CementChemistry #hide
fH2O = 2 * :H + :O
H2O = Species(fH2O)
HSO4 = Species("HSO₄⁻")
CO2 = Species(Dict(:C => 1, :O => 2); symbol="CO₂")
species = [H2O, HSO4, CO2]
canonical_stoich_matrix(species) ;

using PrettyTables #hide
```

## Stochiometric matrix for cement species

A cement species vector can also be expressed in terms of other species on which they depend. Here, the cement species are expressed in terms of the oxides from which they are composed.

```@example
using CementChemistry #hide
C3S = CemSpecies("C3S")
C2S = CemSpecies("C2S")
C3A = CemSpecies("C3A")
C4AF = CemSpecies(Dict(:C=>4, :A=>1, :F=>1); name="C4AF")
cemspecies = [C3S, C2S, C3A, C4AF]
A, indep_comp = canonical_stoich_matrix(cemspecies)

using PrettyTables #hide
```

---

## Stochiometric matrix for species and primary species

The decomposition can also be done according to the primary species previously defined.

```@example
using CementChemistry #hide
fH2O = 2 * :H + :O
H2O = Species(fH2O)
HSO4 = Species("HSO₄⁻")
CO2 = Species(Dict(:C => 1, :O => 2); symbol="CO₂")
species = [H2O, HSO4, CO2]
A, indep_comp, dep_comp = stoich_matrix(species)

using PrettyTables #hide
```

The primary species of Cemdata18 can be listed with the following command:

```julia
using CementChemistry #hide
df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("../../../data/cemdata18-merged.json") #hide
df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat")
```

# Advanced Stoichiometric Matrix

Let's imagine that we now want to form the stochiometric matrix of a list of solid and water species.
As in the previous example, we need to read the database from which these species originate and retrieve the list of primary species from that database.

```julia
using CementChemistry
df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("../../../data/cemdata18-merged.json")
df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat")
```

It is then necessary to identify the list of secondary species likely to appear during the reaction.

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

```@example adv_example
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

A, indep_comp, dep_comp = stoich_matrix(species, candidate_primaries)

using PrettyTables #hide
```

---