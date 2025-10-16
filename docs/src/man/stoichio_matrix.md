# Stoichiometric Matrix

From the definition of species, it is possible to construct a stoichiometric matrix that establishes the relationship between species and chemical elements for species or oxides for cement species. This is called canonical decomposition.

```@setup database_stoichiometry
    using ChemistryLab
    import Pkg; Pkg.add("PrettyTables")
```

## Stochiometric matrix for species

Any species can be described as a linear combination of chemical elements. A species vector can be expressed as a function of the chemical elements on which they depend. This dependence leads to the creation of a stochiometric matrix.


```@example
using ChemistryLab #hide
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
using ChemistryLab #hide
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
using ChemistryLab #hide
fH2O = 2 * :H + :O
H2O = Species(fH2O)
HSO4 = Species("HSO₄⁻")
CO2 = Species(Dict(:C => 1, :O => 2); symbol="CO₂")
species = [H2O, HSO4, CO2]
A, indep_comp, dep_comp = stoich_matrix(species)

using PrettyTables #hide
```

!!! note "Display of the stoichiometric matrix"
    The stochiometric matrix can be displayed with different column and row labels. Simply add the keyword 'label', which can take the following values: *:name*, *:symbol*, *:formula*

    ```julia
    A, indep_comp, dep_comp = stoich_matrix(species, label=:name)
    ```

The primary species of Cemdata18 can be listed with the following command:

```julia
using ChemistryLab #hide
df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("../../../data/cemdata18-merged.json") #hide
df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat")
```
See [`ChemistryLab.parse_cemdata18_thermofun`](@ref) and [`ChemistryLab.extract_primary_species`](@ref)



---