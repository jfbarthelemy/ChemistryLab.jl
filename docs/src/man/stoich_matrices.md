# Stoichiometric Matrix

## Concept

A stoichiometric matrix encodes how each chemical species is composed of elements (or oxide components). Each row corresponds to a species; each column to an element. The entry `A[i, j]` is the number of atoms of element `j` in species `i`.

For example, three species:

| Species | H | O | C |
| ------- | - | - | - |
| H₂O    | 2 | 1 | 0 |
| H⁺     | 1 | 0 | 0 |
| CO₂    | 0 | 2 | 1 |

This 3×3 matrix has rank 3, so all three species are linearly independent and can each serve as a "primary" (independent component). Any additional species in the system (e.g. HCO₃⁻, CO₃²⁻) can be written as a linear combination of these primaries — which is precisely a balanced chemical reaction.

```@example concept
using ChemistryLab

H2O = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
Hp  = Species("H+";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
CO2 = Species("CO2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)

CSM = CanonicalStoichMatrix([H2O, Hp, CO2])
pprint(CSM; label = :symbol)
```

Adding a dependent species:

```@example concept
HCO3 = Species("HCO3-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
CO3  = Species("CO3-2";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)

primaries = [H2O, Hp, CO2]
SM = StoichMatrix([H2O, Hp, CO2, HCO3, CO3], primaries)
pprint(SM; label = :symbol)
```

```@example concept
# Balanced reactions are extracted from the null space of SM
rxns = reactions(SM)
for r in rxns
    println(r.equation)
end
```

---

From the definition of species, it is possible to construct a stoichiometric matrix that establishes the relationship between species and chemical elements for species or oxides for cement species. This is called canonical decomposition.

```@setup database_stoichiometry
    using ChemistryLab
```

## Stoichiometric matrix for species

Any species can be described as a linear combination of chemical elements. A species vector can be expressed as a function of the chemical elements on which they depend. This dependence leads to the creation of a stoichiometric matrix.

```julia
using ChemistryLab
H2O = Species("H₂O", symbol="H₂O@", aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)
HSO4 = Species("HSO₄⁻", symbol="H₂O@", aggregate_state=AS_AQUEOUS, class=SC_COMPONENT)
CO2 = Species(Dict(:C => 1, :O => 2); symbol="CO₂", aggregate_state=AS_GAS, class=SC_GASFLUID)
species = [H2O, HSO4, CO2]
CSM = CanonicalStoichMatrix(species)
```

!!! note "Display of the stoichiometric matrix"
    The stoichiometric matrix can be pretty-printed with different column and row labels using `pprint`. Simply add the keyword `row_label`, `col_label` or `label` for both, which can take the following values: *:name*, *:symbol*, *:formula*
    ```julia
    SM = StoichMatrix(species)
    pprint(SM; label=:name)
    ```

```@example database_stoichiometry
H2O = Species("H₂O", symbol="H₂O@", aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT) #hide
HSO4 = Species("HSO₄⁻", symbol="H₂O@", aggregate_state=AS_AQUEOUS, class=SC_COMPONENT) #hide
CO2 = Species(Dict(:C => 1, :O => 2); symbol="CO₂", aggregate_state=AS_GAS, class=SC_GASFLUID) #hide
species = [H2O, HSO4, CO2] #hide
CSM = CanonicalStoichMatrix(species) #hide

pprint(CSM)
```

---

## Stoichiometric matrix for cement species

A cement species vector can also be expressed in terms of other species on which they depend. Here, the cement species are expressed in terms of the oxides from which they are composed.

```@example database_stoichiometry
C3S = CemSpecies("C3S", symbol="C₃S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C2S = CemSpecies("C2S", symbol="C₂S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C3A = CemSpecies("C3A", symbol="C₃A", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C4AF = CemSpecies(Dict(:C=>4, :A=>1, :F=>1); name="C4AF", symbol="C₄AF", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
cemspecies = [C3S, C2S, C3A, C4AF]
CSM = CanonicalStoichMatrix(cemspecies)

pprint(CSM)
```

---

## Define stoichiometric matrix from primary species

The decomposition of a set of species can also be done according to a base of primary species.

```@example database_stoichiometry
H2O = Species("H₂O")
HSO4 = Species("HSO₄⁻")
CO2 = Species(Dict(:C => 1, :O => 2); symbol="CO₂")
species = [H2O, HSO4, CO2]
H⁺ = Species("H⁺")
SO₄²⁻ = Species("SO₄²⁻")
CO₃²⁻ = Species("CO₃²⁻")
primary_species = [H⁺, SO₄²⁻, CO₃²⁻, H2O]
SM = StoichMatrix(species, primary_species)

pprint(SM)
```

---

## Construct stoichiometric matrix from a database

The recommended workflow loads all species from a JSON file with `build_species`, filters them with `speciation`, then builds a `ChemicalSystem` which internally computes both the canonical stoichiometric matrix (`CSM`) and the primary-based stoichiometric matrix (`SM`).

```julia
using ChemistryLab

# 1. Load all species from the database
all_species = build_species(datapath("cemdata18-merged.json"))
```

For a given set of seed species (e.g. Portlandite and water), `speciation` selects all species whose atomic composition is a subset of the seed atoms:

```julia
# 2. Filter to the relevant chemical space
species = speciation(all_species, split("Portlandite H2O@");
              aggregate_state=[AS_AQUEOUS], exclude_species=split("H2@ O2@"))
```

Primary species candidates for the Cemdata18 database are available via `CEMDATA_PRIMARIES`:

```julia
# 3. Select primaries present in our species subset
dict_species = Dict(symbol(s) => s for s in species)
candidate_primaries = [dict_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_species, s)]
```

A `ChemicalSystem` then computes both stoichiometric matrices at construction time:

```@setup example1
    using ChemistryLab #hide
    all_species = build_species(datapath("cemdata18-merged.json")) #hide
    species = speciation(all_species, split("Portlandite H2O@");
                  aggregate_state=[AS_AQUEOUS], exclude_species=split("H2@ O2@ FeOH+ Fe+2")) #hide
    dict_species = Dict(symbol(s) => s for s in species) #hide
    candidate_primaries = [dict_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_species, s)] #hide
    cs = ChemicalSystem(species, candidate_primaries) #hide
```

```julia
cs = ChemicalSystem(species, candidate_primaries)
```

The stoichiometric matrices are then directly accessible as fields of the system:

```@example example1
cs.CSM   # canonical stoichiometric matrix (elements × species)
pprint(cs.CSM)
```

```@example example1
cs.SM    # stoichiometric matrix w.r.t. primaries
pprint(cs.SM)
```

And independent reactions are reconstructed from `cs.SM`:

```@example example1
list_reactions = reactions(cs.SM)
```

!!! note "ChemicalSystem vs StoichMatrix"
    `ChemicalSystem` is the recommended entry point for equilibrium workflows: it holds the species list, reaction list, index maps and both stoichiometric matrices in one immutable object. `StoichMatrix` and `CanonicalStoichMatrix` remain available for standalone use when equilibrium solving is not needed.

---
