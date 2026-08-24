# Stoichiometric Matrix

Calculating stoichiometric matrices is a prerequisite for equilibrium calculations by minimizing Gibbs energy. The examples below show how they can be constructed from a thermodynamic database.

## Stoichiometric matrix for a subset of species

The recommended workflow is:
1. Load all species from the database with `build_species`.
2. Filter to the relevant chemical space with `speciation`.
3. Build a `ChemicalSystem`, which automatically computes the stoichiometric matrices.

```julia
using ChemistryLab

# 1. Load all species from the database
all_species = build_species(datapath("cemdata18-merged.json"))
```

Identify the secondary species likely to appear during the reactions of interest (here: C₃S, Portlandite, Jennite and water):

```julia
# 2. Filter species compatible with the chosen seeds
species = speciation(all_species, split("C3S Portlandite Jennite H2O@");
              aggregate_state=[AS_AQUEOUS], exclude_species=split("H2@ O2@"))
dict_species = Dict(symbol(s) => s for s in species)
```

Deduce the primary species present in this subset:

```julia
# 3. Select primaries available in the subset
candidate_primaries = [dict_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_species, s)]
```

Build the `ChemicalSystem`. Both stoichiometric matrices (`CSM` and `SM`) are computed once at construction time:

```@setup example1
    using ChemistryLab #hide
    all_species = build_species(datapath("cemdata18-merged.json")) #hide
    species = speciation(all_species, split("C3S Portlandite Jennite H2O@");
                  aggregate_state=[AS_AQUEOUS], exclude_species=split("H2@ O2@")) #hide
    dict_species = Dict(symbol(s) => s for s in species) #hide
    candidate_primaries = [dict_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_species, s)] #hide
```

```@example example1
cs = ChemicalSystem(species, candidate_primaries)
```

All independent reactions are then reconstructed from `cs.SM`:

```@example example1
lr = reactions(cs.SM)
```

---

## Stoichiometric matrix for all aqueous species in the database

To work with the full set of aqueous species, filter by aggregate state directly in `speciation`:

```julia
# Keep only aqueous species
aqueous_species = speciation(all_species, collect(keys(first(all_species).atoms));
                      aggregate_state=[AS_AQUEOUS])
```

A more direct alternative is to filter using all atom symbols present in the database:

```julia
all_atoms = union_atoms(all_species)
aqueous_species = speciation(all_species, all_atoms; aggregate_state=[AS_AQUEOUS])
dict_aqueous_species = Dict(symbol(s) => s for s in aqueous_species)
candidate_primaries = [dict_aqueous_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_aqueous_species, s)]
```

```@setup example1
    aqueous_species = speciation(all_species, union_atoms(all_species); aggregate_state=[AS_AQUEOUS]) #hide
    dict_aqueous_species = Dict(symbol(s) => s for s in aqueous_species) #hide
    candidate_primaries_aq = [dict_aqueous_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_aqueous_species, s)] #hide
```

```@example example1
cs_aq = ChemicalSystem(aqueous_species, candidate_primaries_aq)
lr_aq = reactions(cs_aq.SM)
```

Here only ionic species appear, given the `AS_AQUEOUS` filter.

---

The exercise can also be done on solid species (`AS_CRYSTAL`) or gases (`AS_GAS`). Primaries are still drawn from the aqueous subset in those cases:

```julia
# Solid species
solid_species = speciation(all_species, union_atoms(all_species); aggregate_state=[AS_CRYSTAL])
cs_solid = ChemicalSystem(solid_species, candidate_primaries)
lr_solid = reactions(cs_solid.SM)

# Gas species
gas_species = speciation(all_species, union_atoms(all_species); aggregate_state=[AS_GAS])
cs_gas = ChemicalSystem(gas_species, candidate_primaries)
lr_gas = reactions(cs_gas.SM)
```
