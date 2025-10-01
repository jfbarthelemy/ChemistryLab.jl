

# Stoichiometric Matrix Construction

From the definition of species, it is possible to construct a stoichiometric matrix that establishes the relationship between species and chemical elements for species or cement species:

```@setup database_stoichiometry
    using CementChemistry
```

```julia
fH2O = 2 * :H + :O
H2O = Species(fH2O)
HSO4 = Species("HSO₄⁻")
CO2 = Species(Dict(:C => 1, :O => 2); name="CO₂")
species = [H2O, HSO4, CO2]
A, indep_comp, dep_comp = stoich_matrix(species)
println(A)
println(indep_comp)
println(dep_comp)
```

---

