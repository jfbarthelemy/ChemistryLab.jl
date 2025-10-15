# Bogue Calculation

The way in which species and cementitious species are constructed in CementChemistry and expressed as a linear combination of reference species opens the door to equilibrium calculations. It also makes it quite natural to retrieve Bogue's formulas and use them simply.

Bogue's formulas allow us to find the masses of C3S, C2S, C3A and C4AF as a function of the oxides (CaO, SiO2, Al2O3 and Fe2O3) that are regularly found in manufacturers' cement data sheets. However, using the `stoich_matrix` functions performs a molar decomposition of the species that we wish to decompose as a function of reference species. It is therefore possible to express the anhydrous of the cement as a function of the oxides in a cement data sheet.

```@example Bogue
using CementChemistry #hide
C3S = CemSpecies("C3S")
C2S = CemSpecies("C2S")
C3A = CemSpecies("C3A")
C4AF = CemSpecies("C4AF")
cemspecies = [C3S, C2S, C3A, C4AF]
CaO = CemSpecies("C")
SiO2 = CemSpecies("S")
Al2O3 = CemSpecies("A")
Fe2O3 = CemSpecies("F")
oxides = [CaO, SiO2, Al2O3, Fe2O3]
A, indep_comp = stoich_matrix(cemspecies,oxides)

using PrettyTables #hide
```
Bogue's formulas are thus found by converting species into mass and inverting matrix A.

```@example Bogue
# Molar mass of anhdrous phases
Mw = map(x -> x.molar_mass, cemspecies)
# Molar mass of each oxide
Mwo = map(x -> x.molar_mass, oxides)
Aoa = Mwo .* A .* inv.(Mw)'
pretty_table(
    inv(Aoa),
    row_labels = map(x -> x.name, cemspecies),
    column_labels = map(x -> x.name, oxides),
    style         = TextTableStyle(;
                        row_label = crayon"magenta bold",
                        first_line_column_label = crayon"cyan bold",
                        table_border = crayon"green bold")
)
```
By taking a cement sheet with a classic percentages of oxides (CaO=65.6%; SiO2=21.5%; Al2O3=5.2% and Fe2O3=2.8%), we then obtain the anhydrous masses of the cementitious material. 
```@example Bogue
inv(Aoa) *  [65.6, 21.5, 5.2, 2.8]
```

Another, faster way to do this is to take advantage of the `mass=true` option of the `canonical_stoich_matrix` and `stoich_matrix` functions. This option allows species to be expressed as a linear combination of the reference species in mass rather than in moles.

```@example Bogue
A, indep_comp = canonical_stoich_matrix(cemspecies; mass=true)
```

Bogue's formulas are then immediate.
```@example Bogue
pretty_table(
    inv(A),
    row_labels = map(x -> x.name, cemspecies),
    column_labels = map(x -> x.name, oxides),
    style         = TextTableStyle(;
                        row_label = crayon"magenta bold",
                        first_line_column_label = crayon"cyan bold",
                        table_border = crayon"green bold")
)
```