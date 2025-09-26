# Tutorials



# CementChemistry.jl Step-by-Step Tutorial

This tutorial progressively introduces the main features of [CementChemistry.jl](src/CementChemistry.jl) using practical examples.

---

## 1. Setup

First, load the package and its dependencies:

```julia
using CementChemistry
using Unicode
```

---

## 2. Chemical Formula Manipulation

Create and manipulate chemical formulas using [`Formula`](src/formulas.jl):

```@example
using CementChemistry #hide
using Unicode #hide
fgen = Formula("A1//2B3C0.4")  # Parses a formula with fractional and decimal coefficients
fCO2 = :C + 2 * :O             # Build formula using AtomGroup operations
```

Convert formula coefficients to Float64:

```@example
using CementChemistry #hide
using Unicode #hide
convert(Float64, fgen)
```

---

## 3. Species Creation

Create chemical species for solution or solid phases using [`Species`](src/species.jl):

```julia
fH2O = 2 * :H + :O
H2O = Species(fH2O)
HSO4 = Species("HSO₄⁻")
CO2 = Species(Dict(:C => 1, :O => 2); name="CO₂")
```

---

## 4. Stoichiometric Matrix Construction

Build stoichiometric matrices for a set of species:

```julia
species = [H2O, HSO4, CO2]
A, indep_comp, dep_comp = stoich_matrix(species)
```

---

## 5. Cement Notation and Solid Phases

Work with cement notation and solid phases using [`CemSpecies`](src/species.jl):

```julia
C3S = CemSpecies("C3S")
C2S = CemSpecies("C2S")
C3A = CemSpecies("C3A")
C4AF = CemSpecies(Dict(:C => 4, :A => 1, :F => 1); name = "C4AF")
cemspecies = [C3S, C2S, C3A, C4AF]
A, indep_comp, dep_comp = stoich_matrix(cemspecies)
```

---

## 6. Database Interoperability
Imagine we want to extract the set of aqueous species from a ThermoFun JSON data structure. For that, we use the function
[`CementChemistry.get_aqueous_species`](@ref):

```@example
using CementChemistry

```

Import data from ThermoFun and Cemdata databases:

```julia
df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("data/cemdata18-merged.json")
df_primaries = extract_primary_species("data/CEMDATA18-31-03-2022-phaseVol.dat")
```


---

## 7. Advanced Stoichiometric Matrix (Database Species)

Filter and combine species from the database:

```julia
given_species = filter(row -> row.symbol ∈ split("C3S Portlandite Jennite H2O@"), df_substances)
secondaries = filter(row -> row.aggregate_state == "AS_AQUEOUS" &&
                          all(k -> first(k) ∈ union_atoms(given_species.atoms), row.atoms) &&
                          row.symbol ∉ split("H2@ O2@"),
                          df_substances)
all_species = unique(vcat(given_species, secondaries), :symbol)
species = [Species(f; name = phreeqc_to_unicode(n)) for (f, n) in zip(all_species.formula, all_species.symbol)]
candidate_primaries = [Species(f; name = phreeqc_to_unicode(n)) for (f, n) in zip(df_primaries.formula, df_primaries.symbol)]
A, indep_comp, dep_comp = stoich_matrix(species, candidate_primaries)
```

---

## 8. Symbolic and Numeric CemSpecies

Use symbolic coefficients with [`SymPy`](https://github.com/JuliaPy/SymPy.jl):

```julia
using SymPy
â, b̂, ĝ = symbols("â b̂ ĝ", real = true)
ox = Dict(:C => â, :S => one(Sym), :A => b̂, :H => ĝ)
CSH = CemSpecies(ox)
numCSH = CemSpecies(map(N, map(subs, cemformula(CSH), â => 1.8, b̂ => 1, ĝ => 5)))
floatCSH = Species(convert(Float64, formula(numCSH)))
```

---

## 9. Conversion to Cement Notation

Convert species to cement notation and Unicode:

```julia
H2O = Species("H₂O")
cemH2O = CemSpecies(H2O)
```

---

## 10. Reaction Parsing and Manipulation

Parse and format chemical equations:

```julia
equation = "13H⁺ + NO₃⁻ + CO₃²⁻ + 10e⁻ = 6H₂O@ + HCN@"
reactants, equal_sign, products = parse_equation(equation)
r = Reaction(equation)
format_equation(Dict(symbol(k) => v for (k, v) in r.species_stoich))
```

Create reactions with [`CemReaction`](src/reactions.jl):

```julia
eqC3S = "C₃S + 5.3H = 1.3CH + C₁.₇SH4"
rC3S = CemReaction(eqC3S)
format_equation(Dict(symbol(k) => v for (k, v) in rC3S.species_stoich))
```

---

## 11. Reaction Operations

Build reactions by combining species:

```julia
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
CSH = CemSpecies("C1.7SH4")
r = C3S + 5.3H ↔ 1.3CH + CSH
species_list(r)
stoich_list(r)
```

---

## 12. Automated Reaction Balancing

Balance reactions automatically:

```julia
r = Reaction(CSH, [H, CH, C3S]; equal_sign = "→")
```

---

This tutorial covers the main workflow and features of CementChemistry.jl. For more details, see the [documentation](https://jfbarthelemy.github.io/CementChemistry.jl/dev/) and [reference](docs/src/reference.md).