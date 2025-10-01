This tutorial progressively introduces the main features of [CementChemistry](@ref) using practical examples.

---

# Quickstart

Install CementChemistry in your chosen environment by entering pkg mode by pressing `]` and then:

```julia
pkg> add CementChermistry
```

In order to use CementChemistry, it is then necessary to load the CementChemistry.jl package:

```julia
julia> using CementChemistry
```

---

# Chemical Formula Manipulation

CementChemistry allows you to create and manipulate chemical formulas. It is based on `Formula` which is a structure (`struct`) which contains an expression, a writing of the formula close to those found in the Phreeqc databases, a unicode expression as well as a composition in the form of dictionaries and a charge.

```julia
struct Formula{T<:Number}
    expr::String
    phreeqc::String
    unicode::String
    composition::Dict{Symbol,T}
    charge::Int8
end
```
 Formulas can be constructed:
- by parsing a string containing eventually fractional or decimal coefficients
```@example 1
using CementChemistry #hide
fgen = Formula("A1//2B3C0.4")
```
- from symbols representing atoms 
```@example
using CementChemistry #hide
fCO2 = :C + 2 * :O
```

Coefficient types can be changed *a posteriori*, the type of the `Formula` `struct` being associated with the most complex type of the set of coefficients.

```@example 1
convert(Float64, fgen)
```

---

# Species

`Species` is also a `struct` and is defined by a name, symbol, structure and properties. It creates chemical species for solution or solid phases:

```julia
struct Species{T<:Number} <: AbstractSpecies
    name::String
    symbol::String
    formula::Formula{T}
    properties::Dict{Symbol,Number}
end
```

`Species` can be created from:
- a `Formula`

```@example
using CementChemistry #hide
fH2O = 2 * :H + :O
H2O = Species(fH2O)
```

- a string
```@example
using CementChemistry #hide
HSO4 = Species("HSO₄⁻")
```

- a dictionary
```@example
using CementChemistry #hide
CO2 = Species(Dict(:C => 1, :O => 2); name="CO₂")
```

> **_NOTE:_**  To add a charge when creating species with a dictionary, you must add, after the dictionary, the value of the charge (charge is considered an argument of the structure).
```@example
using CementChemistry #hide
CO2 = Species(Dict(:Si => 1, :O => 3),-2; name="SiO₃²⁻")
```

> **_NOTE:_** You will also have noticed that a calculation of the molar mass of the species is systematically carried out.

---

# Cement Species

The manipulation of chemical formulas can also be done in cement notation. Here are examples of anhydrous phases:

```@setup example_cemspecies
    using CementChemistry
```

```@example example_cemspecies
C3S = CemSpecies("C3S")
C2S = CemSpecies("C2S")
C3A = CemSpecies("C3A")
C4AF = CemSpecies(Dict(:C => 4, :A => 1, :F => 1); name = "C4AF")
```
---

# Symbolic and Numeric CemSpecies

The previous species were constructed from integer values ​​of the number of chemical elements. However, numerical values ​​are possible, as we have seen for formulas, as well as symbolic values. To do this, you need to use the [`SymPy`](https://github.com/JuliaPy/SymPy.jl) library:

```julia
using SymPy
â, b̂, ĝ = symbols("â b̂ ĝ", real = true)
ox = Dict(:C => â, :S => one(Sym), :A => b̂, :H => ĝ)
CSH = CemSpecies(ox)
```

```julia
numCSH = CemSpecies(map(N, map(subs, cemformula(CSH), â => 1.8, b̂ => 1, ĝ => 5)))
```

Conversion of coefficient types can also be done:

```julia
floatCSH = Species(convert(Float64, formula(numCSH)))
```

---

# Conversion to Cement Notation

Convert species to cement notation and Unicode:

```@example example_cemspecies
H2O = Species("H₂O")
cemH2O = CemSpecies(H2O)
```

---


# Database Interoperability

So far, we have looked at the possibility of creating and manipulating any species, whether they exist or not. Creating an H₂O⁺⁴ molecule, for example, is not a problem.

```@setup database_interoperability
    using CementChemistry
```

```@example database_interoperability
HSO4 = Species("H₂O⁺⁴")
```

However, you will admit that it is a little strange...

## Databases

This is why Cement Chemistry relies on existing databases, in particular [Cemdata18](https://www.empa.ch/web/s308/thermodynamic-data) and [PSI-Nagra-12-07](https://www.psi.ch/en/les/thermodynamic-databases). Cemdata18 is a chemical thermodynamic database for hydrated Portland cements and alkali-activated materials. PSI-Nagra is a Chemical Thermodynamic Database. The formalism adopted for these databases is that of [Thermofun](https://thermohub.org/thermofun/thermofun/) which is a universal open-source client that delivers thermodynamic properties of substances and reactions at the temperature and pressure of interest. The information is stored in json files.

With Cementchemistry, you can parse a ThermoFun-like json file and return DataFrames for elements, substances, and reactions.

[`CementChemistry.parse_cemdata18_thermofun`](@ref)
```@example database_interoperability
df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("../../data/cemdata18-merged.json")
show(df_elements, allcols=true, allrows=true)
```

```@example database_interoperability
show(df_substances, allcols=true, allrows=false)
```

```@example database_interoperability
show(df_reactions, allcols=true, allrows=true)
```

It is also possible to retrieve primary species from the Cemdata18 database, primary species being the designation of a subset of species for which any species can be represented as the linear combination of primary species.

[`CementChemistry.extract_primary_species`](@ref)
```@example database_interoperability
df_primaries = extract_primary_species("../../data/CEMDATA18-31-03-2022-phaseVol.dat")
show(df_primaries, allcols=true, allrows=true)
```

---

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

# Advanced Stoichiometric Matrix (Database Species)

The same exercise can be performed between species and primary species defined from the species and cemdata18 database. :

```@example
using CementChemistry #hide
df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("../../data/cemdata18-merged.json") #hide
df_primaries = extract_primary_species("../../data/CEMDATA18-31-03-2022-phaseVol.dat") #hide
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

# Reaction Parsing and Manipulation

Parse and format chemical equations:

```julia
equation = "13H⁺ + NO₃⁻ + CO₃²⁻ + 10e⁻ = 6H₂O@ + HCN@"
reactants, equal_sign, products = parse_equation(equation)
r = Reaction(equation)
format_equation(Dict(symbol(k) => v for (k, v) in r.species_stoich))
```

Create reactions with :

```julia
eqC3S = "C₃S + 5.3H = 1.3CH + C₁.₇SH4"
rC3S = CemReaction(eqC3S)
format_equation(Dict(symbol(k) => v for (k, v) in rC3S.species_stoich))
```

---

# Reaction Operations

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

# Automated Reaction Balancing

Balance reactions automatically:

```julia
r = Reaction(CSH, [H, CH, C3S]; equal_sign = "→")
```

---

This tutorial covers the main workflow and features of CementChemistry.jl. For more details, see the [documentation](https://jfbarthelemy.github.io/CementChemistry.jl/dev/) and .