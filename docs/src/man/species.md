
# Species

`Species` is a `struct` and is defined by a name, symbol, structure and properties. It creates chemical species for solution or solid phases:

```julia
struct Species{T<:Number} <: AbstractSpecies
    name::String
    symbol::String
    formula::Formula{T}
    properties::OrderedDict{Symbol,Number}
end
```

## Species construction

`Species` can be created from:
- a `Formula`

```@example
using ChemistryLab #hide
fH2O = 2 * :H + :O
H2O = Species(fH2O)
```

- a string
```@example
using ChemistryLab #hide
HSO4 = Species("HSO₄⁻")
```

- a dictionary
```@example
using ChemistryLab #hide
CO2 = Species(Dict(:C => 1, :O => 2); name="CO₂")
```

!!! note "Add charge"
    To add a charge when creating species with a dictionary, you must add, after the dictionary, the value of the charge (charge is considered an argument of the structure).

```@example
using ChemistryLab #hide
CO2 = Species(Dict(:Si => 1, :O => 3),-2; name="SiO₃²⁻")
```

!!! tip "Remark"
    You will also have noticed that a calculation of the molar mass of the species is systematically carried out.

---

## Cement Species

The manipulation of chemical formulas can also be done in cement notation. Here are examples of anhydrous phases:

```@setup example_cemspecies
    using ChemistryLab
```

```@example example_cemspecies
C3S = CemSpecies("C3S")
C2S = CemSpecies("C2S")
C3A = CemSpecies("C3A")
C4AF = CemSpecies(Dict(:C => 4, :A => 1, :F => 1); name = "C4AF")
```
---

## Symbolic and Numeric CemSpecies

The previous species were constructed from integer values ​​of the number of chemical elements. However, numerical values ​​are possible (see [species](./databases.md#formulas)), as well as symbolic values. To do this, you need to use the [`SymPy`](https://github.com/JuliaPy/SymPy.jl) library:

```@example sympy1
using ChemistryLab
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

## Conversion to Cement Notation

Convert species to cement notation and Unicode:

```@example example_cemspecies
H2O = Species("H₂O")
cemH2O = CemSpecies(H2O)
```

Species can also be read from a database (Cemdata18 here). Reading the databases is detailed [here](./databases.md#Database-Interoperability). Conversion is then possible:

```@example example_cemspecies
df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("../../../data/cemdata18-merged.json")
df_Jennite = filter(row->row.symbol == "Jennite", df_substances)
Jennite = Species(df_Jennite.formula[1]; name=df_Jennite.name[1], symbol=df_Jennite.symbol[1])
cemJennite = CemSpecies(Jennite)
println(unicode(Jennite), " ≡ ", unicode(cemJennite))
```

---