# Cement Species

The manipulation of chemical formulas can also be done in cement notation. `CemSpecies` is a composite type similar to `Species`.


```julia
struct Species{T<:Number} <: AbstractSpecies
    name::String
    symbol::String
    formula::Formula{T}
    aggregate_state::AggregateState
    class::Class
    properties::OrderedDict{Symbol,PropertyType}
end
```

## CemSpecies construction

Here are several example of `CemSpecies` construction:

```@setup example_cemspecies
    using ChemistryLab
```

```@example example_cemspecies
C3S = CemSpecies("C3S"; name="Alite", symbol="C₃S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C2S = CemSpecies("C₂S"; name="Belite", symbol="C₂S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C3A = CemSpecies("C3A"; name="Aluminate", symbol="C₃A", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C4AF = CemSpecies(Dict(:C => 4, :A => 1, :F => 1); name="Ferrite", symbol="C₄AF", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
```

!!! warning "Warning"
    Not every molecule can be used to build a cement species. It is necessary for this molecule to decompose into a combination of the oxides present in the manufacturers' cement sheet (e.g. $\ce{CaO}$, $\ce{SiO2}$, $\ce{Fe2O3}$, $\ce{Al2O3}$) and water. Thus, the following code will return an error.
    ```julia
    CemSpecies(Species("Ca(OH)"))
    ```


---

## Numeric and Symbolic CemSpecies

The previous species were constructed from integer values ​​of the number of chemical elements. However, other numerical value types ​​are possible (as for [formulas](./formula_manipulation.md#formulas)), such as fraction or Real values.

```@example example_cemspecies
using ChemistryLab
ox = Dict(:C => 1.666667, :S => 1, :H => 2.1)
jennite = CemSpecies(ox, aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
```

Symbolic values are also allowed. In this case, you need to use [`Symbolics`](https://github.com/JuliaSymbolics/Symbolics.jl) (already a dependency of ChemistryLab):

```julia
using ChemistryLab
using Symbolics
@variables a g
ox = Dict(:C => a, :S => one(Num), :H => g)
CSH = CemSpecies(ox)
```

The value of variables can be substituted *a posteriori* to obtain a numeric species:

```julia
jennite = CemSpecies(map(x -> substitute(x, Dict(a => 1.666667, g => 2.1)), cemformula(CSH)))
```

!!! note "Remark"
    Conversion of coefficient types can also be done.
    ```julia
    floatCSH = Species(convert(Float64, formula(numCSH)))
    ```

---

## Conversion of Species to CemSpecies and vice versa

Convert species to cement notation and Unicode. Conversion can be done on simple species:

```@example example_cemspecies
H2O = Species("H₂O", symbol="H₂O@", aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)
cemH2O = CemSpecies(H2O)
```

```@example example_cemspecies
C3S = CemSpecies("C3S"; name="Alite", symbol="C₃S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
spC3S = Species(C3S)
```

Or more complex one:

```@example example_cemspecies
using ChemistryLab #hide
CSH = Species("(SiO2)1(CaO)1.666667(H2O)2.1", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
jennite = CemSpecies(CSH)
```

---
