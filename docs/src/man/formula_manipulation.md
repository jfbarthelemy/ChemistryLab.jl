# Chemical Formula Manipulation

ChemistryLab allows you to create and manipulate chemical formulas. It is based on `Formula` which is a structure (`struct`) which contains an expression, a writing of the formula close to those found in the Phreeqc databases, a unicode expression as well as a composition in the form of dictionaries and a charge.

```julia
struct Formula{T<:Number}
    expr::String
    phreeqc::String
    unicode::String
    colored::String
    composition::OrderedDict{Symbol,T}
    charge::Int8
end
```

## Formula construction

 Formulas can be constructed:
- by parsing a string containing eventually fractional or decimal coefficients
```@example 1
using ChemistryLab #hide
fgen = Formula("C3AFS5//8H4.32")
```

- from symbols representing atoms 
```@example
using ChemistryLab #hide
fCO2 = :C + 2 * :O
```

Charges can also be included during the creation in two different ways:
```@example 1
fHSO₄⁻ = :H+:S+4*:O+:e
```

Or:
```@example 1
fNa⁺ = :Na+:Zz
```

## Type of Formula

The type of the `Formula` `struct` being associated with the most complex type of the set of coefficients.

```@example 1
typeof(Formula("H2O"))
```

```@example 1
typeof(Formula("C3AFS5//8H4.32"))
```

## Change of type

Coefficient types can be converted *a posteriori*.

```@example 1
convert(Float64, Formula("H2O"))
```



---


