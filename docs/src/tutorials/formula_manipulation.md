# Chemical Formula Manipulation

ChemistryLab allows you to create and manipulate chemical formulas. It provides the `Formula` type (holds a formula string, multiple string representations, an atom composition map and a charge) and `AtomGroup` helper values, plus utilities to convert between notations (Phreeqc ↔ Unicode), to format/color formulas, to perform arithmetic on formulas, and to validate atom symbols.

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
f = Formula("C3AFS5//8H4.32")
```

- from a dictionary

```@example 1
using ChemistryLab #hide
fCaCO3 = Formula(Dict(:Ca => 1, :C => 1, :O => 3))
```

- from atom groups

```@example
using ChemistryLab #hide
fCO2 = AtomGroup(:C) + AtomGroup(2,:O)
```

Charges can also be included during the creation in two different ways:

```@example 1
fHSO₄⁻ = AtomGroup(:H)+AtomGroup(:S)+AtomGroup(4,:O)+AtomGroup(:e)
```

Or:

```@example 1
fNa⁺ = AtomGroup(:Na)+AtomGroup(:Zz)
```

## Formula accessors

All fields of a `Formula` are accessible through dedicated functions:

```@example 1
f = Formula("Ca(HSiO3)+")
expr(f)         # canonical string (as given)
phreeqc(f)      # Phreeqc notation: "Ca(HSiO3)+"
unicode(f)      # Unicode notation: "Ca(HSiO₃)⁺"
colored(f)      # colored terminal string (ANSI)
composition(f)  # OrderedDict{Symbol,Int}: composition by atom
charge(f)       # Int8: net charge
```

## Notation conversion

Switch between Phreeqc (plain text) and Unicode (pretty) notations without creating a full `Formula`:

```@example 1
phreeqc_to_unicode("SO4-2")   # "SO₄²⁻"
unicode_to_phreeqc("SO₄²⁻")  # "SO4-2"
phreeqc_to_unicode("Ca+2")    # "Ca²⁺"
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

## Formula arithmetic

Scalar multiplication and division return a new `Formula` with all coefficients scaled:

```@example 1
f = Formula("H2O")

f * 2     # H₄O₂  — all coefficients × 2
f / 2     # H₁O½  — all coefficients ÷ 2 (Float64)
f // 2    # H₁O½  — rational coefficients (Rational)
```

Individual atoms can be appended with `+ AtomGroup(n, :Symbol)`:

```@example 1
fCO2 = Formula("CO2")
fCO2 + AtomGroup(2, :H)   # CO₂H₂  — add 2 hydrogen atoms
```

For element-wise transformations, use `apply`:

```@example 1
apply(x -> x * 0.5, Formula("H2O"))   # H₁O½  — same as f / 2
```

!!! note "Advanced operations"
    For fractional stoichiometry, reaction algebra, and notation conversion see the [Advanced Topics](@ref) page.

---
