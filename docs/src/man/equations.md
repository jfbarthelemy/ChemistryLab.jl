# Chemical reactions

In CementChermistry it is possible to build chemical reactions and manipulate them. A reaction is constructed as a structure, "a composite data type that allows you to store multiple values in a single object". The `struct` is organized as follows:

```julia
struct Reaction{SR<:AbstractSpecies, TR<:Number, SP<:AbstractSpecies, TP<:Number}
    equation::String
    colored::String
    reactants::OrderedDict{SR, TR}
    products::OrderedDict{SP, TP}
    equal_sign::Char
end
```

## Parsing reactions

`Reaction` is a structure which can be build from:
- a string containing [species](./databases.md#species)
```julia
equation = "13H⁺ + NO₃⁻ + CO₃²⁻ + 10e⁻ = 6H₂O@ + HCN@"
reac, prod, equal_sign = parse_equation(equation)
```
- a string containing [cement species](./databases.md#species)
```julia
eqC3S = "C₃S + 5.3H = 1.3CH + C₁.₇SH₄"
rC3S = CemReaction(eqC3S)
```
- an operation on species
```@example
using CementChemistry
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
CSH = CemSpecies("C1.7SH4")
r = C3S + 5.3H ↔ 1.3CH + CSH
typeof(r)
```
- a balance calculation
```@example
using CementChemistry
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
CSH = CemSpecies("C1.7SH4")
r = Reaction([C3S, H, CH, CSH]; equal_sign='→')
```
- a balance calculation with symbolic numbers
```julia
using CementChemistry
using SymPy
â, b̂, ĝ = symbols("â b̂ ĝ", real=true)
CSH = CemSpecies(Dict(:C => â, :S => one(Sym), :H => ĝ))
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
r = Reaction([CSH, C3S, H, CH]; equal_sign='→')
```

```julia
using CementChemistry
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
CSH = CemSpecies("C1.7SH4")
r = map(simplify, Reaction([C3S, H], [CH, CSH]; equal_sign='→'))
A, _, _ = stoich_matrix([C3S], [CSH, H, CH]; involve_all_atoms=true) ;
``` 
