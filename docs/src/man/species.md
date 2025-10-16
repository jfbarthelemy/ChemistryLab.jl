
# Species

`Species` is a composite type (introduced by the keyword `struct`) and is defined by a name, a symbol, a formula, an aggregate state, a class and properties. It creates chemical species for solution or solid phases:

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

!!! info "Advanced description"
    - `aggregate_state` denotes the state of the species (solid, liquid, gas) for which the possible keywords are AS_AQUEOUS, AS_CRYSTAL, AS_GAS and AS_UNDEF
    - `class` defines the role played by the species in the solution. The possible keywords are SC_AQSOLVENT, SC_AQSOLUTE, SC_COMPONENT, SC_GAS_FLUID and SC_UNDEF
    - `properties` refers to the set of properties intrinsic to the species. These properties are detailed below ([Species properties](@ref)). 

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
HSO4⁻ = Species("HSO₄⁻")
```

- a dictionary
```@example
using ChemistryLab #hide
CO2 = Species(Dict(:C => 1, :O => 2))
```

!!! note "Adding charge"
    To add a charge when creating species with a dictionary, you must add, after the dictionary, the value of the charge (charge is considered an argument of the composite type).

```@example
using ChemistryLab #hide
CO2 = Species(Dict(:Si => 1, :O => 3),-2)
```

Keyword arguments such as `name`, `symbol`, `aggregate_state`, `class` can be added during construction.

```@example H2O
using ChemistryLab #hide
fH₂O = 2*:H + :O
H₂O = Species(fH₂O; name="Water", symbol="H₂O@", aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)
```

And `symbol` accept unicode characters.
```@example
using ChemistryLab #hide
CO₂ = Species(Dict(:C=>1, :O=>2); name="Carbon dioxide", symbol="CO₂⤴", aggregate_state=AS_GAS, class=SC_GAS_FLUID)
```

!!! note "Comparison between species"
    Comparison between species are done by comparing atoms, aggregate_state and class. In the example below, vapour is not equal to H₂O since *aggregate_state* and *class* are different despite atoms are identical.
    ```julia
    vapour = Species(2*:H + :O; name="Vapour", symbol="H₂O⤴", aggregate_state=AS_GAS, class=SC_GAS_FLUID)
    vapour == H₂O
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

!!! warning "Warning"
    Not every molecule can be used to build a cement species. It is necessary for this molecule to decompose into a combination of the oxides present in the manufacturers' cement sheet (e.g. $CaO$, $SiO_2$, $Fe_2O_3$, $Al_2O_3$) and water. Thus, the following code will return an error.
    ```julia
    CemSpecies(Species("Ca(OH)"))
    ```


---

## Numeric and Symbolic CemSpecies

The previous species were constructed from integer values ​​of the number of chemical elements. However, other numerical value types ​​are possible (see [species](./databases.md#formulas)), such as fraction or Real values.

```@example
using ChemistryLab
ox = Dict(:C => 1.666667, :S => 1, :H => 2.1)
jennite = CemSpecies(ox)
```

Symbolic values are also allowed. In this case, you need to use the [`SymPy`](https://github.com/JuliaPy/SymPy.jl) library:

```@example sympy1
using ChemistryLab
using SymPy
â, ĝ = symbols("â ĝ", real = true)
ox = Dict(:C => â, :S => one(Sym), :H => ĝ)
CSH = CemSpecies(ox)
```

The value of variables can be defined *a posteriori*.

```julia
jennite = CemSpecies(map(N, map(subs, cemformula(CSH), â => 1.666667, ĝ => 2.1)))
```

!!! note "Remark"
    Conversion of coefficient types can also be done.
    ```julia
    floatCSH = Species(convert(Float64, formula(numCSH)))
    ```

---

## Conversion to Cement Notation

Convert species to cement notation and Unicode. Conversion can be done on simple species:

```@example example_cemspecies
H2O = Species("H₂O")
cemH2O = CemSpecies(H2O)
```

Or more complex one:

```@example CSH
using ChemistryLab #hide
CSH = Species("(SiO2)1(CaO)1.666667(H2O)2.1")
jennite = CemSpecies(CSH)
```

---

## Species properties

Species properties are open and left to the discretion of users. Only the molar mass is systematically calculated and integrated into the species properties, for now. We can of course imagine that these properties could contain thermodynamic properties such as the Gibbs energy of formation, the heat capacity or even the entropy variation, these properties themselves being temperature dependent. These properties must nevertheless respect one of the following types: `Number`, `AbstractVector{<:Number}`, `Function`, `AbstractString`.

Imagine, for example, that we wanted to construct the jennite ($C_{1.67}SH_{2.1}$) molecule with some of its thermodynamic properties. The Gibbs energy of formation of this species is equal to -2480.81 KJ/mol. This property, intrinsic to the species, can be added simply as follows:

```@example CSH
import Unitful: u, g, cm, K, J, mol, Quantity, uconvert, ustrip, unit, uparse
jennite.ΔfG⁰ = -2480.81*uparse("K * J/mol")
```

```@example CSH
function Cₚ(T, a)
    y= a[1] + a[2] * T + a[3] * T^(-2) + a[4] * T^(−0.5)
    return y
end
jennite.Cₚ = f(T) = Cₚ(T, [210.0, 0.120, -3.07e6, 0.0])
```


