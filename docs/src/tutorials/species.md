
# [Species](@id sec-species)

`Species` is a composite type (introduced by the keyword `struct`) and is defined by a human-readable name, a chemical symbol/notation, an underlying `Formula` object holding composition, charge, and string representations, a physical state `aggregate_state`, a species `class`, as well as an extensible map of custom properties (molar mass, thermodynamic data, etc.)


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
    - `aggregate_state` denotes the state of the species (solid, liquid, gas) for which the possible keywords are `AS_AQUEOUS`, `AS_CRYSTAL`, `AS_GAS` and `AS_UNDEF`
    - `class` defines the role played by the species in the solution. The possible keywords are `SC_AQSOLVENT`, `SC_AQSOLUTE`, `SC_COMPONENT`, `SC_GASFLUID` and `SC_UNDEF`
    - `properties` refers to the set of properties intrinsic to the species. These properties are detailed below. 

## Species construction

`Species` can be created from:

- a `Formula`

```@example
using ChemistryLab #hide
fH2O = Formula("H2O")
H2O = Species(fH2O, aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)
```

- a string

```@example
using ChemistryLab #hide
HSO4⁻ = Species("HSO₄⁻", aggregate_state=AS_AQUEOUS, class=SC_COMPONENT)
```

- a dictionary

```@example
using ChemistryLab #hide
CO2 = Species(Dict(:C => 1, :O => 2), aggregate_state=AS_GAS, class=SC_GASFLUID)
```

!!! note "Adding charge"
    To add a charge when creating species with a dictionary, you must add, after the dictionary, the value of the charge (charge is considered an argument of the composite type).

```@example
using ChemistryLab #hide
SiO₃²⁻ = Species(Dict(:Si => 1, :O => 3), -2, aggregate_state=AS_AQUEOUS, class=SC_COMPONENT)
```

Keyword arguments such as `name`, `symbol`, `aggregate_state`, `class` can be added during construction.

```@example H2O
using ChemistryLab #hide
fH₂O = Formula("H2O")
H₂O = Species(fH₂O; name="Water", symbol="H₂O@", aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)
```

And `symbol` accept unicode characters.

```@example CO2
using ChemistryLab #hide
CO₂ = Species(Dict(:C=>1, :O=>2); name="Carbon dioxide", symbol="CO₂⤴", aggregate_state=AS_GAS, class=SC_GASFLUID)
```

!!! note "Comparison between species"
    Comparison between species (or cemspecies) are done by comparing atoms, aggregate_state and class. In the example below, vapor is not equal to H₂O since *aggregate_state* and *class* are different despite atoms are identical.
    ```julia
    vapour = Species("H2O"; name="Vapour", symbol="H₂O⤴", aggregate_state=AS_GAS, class=SC_GASFLUID)
    vapour == H₂O
    ```

!!! tip "Remark"
    You will also have noticed that a calculation of the molar mass of the species is systematically carried out.

---

## Species properties

The molar mass is automatically calculated and stored in the species `properties` dict under the key `:M`. Beyond that, the properties dict is open: any value of type `Number`, `AbstractVector{<:Number}`, `Function`, or `AbstractString` can be added at any time.

**Predefined property:**

| Key | Type | Description |
| :-- | :-- | :-- |
| `:M` | `Quantity` (g/mol) | Molar mass, computed automatically from the formula |

**Common user-added thermodynamic properties (loaded from databases or set manually):**

| Key | Type | Description |
| :-- | :-- | :-- |
| `:Cp⁰` | `SymbolicFunc(T)` | Standard heat capacity (J mol⁻¹ K⁻¹) |
| `:ΔₐH⁰` | `SymbolicFunc(T)` | Standard enthalpy of formation (J mol⁻¹) |
| `:S⁰` | `SymbolicFunc(T)` | Standard entropy (J mol⁻¹ K⁻¹) |
| `:ΔₐG⁰` | `SymbolicFunc(T)` | Standard Gibbs free energy of formation (J mol⁻¹) |
| `:V⁰` | `Quantity` or `Function` | Molar volume (m³ mol⁻¹) |

Properties are accessed and mutated via `[]`:

```@example props_simple
using ChemistryLab
using DynamicQuantities

H2O = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)

# Molar mass is always available
H2O[:M]
```

```@example props_simple
# Add a scalar property
H2O[:V⁰] = 18.07e-6u"m^3/mol"
H2O[:V⁰]
```

```@example props_simple
# Add a string annotation
H2O[:source] = "CRC Handbook 2024"
H2O[:source]
```

### Attaching temperature-dependent thermodynamic functions

The standard workflow uses [`build_thermo_functions`](@ref) to construct callable `SymbolicFunc` objects from reference data and Cp polynomial coefficients, then assigns them to the species:

```@example CO2
using ChemistryLab
using DynamicQuantities

CO₂ = Species("CO2"; name = "Carbon dioxide", aggregate_state = AS_GAS, class = SC_GASFLUID)

# Reference data at 298.15 K (from e.g. thermoddem.brgm.fr)
params_Cp_CO2 = Dict(
    :S⁰   => 213.785u"J/K/mol",
    :ΔₐH⁰ => -393510u"J/mol",
    :ΔₐG⁰ => -394373u"J/mol",
    :a₀   => 33.98u"J/K/mol",       # Cp polynomial: only a₀ and a₁ non-zero here
    :a₁   => 23.88e-3u"J/(mol*K^2)",
    :a₂   => 0.0u"J*K/mol",
    :a₃   => 0.0u"J/(mol*K^0.5)",
    :T    => 298.15u"K",             # reference temperature
)

dtf_CO2 = build_thermo_functions(:cp_ft_equation, params_Cp_CO2)
```

```@example CO2
# Assign each function to the species
for (k, v) in dtf_CO2
    CO₂[k] = v
end

# Evaluate at different temperatures
CO₂[:Cp⁰](T = 298.15)    # J/mol/K at 25 °C
```

```@example CO2
CO₂[:Cp⁰](T = 500.0)     # J/mol/K at 227 °C
```

```@example CO2
CO₂[:ΔₐG⁰](T = 500.0u"K", unit = true)   # Gibbs energy at 227 °C with units
```

!!! note "Cp polynomial"
    The built-in `:cp_ft_equation` model uses a 10-term polynomial:
    $\text{Cp}°(T) = a_0 + a_1 T + a_2 T^{-2} + a_3 T^{-0.5} + \ldots + a_{10} \log T$
    Unused coefficients should be set to `0.0` with appropriate units. Only non-zero terms affect the result. See the [Thermodynamic Functions](@ref sec-thermodynamics) tutorial for the full list of models and parameters.

!!! danger "Reference temperature"
    Always include `:T => reference_temperature` in the parameter dict. The functions are adjusted so that S°(T_ref), ΔₐH°(T_ref), ΔₐG°(T_ref) match the provided reference values exactly.

---

## Requalifying a species: `with_class`

`Species` is an **immutable** struct — none of its fields can be modified after construction.
[`with_class`](@ref) returns a **copy** of the species with only the `class` field changed.
All other fields — name, symbol, formula, aggregate state, properties and thermodynamic
functions — are shared by reference and unchanged.

```@example with_class
using ChemistryLab

# A crystal species with the class that databases typically assign
cal = Species("CaCO3"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)

println("before: ", class(cal))        # SC_COMPONENT

em = with_class(cal, SC_SSENDMEMBER)

println("after:  ", class(em))         # SC_SSENDMEMBER
println("formula preserved: ", formula(em) == formula(cal))
```

!!! note "Solid solutions: no requalification needed"
    [`SolidSolutionPhase`](@ref) automatically promotes end-members to `SC_SSENDMEMBER`
    at construction time. Database species with `SC_COMPONENT` can therefore be passed
    directly — no prior call to `with_class` is required:

    ```julia
    substances = build_species(datapath("cemdata18-thermofun.json"))
    dict = Dict(symbol(s) => s for s in substances)

    # Pass database species directly — SolidSolutionPhase requalifies internally
    cshq = SolidSolutionPhase("CSHQ", [dict["CSHQ-TobD"], dict["CSHQ-TobH"],
                                        dict["CSHQ-JenH"], dict["CSHQ-JenD"]])
    ```

    `with_class` is still useful when you need to track the requalified object
    explicitly (e.g. to inspect its class, or for types other than solid solution
    end-members).
