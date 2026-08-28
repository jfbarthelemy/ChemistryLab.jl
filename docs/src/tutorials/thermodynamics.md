# [Thermodynamic Functions](@id sec-thermodynamics)

ChemistryLab represents temperature-dependent thermodynamic properties (Cp°, ΔₐH°, S°, ΔₐG°, log K°, …) as **callable function objects** rather than plain numbers. This page explains the three concrete types — `SymbolicFunc`, `NumericFunc`, and `ThermoFactory` — and shows how to build, combine, and extend them.

---

## Overview

| Type | Backed by | When to use |
| :-- | :-- | :-- |
| [`SymbolicFunc`](@ref) | Symbolics.jl expression, compiled to a `RuntimeGeneratedFunction` | Any model expressible as a closed-form symbolic formula |
| [`NumericFunc`](@ref) | Plain Julia closure | Complex models with no closed-form (e.g. HKF electrostatic integrals) |
| [`ThermoFactory`](@ref) | A `SymbolicFunc` template with free parameters | Factories that stamp out many `SymbolicFunc` instances from the same expression |

All three types share the same calling convention:

```julia
f(; T = 298.15, P = 1e5)          # raw numeric result (SI units)
f(; T = 350.0u"K", unit = true)   # result as a Quantity with units
f()                                # uses defaults stored in refs
```

---

## Simple function objects

### Constant

```@example sf_basics
using ChemistryLab
using DynamicQuantities

# Numeric constant (unitless)
f0 = SymbolicFunc(42.0)
f0()
```

```@example sf_basics
# Constant with physical unit
Cp_ref = SymbolicFunc(37.14u"J/mol/K")
Cp_ref()
```

```@example sf_basics
Cp_ref(unit = true)    # returns a Quantity
```

### Single variable

```@example sf_basics
# Just T — useful as a building block for arithmetic
f_T = SymbolicFunc(:T)
f_T(T = 350.0)
```

### Expression with parameters

Provide the expression as a Julia `Expr` literal. Any symbol that is not a declared variable (default list: `[:T, :P, :t, :x, :y, :z]`) is treated as a parameter and must be given as a keyword argument.

```@example sf_basics
# Linear: Cp = a₀ + a₁·T
f_lin = SymbolicFunc(:(a₀ + a₁ * T); a₀ = 30.0, a₁ = 5e-3)
f_lin(T = 298.15)
```

```@example sf_basics
# Quadratic: Cp = a + b·T + c/T²
f_quad = SymbolicFunc(:(a + b * T + c / T^2); a = 25.0, b = 8e-3, c = -1.5e5)
f_quad(T = 500.0)
```

### Calling conventions

```@example sf_basics
# Pass T as a plain number (K)
f_lin(T = 400.0)
```

```@example sf_basics
# Pass T as a Quantity — unit-aware evaluation
f_lin(T = 400.0u"K")
```

```@example sf_basics
# Return result with units attached
f_lin(T = 400.0u"K", unit = true)
```

---

## Arithmetic on function objects

`SymbolicFunc` objects support all standard arithmetic operations. The result is always another `SymbolicFunc` that can be evaluated the same way.

```@example sf_arith
using ChemistryLab
using DynamicQuantities

f1 = SymbolicFunc(:(a₀ + a₁ * T); a₀ = 30.0, a₁ = 2e-3)
f2 = SymbolicFunc(:(b₀ + b₁ * T); b₀ = 10.0, b₁ = 1e-3)

f_sum  = f1 + f2      # Cp of a reaction = Σ νᵢ Cpᵢ
f_diff = f1 - f2
f_scal = 2.5 * f1     # scaling by a stoichiometric coefficient

f_sum(T = 298.15)
```

```@example sf_arith
f_diff(T = 298.15)
```

```@example sf_arith
f_scal(T = 500.0)
```

!!! tip "Reaction thermodynamics"
    This is exactly how ChemistryLab computes `r.ΔᵣH⁰`, `r.ΔᵣG⁰`, etc.: each species contributes its `ΔₐH⁰` function multiplied by its stoichiometric coefficient, and the sum is another callable function of T.

---

## Built-in thermodynamic models

ChemistryLab ships three built-in models, registered in [`THERMO_MODELS`](@ref) and pre-compiled in [`THERMO_FACTORIES`](@ref).

### `:cp_ft_equation` — Maier-Kelley polynomial (default)

The standard heat-capacity polynomial for crystalline and gas-phase species:

```
Cp°(T) = a₀ + a₁T + a₂/T² + a₃/√T + a₄T² + a₅T³ + a₆T⁴ + a₇/T³ + a₈/T + a₉√T + a₁₀log(T)
```

S°, ΔₐH°, and ΔₐG° are the analytical integrals, adjusted so that the reference-state values S°(Tᵣ), ΔₐH°(Tᵣ), ΔₐG°(Tᵣ) match the provided data.

**Parameters** (all keyword arguments, in SI):

| Parameter | Unit | Description |
| :-- | :-- | :-- |
| `a₀` … `a₁₀` | see table below | Cp coefficients (unused ones set to 0) |
| `S⁰` | J mol⁻¹ K⁻¹ | Entropy at reference T |
| `ΔₐH⁰` | J mol⁻¹ | Enthalpy of formation at reference T |
| `ΔₐG⁰` | J mol⁻¹ | Gibbs free energy of formation at reference T |
| `T` | K | Reference temperature (usually 298.15 K) |

Coefficient units:

| Parameter | Unit |
| :-- | :-- |
| `a₀`, `a₁₀` | J/(mol·K) |
| `a₁` | J/(mol·K²) |
| `a₂` | J·K/mol |
| `a₃` | J/(mol·K^0.5) |
| `a₄` | J/(mol·K³) |
| `a₅` | J/(mol·K⁴) |
| `a₆` | J/(mol·K⁵) |
| `a₇` | J·K²/mol |
| `a₈` | J/mol |
| `a₉` | J/(mol·K^1.5) |

**Example — CO₂ (gas):**

```@example build_thermo
using ChemistryLab
using DynamicQuantities

params_CO2 = Dict(
    :S⁰   => 213.785u"J/K/mol",
    :ΔₐH⁰ => -393510.0u"J/mol",
    :ΔₐG⁰ => -394373.0u"J/mol",
    :a₀   => 33.98u"J/K/mol",
    :a₁   => 23.88e-3u"J/(mol*K^2)",
    :a₂   => 0.0u"J*K/mol",
    :a₃   => 0.0u"J/(mol*K^0.5)",
    :T    => 298.15u"K",
)

dtf = build_thermo_functions(:cp_ft_equation, params_CO2)
```

The returned `OrderedDict` has four entries:

```@example build_thermo
keys(dtf)
```

Evaluate at different temperatures:

```@example build_thermo
dtf[:Cp⁰](T = 298.15)   # J/mol/K at 25 °C
```

```@example build_thermo
dtf[:Cp⁰](T = 500.0)    # J/mol/K at 227 °C
```

```@example build_thermo
dtf[:ΔₐG⁰](T = 500.0)   # J/mol at 227 °C
```

```@example build_thermo
dtf[:ΔₐG⁰](T = 500.0u"K", unit = true)   # with units
```

### `:logk_fpt_function` — van't Hoff log K fit

Fits the equilibrium constant as a polynomial in T:

```
log₁₀K(T) = A₀ + A₁T + A₂/T + A₃ log T + A₄/T² + A₅T² + A₆/√T
```

**Parameters:** `A₀`…`A₆` (dimensionless / appropriate T-powers), `T` (reference, K).

This model produces a single function `:logKr` (the fit directly) rather than the full {Cp, H, S, G} set. It is used internally when loading database reactions.

### `:solute_hkf88_reaktoro` — HKF aqueous model

The Helgeson-Kirkham-Flowers (1981/1988) model for aqueous solutes. This is a **numeric** model (returns `NumericFunc` objects) because it depends on the equation of state of water and cannot be written as a simple closed-form polynomial.

**Parameters (all in SI):** `a1`–`a4`, `c1`, `c2`, `wref`, `z` (formal charge), `S⁰`, `ΔfH⁰`, `ΔfG⁰`.

---

## Attaching properties to a species

After building the thermodynamic function dictionary, assign each entry to the species `properties` dict:

```@example attach_props
using ChemistryLab
using DynamicQuantities

params_CO2 = Dict(
    :S⁰   => 213.785u"J/K/mol",
    :ΔₐH⁰ => -393510.0u"J/mol",
    :ΔₐG⁰ => -394373.0u"J/mol",
    :a₀   => 33.98u"J/K/mol",
    :a₁   => 23.88e-3u"J/(mol*K^2)",
    :a₂   => 0.0u"J*K/mol",
    :a₃   => 0.0u"J/(mol*K^0.5)",
    :T    => 298.15u"K",
)

CO2 = Species("CO2"; aggregate_state = AS_GAS, class = SC_GASFLUID)
dtf = build_thermo_functions(:cp_ft_equation, params_CO2)
for (k, v) in dtf
    CO2[k] = v
end

CO2[:Cp⁰](T = 400.0)
```

```@example attach_props
CO2[:ΔₐG⁰](T = 400.0u"K", unit = true)
```

---

## Defining a custom model

### From a Cp expression (automatic integration)

The simplest way to add a new model is to provide a Julia expression for Cp(T). ChemistryLab will **automatically integrate** it to obtain H, S, and G.

```@example custom_model
using ChemistryLab
using DynamicQuantities
using SymbolicNumericIntegration  # required for add_thermo_model with Cp expression

# Linear Cp = a + b·T — register with explicit units so build_thermo_functions works
add_thermo_model(:linear_Cp, :(a + b * T), [:T => u"K", :a => u"J/mol/K", :b => u"J/(mol*K^2)"])

# Build functions for a hypothetical species
params = Dict(
    :a    => 28.0u"J/mol/K",
    :b    => 4e-3u"J/(mol*K^2)",
    :S⁰   => 130.0u"J/mol/K",
    :ΔₐH⁰ => -50e3u"J/mol",
    :ΔₐG⁰ => -35e3u"J/mol",
    :T    => 298.15u"K",
)

dtf = build_thermo_functions(:linear_Cp, params)
dtf[:Cp⁰](T = 400.0)
```

```@example custom_model
dtf[:ΔₐG⁰](T = 400.0)
```

!!! note "Integration"
    `add_thermo_model(name, Cpexpr)` uses Symbolics.jl to integrate Cp analytically: H = ∫Cp dT, S = ∫(Cp/T) dT, G = -∫S dT. All three are then adjusted so that the reference-state quantities S°(Tᵣ), ΔₐH°(Tᵣ), ΔₐG°(Tᵣ) are satisfied.

### Constant Cp model

```@example custom_model2
using ChemistryLab
using DynamicQuantities
using SymbolicNumericIntegration  # required for add_thermo_model with Cp expression

# Cp = Cp0 — constant, independent of T
# Use :(Cp0 * T^0) so the expression is an Expr (not a bare Symbol) and T stays in the vars list
add_thermo_model(:const_Cp, :(Cp0 * T^0), [:T => u"K", :Cp0 => u"J/mol/K"])

params = Dict(
    :Cp0  => 75.3u"J/mol/K",       # water
    :S⁰   => 69.95u"J/mol/K",
    :ΔₐH⁰ => -285830.0u"J/mol",
    :ΔₐG⁰ => -237140.0u"J/mol",
    :T    => 298.15u"K",
)

dtf = build_thermo_functions(:const_Cp, params)
# Cp is constant → H is linear → S = Cp·ln(T) → G from H and S
dtf[:Cp⁰](T = 373.15)
```

```@example custom_model2
dtf[:ΔₐG⁰](T = 373.15)
```

### From a full dictionary

For maximum control, supply all four expressions directly:

```julia
using ChemistryLab
using DynamicQuantities

add_thermo_model(:my_model, Dict(
    :Cp => :(a + b * T + c / T^2),
    :H  => :(a * T + b * T^2 / 2 - c / T),
    :S  => :(a * log(T) + b * T - c / (2 * T^2)),
    :G  => :(-(a * T * log(T) - a * T) - b * T^2 / 2 - c / (2 * T)),
    :units => [:a => u"J/mol/K", :b => u"J/(mol*K^2)", :c => u"J*K/mol", :T => u"K"],
))
```

---

## ThermoFactory directly

If you need to evaluate a symbolic expression at many parameter values without going through a registered model, use [`ThermoFactory`](@ref) directly:

```@example factory_direct
using ChemistryLab
using DynamicQuantities

factory = ThermoFactory(:(a + b * T + c / T^2), [:T])

# Stamp out SymbolicFunc instances for different species
f_CO2 = factory(; a = 33.98, b = 23.88e-3, c = 0.0)
f_H2O = factory(; a = 30.54, b = 10.29e-3, c = 0.0)

f_CO2(T = 500.0), f_H2O(T = 500.0)
```

---

## NumericFunc

[`NumericFunc`](@ref) wraps a plain Julia closure. The calling convention is identical to `SymbolicFunc`.

```@example numericfunc
using ChemistryLab
using DynamicQuantities

# Build a NumericFunc manually (e.g. for a piecewise model)
f_piece = NumericFunc(
    (T, P) -> T < 500.0 ? 30.0 + 5e-3 * T : 32.5 + 2e-3 * T,
    (:T, :P),
    (T = 298.15u"K", P = 1e5u"Pa"),
    1.0u"J/mol/K",
)

f_piece(T = 300.0)
```

```@example numericfunc
f_piece(T = 600.0)
```

`NumericFunc` objects support the same arithmetic operations as `SymbolicFunc`. Cross-type arithmetic (e.g. `SymbolicFunc + NumericFunc`) returns a `NumericFunc`.

---

## Summary

| Task | Code |
| :-- | :-- |
| Constant with unit | `SymbolicFunc(37.14u"J/mol/K")` |
| Single variable T | `SymbolicFunc(:T)` |
| Expression with parameters | `SymbolicFunc(:(a + b*T); a=30.0, b=5e-3)` |
| Build from existing model | `build_thermo_functions(:cp_ft_equation, params)` |
| Register new model from Cp | `add_thermo_model(:mymodel, :(a + b*T))` |
| Call with units | `f(T=400.0u"K", unit=true)` |
| Arithmetic | `f1 + f2`, `2*f`, `-f` |
