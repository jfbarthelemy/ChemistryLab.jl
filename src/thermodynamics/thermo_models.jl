# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
# Portions of this file are Julia ports adapted from the Reaktoro C++ library
# (https://github.com/reaktoro/reaktoro), Copyright © 2014-2024 Allan Leal,
# distributed under the LGPL-2.1-or-later.

using DynamicQuantities
using OrderedCollections

"""
    THERMO_MODELS

Dictionary storing raw thermodynamic model expressions and units.
Keys are model names (symbols), values are dictionaries containing:

  - Symbolic expressions for thermodynamic functions (Cp, H, S, G).
  - Units for parameters and variables.
"""
const THERMO_MODELS = Dict(
    :cp_ft_equation => Dict(
        :Cp => :(
            a₀ +
                a₁ * T +
                a₂ / T^2 +
                a₃ / sqrt(T) +
                a₄ * T^2 +
                a₅ * T^3 +
                a₆ * T^4 +
                a₇ / T^3 +
                a₈ / T +
                a₉ * sqrt(T) +
                a₁₀ * log(T)
        ),
        :S => :(
            a₀ * log(T) +
                a₁ * T +
                -(a₂ / 2) / T^2 +
                -2 * a₃ / sqrt(T) +
                (a₄ / 2) * T^2 +
                (a₅ / 3) * T^3 +
                (a₆ / 4) * T^4 +
                -(a₇ / 3) / T^3 +
                -a₈ / T +
                2 * a₉ * sqrt(T) +
                (a₁₀ / 2) * (log(T))^2
        ),
        :H => :(
            a₀ * T +
                a₁ * T^2 / 2 +
                -a₂ / T +
                2 * a₃ * sqrt(T) +
                (a₄ / 3) * T^3 +
                (a₅ / 4) * T^4 +
                (a₆ / 5) * T^5 +
                -(a₇ / 2) / T^2 +
                a₈ * log(T) +
                (2 / 3) * a₉ * T^(3 / 2) +
                a₁₀ * T * log(T) - a₁₀ * T
        ),
        :G => :(
            -a₀ * T * log(T) +
                a₀ * T +
                -(a₁ / 2) * T^2 +
                -(a₂ / 2) / T +
                4 * a₃ * sqrt(T) +
                -(a₄ / 6) * T^3 +
                -(a₅ / 12) * T^4 +
                -(a₆ / 20) * T^5 +
                -(a₇ / 6) / T^2 +
                a₈ * log(T) +
                -(4 / 3) * a₉ * T^(3 / 2) +
                -(a₁₀ / 2) * T * (log(T))^2 +
                a₁₀ * T * log(T) - a₁₀ * T
        ),
        :units => [
            :a₀ => u"J/(mol*K)",
            :a₁ => u"J/(mol*K^2)",
            :a₂ => u"J*K/mol",
            :a₃ => u"J/(mol*K^(1//2))",
            :a₄ => u"J/(mol*K^3)",
            :a₅ => u"J/(mol*K^4)",
            :a₆ => u"J/(mol*K^5)",
            :a₇ => u"J*K^2/mol",
            :a₈ => u"J/mol",
            :a₉ => u"J/(mol*K^(3//2))",
            :a₁₀ => u"J/(mol*K)",
            :T => u"K",
            # :Cp => "J/(mol*K)",
            # :S => "J/(mol*K)",
            # :H => "J/mol",
            # :G => "J/mol",
        ],
    ),
    :logk_fpt_function => Dict(
        :logKr =>
            :(A₀ + A₁ * T + A₂ / T + A₃ * log(T) + A₄ / T^2 + A₅ * T^2 + A₆ * sqrt(T)),
        :units => [
            :A₀ => u"1",
            :A₁ => u"K^(-1)",
            :A₂ => u"K",
            :A₃ => u"1",
            :A₄ => u"K^2",
            :A₅ => u"K^(-2)",
            :A₆ => u"K^(-1//2)",
            :T => u"K",
            # :logKr => "1",
        ],
    ),
)

"""
    THERMO_FACTORIES

Dictionary storing compiled `ThermoFactory` objects for each model.
Used to efficiently generate `SymbolicFunc` instances.
"""
const THERMO_FACTORIES = Dict{Symbol, Dict{Symbol, ThermoFactory}}()

"""
    build_thermo_functions(model_name, params) -> OrderedDict

Build thermodynamic function objects for a specific model and parameters.
Dispatches on `Val(model_name)` — add a new method for each new model.

# Arguments

  - `model_name`: symbol identifying the thermodynamic model (e.g., `:cp_ft_equation`).
  - `params`: dictionary or pair list of parameter values.

# Returns

  - `OrderedDict` containing the constructed thermodynamic functions (`Cp⁰`, `ΔₐH⁰`, `S⁰`, `ΔₐG⁰`).
"""
build_thermo_functions(model_name::Symbol, params) =
    build_thermo_functions(Val(model_name), params)

# Default: use THERMO_FACTORIES (symbolic models)
function build_thermo_functions(::Val{M}, params) where {M}
    dict_factories = THERMO_FACTORIES[M]
    dict_params = Dict(params)

    STref = dict_params[:S⁰]
    HTref = get(
        dict_params, :ΔfH⁰, get(dict_params, :ΔₐH⁰, get(dict_params, :ΔaH⁰, missing))
    )
    GTref = get(
        dict_params, :ΔfG⁰, get(dict_params, :ΔₐG⁰, get(dict_params, :ΔaG⁰, missing))
    )
    Tref = dict_params[:T]

    Cp⁰ = dict_factories[:Cp](; params...)

    H = dict_factories[:H](; params...)
    ΔₐH⁰ = H + (HTref - H(; T = Tref, unit = true))

    S = dict_factories[:S](; params...)
    δS⁰ = STref - S(; T = Tref, unit = true)
    S⁰ = S + δS⁰

    T = SymbolicFunc(:T; units = [:T => u"K"], output_unit = u"K")
    if haskey(dict_factories, :G)
        G = dict_factories[:G](; params...)
        ΔₐG⁰ = (G - T * δS⁰) + (GTref - G(; T = Tref, unit = true) + Tref * δS⁰)
    else
        ΔₐG⁰ = (H - T * S⁰) + (GTref - H(; T = Tref, unit = true) + Tref * STref)
    end

    return OrderedDict(:Cp⁰ => Cp⁰, :ΔₐH⁰ => ΔₐH⁰, :S⁰ => S⁰, :ΔₐG⁰ => ΔₐG⁰)
end

# HKF dispatch
function build_thermo_functions(::Val{:solute_hkf88_reaktoro}, params)
    return _build_hkf_thermo_functions(params)
end

# ============================================================
#  HKF model constants
# ============================================================

const _HKF_Tr = 298.15   # reference temperature (K)
const _HKF_Pr = 1.0e+5   # reference pressure    (Pa)
const _HKF_Zr = -1.278055636e-2  # Born function Z at (Tr, Pr)
const _HKF_Yr = -5.795424563e-5  # Born function Y at (Tr, Pr)
const _HKF_θ = 228.0    # θ constant (K)
const _HKF_Ψ = 2.6e+8   # Ψ constant (Pa)

"""
    _build_hkf_thermo_functions(params) -> OrderedDict

Internal builder for the HKF (Helgeson-Kirkham-Flowers 1981/1988) standard thermodynamic
model for aqueous solutes.

`params` must contain (all in SI units):
  - `:a1`, `:a2`, `:a3`, `:a4`  — equation-of-state coefficients (J·mol⁻¹·Pa⁻¹, etc.)
  - `:c1`, `:c2`                 — heat-capacity coefficients (J·mol⁻¹·K⁻¹, J·K·mol⁻¹)
  - `:wref`                      — reference Born coefficient ω_ref (J/mol)
  - `:z`                         — species charge (dimensionless)
  - `:S⁰` or `:Sr`               — standard entropy at (Tr, Pr) (J·mol⁻¹·K⁻¹)
  - `:ΔₐH⁰` or `:ΔfH⁰`          — standard enthalpy of formation (J/mol)
  - `:ΔₐG⁰` or `:ΔfG⁰`          — standard Gibbs energy of formation (J/mol)
"""
function _build_hkf_thermo_functions(params)
    dp = Dict(params)

    # Extract SI values (strip units if present)
    _strip(x::AbstractQuantity) = ustrip(uexpand(x))
    _strip(x::Real) = Float64(x)

    a1 = _strip(dp[:a1])
    a2 = _strip(dp[:a2])
    a3 = _strip(dp[:a3])
    a4 = _strip(dp[:a4])
    c1 = _strip(dp[:c1])
    c2 = _strip(dp[:c2])
    wref = _strip(dp[:wref])
    z = _strip(dp[:z])

    Sr = _strip(get(dp, :S⁰, get(dp, :Sr, 0.0)))
    Hf = _strip(get(dp, :ΔₐH⁰, get(dp, :ΔfH⁰, 0.0)))
    Gf = _strip(get(dp, :ΔₐG⁰, get(dp, :ΔfG⁰, 0.0)))

    # Build refs as Quantities in SI (consistent with SymbolicFunc.refs)
    T_raw = get(dp, :T, _HKF_Tr)
    P_raw = get(dp, :P, _HKF_Pr)
    refs = (
        T = T_raw isa AbstractQuantity ? force_uconvert(u"K", T_raw) : T_raw * u"K",
        P = P_raw isa AbstractQuantity ? force_uconvert(u"Pa", P_raw) : P_raw * u"Pa",
    )
    vars = (:T, :P)

    Tr = _HKF_Tr
    Pr = _HKF_Pr
    θ = _HKF_θ
    Ψ = _HKF_Ψ
    Zr = _HKF_Zr
    Yr = _HKF_Yr

    # -- Closures (T in K, P in Pa) --

    function _Cp(T::Real, P::Real)
        wtp = water_thermo_props(T, P)
        wep = water_electro_props_jn(T, P, wtp)
        gs = hkf_g_function(T, P, wtp)
        ae = species_electro_props_hkf(gs, z, wref)
        Tth = T - θ
        return c1 + c2 / (Tth * Tth) -
            2 * T / (Tth^3) * (a3 * (P - Pr) + a4 * log((Ψ + P) / (Ψ + Pr))) +
            ae.w * T * wep.bornX + 2 * T * wep.bornY * ae.wT +
            T * (wep.bornZ + 1) * ae.wTT
    end

    function _H(T::Real, P::Real)
        wtp = water_thermo_props(T, P)
        wep = water_electro_props_jn(T, P, wtp)
        gs = hkf_g_function(T, P, wtp)
        ae = species_electro_props_hkf(gs, z, wref)
        Tth = T - θ
        Tth2 = Tth * Tth
        return Hf + c1 * (T - Tr) - c2 * (1 / Tth - 1 / (Tr - θ)) +
            a1 * (P - Pr) + a2 * log((Ψ + P) / (Ψ + Pr)) +
            (2 * T - θ) / Tth2 * (a3 * (P - Pr) + a4 * log((Ψ + P) / (Ψ + Pr))) -
            ae.w * (wep.bornZ + 1) + ae.w * T * wep.bornY +
            T * (wep.bornZ + 1) * ae.wT + wref * (Zr + 1) - wref * Tr * Yr
    end

    function _S(T::Real, P::Real)
        wtp = water_thermo_props(T, P)
        wep = water_electro_props_jn(T, P, wtp)
        gs = hkf_g_function(T, P, wtp)
        ae = species_electro_props_hkf(gs, z, wref)
        Tth = T - θ
        Tth2 = Tth * Tth
        return Sr + c1 * log(T / Tr) -
            c2 / θ * (1 / Tth - 1 / (Tr - θ) + log(Tr / T * Tth / (Tr - θ)) / θ) +
            1 / Tth2 * (a3 * (P - Pr) + a4 * log((Ψ + P) / (Ψ + Pr))) +
            ae.w * wep.bornY + (wep.bornZ + 1) * ae.wT - wref * Yr
    end

    function _G(T::Real, P::Real)
        wtp = water_thermo_props(T, P)
        wep = water_electro_props_jn(T, P, wtp)
        gs = hkf_g_function(T, P, wtp)
        ae = species_electro_props_hkf(gs, z, wref)
        Tth = T - θ
        return Gf - Sr * (T - Tr) - c1 * (T * log(T / Tr) - T + Tr) +
            a1 * (P - Pr) + a2 * log((Ψ + P) / (Ψ + Pr)) -
            c2 * (
            (1 / Tth - 1 / (Tr - θ)) * (θ - T) / θ -
                T / (θ * θ) * log(Tr / T * Tth / (Tr - θ))
        ) +
            1 / Tth * (a3 * (P - Pr) + a4 * log((Ψ + P) / (Ψ + Pr))) -
            ae.w * (wep.bornZ + 1) + wref * (Zr + 1) + wref * Yr * (T - Tr)
    end

    function _V(T::Real, P::Real)
        wtp = water_thermo_props(T, P)
        wep = water_electro_props_jn(T, P, wtp)
        gs = hkf_g_function(T, P, wtp)
        ae = species_electro_props_hkf(gs, z, wref)
        Tth = T - θ
        return a1 + a2 / (Ψ + P) + (a3 + a4 / (Ψ + P)) / Tth -
            ae.w * wep.bornQ - (wep.bornZ + 1) * ae.wP
    end

    return OrderedDict(
        :Cp⁰ => NumericFunc(_Cp, vars, refs, u"J/(mol*K)"),
        :ΔₐH⁰ => NumericFunc(_H, vars, refs, u"J/mol"),
        :S⁰ => NumericFunc(_S, vars, refs, u"J/(mol*K)"),
        :ΔₐG⁰ => NumericFunc(_G, vars, refs, u"J/mol"),
        :V⁰ => NumericFunc(_V, vars, refs, u"m^3/mol"),
    )
end

"""
    _THERMO_OUTPUT_UNITS

Explicit output units for each thermodynamic function key.
Used by `build_thermo_factories` to populate `ThermoFactory.output_unit`
without relying on symbolic unit propagation.
"""
const _THERMO_OUTPUT_UNITS = Dict{Symbol, AbstractQuantity}(
    :Cp => u"J/(mol*K)",
    :S => u"J/(mol*K)",
    :H => u"J/mol",
    :G => u"J/mol",
    :logKr => u"1",
)

"""
    build_thermo_factories(dict_expr) -> Dict

Helper function to build `ThermoFactory` objects from a model dictionary.
"""
function build_thermo_factories(dict_expr)
    return Dict(
        k => ThermoFactory(
            v,
            [:T, :P];
            units = get(dict_expr, :units, nothing),
            output_unit = get(_THERMO_OUTPUT_UNITS, k, u"1"),
        ) for (k, v) in dict_expr if k != :units
    )
end

"""
    add_thermo_model(model_name, dict_model::AbstractDict)

Add a new thermodynamic model to the registry using a dictionary of expressions.

# Arguments

  - `model_name`: unique symbol for the model.
  - `dict_model`: dictionary containing symbolic expressions and units.

See also: `add_thermo_model(model_name, Cpexpr::Expr, units)` — variant that automatically
integrates a Cp expression to obtain H, S and G. Requires loading `SymbolicNumericIntegration`:

```julia
using SymbolicNumericIntegration
add_thermo_model(:my_model, :(a + b*T), [:T => u"K", :a => u"J/mol/K", :b => u"J/(mol*K^2)"])
```
"""
function add_thermo_model(model_name, dict_model::AbstractDict)
    THERMO_MODELS[model_name] = dict_model
    return THERMO_FACTORIES[model_name] = build_thermo_factories(dict_model)
end
