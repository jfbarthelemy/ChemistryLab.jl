# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using LinearAlgebra

# ── AbstractCalorimeter ───────────────────────────────────────────────────────

"""
    abstract type AbstractCalorimeter end

Base type for calorimeter models that can be coupled to a kinetics simulation.

Concrete subtypes:
- [`IsothermalCalorimeter`](@ref): T = constant, tracks Q(t) = ∫q̇dt.
- [`SemiAdiabaticCalorimeter`](@ref): variable-T cell (Lavergne et al. 2018).
"""
abstract type AbstractCalorimeter end

# _ensure_unit is defined in utils/misc.jl (imported at module level)

# ── Heat-rate from kinetic reactions ─────────────────────────────────────────

"""
    heat_rate(kinetic_reactions, rates, T_K) -> Real

Compute the instantaneous heat generation rate [W = J/s]:

```
q̇ = Σᵢ rᵢ(t) × ΔHᵣ,ᵢ(T)
```

where `rᵢ` [mol/s] is the net rate of the i-th kinetic reaction (positive =
dissolution/forward) and `ΔHᵣ,ᵢ(T)` [J/mol] is the enthalpy of reaction.

AD-compatible: the ΔₐH⁰ callables accept `ForwardDiff.Dual` T inputs.
"""
function heat_rate(
        kinetic_reactions::AbstractVector,
        rates::AbstractVector,
        T_K;
        kwargs...,
    )
    T_val = ustrip(T_K)
    q = zero(promote_type(eltype(rates), typeof(T_val)))
    for (kr, r) in zip(kinetic_reactions, rates)
        ΔHr = _reaction_enthalpy(kr, T_val)
        q = q + r * ΔHr
    end
    return q
end

# ── _reaction_enthalpy dispatch hierarchy ────────────────────────────────────
#
# Priority:
#   1. KineticReaction{R, F, Float64}  — explicit heat_per_mol
#   2. KineticReaction{R, F, Nothing}  — delegate to reaction stoichiometry
#   3. AbstractReaction                — stoichiometric sum of ΔₐH⁰

function _reaction_enthalpy(kr::KineticReaction{<:Any, <:Any, Float64}, ::Real)
    return kr.heat_per_mol
end

function _reaction_enthalpy(kr::KineticReaction{<:Any, <:Any, Nothing}, T_K::Real)
    return _reaction_enthalpy(kr.reaction, T_K)
end

function _reaction_enthalpy(reaction::AbstractReaction, T_K::Real)
    # reaction[:ΔᵣH⁰] is a SymbolicFunc or NumericFunc (T-dependent) built lazily
    # by complete_thermo_functions! from species :ΔₐH⁰ properties.
    # Thermodynamic convention: ΔᵣH⁰ < 0 for exothermic.
    # heat_rate requires "heat generated" (positive = exothermic), so we negate.
    if haskey(reaction, :ΔᵣH⁰)
        return -ustrip(reaction[:ΔᵣH⁰](; T = T_K * u"K", unit = true))
    end
    return zero(T_K)
end

# ── Total-enthalpy helper ─────────────────────────────────────────────────────

"""
    _total_enthalpy(n_full, h_fns, T_K) -> Real

Total molar enthalpy `H = Σᵢ nᵢ ΔₐH⁰ᵢ(T)`.
Used by the `DiscreteCallback` in `KineticsOrdinaryDiffEqExt`.
"""
function _total_enthalpy(n_full::AbstractVector, h_fns, T_K::Real)
    H = zero(promote_type(eltype(n_full), typeof(T_K)))
    for (i, hf) in enumerate(h_fns)
        isnothing(hf) && continue
        h_i = ustrip(hf(; T = T_K * u"K", unit = true))
        H += n_full[i] * h_i
    end
    return H
end

# ── IsothermalCalorimeter ─────────────────────────────────────────────────────

"""
    struct IsothermalCalorimeter{T} <: AbstractCalorimeter

Isothermal calorimeter: temperature held constant at `T` [K]; cumulative heat
`Q(t) = ∫₀ᵗ q̇(τ) dτ` [J] integrated as the trailing ODE state, with `q̇` the heat
of the **kinetic** reactions — see the caveat on [`cumulative_heat`](@ref).

# Examples

```julia
cal = IsothermalCalorimeter(298.15u"K")
kp = KineticsProblem(cs, reactions, state0, tspan; calorimeter = cal)
sol = integrate(kp, ks)
t, Q    = cumulative_heat(sol, cal)
t, qdot = heat_flow(sol, cal)
```
"""
struct IsothermalCalorimeter{T} <: AbstractCalorimeter
    T::T
    IsothermalCalorimeter{T}(T_K::T) where {T} = new{T}(T_K)
end

"""
    IsothermalCalorimeter(T) -> IsothermalCalorimeter

Plain `Real` → assumed SI [K]; `Quantity` → converted to K.
"""
function IsothermalCalorimeter(T_K)
    q = _ensure_unit(us"K", T_K)
    return IsothermalCalorimeter{typeof(q)}(q)
end

n_extra_states(::IsothermalCalorimeter) = 1

function extend_u0(u0::AbstractVector, ::IsothermalCalorimeter)
    return vcat(u0, zero(eltype(u0)))
end

function extend_ode!(du, ::Any, p, n_kin::Int, cal::IsothermalCalorimeter)
    T_val = Float64(safe_ustrip(us"K", cal.T))
    qdot = heat_rate(p.kin_rxns, p.rates_buf, T_val)
    du[n_kin + 1] = qdot
    return nothing
end

# ── SemiAdiabaticCalorimeter ──────────────────────────────────────────────────

"""
    struct SemiAdiabaticCalorimeter{C, T, F} <: AbstractCalorimeter

Semi-adiabatic calorimeter following the Lavergne et al. (2018) energy balance:

```math
\\frac{dT}{dt} = \\frac{\\dot{q}(t) - \\varphi(T - T_{\\rm env})}{C_p + \\sum_i n_i C^\\circ_{p,i}(T)}
```

where:
- `q̇(t)` [W] is the instantaneous heat-generation rate,
- `φ(ΔT)` [W] is the heat-loss function (e.g. linear `L·ΔT` or quadratic `a·ΔT + b·ΔT²`),
- `Cp` [J/K] is the fixed calorimeter heat capacity,
- `Σᵢ nᵢ Cp°ᵢ(T)` is the temperature- and mole-dependent sample heat capacity
  (computed from `p.cp_fns` at every ODE step when available).

# Fields

  - `Cp`: heat capacity of everything the ODE does **not** compute for itself —
    the vessel, the flask, an inert filler — in J/K, stored as a `Quantity`.
  - `heat_loss`: callable `φ(ΔT::Real) -> Real [W]`.
  - `T_env`: ambient temperature [K] (stored as `Quantity`).
  - `T0`: initial temperature [K] (stored as `Quantity`).

!!! warning "`Cp` must NOT include the sample"
    The denominator of the balance above is `Cp + Σᵢ nᵢ Cp°ᵢ(T)`, and the second
    term is recomputed from the database at every step. Adding the sample's own
    heat capacity to `Cp` therefore counts it **twice** and understates the
    temperature rise — by a factor of about 1.75 on a cement paste at w/c = 0.4,
    where the paste contributes some 2.5 kJ/K against 0.9 kJ/K for the flask.
    Three scripts shipped with this package made exactly that mistake, and the
    field description above used to invite it by saying "calorimeter + sample".
    Pass the vessel, and let the solver add the paste.

# Examples

```julia
# Linear heat loss — Newton cooling
cal = SemiAdiabaticCalorimeter(; Cp=4000.0u"J/K", T_env=293.15u"K", L=0.5u"W/K", T0=293.15u"K")

# Quadratic heat loss (Lavergne et al. 2018). `Cp` is the vessel alone;
# the paste's own heat capacity is added by the solver from the database.
cal = SemiAdiabaticCalorimeter(;
    Cp        = 900.0u"J/K",
    T_env     = 293.15u"K",
    heat_loss = ΔT -> 0.3*ΔT + 0.003*ΔT^2,
    T0        = 293.15u"K",
)

kp = KineticsProblem(cs, reactions, state0, tspan; calorimeter = cal)
sol = integrate(kp, ks)
t, T_vec = temperature_profile(sol, cal)
t, qdot  = heat_flow(sol, cal)
```

# References

  - Lavergne, F., Ben Fraj, A., Bayane, I. & Barthélémy, J.-F. (2018).
    *Cement and Concrete Research* **104**, 37–60.
"""
struct SemiAdiabaticCalorimeter{C, T, F} <: AbstractCalorimeter
    Cp::C       # heat capacity [J/K]
    heat_loss::F
    T_env::T    # ambient temperature [K]
    T0::T       # initial temperature [K]
end

"""
    SemiAdiabaticCalorimeter(; Cp, T_env, T0, heat_loss=nothing, L=nothing)

Keyword constructor for [`SemiAdiabaticCalorimeter`](@ref).

Exactly one of `heat_loss` or `L` must be provided:
  - `heat_loss`: callable `ΔT -> [W]` (e.g. quadratic `ΔT -> a*ΔT + b*ΔT^2`).
  - `L`: linear Newton cooling coefficient [W/K]. Sets `heat_loss = ΔT -> L * ΔT`.

All scalar fields accept plain `Real` (assumed SI) or `Quantity`:
  - `Cp` → J/K; `T_env` → K; `T0` → K; `L` → W/K.
"""
function SemiAdiabaticCalorimeter(; Cp, T_env, T0, heat_loss = nothing, L = nothing)
    Cp_q = _ensure_unit(us"J/K", Cp)
    T_env_q = _ensure_unit(us"K", T_env)
    T0_q = _ensure_unit(us"K", T0)
    hl = if !isnothing(heat_loss)
        heat_loss
    elseif !isnothing(L)
        L_f = Float64(safe_ustrip(us"W/K", L))
        ΔT -> L_f * ΔT
    else
        throw(
            ArgumentError(
                "SemiAdiabaticCalorimeter requires either `heat_loss` or `L`",
            ),
        )
    end
    return SemiAdiabaticCalorimeter(Cp_q, hl, T_env_q, T0_q)
end

n_extra_states(::SemiAdiabaticCalorimeter) = 1

function extend_u0(u0::AbstractVector, cal::SemiAdiabaticCalorimeter)
    T0_f = Float64(safe_ustrip(us"K", cal.T0))
    return vcat(u0, eltype(u0)(T0_f))
end

"""
    extend_ode!(du, u, p, n_kin, cal::SemiAdiabaticCalorimeter)

Append `dT/dt = (q̇ − φ(ΔT)) / Cp_total(T, n)` to the ODE right-hand side.

`Cp_total = Cp + Σᵢ nᵢ Cp°ᵢ(T)` is recomputed at every ODE step
from `p.cp_fns` and `p.n_full` (Lavergne et al. 2018).
"""
function extend_ode!(du, u, p, n_kin::Int, cal::SemiAdiabaticCalorimeter)
    T_curr = u[n_kin + 1]
    Cp_f = Float64(safe_ustrip(us"J/K", cal.Cp))
    T_env_f = Float64(safe_ustrip(us"K", cal.T_env))
    # Variable total heat capacity: Cp_calorimeter + Σᵢ nᵢ Cp°ᵢ(T)
    Cp_total = Cp_f
    for (i, cp_fn) in enumerate(p.cp_fns)
        isnothing(cp_fn) && continue
        cp_i = cp_fn(; T = T_curr, unit = false)   # J/(mol·K)
        Cp_total = Cp_total + p.n_full[i] * cp_i
    end
    qdot = heat_rate(p.kin_rxns, p.rates_buf, T_curr)
    ΔT = T_curr - T_env_f
    du[n_kin + 1] = (qdot - cal.heat_loss(ΔT)) / Cp_total
    return nothing
end

# ── Result extraction ─────────────────────────────────────────────────────────

"""
    heat_flow(sol, cal::IsothermalCalorimeter) -> (t, qdot)

Instantaneous heat-generation rate `q̇(t)` [W], by differencing
[`cumulative_heat`](@ref) — and carrying the same caveat about partial
equilibrium.
"""
function heat_flow(sol, cal::IsothermalCalorimeter)
    t, Q = cumulative_heat(sol, cal)
    qdot = similar(Q)
    qdot[1] = zero(eltype(Q))
    for i in 2:lastindex(t)
        dt = t[i] - t[i - 1]
        qdot[i] = dt > 0 ? (Q[i] - Q[i - 1]) / dt : zero(eltype(Q))
    end
    return t, qdot
end

"""
    heat_flow(sol, cal::SemiAdiabaticCalorimeter) -> (t, qdot)

Reconstruct q̇(t) [W] from the temperature ODE via the energy balance
`q̇ ≈ Cp × dT/dt + φ(T − T_env)`.

Note: uses the fixed `cal.Cp` (not the variable Cp_total) for this
post-processing reconstruction.
"""
function heat_flow(sol, cal::SemiAdiabaticCalorimeter)
    t = sol.t
    Cp_f = Float64(safe_ustrip(us"J/K", cal.Cp))
    T_env_f = Float64(safe_ustrip(us"K", cal.T_env))
    n_kin = length(sol.u[1]) - n_extra_states(cal)
    T_vec = [u[n_kin + 1] for u in sol.u]
    qdot = similar(T_vec)
    qdot[1] = zero(eltype(T_vec))
    for i in 2:lastindex(t)
        dt = t[i] - t[i - 1]
        dTdt = dt > 0 ? (T_vec[i] - T_vec[i - 1]) / dt : zero(eltype(T_vec))
        ΔT = T_vec[i] - T_env_f
        qdot[i] = Cp_f * dTdt + cal.heat_loss(ΔT)
    end
    return t, qdot
end

"""
    cumulative_heat(sol, cal::IsothermalCalorimeter) -> (t, Q)

Cumulative heat `Q(t) = ∫₀ᵗ q̇(τ) dτ` [J], read off the ODE state the isothermal
calorimeter adds.

!!! warning "This is the heat of the kinetic reactions only"
    `q̇` is [`heat_rate`](@ref), i.e. `Σᵢ rᵢ(−ΔᵣH⁰ᵢ)` over the kinetic reactions.
    That is the heat of hydration when those reactions produce the hydrates — the
    stoichiometric formulation. Under **partial equilibrium** they only dissolve
    the anhydrous phases into ions and the hydrates are precipitated by the Gibbs
    minimization, whose heat this sum cannot see; on an ordinary Portland cement
    the omission is worth hundreds of joules per gram. Use
    [`heat_release`](@ref), which differences the enthalpy of certified
    speciations, whenever `kp.equilibrium_solver !== nothing`.
"""
function cumulative_heat(sol, cal::IsothermalCalorimeter)
    n_kin = length(sol.u[1]) - n_extra_states(cal)
    Q = [u[n_kin + 1] for u in sol.u]
    return sol.t, Q
end

"""
    cumulative_heat(sol, cal::SemiAdiabaticCalorimeter) -> (t, Q)

Integrate the reconstructed heat-flow rate to obtain Q(t) [J].
"""
function cumulative_heat(sol, cal::SemiAdiabaticCalorimeter)
    t, qdot = heat_flow(sol, cal)
    Q = similar(qdot)
    Q[1] = zero(eltype(qdot))
    for i in 2:lastindex(t)
        dt = t[i] - t[i - 1]
        Q[i] = Q[i - 1] + qdot[i] * dt
    end
    return t, Q
end

"""
    temperature_profile(sol, cal::SemiAdiabaticCalorimeter) -> (t, T)

Extract the temperature profile T(t) [K].
"""
function temperature_profile(sol, ::SemiAdiabaticCalorimeter)
    n_kin = length(sol.u[1]) - 1
    T_vec = [u[n_kin + 1] for u in sol.u]
    return sol.t, T_vec
end
