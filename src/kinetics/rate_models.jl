# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
# Portions of this file (the Arrhenius rate constant and the saturation-ratio
# formulation used by the transition-state theory rate model) are Julia ports
# adapted from the Reaktoro C++ library (https://github.com/reaktoro/reaktoro),
# Copyright © 2014-2024 Allan Leal, distributed under the LGPL-2.1-or-later.

using DynamicQuantities
using OrderedCollections

# ── StateView ─────────────────────────────────────────────────────────────────

"""
    StateView{T, I <: AbstractDict}

Thin wrapper giving O(1) named access to a species data vector.

```julia
sv["C3S"] === sv.data[sv.index["C3S"]]
```

The `index` dict is built once at [`KineticsProblem`](@ref) construction;
`data` is a plain vector (mutated in-place or re-wrapped each ODE step) — no
dict allocation in the hot path.

# Examples

```jldoctest
julia> idx = Dict("Ca++" => 1, "C3S" => 2);

julia> sv = StateView([0.5, 1.0], idx);

julia> sv["C3S"]
1.0

julia> haskey(sv, "Ca++")
true
```
"""
struct StateView{T, I <: AbstractDict}
    data::AbstractVector{T}
    index::I
end

Base.getindex(sv::StateView, name::AbstractString) = sv.data[sv.index[name]]
Base.haskey(sv::StateView, name::AbstractString) = haskey(sv.index, name)

# ── KineticFunc ───────────────────────────────────────────────────────────────

"""
    KineticFunc{F, R <: NamedTuple, Q}

Compiled kinetic rate function, analogous to [`NumericFunc`](@ref) for
thermodynamics.

Calling convention (positional, **not** keyword):

```julia
kf(T, P, t, n, lna, n_initial) -> Real   # [mol/s]
```

where:
  - `T` [K], `P` [Pa]: temperature and pressure (plain `Real` or `ForwardDiff.Dual`).
  - `t` [s]: current time.
  - `n::StateView`: moles of all species (named access: `n["C3S"]`).
  - `lna::StateView`: log-activities of all species.
  - `n_initial::StateView`: initial moles (always `Float64`).
  - return: net dissolution rate [mol/s], positive = dissolution.

AD-compatible when the compiled closure is AD-compatible.

# Examples

```jldoctest
julia> idx = Dict("C3S" => 1);

julia> n_sv  = StateView([1.0], idx);

julia> lna_sv = StateView([0.0], idx);

julia> pk = parrot_killoh(PK_PARAMS_C3S, "C3S");

julia> pk(293.15, 1e5, 0.0, n_sv, lna_sv, n_sv) > 0
true
```
"""
struct KineticFunc{F, R <: NamedTuple, Q} <: Function
    compiled::F
    vars::NTuple{6, Symbol}   # (T, P, t, n, lna, n_initial) — positional arg names
    refs::R                   # default variable values as Quantity (for documentation/display)
    unit::Q                   # output unit, always u"mol/s"
end

const _KF_VARS = (:T, :P, :t, :n, :lna, :n_initial)
const _KF_DEFAULT_REFS = (T = 298.15u"K", P = 1.0e5u"Pa")

# Convenience constructor — vars defaults to the standard 6-argument names.
KineticFunc(compiled, refs::NamedTuple, unit) = KineticFunc(compiled, _KF_VARS, refs, unit)

# Positional call — hot path for ODE integration (no allocation, no unit handling).
(kf::KineticFunc)(T, P, t, n, lna, n_initial) = kf.compiled(T, P, t, n, lna, n_initial)

# Keyword call — user convenience (REPL, scripts), mirroring NumericFunc/SymbolicFunc.
# T, P, t are ustripped (accept Quantity or plain Real); n, lna, n_initial pass through.
@inline function (kf::KineticFunc)(; kwargs...)
    T_raw = haskey(kwargs, :T) ? kwargs[:T] :
        get(kf.refs, :T, get(_KF_DEFAULT_REFS, :T, nothing))
    P_raw = haskey(kwargs, :P) ? kwargs[:P] :
        get(kf.refs, :P, get(_KF_DEFAULT_REFS, :P, nothing))
    t_raw = haskey(kwargs, :t) ? kwargs[:t] : get(kf.refs, :t, 0.0)
    n_val = get(kwargs, :n, nothing)
    lna_val = get(kwargs, :lna, nothing)
    n0_val = get(kwargs, :n_initial, nothing)
    val = kf.compiled(ustrip(T_raw), ustrip(P_raw), ustrip(t_raw), n_val, lna_val, n0_val)
    return get(kwargs, :unit, false) ? val * kf.unit : val
end

function Base.show(io::IO, kf::KineticFunc)
    print(io, "KineticFunc [", dimension(kf.unit), "]")
    if !isempty(kf.vars)
        print(io, " ◆ vars=(", join(kf.vars, ", "), ")")
    end
    return if !isempty(kf.refs)
        print(io, " ◆ ", join(["$k=$v" for (k, v) in pairs(kf.refs)], ", "))
    end
end

function Base.show(io::IO, ::MIME"text/plain", kf::KineticFunc)
    println(io, "KineticFunc:")
    print(io, "  Unit: [", dimension(kf.unit), "]")
    print(io, "\n  Variables: ", join(kf.vars, ", "))
    return if !isempty(kf.refs)
        print(io, "\n  References: ", join(["$k=$v" for (k, v) in pairs(kf.refs)], ", "))
    end
end

# ── KINETICS_RATE_MODELS / KINETICS_RATE_FACTORIES ────────────────────────────

"""
    KINETICS_RATE_MODELS

Dictionary of raw kinetic rate-constant model expressions, analogous to
[`THERMO_MODELS`](@ref).

Each entry maps a model name (`:arrhenius`, …) to a `Dict` containing:
  - `:k` — symbolic `Expr` for the rate constant as a function of variables.
  - `:vars` — list of variable symbols (e.g. `[:T]`).
  - `:units` — list of `Symbol => Quantity` pairs for parameters and variables.
  - `:output_unit` — `Quantity` representing the output unit.

At package initialization, every entry is compiled into a `ThermoFactory`
stored in [`KINETICS_RATE_FACTORIES`](@ref).

# Example

```julia
k_acid = KINETICS_RATE_FACTORIES[:arrhenius](;
    k₀    = 5.012e-1,   # mol/(m² s) at T_ref
    Ea    = 14400.0,    # J/mol
    T_ref = 298.15,     # K
)
k_acid(; T = 310.0)   # → Float64 rate constant
```
"""
const KINETICS_RATE_MODELS = Dict{Symbol, Dict}(
    :arrhenius => Dict(
        :k => :(k₀ * exp(-Ea / R_gas * (1 / T - 1 / T_ref))),
        :vars => [:T],
        :units => [
            :T => u"K",
            :T_ref => u"K",
            :k₀ => u"mol/(m^2*s)",
            :Ea => u"J/mol",
            :R_gas => u"J/(mol*K)",
        ],
        :output_unit => u"mol/(m^2*s)",
    ),
)

"""
    KINETICS_RATE_FACTORIES

Compiled `ThermoFactory` objects for each kinetic rate model.
Populated by `__init__()` from [`KINETICS_RATE_MODELS`](@ref).

Keys are model name symbols (e.g. `:arrhenius`).
Values are `ThermoFactory` callables that return `SymbolicFunc{1}` instances.

# Usage

```julia
factory = KINETICS_RATE_FACTORIES[:arrhenius]
k = factory(; k₀=1e-5, Ea=50000.0, T_ref=298.15, R_gas=8.31446)
k(; T = 298.15)   # → 1e-5  (rate constant at reference temperature)
```
"""
const KINETICS_RATE_FACTORIES = Dict{Symbol, ThermoFactory}()

"""
    add_kinetics_rate_model(name::Symbol, dict_model::Dict)

Register a new kinetic rate-constant model in [`KINETICS_RATE_MODELS`](@ref)
and compile it into [`KINETICS_RATE_FACTORIES`](@ref).

`dict_model` must contain at minimum `:k` (expression), `:vars` (variable list),
`:units` (parameter units), and `:output_unit`.

# Example

```julia
add_kinetics_rate_model(:power_law, Dict(
    :k    => :(k₀ * (T / T_ref)^n),
    :vars => [:T],
    :units => [:T => u"K", :T_ref => u"K", :k₀ => u"mol/(m^2*s)", :n => u"1"],
    :output_unit => u"mol/(m^2*s)",
))
```
"""
function add_kinetics_rate_model(name::Symbol, dict_model::Dict)
    KINETICS_RATE_MODELS[name] = dict_model
    KINETICS_RATE_FACTORIES[name] = _build_kinetics_rate_factory(dict_model)
    return nothing
end

# Internal: compile one KINETICS_RATE_MODELS entry → ThermoFactory
function _build_kinetics_rate_factory(d::Dict)
    return ThermoFactory(
        d[:k],
        get(d, :vars, [:T]);
        units = get(d, :units, nothing),
        output_unit = get(d, :output_unit, u"1"),
    )
end

# ── arrhenius_rate_constant ────────────────────────────────────────────────────

"""
    arrhenius_rate_constant(k₀, Ea; T_ref=298.15, R_gas=8.31446261815324) -> NumericFunc

Build a temperature-dependent Arrhenius rate constant as a [`NumericFunc`](@ref):

```
k(T) = k₀ × exp(-Eₐ / R × (1/T - 1/T_ref))
```

The returned object is callable as `k(; T=...)` and fully AD-compatible
(ForwardDiff-safe: the closure captures `k₀`, `Ea`, `T_ref`, `R_gas` directly,
so dual numbers propagate correctly through all parameters).

Arithmetic between `SymbolicFunc`/`NumericFunc` objects is supported, so rate
constants can be composed with activity or surface-area functions.

# Arguments

  - `k₀`: pre-exponential factor at `T_ref`. Plain `Real` → SI [mol/(m² s)];
    `Quantity` → automatically converted (e.g. `5e-4u"mol/(m^2*s)"`).
  - `Ea`: activation energy. Plain `Real` → SI [J/mol]; `Quantity` → converted
    (e.g. `62.0u"kJ/mol"`).
  - `T_ref`: reference temperature. Plain `Real` → SI [K]; `Quantity` → converted
    (e.g. `298.15u"K"`). Default `298.15`.
  - `R_gas`: gas constant [J/(mol K)] (plain `Real` only; default `8.31446261815324`).

# Returns

A `NumericFunc` with variable `T` (in K) and `refs = (T = T_ref * u"K",)`.

# Examples

```jldoctest
julia> k = arrhenius_rate_constant(5.0e-4, 62000.0);

julia> isapprox(k(; T = 298.15), 5.0e-4; rtol = 1e-10)
true

julia> k(; T = 350.0) > k(; T = 298.15)   # higher T → higher k
true
```

Unit-aware: `k₀` in mmol/(m²·s), `Ea` in kJ/mol, `T_ref` in K — all converted to SI:
```julia
k = arrhenius_rate_constant(0.5u"mmol/(m^2*s)", 62.0u"kJ/mol"; T_ref = 298.15u"K")
```

AD-compatible through all parameters:
```julia
ForwardDiff.derivative(T  -> arrhenius_rate_constant(5e-4, 62000.0)(; T = T),  298.15)
ForwardDiff.derivative(Ea -> arrhenius_rate_constant(5e-4, Ea)(; T = 350.0),   62000.0)
ForwardDiff.derivative(k₀ -> arrhenius_rate_constant(k₀,   62000.0)(; T = 298.15), 5e-4)
```
"""
function arrhenius_rate_constant(
        k₀,
        Ea;
        T_ref = 298.15,
        R_gas::Real = 8.31446261815324,
    )
    k₀_si = safe_ustrip(us"mol/(m^2*s)", k₀)
    Ea_si = safe_ustrip(us"J/mol", Ea)
    T_ref_si = safe_ustrip(us"K", T_ref)
    # Closure captures SI values; no Float64 cast → ForwardDiff.Dual propagates correctly
    # through k₀, Ea, or T_ref when differentiating through construction.
    f = (T) -> k₀_si * exp(-Ea_si / R_gas * (1 / T - 1 / T_ref_si))
    # refs is metadata for default call values — always stored as plain Float64
    refs = (T = Float64(_primal(T_ref_si)) * u"K",)
    return NumericFunc(f, (:T,), refs, u"mol/(m^2*s)")
end

# ── Saturation ratio ───────────────────────────────────────────────────────────

"""
    saturation_ratio(stoich::AbstractVector, lna::AbstractVector,
                     ΔₐG⁰overT::AbstractVector; ϵ=1e-16) -> Real

Compute the saturation ratio Ω = IAP / K for a kinetic reaction.

```
ln Ω = Σᵢ νᵢ ln aᵢ + ln K
     = Σᵢ νᵢ ln aᵢ + ΔᵣG⁰/(RT)   (note: ΔᵣG⁰/RT = Σᵢ νᵢ ΔₐG⁰ᵢ/RT for reactants→products)
```

where `stoich[i]` is the stoichiometric coefficient (positive for products,
negative for reactants), `lna[i]` is the log-activity of species `i`,
and `ΔₐG⁰overT[i]` is the dimensionless standard Gibbs energy of formation
`ΔₐG⁰ᵢ / RT` for species `i`.

# Arguments

  - `stoich`: stoichiometric coefficient vector for this reaction (length = number of species).
  - `lna`: log-activity vector (same indexing as species in system).
  - `ΔₐG⁰overT`: dimensionless standard Gibbs energies `ΔₐG⁰ᵢ/RT`.
  - `ϵ`: floor to avoid `exp` overflow when Ω → ∞.

# Returns

`Ω = exp(ln_IAP - ln_K)` where `ln_K = -ΔᵣG⁰/RT`.

AD-compatible (ForwardDiff-safe).
"""
function saturation_ratio(
        stoich::AbstractVector,
        lna::AbstractVector,
        ΔₐG⁰overT::AbstractVector;
        ϵ::Real = 1.0e-16,
    )
    # ln IAP = Σᵢ νᵢ ln aᵢ
    ln_iap = sum(stoich[i] * lna[i] for i in eachindex(stoich))
    # ln K = -ΔᵣG⁰/RT = -Σᵢ νᵢ ΔₐG⁰ᵢ/RT
    ln_K = -sum(stoich[i] * ΔₐG⁰overT[i] for i in eachindex(stoich))
    return exp(ln_iap - ln_K)
end

# ── RateModelCatalyst ─────────────────────────────────────────────────────────

"""
    struct RateModelCatalyst{T<:Real}

Describes the contribution of a catalyst species to a reaction mechanism rate.

The catalyst multiplies the base rate by `exp(n * ln aᵢ) = aᵢ^n`, where
`aᵢ` is the activity of the catalyst species.

# Fields

  - `species`: PHREEQC-format formula string of the catalyst species (e.g. `"H+"`, `"OH-"`).
  - `n`: power exponent (dimensionless).

# Examples

```julia
acid_catalyst   = RateModelCatalyst("H+",  0.5)    # ∝ a(H+)^0.5
base_catalyst   = RateModelCatalyst("OH-", 0.5)    # ∝ a(OH-)^0.5
co2_catalyst    = RateModelCatalyst("CO2", 1.0)    # ∝ a(CO2)
```
"""
struct RateModelCatalyst{T <: Real}
    species::String
    n::T
end

# ── RateMechanism ─────────────────────────────────────────────────────────────

"""
    struct RateMechanism{F<:AbstractFunc, T<:Real}

A single kinetic mechanism (acid/neutral/base/…) contributing to the overall
mineral dissolution or precipitation rate.

The mechanism rate is:
```
r_mech = k(T) × [Π_catalysts aᵢ^nᵢ] × sign(1 - Ω) × |1 - Ω^p|^q
```

# Fields

  - `k`: rate constant as `AbstractFunc` (typically `SymbolicFunc{1}` from
    [`arrhenius_rate_constant`](@ref)). Called as `k(; T=...)`.
  - `p`: saturation exponent `p` in `(1 - Ω^p)^q`. Default 1.0.
  - `q`: outer exponent `q`. Default 1.0.
  - `catalysts`: vector of [`RateModelCatalyst`](@ref) (may be empty).

# Examples

```julia
k_acid = arrhenius_rate_constant(5.012e-1, 14400.0)
mech   = RateMechanism(k_acid, 1.0, 1.0, [RateModelCatalyst("H+", 1.0)])
```
"""
struct RateMechanism{F <: AbstractFunc, T <: Real}
    k::F
    p::T
    q::T
    catalysts::Vector{RateModelCatalyst{T}}
end

"""
    RateMechanism(k::AbstractFunc, p::Real, q::Real) -> RateMechanism

Construct a [`RateMechanism`](@ref) with no catalyst contributions.
"""
function RateMechanism(k::AbstractFunc, p::Real, q::Real)
    T = typeof(promote(p, q)[1])
    return RateMechanism{typeof(k), T}(k, T(p), T(q), RateModelCatalyst{T}[])
end

# ── parrot_killoh factory ──────────────────────────────────────────────────────

"""
    parrot_killoh(params::NamedTuple, mineral_name::AbstractString; α_max=1.0) -> KineticFunc

Build the Parrot & Killoh (1984) cement clinker hydration rate as a
[`KineticFunc`](@ref).

`params` must be a `NamedTuple` with keys `K₁`, `N₁`, `K₂`, `N₂`, `K₃`, `N₃`,
`B`, `Ea`, `T_ref`. All dimensional values accept plain `Real` (SI) or
`DynamicQuantities.Quantity`.

`mineral_name` is the PHREEQC formula string (e.g. `"C3S"`) used to look up the
mineral moles in the `n` and `n_initial` [`StateView`](@ref)s.

Three competing mechanisms determine the rate (Parrot & Killoh 1984):

| Mechanism | Formula |
|-----------|---------|
| Nucleation–growth | `r_NG = (K₁/N₁)(1-ξ)^N₁ / (1 + B·ξ^N₃)` |
| Interaction | `r_I = K₂(1-ξ)^N₂` |
| Diffusion | `r_D = 3K₃(1-ξ)^(2/3) / (N₃·(1-(1-ξ)^(1/3)))` |

The rate [mol/s] is `n_initial × Aₜ × min(max(r_NG, r_I), r_D)` where
`ξ = α / α_max` is the normalized degree of hydration and
`Aₜ = exp(-Ea/R × (1/T - 1/T_ref))` is the Arrhenius factor.

`α_max` can be set to apply the Powers (1948) water/cement ratio limit:
`α_max = min(1.0, w_c / 0.42)`.

# Returns

A [`KineticFunc`](@ref) — callable as
`pk(T, P, t, n::StateView, lna::StateView, n_initial::StateView) -> Real [mol/s]`.
AD-compatible (ForwardDiff-safe): no `Float64` casts in the evaluation path.

# Examples

```jldoctest
julia> pk = parrot_killoh(PK_PARAMS_C3S, "C3S");

julia> idx = Dict("C3S" => 1);

julia> n0  = StateView([1.0], idx);

julia> lna = StateView([0.0], idx);

julia> pk(293.15, 1e5, 0.0, n0, lna, n0) > 0
true
```

See also: [`PK_PARAMS_C3S`](@ref), [`PK_PARAMS_C2S`](@ref),
[`PK_PARAMS_C3A`](@ref), [`PK_PARAMS_C4AF`](@ref).
"""
function parrot_killoh(params::NamedTuple, mineral_name::AbstractString; α_max::Real = 1.0)
    K₁ = safe_ustrip(us"1/s", params.K₁)
    N₁ = float(params.N₁)
    K₂ = safe_ustrip(us"1/s", params.K₂)
    N₂ = float(params.N₂)
    K₃ = safe_ustrip(us"1/s", params.K₃)
    N₃ = float(params.N₃)
    B = float(params.B)
    Ea = safe_ustrip(us"J/mol", params.Ea)
    T_ref = safe_ustrip(us"K", params.T_ref)
    α_max_f = float(α_max)
    R_gas = 8.31446261815324

    f = (T, _P, _t, n, _lna, n_initial) -> begin
        n_m = n[mineral_name]
        n_init = max(n_initial[mineral_name], oneunit(n_m) * 1.0e-30)
        # degree of hydration α ∈ [0, α_max)
        α = min(max(one(T) - n_m / n_init, zero(T)), α_max_f - oftype(T, 1.0e-10))
        ξ = α / α_max_f
        # Arrhenius temperature correction
        Aₜ = exp(-Ea / R_gas * (one(T) / T - one(T) / T_ref))
        one_m_ξ = one(ξ) - ξ
        # r_NG: nucleation–growth [s⁻¹]
        r_NG = (K₁ / N₁) * one_m_ξ^N₁ / (one(ξ) + B * ξ^N₃)
        # r_I: interaction [s⁻¹]
        r_I = K₂ * one_m_ξ^N₂
        # r_D: diffusion [s⁻¹] (denominator clamped to avoid 0/0 at α=0)
        denom_D = max(one(ξ) - one_m_ξ^(one(ξ) / 3), oftype(ξ, 1.0e-10))
        r_D = 3 * K₃ * one_m_ξ^(2 * one(ξ) / 3) / (N₃ * denom_D)
        return n_init * Aₜ * min(max(r_NG, r_I), r_D)
    end

    refs = (T = Float64(_primal(T_ref)) * u"K", P = 1.0e5u"Pa")
    return KineticFunc(f, refs, u"mol/s")
end

# ── Predefined Parrot & Killoh (1984) parameters ─────────────────────────────

"""
    PK_PARAMS_C3S :: NamedTuple

Parrot & Killoh (1984) parameters for alite (C₃S = Ca₃SiO₅).

Original paper values (K₁=1.5, K₂=0.018, K₃=0.0024 d⁻¹).
Activation energy from Schindler & Folliard (2005).
Reference temperature: 293.15 K (20 °C).

Pass to [`parrot_killoh`](@ref) to build a [`KineticFunc`](@ref):

```julia
pk = parrot_killoh(PK_PARAMS_C3S, "C3S")
# or with α_max limit (Powers 1948):
pk = parrot_killoh(PK_PARAMS_C3S, "C3S"; α_max = min(1.0, w_c / 0.42))
```
"""
const PK_PARAMS_C3S = (
    K₁ = 1.5u"1/d",
    N₁ = 3.3,
    K₂ = 0.018u"1/d",
    N₂ = 2.5,
    K₃ = 0.0024u"1/d",
    N₃ = 4.0,
    B = 0.5,
    Ea = 41_570.0u"J/mol",
    T_ref = 293.15u"K",
)

"""
    PK_PARAMS_C2S :: NamedTuple

Parrot & Killoh (1984) parameters for belite (C₂S = Ca₂SiO₄).

Original paper values (K₁=0.95, K₂=0.0005, K₃=0.0024 d⁻¹).
Activation energy from Schindler & Folliard (2005).
Reference temperature: 293.15 K (20 °C).
"""
const PK_PARAMS_C2S = (
    K₁ = 0.95u"1/d",
    N₁ = 0.5,
    K₂ = 0.0005u"1/d",
    N₂ = 2.5,
    K₃ = 0.0024u"1/d",
    N₃ = 4.0,
    B = 0.2,
    Ea = 43_670.0u"J/mol",
    T_ref = 293.15u"K",
)

"""
    PK_PARAMS_C3A :: NamedTuple

Parrot & Killoh (1984) parameters for tricalcium aluminate (C₃A = Ca₃Al₂O₆)
in the presence of sulfate (gypsum), corresponding to ettringite formation.

Original paper values (K₁=0.082, K₂=0.00024, K₃=0.0024 d⁻¹).
Activation energy from Schindler & Folliard (2005).
Reference temperature: 293.15 K (20 °C).
"""
const PK_PARAMS_C3A = (
    K₁ = 0.082u"1/d",
    N₁ = 0.87,
    K₂ = 0.00024u"1/d",
    N₂ = 2.0,
    K₃ = 0.0024u"1/d",
    N₃ = 4.0,
    B = 0.04,
    Ea = 54_040.0u"J/mol",
    T_ref = 293.15u"K",
)

"""
    PK_PARAMS_C4AF :: NamedTuple

Parrot & Killoh (1984) parameters for tetracalcium aluminoferrite
(C₄AF = Ca₄Al₂Fe₂O₁₀).

Original paper values (K₁=0.165, K₂=0.0015, K₃=0.0024 d⁻¹).
Activation energy from Schindler & Folliard (2005).
Reference temperature: 293.15 K (20 °C).
"""
const PK_PARAMS_C4AF = (
    K₁ = 0.165u"1/d",
    N₁ = 3.7,
    K₂ = 0.0015u"1/d",
    N₂ = 2.5,
    K₃ = 0.0024u"1/d",
    N₃ = 4.0,
    B = 0.5,
    Ea = 34_420.0u"J/mol",
    T_ref = 293.15u"K",
)

# ── parrot_killoh_avrami — the canonical 1984 formulation ────────────────────

"""
    PK_AVRAMI_SEED

Lower bound imposed on the normalized degree of hydration inside the Avrami
branch of [`parrot_killoh_avrami`](@ref), so that the rate is strictly positive
at `α = 0` and the ODE leaves its degenerate initial point. See the discussion
in [`parrot_killoh_avrami`](@ref).

One visible consequence: a phase whose rate is governed by the Avrami branch
near `α = 0` (C₃S, C₃A, C₄AF) starts more slowly than one governed by the power
law (C₂S), so alite only overtakes belite after a few minutes. The crossover
falls inside the induction period, which this model does not describe anyway.
"""
const PK_AVRAMI_SEED = 1.0e-6

"""
    parrot_killoh_avrami(params::NamedTuple, mineral_name::AbstractString;
                         α_max = 1.0, blaine = nothing, humidity = nothing) -> KineticFunc

Build the Parrot & Killoh (1984) clinker hydration rate in its **canonical
formulation**, as reported by Lothenbach et al. (2008) and used by Lavergne
et al. (2018).

`params` must be a `NamedTuple` with keys `k₁`, `n₁`, `k₂`, `k₃`, `n₃`, `Ea`,
`T_ref` — see [`PK84_PARAMS_C3S`](@ref) and siblings. Dimensional values accept
plain `Real` (SI) or `DynamicQuantities.Quantity`.

Three competing mechanisms limit the rate, and the **slowest one wins**:

| Mechanism | Formula |
|-----------|---------|
| Nucleation–growth (Avrami) | `α̇₁ = (k₁/n₁)(1-ξ)(-ln(1-ξ))^(1-n₁)` |
| Diffusion (Jander) | `α̇₂ = k₂(1-ξ)^(2/3) / (1-(1-ξ)^(1/3))` |
| Shell formation (power law) | `α̇₃ = k₃(1-ξ)^n₃` |

so that `α̇ = min(α̇₁, α̇₂, α̇₃)`, with `ξ = α/α_max` the normalized degree of
hydration. The returned rate [mol/s] is `n_initial × Aₜ × β_B × β_h × α̇`, where
`Aₜ = exp(-Ea/R × (1/T - 1/T_ref))` is the Arrhenius factor, `β_B` the Blaine
fineness factor ([`blaine_factor`](@ref)) and `β_h` the relative-humidity
reduction ([`humidity_factor`](@ref)). Both default to 1 when their keyword is
`nothing`.

!!! note "Two Parrot–Killoh variants ship with ChemistryLab"
    [`parrot_killoh`](@ref) implements a *different*, smoothed variant
    (`min(max(r_NG, r_I), r_D)` with a `B`-damped nucleation term) together with
    the parameter set of [`PK_PARAMS_C3S`](@ref) and siblings. The two are not
    interchangeable: their parameters are **not** transferable, and only
    `parrot_killoh_avrami` with [`PK84_PARAMS_C3S`](@ref) reproduces the α(t)
    curves published in the cement literature cited above.

With the canonical parameters, C₂S has no nucleation–growth stage and C₃S has no
diffusion-controlled stage — an artifact of the 1984 fit that the original
authors acknowledged, and a convenient signature to check an implementation
against.

# Keyword arguments

  - `α_max`: Powers (1948) water availability cap — see [`powers_alpha_max`](@ref).
  - `blaine`: Blaine fineness of the binder, as a `Quantity` or a plain `Real`
    in m²/kg. `nothing` (default) means no correction.
  - `humidity`: internal relative humidity, either a constant in `[0, 1]` or a
    callable `t -> h(t)`. `nothing` (default) means no correction.

# Returns

A [`KineticFunc`](@ref) — callable as
`pk(T, P, t, n::StateView, lna::StateView, n_initial::StateView) -> Real [mol/s]`.
AD-compatible (ForwardDiff-safe): no `Float64` casts in the evaluation path.

# Examples

```jldoctest
julia> pk = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S"; blaine = 380u"m^2/kg");

julia> idx = Dict("C3S" => 1);

julia> n0 = StateView([1.0], idx);

julia> lna = StateView([0.0], idx);

julia> pk(293.15, 1e5, 3600.0, StateView([0.9], idx), lna, n0) > 0
true
```

See also: [`PK84_PARAMS_C3S`](@ref), [`waller`](@ref), [`blaine_factor`](@ref),
[`humidity_factor`](@ref), [`powers_alpha_max`](@ref).
"""
function parrot_killoh_avrami(
        params::NamedTuple, mineral_name::AbstractString;
        α_max::Real = 1.0, blaine = nothing, humidity = nothing
    )
    k₁ = safe_ustrip(us"1/s", params.k₁)
    n₁ = float(params.n₁)
    k₂ = safe_ustrip(us"1/s", params.k₂)
    k₃ = safe_ustrip(us"1/s", params.k₃)
    n₃ = float(params.n₃)
    Ea = safe_ustrip(us"J/mol", params.Ea)
    T_ref = safe_ustrip(us"K", params.T_ref)
    α_max_f = float(α_max)
    β_B = blaine === nothing ? 1.0 : blaine_factor(blaine)
    R_gas = 8.31446261815324

    f = (T, _P, t, n, _lna, n_initial) -> begin
        n_m = n[mineral_name]
        n_init = max(n_initial[mineral_name], oneunit(n_m) * 1.0e-30)
        α = min(max(one(T) - n_m / n_init, zero(T)), α_max_f - oftype(T, 1.0e-10))
        ξ = α / α_max_f
        Aₜ = exp(-Ea / R_gas * (one(T) / T - one(T) / T_ref))
        β_h = humidity === nothing ? one(ξ) : humidity_factor(_humidity_at(humidity, t))
        one_m_ξ = max(one(ξ) - ξ, oftype(ξ, 1.0e-12))
        # α̇₁ — Avrami nucleation and growth. For n₁ < 1 the (-ln(1-ξ))^(1-n₁)
        # factor VANISHES at ξ = 0, so α̇ = 0 and α ≡ 0 solves the ODE: hydration
        # would never start. Parrot & Killoh's own discrete scheme escapes this by
        # evaluating the integrated Avrami law over the first finite time step;
        # a continuous solver cannot, so the argument is floored at ξ_seed. The
        # seed only sets how fast the solution leaves the degenerate point — by
        # ξ ≈ 1e-3 the Avrami branch is already three orders of magnitude above it.
        ξ_avrami = max(ξ, oftype(ξ, PK_AVRAMI_SEED))
        r₁ = (k₁ / n₁) * one_m_ξ * (-log(one(ξ) - ξ_avrami))^(one(ξ) - n₁)
        # α̇₂ — Jander diffusion through the hydrate layer. Denominator → 0 as ξ → 0.
        denom = max(one(ξ) - one_m_ξ^(one(ξ) / 3), oftype(ξ, 1.0e-12))
        r₂ = k₂ * one_m_ξ^(2 * one(ξ) / 3) / denom
        # α̇₃ — power law, thick shell around the unreacted grain.
        r₃ = k₃ * one_m_ξ^n₃
        return n_init * Aₜ * β_B * β_h * min(r₁, r₂, r₃)
    end

    refs = (T = Float64(_primal(T_ref)) * u"K", P = 1.0e5u"Pa")
    return KineticFunc(f, refs, u"mol/s")
end

# Internal: a humidity keyword is either a constant or a function of time.
@inline _humidity_at(h::Real, _t) = h
@inline _humidity_at(h, t) = h(t)

# ── Canonical Parrot & Killoh (1984) parameters ──────────────────────────────
#
# Table 3 of Lavergne et al. (2018), themselves quoting Parrot & Killoh (1984)
# as reported by Lothenbach et al. (2008); activation energies from Table 4
# (Maekawa et al.), which differ markedly from the uniform values often assumed:
# the minerals that hydrate later have the *lower* apparent activation energy,
# E_C3A > E_C3S > E_C4AF > E_C2S.

"""
    PK84_PARAMS_C3S :: NamedTuple

Canonical Parrot & Killoh (1984) parameters for alite (C₃S = Ca₃SiO₅), valid for
a Blaine fineness of 385 m²/kg and a reference temperature of 20 °C.

`k₁ = 1.5 d⁻¹`, `n₁ = 0.7`, `k₂ = 0.05 d⁻¹`, `k₃ = 1.1 d⁻¹`, `n₃ = 3.3`,
`Ea = 42 kJ/mol`.

Pass to [`parrot_killoh_avrami`](@ref), **not** to [`parrot_killoh`](@ref) —
the two use different functional forms and their parameters are not transferable.
"""
const PK84_PARAMS_C3S = (
    k₁ = 1.5u"1/d", n₁ = 0.7, k₂ = 0.05u"1/d", k₃ = 1.1u"1/d", n₃ = 3.3,
    Ea = 42_000.0u"J/mol", T_ref = 293.15u"K",
)

"""
    PK84_PARAMS_C2S :: NamedTuple

Canonical Parrot & Killoh (1984) parameters for belite (C₂S = Ca₂SiO₄).

`k₁ = 0.5 d⁻¹`, `n₁ = 1.0`, `k₂ = 0.006 d⁻¹`, `k₃ = 0.2 d⁻¹`, `n₃ = 5.0`,
`Ea = 21 kJ/mol`.

With `n₁ = 1` the Avrami branch reduces to `k₁(1-ξ)`, which never limits the
rate: belite hydration is governed by the power law throughout.

See [`PK84_PARAMS_C3S`](@ref).
"""
const PK84_PARAMS_C2S = (
    k₁ = 0.5u"1/d", n₁ = 1.0, k₂ = 0.006u"1/d", k₃ = 0.2u"1/d", n₃ = 5.0,
    Ea = 21_000.0u"J/mol", T_ref = 293.15u"K",
)

"""
    PK84_PARAMS_C3A :: NamedTuple

Canonical Parrot & Killoh (1984) parameters for tricalcium aluminate
(C₃A = Ca₃Al₂O₆).

`k₁ = 1.0 d⁻¹`, `n₁ = 0.85`, `k₂ = 0.04 d⁻¹`, `k₃ = 1.0 d⁻¹`, `n₃ = 3.2`,
`Ea = 54 kJ/mol`.

See [`PK84_PARAMS_C3S`](@ref).
"""
const PK84_PARAMS_C3A = (
    k₁ = 1.0u"1/d", n₁ = 0.85, k₂ = 0.04u"1/d", k₃ = 1.0u"1/d", n₃ = 3.2,
    Ea = 54_000.0u"J/mol", T_ref = 293.15u"K",
)

"""
    PK84_PARAMS_C4AF :: NamedTuple

Canonical Parrot & Killoh (1984) parameters for tetracalcium aluminoferrite
(C₄AF = Ca₄Al₂Fe₂O₁₀).

`k₁ = 0.37 d⁻¹`, `n₁ = 0.7`, `k₂ = 0.015 d⁻¹`, `k₃ = 0.4 d⁻¹`, `n₃ = 3.7`,
`Ea = 32 kJ/mol`.

See [`PK84_PARAMS_C3S`](@ref).
"""
const PK84_PARAMS_C4AF = (
    k₁ = 0.37u"1/d", n₁ = 0.7, k₂ = 0.015u"1/d", k₃ = 0.4u"1/d", n₃ = 3.7,
    Ea = 32_000.0u"J/mol", T_ref = 293.15u"K",
)

# ── waller — pozzolanic and latent-hydraulic additions ───────────────────────

"""
    waller(params::NamedTuple, mineral_name::AbstractString;
           α_max = 1.0, blaine = nothing, humidity = nothing) -> KineticFunc

Build the Waller (1999) reaction rate of a pozzolanic or latent-hydraulic
addition (fly ash, silica fume, ground granulated slag) as a [`KineticFunc`](@ref).

The degree of reaction follows a sigmoid in log-time,

```math
α(t) = \\frac{1}{1 + (τ/t)^n}
```

whose rate, written as a function of the current degree so that it composes with
temperature, fineness and humidity corrections, is

```math
α̇ = \\frac{n}{τ} (1 - α)^{1 + 1/n} α^{1 - 1/n}.
```

`params` must be a `NamedTuple` with keys `τ`, `n`, `Ea`, `T_ref` and, optionally,
`blaine_ref` — see [`WALLER_PARAMS_FLY_ASH`](@ref).

Pozzolanic reactions are markedly more temperature-sensitive than the hydraulic
reactions of clinker: the shipped activation energy is 83.14 kJ/mol against
21–54 kJ/mol for the clinker phases.

Silica fume is far finer than the cement (about 20 000 m²/kg by BET). Its
reactivity is nonetheless represented here through the *Blaine* fineness, for
which an effective 2000 m²/kg is the recommended default
([`WALLER_PARAMS_SILICA_FUME`](@ref)) — the two measurements probe different
physical phenomena and are not interchangeable.

# Keyword arguments

Identical to [`parrot_killoh_avrami`](@ref). The Blaine correction is taken
relative to `params.blaine_ref` (400 m²/kg for fly ash), not to the clinker
reference of 385 m²/kg.

# Examples

```jldoctest
julia> fa = waller(WALLER_PARAMS_FLY_ASH, "FlyAsh");

julia> idx = Dict("FlyAsh" => 1);

julia> n0 = StateView([1.0], idx);

julia> lna = StateView([0.0], idx);

julia> fa(293.15, 1e5, 86400.0, StateView([0.9], idx), lna, n0) > 0
true
```

See also: [`WALLER_PARAMS_FLY_ASH`](@ref), [`WALLER_PARAMS_SILICA_FUME`](@ref),
[`WALLER_PARAMS_SLAG`](@ref), [`parrot_killoh_avrami`](@ref).
"""
function waller(
        params::NamedTuple, mineral_name::AbstractString;
        α_max::Real = 1.0, blaine = nothing, humidity = nothing
    )
    τ = safe_ustrip(us"s", params.τ)
    n_w = float(params.n)
    Ea = safe_ustrip(us"J/mol", params.Ea)
    T_ref = safe_ustrip(us"K", params.T_ref)
    α_max_f = float(α_max)
    blaine_ref = hasproperty(params, :blaine_ref) ? params.blaine_ref : 400.0u"m^2/kg"
    β_B = blaine === nothing ? 1.0 : blaine_factor(blaine; blaine_ref = blaine_ref)
    R_gas = 8.31446261815324

    f = (T, _P, t, n, _lna, n_initial) -> begin
        n_m = n[mineral_name]
        n_init = max(n_initial[mineral_name], oneunit(n_m) * 1.0e-30)
        α = min(max(one(T) - n_m / n_init, zero(T)), α_max_f - oftype(T, 1.0e-10))
        ξ = α / α_max_f
        Aₜ = exp(-Ea / R_gas * (one(T) / T - one(T) / T_ref))
        β_h = humidity === nothing ? one(ξ) : humidity_factor(_humidity_at(humidity, t))
        # At ξ = 0 the closed form α̇(α) is singular (α^(1-1/n) → ∞ for n < 1).
        # Fall back to the explicit α̇(t) of the sigmoid, which is finite for t > 0
        # and vanishes as t → 0 — the induction period the sigmoid encodes.
        if _primal(ξ) < 1.0e-12
            t_pos = max(t, oftype(ξ, 1.0e-6))
            x = (τ / t_pos)^n_w
            r = (n_w / t_pos) * x / (one(ξ) + x)^2
        else
            one_m_ξ = max(one(ξ) - ξ, oftype(ξ, 1.0e-12))
            r = (n_w / τ) * one_m_ξ^(one(ξ) + one(ξ) / n_w) * ξ^(one(ξ) - one(ξ) / n_w)
        end
        return n_init * Aₜ * β_B * β_h * r
    end

    refs = (T = Float64(_primal(T_ref)) * u"K", P = 1.0e5u"Pa")
    return KineticFunc(f, refs, u"mol/s")
end

"""
    WALLER_PARAMS_FLY_ASH :: NamedTuple

Waller (1999) parameters for class-F fly ash: `τ = 80 d`, `n = 0.7`,
`blaine_ref = 400 m²/kg`, `Ea = 83.14 kJ/mol`.

Adjusted on SEM image analysis assuming a fly ash of 60 % pozzolanic activity.

Pass to [`waller`](@ref).
"""
const WALLER_PARAMS_FLY_ASH = (
    τ = 80.0u"d", n = 0.7, blaine_ref = 400.0u"m^2/kg",
    Ea = 83_140.0u"J/mol", T_ref = 293.15u"K",
)

"""
    WALLER_PARAMS_SILICA_FUME :: NamedTuple

Waller (1999) parameters applied to silica fume — identical kinetics to
[`WALLER_PARAMS_FLY_ASH`](@ref), the higher reactivity being carried by the
fineness. Pass `blaine = 2000u"m^2/kg"` to [`waller`](@ref), the effective value
recommended by Lavergne et al. (2018); the BET surface of silica fume (about
20 000 m²/kg) is **not** a Blaine fineness and must not be used here.
"""
const WALLER_PARAMS_SILICA_FUME = (
    τ = 80.0u"d", n = 0.7, blaine_ref = 400.0u"m^2/kg",
    Ea = 83_140.0u"J/mol", T_ref = 293.15u"K",
)

"""
    WALLER_PARAMS_SLAG :: NamedTuple

Waller (1999) parameters for ground granulated blast-furnace slag:
`τ = 100 d`, `n = 0.7`, `blaine_ref = 400 m²/kg`, `Ea = 83.14 kJ/mol`.

Slag is latent-hydraulic rather than pozzolanic; the longer characteristic time
reflects its slower long-term reaction. Combine with an `α_max` below 1 (0.9 is
customary) to account for the unreactive crystalline fraction.
"""
const WALLER_PARAMS_SLAG = (
    τ = 100.0u"d", n = 0.7, blaine_ref = 400.0u"m^2/kg",
    Ea = 83_140.0u"J/mol", T_ref = 293.15u"K",
)

# ── Correction factors ───────────────────────────────────────────────────────

"""
    blaine_factor(blaine; blaine_ref = 385u"m^2/kg") -> Real

Fineness correction of the hydration rate: the Parrot & Killoh parameters were
adjusted for a cement of Blaine fineness `blaine_ref`, and the rate scales as
`blaine / blaine_ref`.

Both arguments accept a `DynamicQuantities.Quantity` or a plain `Real` in m²/kg.
The default reference is 385 m²/kg for clinker phases; pass
`blaine_ref = 400u"m^2/kg"` for the Waller kinetics of additions.

# Examples

```jldoctest
julia> blaine_factor(385u"m^2/kg")
1.0

julia> round(blaine_factor(462u"m^2/kg"); digits = 4)
1.2
```
"""
function blaine_factor(blaine; blaine_ref = 385.0u"m^2/kg")
    B = safe_ustrip(us"m^2/kg", blaine)
    B₀ = safe_ustrip(us"m^2/kg", blaine_ref)
    return B / B₀
end

"""
    humidity_factor(h) -> Real

Reduction coefficient `β_h` applied to the hydration rate at internal relative
humidity `h ∈ [0, 1]` (Parrot et al., as used by van Breugel):

```math
β_h = \\left(\\frac{h - 0.55}{0.45}\\right)^4 \\ \\text{if } h > 0.80,
\\qquad β_h = 0 \\ \\text{otherwise.}
```

Hydration is taken to stop below 80 % relative humidity, on thermodynamic
grounds. The cut is a genuine discontinuity: the one-sided limit from above is
`β_h(0.80⁺) ≈ 0.0953` while the value at and below 0.80 is exactly 0. The jump
is mild in practice, the rate having already fallen by an order of magnitude
from `β_h(0.99) ≈ 0.914`.

# Examples

```jldoctest
julia> humidity_factor(0.75)
0.0

julia> humidity_factor(0.80)
0.0

julia> round(humidity_factor(0.801); digits = 4)
0.0968
```
"""
function humidity_factor(h::Real)
    return _primal(h) > 0.8 ? ((h - oftype(h, 0.55)) / oftype(h, 0.45))^4 : zero(h)
end

"""
    powers_alpha_max(w_c) -> Real

Powers (1948) upper bound on the degree of hydration set by the availability of
water, `α_max = min(1, w/c / 0.42)`: hydrating one gram of cement binds about
0.42 g of water, so a paste below `w/c = 0.42` cannot hydrate completely.

Pass the result as the `α_max` keyword of [`parrot_killoh`](@ref),
[`parrot_killoh_avrami`](@ref) or [`waller`](@ref).

# Examples

```jldoctest
julia> powers_alpha_max(0.5)
1.0

julia> round(powers_alpha_max(0.32); digits = 4)
0.7619
```
"""
powers_alpha_max(w_c::Real) = min(one(w_c), w_c / oftype(w_c, 0.42))
