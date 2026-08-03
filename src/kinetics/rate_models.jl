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
