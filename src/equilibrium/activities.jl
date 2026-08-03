# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using LinearAlgebra

# ── Abstract activity model ───────────────────────────────────────────────────

"""
    abstract type AbstractActivityModel end

Base type for all activity models. Each concrete subtype must implement
`activity_model(cs::ChemicalSystem, model::AbstractActivityModel)`
which returns a closure `a(n, p) -> Vector{Float64}` of log-activities.
"""
abstract type AbstractActivityModel end

# ── Concrete models ───────────────────────────────────────────────────────────

"""
    struct DiluteSolutionModel <: AbstractActivityModel

Ideal dilute solution model:
- Solvent:  Raoult's law  — `ln a = ln(x_solvent)`
- Solutes:  Henry's law   — `ln a = ln(c_i / c°)` where `c° = 1 mol/L`
- Crystals: pure solid    — `ln a = 0`
- Gas:      ideal mixture — `ln a = ln(x_i)`
- SS end-members: ideal mixing — `ln a = ln(xᵢ)` within the solid-solution phase
"""
struct DiluteSolutionModel <: AbstractActivityModel end

# ── Activity model factory ────────────────────────────────────────────────────

"""
    activity_model(cs::ChemicalSystem, ::DiluteSolutionModel) -> Function

Return a closure `lna(n, p) -> Vector{Float64}` computing the vector of
log-activities for the dilute ideal solution model.

The returned function has signature `lna(n, p)` where:
- `n`: dimensionless mole vector (same indexing as `cs.species`)
- `p`: `NamedTuple` containing at least `ϵ` (floor value to avoid log(0))

Solid-solution end-members (class `SC_SSENDMEMBER`) receive `ln aᵢ = ln xᵢ`
where `xᵢ = nᵢ / Σnⱼ` within the same solid-solution phase.

All quantities are dimensionless — units are stripped at construction time.
"""
function activity_model(cs::ChemicalSystem, ::DiluteSolutionModel)

    has_aqueous = !isempty(cs.idx_solvent)
    idx_solvent = has_aqueous ? only(cs.idx_solvent) : 0
    idx_solutes = cs.idx_solutes
    idx_gas = cs.idx_gas

    ln_c_solvent = if has_aqueous
        M_solvent = ustrip(us"kg/mol", cs.species[idx_solvent][:M])
        log(1.0 / M_solvent)    # c° = ρ/M ≈ 1/M (ρ ≈ 1 kg/L)
    else
        0.0
    end

    ss_groups = cs.ss_groups
    has_ss = !isempty(ss_groups)
    has_gas = !isempty(idx_gas)
    ss_models = has_ss ? map(ss -> ss.model, cs.solid_solutions) : nothing

    function lna(n::AbstractVector, p)
        ϵ = p.ϵ
        _n = max.(n, ϵ)     # ϵ::Float64 — promotion vers Dual automatique si n est Dual

        out = zeros(eltype(_n), length(_n))

        if has_aqueous
            # n_aqueous ≥ ϵ > 0 always (because _n[i] ≥ ϵ), so no iszero guard needed
            n_aqueous = _n[idx_solvent] + sum((_n[i] for i in idx_solutes); init = zero(eltype(_n)))
            out[idx_solvent] = log(_n[idx_solvent] / n_aqueous)
            @inbounds for i in idx_solutes
                out[i] = log(_n[i] / _n[idx_solvent]) + ln_c_solvent
            end
        end

        if has_gas
            n_gas = sum((_n[i] for i in idx_gas); init = zero(eltype(_n)))
            @inbounds for i in idx_gas
                out[i] = log(_n[i] / n_gas)
            end
        end

        if has_ss
            T_val = hasproperty(p, :T) ? p.T : 298.15
            _solid_solution_lna!(out, _n, ss_groups, ss_models, T_val, ϵ)
        end

        return out
    end

    return lna
end

# ── Potential builder ─────────────────────────────────────────────────────────

"""
    build_potentials(cs::ChemicalSystem, model::AbstractActivityModel) -> Function

Return a closure `μ(n, p) -> Vector{Float64}` computing dimensionless chemical
potentials `μ_i / RT` for all species.

``\\mu_i / RT = \\Delta_a G_i^0 / RT + \\ln a_i``

The returned function is compatible with SciML solvers:
- `n`: dimensionless mole vector
- `p`: `NamedTuple` containing:
  - `ΔₐG⁰overT`: vector of standard Gibbs energies of formation divided by RT
  - `ϵ`: regularization floor (e.g. `1e-30`)

All quantities are dimensionless — caller is responsible for stripping units
from `ΔₐG⁰overT` before passing them in `p`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("Na+"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
       ]);

julia> μ = build_potentials(cs, DiluteSolutionModel());

julia> n = [55.5, 0.1];

julia> p = (ΔₐG⁰overT = [-95.6, -105.6], ϵ = 1e-30);

julia> length(μ(n, p)) == 2
true
```
"""
function build_potentials(cs::ChemicalSystem, model::AbstractActivityModel)

    # Build the activity closure once — captures precomputed indices and constants
    lna = activity_model(cs, model)

    function μ(n::AbstractVector, p)
        return p.ΔₐG⁰overT .+ lna(n, p)       # μ_i/RT = ΔₐG⁰_i/RT + ln(a_i)
    end

    return μ
end

# ── HKF / Debye-Hückel B-dot activity model ──────────────────────────────────

"""
    REJ_HKF::Dict{String,Float64}

Effective electrostatic radii åᵢ [Å] for aqueous ions from
Helgeson, Kirkham & Flowers (1981), *Am. J. Sci.* **281**, Table 3.

Keys are PHREEQC-format formula strings (e.g. `"Na+"`, `"Ca+2"`, `"SO4-2"`).
Used by [`HKFActivityModel`](@ref) with priority 2 in the radius lookup chain:
`sp[:å]` > `REJ_HKF` > [`REJ_CHARGE_DEFAULT`](@ref) > `model.å_default`.

See also: [`REJ_CHARGE_DEFAULT`](@ref), [`HKFActivityModel`](@ref).
"""
const REJ_HKF = Dict{String, Float64}(
    "H+" => 3.08, "Li+" => 1.64, "Na+" => 1.91, "K+" => 2.27,
    "Rb+" => 2.41, "Cs+" => 2.61, "NH4+" => 2.31, "Ag+" => 2.2,
    "Mg+2" => 2.54, "Ca+2" => 2.87, "Sr+2" => 3.0, "Ba+2" => 3.22,
    "Fe+2" => 2.62, "Al+3" => 3.33, "Fe+3" => 3.46, "La+3" => 3.96,
    "F-" => 1.33, "Cl-" => 1.81, "Br-" => 1.96, "I-" => 2.2,
    "OH-" => 1.4, "HS-" => 1.84, "NO3-" => 2.81, "HCO3-" => 2.1,
    "HSO4-" => 2.37, "SO4-2" => 3.15, "CO3-2" => 2.81,
)

"""
    REJ_CHARGE_DEFAULT::Dict{Int,Float64}

Fallback effective electrostatic radii åᵢ [Å] indexed by formal charge,
from ToughReact V2 (Xu et al. 2011, Table A2; after Helgeson et al. 1981).

Used by [`HKFActivityModel`](@ref) with priority 3 in the radius lookup chain:
`sp[:å]` > [`REJ_HKF`](@ref) > `REJ_CHARGE_DEFAULT` > `model.å_default`.

See also: [`REJ_HKF`](@ref), [`HKFActivityModel`](@ref).
"""
const REJ_CHARGE_DEFAULT = Dict{Int, Float64}(
    -3 => 4.2, -2 => 3.0, -1 => 1.81,
    1 => 2.31, 2 => 2.8, 3 => 3.6, 4 => 4.5,
)

# ── Internal helpers ──────────────────────────────────────────────────────────

"""
    _hkf_sigma(x) -> Real

Compute the σ function used in the osmotic coefficient formula
(Helgeson et al. 1981, Eq. 132–137):

```
σ(x) = (3/x³)(x − 2 ln(1+x) − 1/(1+x) + 1)
```

For `|x| < 1e-3`, a Taylor series `1 − (3/2)x + (9/5)x²` is used to avoid
catastrophic cancellation. Both branches agree to O(x³) at the threshold,
so the gradient is continuous.

AD-compatible: branching is on `ForwardDiff.value(x)`, not on the Dual itself.
"""
function _hkf_sigma(x::T) where {T <: Real}
    if abs(_primal(x)) < 1.0e-3
        return one(T) - (3 // 2) * x + (9 // 5) * x^2
    else
        return (3 / x^3) * (x - 2 * log1p(x) - one(T) / (one(T) + x) + one(T))
    end
end

"""
    hkf_debye_huckel_params(T_K, P_Pa) -> NamedTuple{(:A, :B)}

Compute the Debye-Hückel A and B parameters from the water density ρ [g/cm³]
and dielectric constant εᵣ at temperature `T_K` (K) and pressure `P_Pa` (Pa).

Formulas (Helgeson et al. 1981):
```
A(T,P) = 1.824829238×10⁶ × √ρ / (εᵣ T)^(3/2)    [(kg/mol)^(1/2)]
B(T,P) = 50.29158649      × √ρ / √(εᵣ T)          [Å⁻¹ (kg/mol)^(1/2)]
```
where ρ is in g/cm³.

Uses `water_thermo_props` (HGK equation of state) and
`water_electro_props_jn` (Johnson-Norton dielectric constant).

AD-compatible (ForwardDiff-safe). Returns a `NamedTuple` `(A=..., B=...)`.

# Examples
```jldoctest
julia> p = hkf_debye_huckel_params(298.15, 1e5);

julia> isapprox(p.A, 0.5114; rtol=1e-3)
true

julia> isapprox(p.B, 0.3288; rtol=1e-3)
true
```
"""
function hkf_debye_huckel_params(T_K, P_Pa)
    T_K, P_Pa = promote(T_K, P_Pa)
    wtp = water_thermo_props(T_K, P_Pa)
    wep = water_electro_props_jn(T_K, P_Pa, wtp)
    ρ_gcm3 = wtp.D / 1000                       # kg/m³ → g/cm³
    εT = wep.epsilon * T_K                   # dimensionless × K
    A = 1.824829238e6 * sqrt(ρ_gcm3) / εT^(3 // 2)
    B = 50.29158649 * sqrt(ρ_gcm3) / sqrt(εT)
    return (A = A, B = B)
end

# Internal: four-level ionic radius priority lookup.
function _hkf_lookup_å(sp::AbstractSpecies, model)
    if haskey(properties(sp), :å)
        v = sp[:å]
        return v isa Number ? float(v) : float(safe_ustrip(1.0u"Å", v))
    end
    pf = phreeqc(formula(sp))
    haskey(REJ_HKF, pf)           && return REJ_HKF[pf]
    z = Int(charge(sp))
    haskey(REJ_CHARGE_DEFAULT, z) && return REJ_CHARGE_DEFAULT[z]
    return float(model.å_default)
end

# ── HKFActivityModel ──────────────────────────────────────────────────────────

"""
    struct HKFActivityModel{T<:Real} <: AbstractActivityModel

Extended Debye-Hückel activity model with B-dot term (Helgeson 1969 /
Helgeson, Kirkham & Flowers 1981). This is the model used by PHREEQC and EQ3/6.

# Activity coefficient formulas

For ionic species (charge `z ≠ 0`):
```
log₁₀ γᵢ = −A zᵢ² √I / (1 + B åᵢ √I)  +  Ḃ I
```

For neutral aqueous species (charge `z = 0`):
```
log₁₀ γᵢ = Kₙ I
```

Log-activity (molality convention, standard state `m° = 1 mol/kg`):
```
ln aᵢ = ln(10) × log₁₀(γᵢ) + ln(mᵢ)
```

Water activity is computed from the osmotic coefficient φ (Gibbs-Duhem):
```
ln a_w = −Mw × Σmⱼ × φ
```

Ionic strength: `I = ½ Σ mⱼ zⱼ²`

# Fields

  - `A`: Debye-Hückel A parameter [(kg/mol)^(1/2)]. Default 0.5114 at 25 °C/1 bar.
  - `B`: Debye-Hückel B parameter [Å⁻¹(kg/mol)^(1/2)]. Default 0.3288.
  - `Ḃ`: B-dot extended term [kg/mol]. Default 0.041.
  - `Kₙ`: salting-out coefficient for neutral species [kg/mol]. Default 0.1.
  - `å_default`: global fallback effective ionic radius [Å]. Default 3.72.
  - `temperature_dependent`: if `true`, recompute A and B from `p.T`, `p.P`
    at each call to the activity closure (requires `T` and `P` in `p`).
    Default `false`.

# Ionic radius lookup

The effective radius åᵢ is resolved in order:
1. `sp[:å]` — explicit value in the species properties dict.
2. [`REJ_HKF`](@ref) — Helgeson et al. (1981) Table 3, keyed by PHREEQC formula.
3. [`REJ_CHARGE_DEFAULT`](@ref) — fallback by formal charge.
4. `model.å_default` — global fallback.

# Valid range

`I ≲ 1 mol/kg`. Beyond ~2 mol/kg, the Pitzer model is recommended.

# References

  - Helgeson, H.C. (1969). Am. J. Sci. **267**, 729–804.
  - Helgeson, H.C., Kirkham, D.H. & Flowers, G.C. (1981). Am. J. Sci. **281**, 1249–1516.
  - Parkhurst, D.L. & Appelo, C.A.J. (2013). USGS Techniques Methods, Book 6, ch. A43.

# Examples
```julia
# Default model at 25 °C / 1 bar (fixed A, B)
model = HKFActivityModel()

# Temperature-dependent A and B (recomputed at each solve)
model_tdep = HKFActivityModel(temperature_dependent=true)

state_eq = equilibrate(state; model=HKFActivityModel())
```
"""
struct HKFActivityModel{T <: Real} <: AbstractActivityModel
    A::T
    B::T
    Ḃ::T
    Kₙ::T
    å_default::T
    temperature_dependent::Bool
end

"""
    HKFActivityModel(; A=0.5114, B=0.3288, Ḃ=0.041, Kₙ=0.1, å_default=3.72,
                       temperature_dependent=false) -> HKFActivityModel

Construct an [`HKFActivityModel`](@ref) with the given parameters.

Default values are from Helgeson et al. (1981), Table 1, at 25 °C / 1 bar.
"""
function HKFActivityModel(;
        A::Real = 0.5114,
        B::Real = 0.3288,
        Ḃ::Real = 0.041,
        Kₙ::Real = 0.1,
        å_default::Real = 3.72,
        temperature_dependent::Bool = false,
    )
    vals = promote(A, B, Ḃ, Kₙ, å_default)
    return HKFActivityModel{eltype(vals)}(vals..., temperature_dependent)
end

"""
    activity_model(cs::ChemicalSystem, model::HKFActivityModel) -> Function

Return a closure `lna(n, p) -> Vector` computing log-activities for the
extended Debye-Hückel (B-dot) model of Helgeson (1969).

The closure captures all species indices and ionic radii at construction time.
Inside `lna`:
- Solutes: molality convention, B-dot formula for ions, salting-out for neutrals.
- Solvent: osmotic coefficient from Gibbs-Duhem (σ-function).
- Crystals: `ln a = 0` (pure solid).
- Gas: ideal mixture `ln a = ln(xᵢ)`.

If `model.temperature_dependent=true`, `p` must contain `T` (K) and `P` (Pa)
— both are provided automatically by `_build_params`.

AD-compatible: all closure computations accept `ForwardDiff.Dual` inputs.
"""
function activity_model(cs::ChemicalSystem, model::HKFActivityModel)

    # ── Precompute at closure-construction time ────────────────────────────
    idx_solvent = only(cs.idx_solvent)
    idx_solutes = cs.idx_solutes
    idx_gas = cs.idx_gas

    ss_groups = cs.ss_groups
    has_ss = !isempty(ss_groups)
    has_gas = !isempty(idx_gas)
    ss_models = has_ss ? map(ss -> ss.model, cs.solid_solutions) : nothing

    M_w = ustrip(us"kg/mol", cs.species[idx_solvent][:M])   # kg/mol, e.g. 0.018015

    A_fixed = model.A
    B_fixed = model.B
    Ḃ = model.Ḃ
    Kₙ = model.Kₙ
    temp_dep = model.temperature_dependent

    # Per-species data (Float64 — not differentiated).
    zv = Int8[charge(sp) for sp in cs.species]
    åv = Float64[
        iszero(zv[i]) ? 0.0 : _hkf_lookup_å(cs.species[i], model)
            for i in eachindex(zv)
    ]
    n_sp = lastindex(zv)

    idx_ions = [i for i in idx_solutes if !iszero(zv[i])]
    idx_neutrals = [i for i in idx_solutes if  iszero(zv[i])]

    ln10 = log(10.0)

    function lna(n::AbstractVector, p)
        ϵ = p.ϵ
        _n = max.(n, ϵ)

        # ── A and B (fixed or T,P-dependent) ──────────────────────────────
        if temp_dep && hasproperty(p, :T) && hasproperty(p, :P)
            AB = hkf_debye_huckel_params(p.T, p.P)
            A, B = AB.A, AB.B
        else
            A, B = A_fixed, B_fixed
        end

        out = zeros(eltype(_n), n_sp)

        # ── Molality: mᵢ = nᵢ / (n_w × M_w) [mol/kg] ─────────────────────
        n_w = _n[idx_solvent]
        denom_mol = n_w * M_w             # kg of solvent

        # ── Ionic strength I = ½ Σ mⱼ zⱼ² ────────────────────────────────
        I = zero(eltype(_n))
        @inbounds for i in idx_solutes
            mᵢ = _n[i] / denom_mol
            I = I + mᵢ * zv[i]^2
        end
        I = I / 2
        sqrtI = sqrt(I + ϵ)              # regularized to avoid Dual NaN at I=0

        # ── Effective radius å_eff for osmotic coefficient ─────────────────
        sum_mz2a = zero(eltype(_n))
        @inbounds for i in idx_ions
            sum_mz2a = sum_mz2a + (_n[i] / denom_mol) * zv[i]^2 * åv[i]
        end
        sum_mz2 = 2 * I                 # Σ mⱼ zⱼ² = 2I by definition
        # Smooth blend: avoids branching on Dual values at ionic-strength ≈ 0
        å_eff = (sum_mz2a + model.å_default * ϵ) / (sum_mz2 + ϵ)

        # ── Ion log-activity coefficients ──────────────────────────────────
        @inbounds for i in idx_ions
            denom_dh = 1 + B * åv[i] * sqrtI
            log10γᵢ = -A * zv[i]^2 * sqrtI / denom_dh + Ḃ * I
            mᵢ = _n[i] / denom_mol
            out[i] = ln10 * log10γᵢ + log(mᵢ + ϵ)
        end

        # ── Neutral solute log-activities ──────────────────────────────────
        @inbounds for i in idx_neutrals
            log10γᵢ = Kₙ * I
            mᵢ = _n[i] / denom_mol
            out[i] = ln10 * log10γᵢ + log(mᵢ + ϵ)
        end

        # ── Water activity via osmotic coefficient (Gibbs-Duhem) ───────────
        sum_m = zero(eltype(_n))
        @inbounds for i in idx_solutes
            sum_m = sum_m + _n[i] / denom_mol
        end
        x_arg = B * å_eff * sqrtI
        σ = _hkf_sigma(x_arg)
        φ = 1 - (A * ln10 / 3) * (sum_mz2 / (sum_m + ϵ)) * sqrtI * σ +
            (Ḃ * ln10 / 2) * I
        out[idx_solvent] = -M_w * sum_m * φ

        # ── Gas: ideal mixture ─────────────────────────────────────────────
        if has_gas
            n_gas = sum((_n[i] for i in idx_gas); init = zero(eltype(_n)))
            @inbounds for i in idx_gas
                out[i] = log(_n[i] / n_gas)
            end
        end

        # ── Solid solutions ────────────────────────────────────────────────
        if has_ss
            T_val = hasproperty(p, :T) ? p.T : 298.15
            _solid_solution_lna!(out, _n, ss_groups, ss_models, T_val, ϵ)
        end

        return out
    end

    return lna
end

# ── DaviesActivityModel ───────────────────────────────────────────────────────

"""
    struct DaviesActivityModel{T<:Real} <: AbstractActivityModel

Davies (1962) activity model — a simplified Debye-Hückel without ionic radii.

# Activity coefficient formula

For ionic species (charge `z ≠ 0`):
```
log₁₀ γᵢ = −A zᵢ² (√I / (1 + √I)  −  b I)
```

For neutral aqueous species (charge `z = 0`):
```
log₁₀ γᵢ = bₙ I
```

Log-activity (molality convention): `ln aᵢ = ln(10) × log₁₀(γᵢ) + ln(mᵢ)`

Water activity: Raoult approximation `ln a_w = ln(x_w)` (V1; accurate for I < 0.1).

# Fields

  - `A`: Debye-Hückel A parameter [(kg/mol)^(1/2)]. Default 0.5114 at 25 °C/1 bar.
  - `b`: Davies empirical constant (default 0.3).
  - `bₙ`: salting-out coefficient for neutral species (default 0.1).
  - `temperature_dependent`: if `true`, recompute A from `p.T`, `p.P` at each call.

# References

  - Davies, C.W. (1962). *Ion Association*. Butterworths, London.
"""
struct DaviesActivityModel{T <: Real} <: AbstractActivityModel
    A::T
    b::T
    bₙ::T
    temperature_dependent::Bool
end

"""
    DaviesActivityModel(; A=0.5114, b=0.3, bₙ=0.1, temperature_dependent=false)

Construct a [`DaviesActivityModel`](@ref).
"""
function DaviesActivityModel(;
        A::Real = 0.5114,
        b::Real = 0.3,
        bₙ::Real = 0.1,
        temperature_dependent::Bool = false,
    )
    vals = promote(A, b, bₙ)
    return DaviesActivityModel{eltype(vals)}(vals..., temperature_dependent)
end

"""
    activity_model(cs::ChemicalSystem, model::DaviesActivityModel) -> Function

Return a closure `lna(n, p) -> Vector` computing log-activities for the
Davies (1962) model. No species-specific ionic radii are required.

AD-compatible: all closure computations accept `ForwardDiff.Dual` inputs.
"""
function activity_model(cs::ChemicalSystem, model::DaviesActivityModel)

    idx_solvent = only(cs.idx_solvent)
    idx_solutes = cs.idx_solutes
    idx_gas = cs.idx_gas

    ss_groups = cs.ss_groups
    has_ss = !isempty(ss_groups)
    has_gas = !isempty(idx_gas)
    ss_models = has_ss ? map(ss -> ss.model, cs.solid_solutions) : nothing

    M_w = ustrip(us"kg/mol", cs.species[idx_solvent][:M])

    A_fixed = model.A
    b = model.b
    bₙ = model.bₙ
    temp_dep = model.temperature_dependent

    zv = Int8[charge(sp) for sp in cs.species]
    n_sp = lastindex(zv)
    idx_ions = [i for i in idx_solutes if !iszero(zv[i])]
    idx_neutrals = [i for i in idx_solutes if  iszero(zv[i])]

    ln10 = log(10.0)

    function lna(n::AbstractVector, p)
        ϵ = p.ϵ
        _n = max.(n, ϵ)

        A = if temp_dep && hasproperty(p, :T) && hasproperty(p, :P)
            hkf_debye_huckel_params(p.T, p.P).A
        else
            A_fixed
        end

        out = zeros(eltype(_n), n_sp)

        n_w = _n[idx_solvent]
        denom_mol = n_w * M_w

        # Ionic strength
        I = zero(eltype(_n))
        @inbounds for i in idx_solutes
            I = I + (_n[i] / denom_mol) * zv[i]^2
        end
        I = I / 2
        sqrtI = sqrt(I + ϵ)
        dI = sqrtI / (1 + sqrtI)     # √I / (1 + √I)

        # Ions
        @inbounds for i in idx_ions
            log10γᵢ = -A * zv[i]^2 * (dI - b * I)
            mᵢ = _n[i] / denom_mol
            out[i] = ln10 * log10γᵢ + log(mᵢ + ϵ)
        end

        # Neutral solutes
        @inbounds for i in idx_neutrals
            mᵢ = _n[i] / denom_mol
            out[i] = ln10 * bₙ * I + log(mᵢ + ϵ)
        end

        # Water activity — Raoult (mole fraction) approximation
        # n_aqueous ≥ ϵ > 0 always (because _n[i] ≥ ϵ), so no iszero guard needed
        n_aqueous = n_w + sum((_n[i] for i in idx_solutes); init = zero(eltype(_n)))
        out[idx_solvent] = log(n_w / n_aqueous)

        # Gas: ideal mixture
        if has_gas
            n_gas = sum((_n[i] for i in idx_gas); init = zero(eltype(_n)))
            @inbounds for i in idx_gas
                out[i] = log(_n[i] / n_gas)
            end
        end

        # Solid solutions
        if has_ss
            T_val = hasproperty(p, :T) ? p.T : 298.15
            _solid_solution_lna!(out, _n, ss_groups, ss_models, T_val, ϵ)
        end

        return out
    end

    return lna
end

# ── Solid solution activity helpers ───────────────────────────────────────────

"""
    _excess_ln_gamma(model, k, x, T) -> Real

Return the excess log-activity coefficient `ln γₖ` for end-member `k` (1-based index)
of a solid solution with mole-fraction vector `x` at temperature `T` (K).

AD-compatible: all branches preserve `ForwardDiff.Dual` through computations on `x`.

Methods:
- [`IdealSolidSolutionModel`](@ref): returns `zero(eltype(x))`.
- [`RedlichKisterModel`](@ref): binary Redlich-Kister formula (requires `length(x) == 2`).
"""
_excess_ln_gamma(::IdealSolidSolutionModel, k::Int, x::AbstractVector, T::Real) =
    zero(eltype(x))

function _excess_ln_gamma(m::RedlichKisterModel, k::Int, x::AbstractVector, T::Real)
    x1, x2 = x[1], x[2]
    RT = 8.31446261815324 * T   # J/mol
    a0 = m.a0 / RT
    a1 = m.a1 / RT
    a2 = m.a2 / RT
    if k == 1
        return x2^2 * (a0 + a1 * (3 * x1 - x2) + a2 * (x1 - x2) * (5 * x1 - x2))
    else
        return x1^2 * (a0 - a1 * (3 * x2 - x1) + a2 * (x2 - x1) * (5 * x2 - x1))
    end
end

"""
    _solid_solution_lna!(out, _n, ss_groups, ss_models, T, ϵ)

Fill `out[i]` with `ln aᵢ = ln xᵢ + ln γᵢ` for all solid-solution end-members.

`ss_groups[k]` and `ss_models[k]` describe the k-th solid-solution phase.
`T` is the temperature in K (only relevant for non-ideal models).
`ϵ` is a regularization floor to avoid `log(0)`.

ForwardDiff-compatible.
"""
function _solid_solution_lna!(
        out::AbstractVector, _n::AbstractVector{ET},
        ss_groups::Vector{Vector{Int}}, ss_models, T, ϵ
    ) where {ET}
    for (grp, mdl) in zip(ss_groups, ss_models)
        n_total = sum(_n[i] for i in grp) + ϵ
        # ET follows eltype(_n) — AD-compatible (Dual when differentiating)
        x = Vector{ET}(undef, length(grp))
        @inbounds for (j, i) in enumerate(grp)
            x[j] = _n[i] / n_total
        end
        @inbounds for (k, i) in enumerate(grp)
            out[i] = log(x[k] + ϵ) + _excess_ln_gamma(mdl, k, x, T)
        end
    end
    return out
end
