# =============================================================================
#  hydration_calibration.jl — calibrating OPC hydration kinetics against
#                             measured isothermal calorimetry
#
#  `ionic_hydration.jl` runs a complete CEM I forward and reads its calorimetry
#  off a certified replay, with the *published* Parrott & Killoh (1984)
#  parameters. This script closes the loop the other way: given a measured heat
#  curve, which kinetic parameters reproduce it?
#
#  Only KINETIC parameters are fitted. The thermodynamics — ΔₐG⁰, ΔₐH⁰, Cp⁰,
#  log K — come from CEMDATA18 and are not adjustable: they are measured
#  quantities, and moving them to improve a fit would be fabricating chemistry.
#  Blaine fineness, w/b ratio and temperature are not fitted either, because the
#  experiment reports them.
#
#  The data are the CC-BY-4.0 deposit of Šmilauer & Reiterman (2025); see
#  `data/experimental/README.md` and `data/experimental/LICENSE`.
#
#  Usage:
#      julia --project=scripts scripts/hydration_calibration.jl
#      or, from the REPL:  include("scripts/hydration_calibration.jl")
#
#  This script deliberately does NOT call `Pkg.activate`:
#  `docs/src/examples/hydration_calibration.md` and
#  `test/kinetics/test_calibration.jl` include this file, and the active
#  project is global process state — activating an environment halfway through a
#  documentation build changes it for every `@example` block that follows.
#
#  The coupled stage is expensive — one Gibbs minimization per accepted step,
#  plus a certified replay per output instant. `main()` reports its own cost.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using LinearAlgebra
using OptimaSolver
using Optimization
using OptimizationOptimJL
# `NelderMead` explicitly, not through the `using` above: OptimizationOptimJL
# 0.4.19 stopped re-exporting Optim's algorithm names. The binding still exists
# in that module, so an explicit import works on either version — whereas
# relying on the re-export made this script depend on which version the
# environment happened to resolve. It failed in the test environment, which
# resolves fresh from the registry, while a pinned `scripts/Manifest.toml`
# still carrying 0.4.18 kept working.
import OptimizationOptimJL: NelderMead
using OrderedCollections
using OrdinaryDiffEq
using Printf
using SciMLBase
using Statistics

# The forward model, the chemical system and the dissolution reactions are the
# ones of the companion script; nothing about the physics is restated here.
isdefined(Main, :run_ionic_hydration) || include(joinpath(@__DIR__, "ionic_hydration.jl"))

# ── the measured data ─────────────────────────────────────────────────────────

"""
    CALORIMETRY_DIR

Directory holding the vendored measured records. Its `README.md` documents the
format and its `LICENSE` the terms — the data are CC-BY-4.0 and **not** covered
by the package license.
"""
const CALORIMETRY_DIR = datapath("experimental")

"""
    CalorimetryData

One isothermal-calorimetry record.

  - `t` — hydration time [s], strictly increasing, **not** starting at zero:
    every calorimeter record begins after thermal equilibration.
  - `qdot` — heat flow [W/g of binder].
  - `Q` — released heat [J/g of binder], counted from `t[1]`, not from mixing.
  - `Qref` — **not a measurement**: the depositors' own fitted curve at the same
    instants, or empty when the file carries none. It is there so a fit can be
    judged against somebody else's on the very same record.
  - `meta` — what the experiment reported and a fit must therefore not touch:
    `cement`, `blaine`, `wb`, `T`, plus `q_before_start` (the heat already gone
    when the record opens) and the provenance fields `source`, `doi`, `license`.

!!! note "Q is counted from the first sample"
    Both the record and the model are referred to `t[1]`. Aligning them there
    rather than at mixing avoids inventing the missing first hour: the induction
    period is exactly what neither the instrument nor a Parrott–Killoh rate law
    describes. `meta.q_before_start` is carried so the *total* heat can still be
    reported, and it is the depositors' own estimate, not a measurement.
"""
struct CalorimetryData{M}
    t::Vector{Float64}
    qdot::Vector{Float64}
    Q::Vector{Float64}
    Qref::Vector{Float64}
    meta::M
end

Base.length(d::CalorimetryData) = length(d.t)

function Base.show(io::IO, d::CalorimetryData)
    return @printf(
        io, "CalorimetryData(%s, w/b %.2f, Blaine %.0f m²/kg, %.1f–%.1f h, %d points)",
        d.meta.cement, d.meta.wb, d.meta.blaine, d.t[1] / 3600, d.t[end] / 3600, length(d.t)
    )
end

"""
    _meta_field(comments, key) -> Union{String, Nothing}

Value of a `# key: value` provenance line, or `nothing`.
"""
function _meta_field(comments, key)
    for line in comments
        clean = strip(lstrip(line, ['#', ' ']))
        startswith(clean, key * ":") || continue
        return strip(clean[(length(key) + 2):end])
    end
    return nothing
end

"""
    _leading_float(s) -> Float64

First number in a string, so that `"397 m2/kg"` and `"20C"` both parse. Returns
`NaN` when there is none.
"""
function _leading_float(s)
    m = match(r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?", s)
    return m === nothing ? NaN : parse(Float64, m.match)
end

"""
    read_calorimetry(path) -> CalorimetryData

Read one record: `#` comment lines carrying the provenance and the experimental
conditions, then a comma-separated table headed
`time_h,heat_flow_W_per_g,heat_J_per_g`.

**This is the seam for your own data.** Any file in that shape works. The
comment lines the reader looks for are `cement`, `blaine`, `wb`, `temperature`
and `Released heat up to`; a missing one becomes `NaN` or an empty string rather
than an error, but `blaine` and `wb` are needed by the forward model and their
absence will surface there. A fourth column, if present, is read as somebody
else's fitted curve and never as data.
"""
function read_calorimetry(path)
    comments = String[]
    t_h = Float64[]
    qdot = Float64[]
    Q = Float64[]
    Qref = Float64[]
    for line in eachline(path)
        s = strip(line)
        isempty(s) && continue
        if startswith(s, "#")
            push!(comments, s)
        elseif startswith(s, "time_h")
            continue
        else
            f = split(s, ',')
            length(f) >= 3 || error("malformed data line in $path: $s")
            push!(t_h, parse(Float64, f[1]))
            push!(qdot, parse(Float64, f[2]))
            push!(Q, parse(Float64, f[3]))
            length(f) >= 4 && push!(Qref, parse(Float64, f[4]))
        end
    end
    isempty(t_h) && error("no data rows in $path")

    released = _meta_field(comments, "Released heat up to 45 minutes (J/g of binder)")
    meta = (
        cement = something(_meta_field(comments, "cement"), "unknown"),
        blaine = _leading_float(something(_meta_field(comments, "blaine"), "")),
        wb = _leading_float(something(_meta_field(comments, "wb"), "")),
        T = _leading_float(something(_meta_field(comments, "temperature"), "")) + 273.15,
        q_before_start = released === nothing ? NaN : _leading_float(released),
        source = something(_meta_field(comments, "source file"), basename(path)),
        doi = something(_meta_field(comments, "doi"), ""),
        license = something(_meta_field(comments, "license"), ""),
    )
    return CalorimetryData(t_h .* 3600, qdot, Q, Qref, meta)
end

"""
    CEM_I_TARGET, CEM_I_HOLDOUT

The two vendored records: the CEM I 52.5 R of Čížkovice at w/b = 0.50 and Blaine
397 m²/kg, fitted; and the CEM I 52.5 R of Ladce at w/b = 0.45 and Blaine
415 m²/kg, predicted without refitting.

The holdout is the only honest test available here of whether
[`powers_alpha_max`](@ref) and [`blaine_factor`](@ref) carry a calibration
across mixes, which is the claim those two corrections make.
"""
const CEM_I_TARGET = read_calorimetry(
    joinpath(CALORIMETRY_DIR, "smilauer2025-122-cemI-52.5R-cizkovice.csv")
)

"""
    CEM_I_HOLDOUT

The CEM I 52.5 R of Ladce at w/b = 0.45 and Blaine 415 m²/kg, over 26 days —
predicted with the parameters fitted on [`CEM_I_TARGET`](@ref), never fitted
itself. See [`CEM_I_TARGET`](@ref).
"""
const CEM_I_HOLDOUT = read_calorimetry(
    joinpath(CALORIMETRY_DIR, "smilauer2025-116-cemI-52.5R-ladce.csv")
)

"""
    resample_log(data, n; tmax = Inf) -> CalorimetryData

Subset of `data` on `n` log-spaced times, by nearest sample.

A calorimeter samples uniformly in time, which weights a ten-day plateau some
hundred times more heavily than the acceleration peak. Hydration is a decades-of-
time process and the residuals belong on a log grid; `n ≈ 60` carries every
feature of these curves at a hundredth of the forward-solve count.
"""
function resample_log(data::CalorimetryData, n; tmax = Inf)
    hi = something(findlast(<=(tmax), data.t), lastindex(data.t))
    targets = 10 .^ range(log10(data.t[1]), log10(data.t[hi]); length = n)
    idx = sort!(unique!([argmin(abs.(data.t .- τ)) for τ in targets]))
    Qref = isempty(data.Qref) ? Float64[] : data.Qref[idx]
    return CalorimetryData(data.t[idx], data.qdot[idx], data.Q[idx], Qref, data.meta)
end

# ── the clinker, which is assumed and not fitted ──────────────────────────────

"""
    CALIB_CLINKER, CALIB_GYPSUM, CALIB_FILLER

The formulation held fixed through every fit: the CEM I 52.5 N of Lavergne et al.
(2018), Table 9 — C₃S 65 / C₂S 11 / C₃A 11 / C₄AF 8 by mass of clinker, with
4.6 % gypsum and 3.5 % limestone.

!!! warning "The deposit does not report a phase composition"
    It reports the cement type, the Blaine fineness, the w/b ratio and the
    temperature — not the clinker mineralogy. So a composition has to be
    supplied, and the consequence is exact rather than vague: a heat curve
    constrains the *product* of a phase fraction and its reaction rate, never
    the two apart. What is calibrated below is a rate **given this composition**,
    and §"Composition sensitivity" of the documentation page measures what that
    assumption costs.

    This particular composition was chosen because the target record's w/b of
    0.50 and Blaine of 397 m²/kg sit within a few percent of the mix it was
    published for, and because `ionic_hydration.jl` already runs it — so the two
    scripts describe one cement, not two.
"""
const CALIB_CLINKER = (C3S = 0.65, C2S = 0.11, C3A = 0.11, C4AF = 0.08)

"""
    CALIB_GYPSUM

Gypsum, 4.6 % of the binder mass. Fixed; see [`CALIB_CLINKER`](@ref).
"""
const CALIB_GYPSUM = 0.046

"""
    CALIB_FILLER

Limestone filler, 3.5 % of the binder mass. Fixed; see [`CALIB_CLINKER`](@ref).
"""
const CALIB_FILLER = 0.035

# ── the induction period, which Parrott & Killoh does not have ────────────────

"""
    INDUCTION_PHASES

Phases whose dissolution the induction factor slows: the two silicates.

The dormant period is a silicate phenomenon — dissolution inhibited at the alite
surface until the etch-pit mechanism opens up, in the reading of
[Scrivener2019](@cite). The aluminate reacts within minutes of wetting, and its
first peak falls before these records even open, so slowing it would be
inventing a delay nobody measured here.
"""
const INDUCTION_PHASES = ("C3S", "C2S")

"""
    induction_factor(τ, m) -> Function

`β(t) = 1 - exp(-(t/τ)^m)`, a dimensionless factor in `[0, 1)` multiplying a
dissolution rate: an induction period of characteristic length `τ` [s] and
sharpness `m`.

# Why this is needed at all

`parrot_killoh_avrami` has **no induction period**. Its Avrami branch is floored
at `PK_AVRAMI_SEED` precisely so that the ODE leaves `ξ = 0` at all, which means
hydration starts the instant the water does. Real cement does not: it has a
dormant hour or several. Measured against the record used here, the unmodified
model has released 23.8 J/g at 2.7 h where the calorimeter saw 4.7, and 63.3 J/g
at 5.7 h against 23.2 — and having spent that heat early it then runs *behind*,
by some 60 J/g at two days, and never recovers.

No rate constant repairs that, because the two errors have opposite signs: a
slower rate improves the early residual and worsens the late one. It is a missing
term, not a mis-set parameter. `β(t)` is that term, in the standard
nucleation-with-induction form, and it is what makes a mechanistic model
comparable with the empirical affinity functions whose second parameter does the
same job implicitly.

`τ` and `m` are kinetic parameters like any other and are fitted.

!!! note "Why `m` is boxed narrowly"
    A large `m` makes `β` a near-step function. That is not stiffness — an
    L-stable method handles stiffness by construction — but it is a rapid,
    localized variation of the right-hand side in *time*, which any adaptive
    integrator has to find by rejecting steps unless it is told where to look.
    [`induction_tstops`](@ref) tells it. `[0.8, 4.0]` is then a statement about
    the physics rather than about the solver: a dormant period that ends over an
    hour or over six, not one that ends instantaneously.
"""
induction_factor(τ, m) = t -> -expm1(-(max(t, zero(t)) / τ)^m)

"""
    induction_tstops(τ, m; n = 9) -> Vector{Float64}

Times spanning the induction transition, for `tstops`.

**Not used by default, and the reason is measured.** `β(t) = 1 - exp(-(t/τ)^m)`
climbs from ≈0 to ≈1 over roughly `τ/4` to `4τ`, a rapid but purely *explicit*
variation of the right-hand side. Forcing steps there looked like the obvious
thing to do — it is the SciML idiom for a localized feature in time — and on this
problem it makes matters worse: 117 steps and 26 ms with these `tstops` against
113 steps and 12 ms without them, on the same trajectory to the same answer. The
adaptive controller locates a three-hour transition inside a 262-hour integration
without help, and the forced steps are simply extra.

Kept because it is the right tool for a *sharper* transition than this one, and
because the negative result is worth being able to reproduce: pass it as `tstops`
to [`stoichiometric_run`](@ref) or [`run_ionic_hydration`](@ref).
"""
induction_tstops(τ, m; n = 9) = collect(τ .* 10 .^ range(-0.6, 0.6; length = n))

"""
    damped_rate(base::KineticFunc, β) -> KineticFunc

`base` multiplied by `β(t)`, or `base` unchanged when `β === nothing`.

The same closure-wrapping pattern `ionic_hydration.jl` uses for its per-phase
calibration factors, so the rate laws themselves are never edited.
"""
function damped_rate(base, β)
    β === nothing && return base
    return KineticFunc(
        (T, P, t, n, lna, n0) -> β(t) * base(T, P, t, n, lna, n0),
        (T = 293.15u"K", P = 1.0e5u"Pa"), u"mol/s",
    )
end

"""
    INDUCTION_KEY

Pseudo-phase marking a [`CalibParameter`](@ref) that configures
[`induction_factor`](@ref) rather than a Parrott–Killoh field. Its `field` is
`:τ` (seconds) or `:m`.
"""
const INDUCTION_KEY = "induction"

# ── what is calibrated, and what is deliberately not ──────────────────────────

"""
    CalibParameter

One free parameter: its `name`, the `phase` and Parrott–Killoh field it writes
(`field`), the published `prior`, a `(lo, hi)` box, and `logscale` telling the
optimizer to work in `log` of it.

Rates span decades and must be searched multiplicatively; exponents are O(1) and
must not.
"""
struct CalibParameter
    name::Symbol
    phase::String
    field::Symbol
    prior::Float64
    lo::Float64
    hi::Float64
    logscale::Bool
end

"""
    CALIB_SPEC_FULL :: Vector{CalibParameter}

The six **candidate** parameters — one per mechanism with a distinct signature on
the heat curve — chosen so that none is redundant with another *by construction*.

| θ | phase | mechanism it governs | prior |
|---|---|---|---|
| `k₁_C3S` | alite | Avrami nucleation and growth — the acceleration | 1.5 d⁻¹ |
| `n₁_C3S` | alite | Avrami exponent — the shape of the main peak | 0.7 |
| `k₃_C3S` | alite | power-law shell — where deceleration sets in | 1.1 d⁻¹ |
| `n₃_C3S` | alite | power-law exponent — how hard it decelerates | 3.3 |
| `k₁_C3A` | aluminate | Avrami rate — the sulfate-depletion shoulder | 1.0 d⁻¹ |
| `k₃_C2S` | belite | power-law shell — the multi-day plateau | 0.2 d⁻¹ |

Six is what the *model* offers. [`CALIB_SPEC`](@ref) is what the *measurement*
supports, which is fewer; this list is kept for the analysis that establishes
that, and for anyone whose data support more.

# Excluded before any fitting, and why

  - **Every `Ea`.** These records were measured at one temperature, 20 °C. An
    activation energy is a *temperature* sensitivity; at one temperature it is not
    identifiable, and fitting it would report a number the data cannot contain.
  - **A uniform per-phase multiplier**, of the kind `IONIC_CALIBRATION` offers.
    For alite the rate is `min(α̇₁, α̇₃)`, so `f · min(α̇₁, α̇₃) = min(f·α̇₁, f·α̇₃)`:
    scaling the phase is *exactly* scaling `k₁` and `k₃` together. Fitting both
    would make the problem rank-deficient by construction rather than by accident.
    Use `IONIC_CALIBRATION` **or** the rate-law fields, never both.
  - **`k₂` (Jander diffusion) for alite**, but *not* for the reason usually given.
    Parrott & Killoh reported no diffusion-controlled stage for C₃S, and
    [`parrot_killoh_avrami`](@ref)'s docstring repeats it — so the expected
    sensitivity is zero. Measured, it is not: `α̇₂ = k₂(1-ξ)^{2/3}/(1-(1-ξ)^{1/3})`
    falls as ξ grows, so past a high degree of hydration the Jander branch does
    become the minimum, and `k₂` moves the released heat by up to about an eighth
    of what `k₁` does, entirely after the first day. It is left out because that is
    a minor effect competing for a place in a three-dimensional identifiable space,
    not because it is absent. The documentation page measures it rather than
    asserting either version.
  - **Gypsum and calcite dissolution.** [`ionic_reactions`](@ref) deliberately
    makes them fast so sulfate is available to the minimization from the start
    rather than rate-limiting it. A step that does not limit cannot be calibrated.
  - **The whole thermodynamic database.** Not a parameter. See the file header.
"""
const CALIB_SPEC_FULL = [
    CalibParameter(:k₁_C3S, "C3S", :k₁, 1.5, 1.5 / 4, 1.5 * 4, true),
    CalibParameter(:n₁_C3S, "C3S", :n₁, 0.7, 0.3, 1.3, false),
    CalibParameter(:k₃_C3S, "C3S", :k₃, 1.1, 1.1 / 6, 1.1 * 12, true),
    CalibParameter(:n₃_C3S, "C3S", :n₃, 3.3, 1.5, 6.0, false),
    CalibParameter(:k₁_C3A, "C3A", :k₁, 1.0, 1.0 / 4, 1.0 * 4, true),
    CalibParameter(:k₃_C2S, "C2S", :k₃, 0.2, 0.2 / 4, 0.2 * 4, true),
    CalibParameter(:τ_ind, INDUCTION_KEY, :τ, 4.0 * 3600, 0.2 * 3600, 24.0 * 3600, true),
    CalibParameter(:m_ind, INDUCTION_KEY, :m, 2.0, 0.8, 4.0, false),
]

"""
    CALIB_SPEC :: Vector{CalibParameter}

The parameters actually fitted: the three of [`CALIB_SPEC_FULL`](@ref) that a
single isothermal heat curve can carry.

Three, not six, and the reason is measured rather than assumed. Over the six
candidates the singular values of `∂Q/∂log θ` come out at roughly
`[420, 100, 60, 6.3, 1.4, 0.20]` on this record — a factor of nine between the
third and the fourth — so the measurement determines three *combinations* of the
six, not six numbers.

Which three is decided by the correlation matrix, and it is unusually clean. Two
pairs are almost perfectly collinear:

  - `k₁_C3S` with `n₁_C3S`, correlation −0.985;
  - `k₃_C3S` with `n₃_C3S`, correlation −0.977 — between them they are essentially
    the whole leading singular direction.

Within each pair the data see a combination and not its members, so one per pair
is all that can be fitted. The **rate constant** is the one kept in each: an
Avrami or power-law exponent is a shape parameter of an empirical 1984 fit,
whereas a rate constant is the quantity a different clinker plausibly changes.
`n₁_C3S` and `n₃_C3S` therefore stay at their published values.

The third is `k₁_C3A`, which is separately visible — very nearly the fourth
singular direction on its own — if weakly. `k₃_C2S` is dropped outright: over
eleven days belite carries a weight of 0.001 in the three directions the data see,
which is a quantitative way of saying invisible.

!!! note "A parameter at a bound is not a result"
    Each box spans a factor of four either side of the published value. A fit that
    pins a parameter against a bound is not reporting that value — it is reporting
    that the model cannot make this curve and is saturating whatever it is handed.
    Fitting all six does exactly that, on four of them at once.
    [`box_position`](@ref) is printed with every fit for that reason.
"""
const CALIB_SPEC = CALIB_SPEC_FULL[[1, 3, 5, 7, 8]]

"""
    PK84_PRIOR :: Dict{String, NamedTuple}

The published parameter sets, by phase symbol — the starting point of every fit
and the reference every result is reported against. Aliases `IONIC_PK84` of the
companion script, so the two describe one cement and not two.
"""
const PK84_PRIOR = IONIC_PK84

"""
    to_native(z, spec = CALIB_SPEC) -> Vector{Float64}
    to_search(θ, spec = CALIB_SPEC) -> Vector{Float64}

Map between the unbounded space `z ∈ ℝⁿ` the optimizer searches and native
parameter values inside their boxes.

Two changes of variable, composed. A logistic squashes `ℝ` onto `(0, 1)`, and
`(0,1)` is mapped affinely onto `[lo, hi]` — in `log` of the parameter when it is
a rate, linearly when it is an exponent. Rates span decades and must be searched
multiplicatively; exponents are O(1) and must not.

The point of folding the box into the transform rather than declaring it to the
solver is that *every* optimizer then works unconstrained, derivative-free ones
included, and no iterate can ever leave the physically admissible region. The
price is that a parameter running to a bound shows up as a large `|z|` instead of
an active constraint, which is why [`box_position`](@ref) exists.
"""
function to_native(z, spec = CALIB_SPEC)
    return map(zip(spec, z)) do (p, zi)
        u = 1 / (1 + exp(-zi))
        lo, hi = p.logscale ? (log(p.lo), log(p.hi)) : (p.lo, p.hi)
        w = lo + (hi - lo) * u
        # `clamp` is not decoration. At large |z| the logistic saturates to exactly
        # 1.0, and `exp(lo + (hi - lo))` is then not bit-identical to `p.hi`: the
        # box top came out 4e-15 high on `k₃_C3S`. Harmless as a rate, but the
        # invariant this transform exists to guarantee is that no iterate ever
        # leaves the admissible region, and an invariant that holds to 4e-15 is
        # not an invariant. The test suite caught it.
        return p.logscale ? clamp(exp(w), p.lo, p.hi) : clamp(w, p.lo, p.hi)
    end
end

@doc (@doc to_native)
function to_search(θ, spec = CALIB_SPEC)
    return map(zip(spec, θ)) do (p, v)
        lo, hi = p.logscale ? (log(p.lo), log(p.hi)) : (p.lo, p.hi)
        w = p.logscale ? log(v) : v
        u = clamp((w - lo) / (hi - lo), 1.0e-6, 1 - 1.0e-6)
        return log(u / (1 - u))
    end
end

"""
    box_position(θ, spec = CALIB_SPEC) -> Vector{Float64}

Where each parameter sits in its box, as a fraction in `[0, 1]`.

A value that has run to 0 or 1 is not a result: it says the data pushed the
parameter as far as it was allowed to go, and the box — not the measurement —
decided the answer. Reported alongside every fit for that reason.
"""
function box_position(θ, spec = CALIB_SPEC)
    return map(zip(spec, θ)) do (p, v)
        lo, hi = p.logscale ? (log(p.lo), log(p.hi)) : (p.lo, p.hi)
        return ((p.logscale ? log(v) : v) - lo) / (hi - lo)
    end
end

"""
    prior_vector(spec = CALIB_SPEC) -> Vector{Float64}

The published values, in native units, in the order of `spec`.
"""
prior_vector(spec = CALIB_SPEC) = [p.prior for p in spec]

"""
    apply_parameters(θ, spec = CALIB_SPEC) -> Dict{String, NamedTuple}

The Parrott–Killoh parameter sets with `θ` written into them, everything else
left at its published value.

Rate fields are reunited with their day⁻¹ unit here, which is where they came
from: the search vector carries bare numbers so the optimizer never sees a
`Quantity`.
"""
function apply_parameters(θ, spec = CALIB_SPEC)
    out = Dict{String, NamedTuple}(k => v for (k, v) in PK84_PRIOR)
    for (p, v) in zip(spec, θ)
        p.phase == INDUCTION_KEY && continue
        unitful = startswith(String(p.field), "k") ? v * u"1/d" : v
        out[p.phase] = merge(out[p.phase], NamedTuple{(p.field,)}((unitful,)))
    end
    return out
end

"""
    induction_of(θ, spec = CALIB_SPEC) -> Union{Function, Nothing}

The induction factor `θ` asks for, or `nothing` when `spec` holds no induction
parameters — in which case the model has no induction period, as published.
"""
function induction_of(θ, spec = CALIB_SPEC)
    τ, m = _induction_params(θ, spec)
    τ === nothing && return nothing
    return induction_factor(τ, something(m, 2.0))
end

"""
    induction_tstops_of(θ, spec = CALIB_SPEC) -> Vector{Float64}

The `tstops` of [`induction_tstops`](@ref) for the induction period `θ` asks for,
or empty when it asks for none.
"""
function induction_tstops_of(θ, spec = CALIB_SPEC)
    τ, m = _induction_params(θ, spec)
    τ === nothing && return Float64[]
    return induction_tstops(τ, something(m, 2.0))
end

# The (τ, m) pair a parameter vector carries, or (nothing, nothing).
function _induction_params(θ, spec)
    τ = nothing
    m = nothing
    for (p, v) in zip(spec, θ)
        p.phase == INDUCTION_KEY || continue
        p.field === :τ && (τ = v)
        p.field === :m && (m = v)
    end
    return τ, m
end

"""
    SURROGATE_SOLVER

The integrator used by [`stoichiometric_run`](@ref), chosen by measurement rather
than inherited — the documentation page benchmarks the alternatives.

Measured on the target record, all seven candidates agree on `Q(t_end)` to better
than 0.01 J/g, and the surrogate turns out **not to be stiff**:

| solver | no induction | with induction |
|---|---|---|
| `Tsit5` (explicit) | 0.023 s | 0.023 s |
| `FBDF` (multistep) | 0.023 s | 0.025 s |
| `AutoTsit5(Rosenbrock23)` | 0.023 s | 0.023 s |
| `AutoVern7(Rodas5P)` | 0.025 s | 0.027 s |
| `Rodas5P` | 0.036 s | 0.040 s |
| `Rosenbrock23` | 0.061 s | 0.369 s |

`Rodas5P`, which the other scripts in this directory use, costs 1.6× here for
nothing, and `Rosenbrock23` collapses by a factor of ten once the induction
transition is present. `AutoTsit5(Rosenbrock23())` matches the fastest: SciML's
automatic-switching wrapper runs the explicit `Tsit5` while the problem is
non-stiff and hands over to the stiff `Rosenbrock23` if the detector fires. It
costs nothing when the problem is easy and protects against a parameter set that
makes it hard — which, in an optimizer that will visit hundreds of parameter sets,
is the property that matters.
"""
const SURROGATE_SOLVER = AutoTsit5(Rosenbrock23())

"""
    COUPLED_SOLVER

The integrator used by the coupled mode.

!!! warning "This one is inherited, not measured"
    The surrogate table above is measured; this is not. The coupled mode is a
    different regime — a Gibbs minimization runs on every accepted step, so the
    cost per step dwarfs the cost of taking it, and the right choice is whatever
    needs the fewest steps, which is not what the surrogate benchmark answers.
    `Rodas5P()` is kept because it is what the rest of this repository uses and is
    known to work here, not because it was shown to be best. Measuring it costs one
    coupled solve per candidate, a few minutes, and is worth doing; `forward_Q`
    takes an `ode_solver` keyword precisely so that it can be.

See [`SURROGATE_SOLVER`](@ref).
"""
const COUPLED_SOLVER = Rodas5P()

# ── the forward map ───────────────────────────────────────────────────────────

"""
    STOICH_SPECIES

Species of the fast surrogate model of [`stoichiometric_run`](@ref).
"""
const STOICH_SPECIES = split(
    "C3S C2S C3A C4AF Gp Cal " *
        "Portlandite Jennite ettringite C3AH6 C3FH6 H2O@"
)

const _STOICH_SYSTEM = Ref{Any}(nothing)

"""
    stoichiometric_system() -> ChemicalSystem

The surrogate's chemical system, built once and cached — an optimizer that
rebuilds it per iteration spends most of its time reading JSON.
"""
function stoichiometric_system()
    if _STOICH_SYSTEM[] === nothing
        subs = build_species(datapath("cemdata18-thermofun.json"))
        sp = speciation(subs, STOICH_SPECIES; aggregate_state = [AS_AQUEOUS])
        _STOICH_SYSTEM[] = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    end
    return _STOICH_SYSTEM[]
end

"""
    stoichiometric_run(pk; wb, blaine, clinker, gypsum, filler, tend, ...) -> NamedTuple

The **surrogate**: the same Parrott–Killoh rates, but with the hydrate assemblage
written down by hand instead of chosen by a Gibbs minimization, and therefore no
equilibrium solve at all.

    C₃S  + 103/30 H₂O             → Jennite + 4/3 Portlandite
    C₂S  +  73/30 H₂O             → Jennite + 1/3 Portlandite
    C₃A  + 6 H₂O                  → C₃AH₆
    C₄AF + 2 Portlandite + 10 H₂O → C₃AH₆ + C₃FH₆

Gypsum and limestone are inert here.

!!! warning "This is a surrogate, and its absolute heat is wrong"
    Imposing the assemblage is exactly what `ionic_hydration.jl` exists to avoid,
    and the price shows up immediately in the released heat. Two reasons, both
    structural rather than a matter of tuning:

    1. **The aluminate cannot go to ettringite.** Writing `C₃A + 3 Gp + 26 H₂O →
       ettringite` is stoichiometrically impossible here: 11 % C₃A needs about
       1.1 mol of gypsum per kilogram of binder and 4.6 % gypsum supplies 0.27,
       so a rate law that reads only C₃A drives the gypsum **negative** and the
       enthalpy with it. The aluminate is therefore sent to hydrogarnet, which
       real cement only reaches after the sulfate is gone.
    2. **Jennite plus portlandite is not what a cement precipitates.** Its
       enthalpy drop per mole of alite is larger than that of the assemblage a
       Gibbs minimization selects, so the surrogate overshoots the measured heat
       several-fold. That gap is measured in the documentation page, and it is the
       quantitative argument for the coupled model.

    So the surrogate is **not** used to calibrate. It is used where it is honest
    and where being 700 times cheaper matters: exploring the shape of the loss and
    the sensitivity structure at a cost the coupled model cannot match.

Returns `(; cs, state0, kp, sol)`.
"""
function stoichiometric_run(
        pk;
        wb, blaine, clinker = CALIB_CLINKER,
        gypsum = CALIB_GYPSUM, filler = CALIB_FILLER,
        tend, T = 293.15, binder_mass = 1.0,
        reltol = 1.0e-7, abstol = 1.0e-10, induction = nothing,
        ode_solver = SURROGATE_SOLVER, tstops = Float64[], solver_kwargs...,
    )
    cs = stoichiometric_system()
    α_max = powers_alpha_max(wb)
    f_clinker = 1.0 - gypsum - filler

    state0 = ChemicalState(cs; T = T * u"K")
    for (nm, w) in pairs(clinker)
        set_quantity!(state0, string(nm), (binder_mass * f_clinker * w)u"kg")
    end
    set_quantity!(state0, "Gp", (binder_mass * gypsum)u"kg")
    set_quantity!(state0, "Cal", (binder_mass * filler)u"kg")
    set_quantity!(state0, "H2O@", (binder_mass * wb)u"kg")

    s(name) = cs[name]
    rates = Dict(
        ph => damped_rate(
                parrot_killoh_avrami(pk[ph], ph; α_max, blaine = blaine * u"m^2/kg"),
                ph in INDUCTION_PHASES ? induction : nothing,
            ) for ph in ("C3S", "C2S", "C3A", "C4AF")
    )

    specs = (
        (
            "C3S", OrderedDict(s("C3S") => 1.0, s("H2O@") => 103 / 30),
            OrderedDict(s("Jennite") => 1.0, s("Portlandite") => 4 / 3), "C₃S hydration",
        ),
        (
            "C2S", OrderedDict(s("C2S") => 1.0, s("H2O@") => 73 / 30),
            OrderedDict(s("Jennite") => 1.0, s("Portlandite") => 1 / 3), "C₂S hydration",
        ),
        (
            "C3A", OrderedDict(s("C3A") => 1.0, s("H2O@") => 6.0),
            OrderedDict(s("C3AH6") => 1.0), "C₃A hydration to hydrogarnet",
        ),
        (
            "C4AF", OrderedDict(s("C4AF") => 1.0, s("Portlandite") => 2.0, s("H2O@") => 10.0),
            OrderedDict(s("C3AH6") => 1.0, s("C3FH6") => 1.0), "C₄AF hydration",
        ),
    )

    rxns = map(specs) do (ph, reactants, products, label)
        r = Reaction(reactants, products; symbol = label)
        r[:rate] = rates[ph]
        return r
    end

    kp = KineticsProblem(
        cs, collect(rxns), state0, (0.0, tend); equilibrium_solver = nothing,
    )
    ks = KineticsSolver(; ode_solver, reltol, abstol, tstops, solver_kwargs...)
    return (; cs, state0, kp, sol = integrate(kp, ks))
end

"""
    forward_Q(θ, data; mode = :surrogate, spec = CALIB_SPEC, kwargs...)
        -> Vector{Float64}

Released heat [J/g of binder] at `data.t`, referred to `data.t[1]` exactly as the
record is.

`mode` is `:surrogate` (fast, imposed assemblage) or `:coupled` (the model of
`ionic_hydration.jl`: dissolution plus one Gibbs minimization per accepted step,
read back through a certified replay). Both go through
[`heat_release`](@ref) — never `cumulative_heat` — because under partial
equilibrium the hydrates are precipitated by the minimization and only the
enthalpy difference between certified states sees their heat.

Returns a vector of `NaN` if the solve fails, which [`calorimetry_loss`](@ref)
turns into a penalty rather than an exception.
"""
function forward_Q(
        θ, data; mode::Symbol = :surrogate, spec = CALIB_SPEC,
        ode_solver = nothing, kwargs...,
    )
    pk = apply_parameters(θ, spec)
    induction = induction_of(θ, spec)
    solver = something(
        ode_solver, mode === :surrogate ? SURROGATE_SOLVER : COUPLED_SOLVER
    )
    tend = data.t[end]
    run = try
        if mode === :surrogate
            stoichiometric_run(
                pk; wb = data.meta.wb, blaine = data.meta.blaine,
                tend, T = data.meta.T, induction,
                ode_solver = solver, kwargs...
            )
        elseif mode === :coupled
            run_ionic_hydration(;
                wb = data.meta.wb, clinker = CALIB_CLINKER,
                gypsum = CALIB_GYPSUM, filler = CALIB_FILLER,
                blaine = data.meta.blaine * u"m^2/kg", tend,
                pk_params = pk, induction,
                induction_phases = INDUCTION_PHASES,
                ode_solver = solver, kwargs...
            )
        else
            throw(ArgumentError("mode must be :surrogate or :coupled, got :$mode"))
        end
    catch err
        err isa ArgumentError && rethrow()
        return fill(NaN, length(data.t))
    end

    SciMLBase.successful_retcode(run.sol) || return fill(NaN, length(data.t))
    _, Q, _ = heat_release(run.sol, run.kp; times = data.t)
    # heat_release works on the 1 kg binder basis of the state; the records are
    # normalized per gram.
    return Q ./ 1000
end

"""
    forward_curves(θ, data; mode, spec, kwargs...) -> (Q, qdot)

Released heat [J/g] **and** heat flow [W/g] at `data.t`, from one forward solve.

`heat_release` returns both — `(times, Q, q̇)`.

!!! warning "Adding q̇ to the loss was tried, and it makes things worse"
    The reasoning looked sound: cumulative heat integrates the position of the
    acceleration peak away, so a loss on `Q` alone discards the timing information
    a kinetic calibration most needs. Measured on the target record it is false.
    Fitting `Q + q̇` rather than `Q` alone took the residual on `Q` from 24.4 to
    27.6 J/g and the `τ_ind`–`k₁_C3S` correlation from 0.982 to **0.994**.

    The cause is a property of `q̇` here: it is a **centered finite difference of
    `Q` over the output grid**, and that grid is log-spaced over two and a half
    decades, so its spacing near the peak is coarse. What the loss gained was a
    smeared, low-information version of the heat flow — noise rather than signal.

    **Repaired, and it still does not help — because the cause was elsewhere.**
    `calorimetry_loss(...; w_flow > 0)` now differences model and measurement with
    the same operator on the same grid ([`grid_slope`](@ref)), so the discretization
    cancels instead of being fitted. Measured: 27.6 J/g and a correlation of 0.994,
    unchanged.

    Two further experiments settle why. First, `grid_slope(Q, t)` is a **linear map
    on the same numbers**, so a residual on the differenced curve is a linear
    combination of the residuals on `Q` — it re-weights information rather than
    adding any. Second, the genuinely independent quantity, the *instrument's* `q̇`
    at 77 s sampling compared on a grid four times finer near the peak (0.24 h
    against 0.91 h), behaves identically: 27.8 J/g and 0.994.

    So the `τ_ind`–`k₁_C3S` degeneracy is **structural in the rate law**, not an
    artifact of the observable or its discretization. Both parameters act on the
    same feature — when, and how fast, the acceleration happens — and calorimetry
    cannot separate them.

    What the flow residual *does* buy is the dormancy timescale in absolute terms:
    every flow-informed fit puts `τ` at 5.2-5.4 h against the 3.1 h of a `Q`-only
    surrogate fit, and the coupled `Q`-only fit independently gives 5.6 h. Three
    routes agreeing, at the price of the `Q` residual. That is why `w_flow` defaults
    to 0 for fitting, and why the flow is worth carrying anyway.

Returns `(NaN…, NaN…)` on a failed solve, which [`calorimetry_loss`](@ref) turns
into a penalty.
"""
function forward_curves(
        θ, data; mode::Symbol = :surrogate, spec = CALIB_SPEC,
        ode_solver = nothing, kwargs...,
    )
    pk = apply_parameters(θ, spec)
    induction = induction_of(θ, spec)
    solver = something(
        ode_solver, mode === :surrogate ? SURROGATE_SOLVER : COUPLED_SOLVER
    )
    nan = fill(NaN, length(data.t))
    run = try
        if mode === :surrogate
            stoichiometric_run(
                pk; wb = data.meta.wb, blaine = data.meta.blaine,
                tend = data.t[end], T = data.meta.T, induction,
                ode_solver = solver, kwargs...,
            )
        elseif mode === :coupled
            run_ionic_hydration(;
                wb = data.meta.wb, clinker = CALIB_CLINKER,
                gypsum = CALIB_GYPSUM, filler = CALIB_FILLER,
                blaine = data.meta.blaine * u"m^2/kg", tend = data.t[end],
                pk_params = pk, induction,
                induction_phases = INDUCTION_PHASES,
                ode_solver = solver, kwargs...,
            )
        else
            throw(ArgumentError("mode must be :surrogate or :coupled, got :$mode"))
        end
    catch err
        err isa ArgumentError && rethrow()
        return (nan, nan)
    end
    SciMLBase.successful_retcode(run.sol) || return (nan, nan)
    _, Q, qdot = heat_release(run.sol, run.kp; times = data.t)
    return (Q ./ 1000, qdot ./ 1000)
end

"""
    grid_slope(y, t) -> Vector{Float64}

Centered difference `dy/dt` at the interior points of `t`, with one-sided
differences at the ends. Length `length(t)`.

# Why the loss differences the measurement too

The first attempt at a heat-flow residual compared the model's `q̇` — which
[`heat_release`](@ref) produces by centered-differencing `Q` over the **output
grid** — against `data.qdot`, the instrument's own reading sampled every 77 s. On a
log grid spanning two and a half decades, consecutive points near the 10-hour peak
are *hours* apart, so the model's `q̇` there is a smeared average of a feature only
a few hours wide, while the measurement resolves it. The residual between the two
was then dominated by the grid spacing rather than by the parameters, and the
optimizer chased discretization error: the fit on `Q` got worse (24.4 → 27.6 J/g)
and the `τ_ind`–`k₁_C3S` correlation tightened (0.982 → 0.994).

Making the model's grid dense enough to resolve the peak would fix it and multiply
the coupled cost by the refinement factor, which is not affordable. Applying the
**same operator to both sides** fixes it for nothing: difference the measured `Q` on
the same grid, and the leading discretization error cancels. What is left is a
genuine constraint on local slope — which is the timing information the cumulative
heat integrates away — expressed in a quantity both sides possess exactly.

This is the observation-operator principle: push the model through whatever the
comparison does to it, rather than trying to undo it on the data.
"""
function grid_slope(y, t)
    n = length(t)
    s = similar(float.(y))
    n < 2 && return fill!(s, zero(eltype(s)))
    s[1] = (y[2] - y[1]) / (t[2] - t[1])
    s[n] = (y[n] - y[n - 1]) / (t[n] - t[n - 1])
    for i in 2:(n - 1)
        s[i] = (y[i + 1] - y[i - 1]) / (t[i + 1] - t[i - 1])
    end
    return s
end

"""
    peak_time(data) -> Float64

Time [s] of the maximum measured heat flow, after the aluminate shoulder.

The search starts at 2 h so the sulfate shoulder, which is a different peak, cannot
be mistaken for the silicate one. On the target record it returns 10.6 h.

!!! warning "Pinning `τ_ind` from this was tried, and it fails"
    The hope was to break the 0.994 correlation between `τ_ind` and `k₁_C3S` by
    determining `τ` from the peak instead of fitting it. Taking
    `τ = 0.55 × t_peak` = 5.8 h and fitting the remaining four parameters gave a
    residual of 28.9 J/g against 24.4 for the five-parameter fit, with `m_ind`
    driven to its floor, `k₃_C3S` to its ceiling, and a condition number of 1e18 —
    numerically singular. The model cannot reach the data with `τ` held there, so
    the other parameters saturate trying.

    Two things were wrong. The factor 0.55 was chosen by comparing against the `τ`
    of the **coupled** fit while the surrogate's own optimum wants 3–5 h, so the
    apparent agreement was between two different models and established nothing.
    And a defensible version would use no factor at all: it would compute where the
    model's own heat-flow peak falls as a function of `τ` and invert that — a
    one-parameter numerical inversion, which has not been done.

    Kept because the peak time is a useful diagnostic to *report*, and because the
    negative result should stay reproducible.
"""
function peak_time(data)
    ok = findall(>=(2 * 3600), data.t)
    isempty(ok) && return data.t[argmax(data.qdot)]
    return data.t[ok[argmax(data.qdot[ok])]]
end

# ── the loss ──────────────────────────────────────────────────────────────────

"""
    LOSS_PENALTY

Value returned for a failed forward solve. Large enough that no optimizer walks
towards it, finite so none of them dies on it.
"""
const LOSS_PENALTY = 1.0e4

"""
    calorimetry_loss(θ, data; mode, spec, shape = false) -> Float64

Root-mean-square residual on the released heat, in **J/g of binder**.

Absolute and unweighted, on purpose. Isothermal calorimetry has an approximately
absolute error — the depositors of these data quote a database-wide 4.11 J/g
standard deviation for their own four-parameter fit, and the international round
robin of [Wadso2016](@cite) puts inter-laboratory scatter at a comparable
scale — so an absolute residual matches the instrument, and the loss can be read
directly against those figures instead of being a dimensionless score. The time
weighting lives in the log-spaced grid of [`resample_log`](@ref), not in the norm.

# `shape = true`

Compare `Q(t)/Q(t_end)` instead of `Q(t)`, rescaled by the measured `Q(t_end)` so
the number stays in J/g and stays comparable.

This exists for the surrogate, and only for it. The surrogate's assemblage is
imposed, so its *level* is wrong by construction — some 9 % low on this record —
and a fit on the absolute heat would spend its parameters compensating for a bias
no rate constant is responsible for. On the normalized curve it is judged on what
it can actually get right, the timing, which is exactly what makes it a useful
starting point for the coupled search.
"""
function calorimetry_loss(
        θ, data; mode::Symbol = :surrogate, spec = CALIB_SPEC, shape::Bool = false,
        w_flow::Real = 0.0,
    )
    if w_flow > 0
        Qm = forward_Q(θ, data; mode, spec)
        (any(isnan, Qm) || iszero(Qm[end])) && return LOSS_PENALTY
        rQ = mean(abs2, Qm .- data.Q)
        # The slope block compares like with like: the SAME centered difference is
        # applied to the model and to the measurement, on the SAME grid. See
        # `grid_slope` for why that matters and why comparing against `data.qdot`
        # instead is the trap.
        sm = grid_slope(Qm, data.t)
        sd = grid_slope(data.Q, data.t)
        scale = maximum(abs, sd)
        iszero(scale) && return sqrt(rQ)
        rq = mean(abs2, (sm .- sd) ./ scale) * mean(abs2, data.Q)
        return sqrt((rQ + w_flow * rq) / (1 + w_flow))
    end
    Qm = forward_Q(θ, data; mode, spec)
    any(isnan, Qm) && return LOSS_PENALTY
    iszero(Qm[end]) && return LOSS_PENALTY
    shape || return sqrt(mean(abs2, Qm .- data.Q))
    return data.Q[end] * sqrt(mean(abs2, Qm ./ Qm[end] .- data.Q ./ data.Q[end]))
end

# ── calibration ───────────────────────────────────────────────────────────────

"""
    CalibrationResult

What a fit returns: the fitted parameters `θ` in native units, the `loss` reached
(RMSE in J/g of binder), the `loss0` at the starting point, the `trace` of every
objective evaluation, the `box` positions of [`box_position`](@ref), the number
of forward solves `nevals`, and the wall-clock `seconds`.

`trace` is kept because it is the only honest evidence that a derivative-free
search converged rather than stopped: the documentation page plots it.
"""
struct CalibrationResult
    θ::Vector{Float64}
    loss::Float64
    loss0::Float64
    trace::Vector{Float64}
    box::Vector{Float64}
    nevals::Int
    seconds::Float64
    mode::Symbol
    spec::Vector{CalibParameter}
end

function Base.show(io::IO, r::CalibrationResult)
    println(
        io, "CalibrationResult(:", r.mode, ")  RMSE ", round(r.loss; digits = 3),
        " J/g  (from ", round(r.loss0; digits = 3), ")"
    )
    println(io, "  ", r.nevals, " forward solves in ", round(r.seconds; digits = 1), " s")
    for (p, v, b) in zip(r.spec, r.θ, r.box)
        flag = (b < 0.02 || b > 0.98) ? "  ← AT A BOUND, not identified" : ""
        @printf(
            io, "  %-8s %10.4g   (prior %8.4g, %5.1f %% of box)%s\n",
            p.name, v, p.prior, 100b, flag
        )
    end
    return nothing
end

"""
    calibrate(data; mode, θ0, spec, optimizer, maxiters) -> CalibrationResult

Fit `spec` to `data` by minimizing [`calorimetry_loss`](@ref) through
`Optimization.jl`.

The search is **unconstrained** in the squashed coordinates of
[`to_native`](@ref), so a derivative-free method needs no box support and cannot
leave the admissible region either way.

# Why not a gradient

`Optimization.jl` would happily take `AutoForwardDiff()` here and it would not
work: the kinetics core is AD-clean in the ODE state and in time — deliberately,
because `Rodas5P` needs the time gradient — but not in the parameters.
`build_u0` returns a `Vector{Float64}`; `build_kinetics_params` casts the
temperature, the initial amounts, the stoichiometry and the calorimeter constants
to `Float64`; and under partial equilibrium `respeciate!` is `Float64`-only by
construction. So a dual number baked into a rate closure cannot reach the
integrator. `AutoFiniteDiff()` is the honest backend until that changes, and with
six parameters it costs seven forward solves per gradient — which is why
`NelderMead` is the default for the expensive mode.
"""
function calibrate(
        data;
        mode::Symbol = :surrogate,
        θ0 = prior_vector(),
        spec = CALIB_SPEC,
        optimizer = NelderMead(),
        maxiters = 400,
        shape::Bool = false,
        w_flow::Real = 0.0,
    )
    trace = Float64[]
    z0 = to_search(θ0, spec)
    loss0 = calorimetry_loss(θ0, data; mode, spec, shape, w_flow)

    objective = function (z, _p)
        L = calorimetry_loss(to_native(z, spec), data; mode, spec, shape, w_flow)
        push!(trace, L)
        return L
    end

    t0 = time()
    f = OptimizationFunction(objective, AutoFiniteDiff())
    sol = solve(OptimizationProblem(f, z0, nothing), optimizer; maxiters)
    seconds = time() - t0

    θ̂ = to_native(sol.u, spec)
    return CalibrationResult(
        θ̂, calorimetry_loss(θ̂, data; mode, spec, shape, w_flow), loss0,
        trace, box_position(θ̂, spec), length(trace), seconds, mode, collect(spec),
    )
end

"""
    CEM_I_DEPOSIT_RECORDS

The fourteen plain CEM I records of the Šmilauer & Reiterman deposit, as the file
stems `data/experimental/regenerate.jl` would need in its `RECORDS` list.

Only two are vendored — the target and the holdout. The rest are one edit away, and
they are the raw material for the joint fit of
[`multirecord_loss`](@ref).

| file stem | class | w/b | Blaine (m²/kg) | duration (h) |
|---|---|---|---|---|
| `100-CEM I 32.5 R Ozarow-330` | 32.5 R | 0.45 | 330 | — |
| `103-CEM I 42.5 Rsc Mokra-256` | 42.5 Rsc | 0.45 | 256 | 357 |
| `105-CEM I 42.5 Rsc Mokra-306` | 42.5 Rsc | 0.45 | 306 | 307 |
| `106-CEM I 42.5 N Mokra-264` | 42.5 N | 0.45 | 264 | 307 |
| `109-CEM I 42.5 Rsc Mokra-320` | 42.5 Rsc | 0.45 | 320 | 378 |
| `114-CEM I 32.5 R Ladce-250` | 32.5 R | 0.45 | 250 | 617 |
| `115-CEM I 42.5 R Ladce-339` | 42.5 R | 0.45 | 339 | 617 |
| `116-CEM I 52.5 R Ladce-415` | 52.5 R | 0.45 | 415 | 617 |
| `122-CEM I 52.5 R Cizkovice-397` | 52.5 R | 0.50 | 397 | 262 |
| `158-CEM I 42.5 Rsc Mokra-312` | 42.5 Rsc | 0.40 | 312 | 407 |
| `163-CEM I 42.5 Rsc Hranice-340` | 42.5 Rsc | 0.45 | 440 | 384 |
| `172-CEM I 42.5 R Rohoznik-390` | 42.5 R | 0.40 | 390 | 423 |
| `179-CEM I 42.5 Rsc Mokra-263` | 42.5 Rsc | 0.40 | 263 | 359 |
| `202-CEM I 42.5 R Mokra-409` | 42.5 R | 0.40 | 409 | 862 |

!!! warning "They are not one cement measured fourteen times"
    Five plants, four strength classes, w/b from 0.40 to 0.50 and fineness from
    250 to 440 m²/kg. `blaine_factor` and `powers_alpha_max` are supposed to carry
    a rate law across the last two; the first two it cannot carry, because the
    clinker composition differs and the deposit does not report it. A joint fit
    therefore estimates **one** kinetic parameter set for **fourteen different
    clinkers**, and whatever composition scatter exists gets absorbed into the
    kinetics. That is a weaker assumption than it sounds only because the
    alternative — one parameter set per record — is what already failed to
    transfer.
"""
const CEM_I_DEPOSIT_RECORDS = [
    "100-CEM I 32.5 R Ozarow-330", "103-CEM I 42.5 Rsc Mokra-256",
    "105-CEM I 42.5 Rsc Mokra-306", "106-CEM I 42.5 N Mokra-264",
    "109-CEM I 42.5 Rsc Mokra-320", "114-CEM I 32.5 R Ladce-250",
    "115-CEM I 42.5 R Ladce-339", "116-CEM I 52.5 R Ladce-415",
    "122-CEM I 52.5 R Cizkovice-397", "158-CEM I 42.5 Rsc Mokra-312",
    "163-CEM I 42.5 Rsc Hranice-340", "172-CEM I 42.5 R Rohoznik-390",
    "179-CEM I 42.5 Rsc Mokra-263", "202-CEM I 42.5 R Mokra-409",
]

"""
    multirecord_loss(θ, records; mode, spec, weights) -> Float64

Root-mean-square residual over **several** records at once, in J/g of binder.

One parameter set, many cements. Each record contributes its own forward solve at
its own w/b, fineness and temperature — all three read from its own metadata — and
the residuals are pooled. `weights` defaults to one per record; pass the reciprocal
of a record's variance if some are trusted more than others.

# Why this is the fix for the failure of §7 of the documentation page

Fitting one record gave five parameters a residual to absorb that was specific to
that cement, and the absorbed part did not transfer: the calibrated set was worse
than the published one on a record it had not seen. A joint objective removes the
opportunity. The parameter count stays at five while the data multiplies, so the
`τ_ind`–`k₁_C3S` degeneracy is constrained by fourteen different peak positions
rather than one, and transferability is what is being *optimized* rather than what
is being *tested afterwards*.

**This has not been run.** See [`multirecord_recipe`](@ref) for what it costs and
how to launch it.
"""
function multirecord_loss(
        θ, records; mode::Symbol = :coupled, spec = CALIB_SPEC,
        weights = ones(length(records)),
    )
    total = 0.0
    wsum = 0.0
    for (w, d) in zip(weights, records)
        Qm = forward_Q(θ, d; mode, spec)
        (any(isnan, Qm) || iszero(Qm[end])) && return LOSS_PENALTY
        total += w * mean(abs2, Qm .- d.Q)
        wsum += w
    end
    return sqrt(total / wsum)
end

"""
    calibrate_multirecord(records; mode, θ0, spec, optimizer, maxiters, weights)
        -> CalibrationResult

[`calibrate`](@ref) against [`multirecord_loss`](@ref) instead of a single record.

**Not run in this repository.** [`multirecord_recipe`](@ref) prints what it would
cost and the four lines that launch it.
"""
function calibrate_multirecord(
        records; mode::Symbol = :coupled, θ0 = prior_vector(), spec = CALIB_SPEC,
        optimizer = NelderMead(), maxiters = 120, weights = ones(length(records)),
    )
    trace = Float64[]
    z0 = to_search(θ0, spec)
    loss0 = multirecord_loss(θ0, records; mode, spec, weights)
    objective = function (z, _p)
        L = multirecord_loss(to_native(z, spec), records; mode, spec, weights)
        push!(trace, L)
        return L
    end
    t0 = time()
    f = OptimizationFunction(objective, AutoFiniteDiff())
    sol = solve(OptimizationProblem(f, z0, nothing), optimizer; maxiters)
    seconds = time() - t0
    θ̂ = to_native(sol.u, spec)
    return CalibrationResult(
        θ̂, multirecord_loss(θ̂, records; mode, spec, weights), loss0,
        trace, box_position(θ̂, spec), length(trace), seconds, mode, collect(spec),
    )
end

"""
    multirecord_recipe(; n_records = 14, n_instants = 20, maxiters = 120)

Print, without running anything, what the joint fit costs and exactly how to
launch it.

The estimate is anchored on a measurement rather than a guess: on the machine used
here — two cores — one coupled forward solve on 20 log-spaced instants took about
35 s, and the single-record fit of [`CALIBRATED_THETA`](@ref) needed 152 objective
evaluations for five parameters.
"""
function multirecord_recipe(; n_records = 14, n_instants = 20, maxiters = 120)
    per_solve = 35.0
    evals = 1.7 * maxiters                      # NelderMead, measured ratio
    hours = n_records * per_solve * evals / 3600
    println("Joint fit over $n_records records, $n_instants instants, $maxiters iterations")
    @printf("  one solve      ≈ %.0f s   (measured, two cores)\n", per_solve)
    @printf("  one evaluation ≈ %.0f s   (%d records)\n", n_records * per_solve, n_records)
    @printf("  ≈%.0f evaluations  ⇒  ≈ %.1f h wall clock\n", evals, hours)
    println()
    println("Step 1 — vendor the records. `data/experimental/regenerate.jl` holds a")
    println("  `RECORDS` list of (source stem, output stem, role) triples; extend it")
    println("  with the stems of `CEM_I_DEPOSIT_RECORDS` and run:")
    println("      julia data/experimental/regenerate.jl")
    println("  It downloads the CC-BY-4.0 deposit once and writes one CSV per record.")
    println()
    println("Step 2 — load them and launch. Four lines:")
    println("      records = [resample_log(read_calorimetry(p), N_RESIDUALS_COUPLED)")
    println("                 for p in readdir(CALORIMETRY_DIR; join = true)")
    println("                 if endswith(p, \".csv\")]")
    println("      r = calibrate_multirecord(records; maxiters = $maxiters)")
    println()
    println("Step 3 — read the result the way §5 of the page reads a fit:")
    println("      report_fit(r, records[1]; label = \"joint\")")
    println("      id = local_identifiability(r.θ, records[1]; mode = :coupled)")
    println("  The number to watch is not the residual. It is whether")
    println("  corr(τ_ind, k₁_C3S) has come down from 0.994 — that is what fourteen")
    println("  different peak positions are supposed to buy, and if it has not, the")
    println("  degeneracy is structural and a second observable is the only way out.")
    println()
    println("Run it detached; it outlives a terminal session:")
    println("      nohup julia --project=scripts -e 'include(\"scripts/hydration_calibration.jl\");")
    println("            multirecord_recipe()' > joint.log 2>&1 &")
    return nothing
end

# ── identifiability ───────────────────────────────────────────────────────────

"""
    sensitivity_matrix(θ, data; mode, spec, relstep) -> Matrix{Float64}

`∂Q/∂log θⱼ` at each instant of `data`, by central differences — the columns are
in units of J/g per *relative* change in the parameter, so they are comparable
across a rate and an exponent.

Costs `2n` forward solves. `spec` may be **wider** than the one that was fitted,
provided `θ` is extended to match — that is how a parameter's sensitivity is
shown to be negligible without spending an optimization on it.
"""
function sensitivity_matrix(θ, data; mode::Symbol = :surrogate, spec = CALIB_SPEC, relstep = 0.05)
    J = Matrix{Float64}(undef, length(data.t), length(spec))
    for j in eachindex(spec)
        θp = copy(collect(float.(θ)));  θp[j] *= (1 + relstep)
        θm = copy(collect(float.(θ)));  θm[j] *= (1 - relstep)
        Qp = forward_Q(θp, data; mode, spec)
        Qm = forward_Q(θm, data; mode, spec)
        J[:, j] = (Qp .- Qm) ./ (2 * relstep)
    end
    return J
end

"""
    local_identifiability(θ, data; mode, spec, relstep)
        -> (; J, U, S, V, cond, correlation, stderr, rmse)

Local identifiability of `θ` from `data`.

`S` holds the singular values of `∂Q/∂log θ`. A ratio `S[1]/S[end]` of a few
means every direction is constrained; several orders of magnitude means the data
determine a *combination* of parameters and not the parameters, and `V[:, end]`
names which combination. `correlation` is the parameter correlation matrix from
`(JᵀJ)⁻¹`, and `stderr` the approximate relative standard errors
`σ√diag((JᵀJ)⁻¹)` with `σ` the residual RMSE.

These are the *linearized* errors at one point, which is all six parameters and a
few hundred forward solves can buy. They are reported to say which numbers in a
fit deserve to be quoted, not as confidence intervals.
"""
function local_identifiability(
        θ, data; mode::Symbol = :surrogate, spec = CALIB_SPEC, relstep = 0.05,
    )
    J = sensitivity_matrix(θ, data; mode, spec, relstep)
    Qm = forward_Q(θ, data; mode, spec)
    rmse = sqrt(mean(abs2, Qm .- data.Q))
    F = svd(J)
    JtJ = J' * J
    C = try
        inv(JtJ)
    catch
        pinv(JtJ)
    end
    d = sqrt.(abs.(diag(C)))
    correlation = C ./ (d * d')
    return (;
        J, U = F.U, S = F.S, V = F.V,
        cond = F.S[1] / max(F.S[end], eps()),
        correlation, stderr = rmse .* d, rmse,
    )
end

# ── driver ────────────────────────────────────────────────────────────────────

"""
    N_RESIDUALS_COUPLED

Number of log-spaced instants the coupled loss is evaluated on.

The coupled forward solve is dominated by the certified replay and the enthalpy
sums, both linear in this number, so the grid is the main cost lever. Twenty
instants over two and a half decades is eight per decade, which resolves every
feature of these curves.
"""
const N_RESIDUALS_COUPLED = 20

"""
    N_RESIDUALS_SURROGATE

As [`N_RESIDUALS_COUPLED`](@ref), for the surrogate, which is cheap enough not to
care.
"""
const N_RESIDUALS_SURROGATE = 60

"""
    COUPLED_ITERS

Iteration budget for the coupled search. `main` measures and prints the cost of a
single solve before spending this many of them.
"""
const COUPLED_ITERS = 60

"""
    CALIBRATED_THETA :: Vector{Float64}

The result of the coupled calibration on [`CEM_I_TARGET`](@ref), in the order of
[`CALIB_SPEC`](@ref).

Stored rather than recomputed because the search costs of the order of two hours,
which does not belong in a documentation build or a test suite. `main()`
reproduces it.

# The run that produced it

152 objective evaluations of `NelderMead` on 20 log-spaced instants, coupled
mode, from the published Parrott–Killoh values with `τ = 4 h`, `m = 2` as the
starting point. The residual went `17.48 → 13.04 J/g`, and every parameter
finished **interior** to its box — at 61, 50, 37, 70 and 86 % of it — so the box
did not choose any of these numbers.

| θ | fitted | published | ratio |
|---|---|---|---|
| `k₁_C3S` | 2.043 d⁻¹ | 1.5 | 1.36 |
| `k₃_C3S` | 1.583 d⁻¹ | 1.1 | 1.44 |
| `k₁_C3A` | 0.687 d⁻¹ | 1.0 | 0.69 |
| `τ_ind` | 20 228 s = 5.6 h | (none published) | — |
| `m_ind` | 3.56 | (none published) | — |

!!! warning "These numbers fit one record and do not transfer to another"
    On the holdout — a nominally same-class CEM I at a different w/b and
    fineness, predicted without refitting — the calibrated parameters are
    **worse** than the published ones: 47.0 against 40.7 J/g. The calibration
    improved the record it was fitted to by 25 % and degraded the one it had not
    seen by 16 %.

    That is not a surprise once the identifiability is read. At the optimum
    `k₁_C3S` and `τ_ind` are correlated at **0.994** — a longer dormant period
    followed by a faster rate makes very nearly the same curve — and the
    approximate relative standard errors are 981 %, 25 %, 123 %, 389 % and 196 %.
    Only `k₃_C3S` is determined to better than a factor of a few. What the fit
    found is one or two combinations plus a target-specific residual, and the
    target-specific part is precisely what does not generalize.

    The remedy is not a better optimizer. It is either fewer parameters — fix
    `τ_ind` from the measured position of the main peak instead of fitting it
    against `k₁_C3S` — or a second observable, which is what §10 of the
    documentation page argues.
"""
const CALIBRATED_THETA = [
    2.0430166816339432, 1.5830829931535122, 0.6872659858977794,
    20228.44148315133, 3.55870720991354,
]

"""
    MEASURED_IDENTIFIABILITY

The two identifiability analyses of the documentation page, stored rather than
recomputed.

Each costs `2n + 1` coupled forward solves — thirteen for the six candidates,
eleven at the optimum — and at tens of seconds apiece that is a quarter of an hour
that has no business inside a documentation build. `main()` recomputes both, and
[`local_identifiability`](@ref) is the one function involved; what is stored here
is only its output.

  - `candidates` — at the published parameters, over all six of
    [`CALIB_SPEC_FULL`](@ref), on 30 log-spaced instants. This is the analysis that
    decided how many parameters to fit.
  - `optimum` — at [`CALIBRATED_THETA`](@ref), over the five of
    [`CALIB_SPEC`](@ref), on 20 instants. This is the analysis that explains why
    the fit does not transfer.

`S` are the singular values of `∂Q/∂log θ`, `V` its right singular vectors as rows
per parameter, `correlation` the parameter correlation matrix from `(JᵀJ)⁻¹`, and
`stderr` the approximate **relative** standard errors.
"""
const MEASURED_IDENTIFIABILITY = (
    candidates = (
        # Spelled out, NOT derived from `CALIB_SPEC_FULL`. Deriving them is how this
        # went wrong once: the stored matrices are 6x6, from when that list held six
        # entries, and adding `τ_ind` and `m_ind` to it later left the names 8 long
        # against 6x6 data — which the documentation build caught as a plotting
        # error rather than as the desynchronization it was. Stored numbers must
        # carry their own labels.
        names = [:k₁_C3S, :n₁_C3S, :k₃_C3S, :n₃_C3S, :k₁_C3A, :k₃_C2S],
        rmse = 26.081, cond = 2101.0,
        S = [421.8, 101.0, 60.1, 6.342, 1.438, 0.2008],
        column_norms = [48.49, 99.92, 207.0, 368.3, 28.89, 12.64],
        stderr = [17.538, 11.706, 3.069, 3.932, 4.394, 129.361],
        V = [
            -0.058  0.335 -0.41 -0.077 -0.84 -0.067
            -0.193  0.436 -0.628 -0.325  0.519  0.053
            0.45  0.765  0.456 -0.034  0.053  0.022
            -0.869  0.271  0.412  0.025 -0.033 -0.03
            -0.032  0.198 -0.245  0.941  0.114  0.0
            0.03  0.01 -0.008 -0.013  0.086 -0.996
        ],
        correlation = [
            1.0 -0.985 -0.738  0.619 -0.433  0.484
            -0.985  1.0  0.806 -0.701  0.288 -0.573
            -0.738  0.806  1.0 -0.977  0.12 -0.941
            0.619 -0.701 -0.977  1.0 -0.063  0.985
            -0.433  0.288  0.12 -0.063  1.0 -0.009
            0.484 -0.573 -0.941  0.985 -0.009  1.0
        ],
    ),
    optimum = (
        names = [:k₁_C3S, :k₃_C3S, :k₁_C3A, :τ_ind, :m_ind],
        rmse = 13.0376, cond = 97.346,
        S = [120.08, 64.597, 13.531, 6.5475, 1.2335],
        stderr = [9.8053, 0.2476, 1.2279, 3.8862, 1.9604],
        V = [
            -0.3135  0.0665 -0.1765  0.0775 -0.9274
            -0.2454 -0.9626  0.1133 -0.0188 -0.0092
            -0.3494  0.1961  0.848 -0.342 -0.0578
            0.8433 -0.1694  0.3547 -0.0174 -0.3661
            0.0909 -0.0437 -0.3332 -0.9362 -0.0487
        ],
        correlation = [
            1.0 0.3839 0.4773 0.9937 0.2504
            0.3839 1.0 0.5516 0.4376 0.178
            0.4773 0.5516 1.0 0.5579 0.5487
            0.9937 0.4376 0.5579 1.0 0.2557
            0.2504 0.178 0.5487 0.2557 1.0
        ],
    ),
)

"""
    report_stored_identifiability(m, label)

Print a stored [`MEASURED_IDENTIFIABILITY`](@ref) entry in the same shape
[`report_identifiability`](@ref) prints a freshly computed one.
"""
function report_stored_identifiability(m, label)
    @printf(
        "%s — %d parameters, residual RMSE %.2f J/g, condition number %.4g\n",
        label, length(m.names), m.rmse, m.cond
    )
    println("singular values: ", join((@sprintf("%.4g", v) for v in m.S), "  "))
    println("approximate relative standard errors:")
    for (n, e) in zip(m.names, m.stderr)
        @printf("   %-10s %8.0f %%\n", n, 100e)
    end
    return nothing
end

"""
    rule(title)

A section heading in the style of the other scripts in this directory.
"""
rule(title) = println("\n── ", title, " ", "─"^max(0, 74 - length(title)))

"""
    report_fit(r, data; label)

Print a fit: the residual against the depositors' own 4.11 J/g yardstick, then
each parameter with its prior, its ratio to the prior, and its position in its
box.
"""
function report_fit(r::CalibrationResult, data; label = "fit")
    @printf(
        "%s (%s): RMSE %.2f → %.2f J/g over %d instants, %d solves, %.1f s\n",
        label, r.mode, r.loss0, r.loss, length(data.t), r.nevals, r.seconds
    )
    println("   parameter      fitted      prior    ratio    box")
    for (p, v, b) in zip(r.spec, r.θ, r.box)
        flag = (b < 0.02 || b > 0.98) ? "  ← at a bound, not identified" : ""
        @printf(
            "   %-10s %10.4g %10.4g %8.2f %5.0f %%%s\n",
            p.name, v, p.prior, v / p.prior, 100b, flag
        )
    end
    return nothing
end

"""
    report_identifiability(id, spec = CALIB_SPEC)

Print the singular spectrum of `∂Q/∂log θ`, the least-constrained direction and
the approximate relative standard errors.
"""
function report_identifiability(id, spec = CALIB_SPEC)
    @printf(
        "residual RMSE %.2f J/g;  condition number of ∂Q/∂log θ = %.3g\n", id.rmse, id.cond
    )
    println("singular values: ", join((@sprintf("%.3g", v) for v in id.S), "  "))
    println("least-constrained direction (right singular vector of the smallest value):")
    for (p, c) in zip(spec, id.V[:, end])
        @printf("   %-10s %+7.3f\n", p.name, c)
    end
    println("approximate relative standard errors:")
    for (p, e) in zip(spec, id.stderr)
        @printf("   %-10s %8.1f %%\n", p.name, 100e)
    end
    return nothing
end

"""
    main()

The whole exercise, in the order the argument runs: look at the data, try the
published parameters, show why the hydrate assemblage cannot be imposed, find out
what the measurement can constrain, calibrate the coupled model, then predict a
record it never saw.

Expect an hour or so. The calibration is a derivative-free search over a model
that costs tens of seconds per evaluation, and `main` prints that cost before
spending it.
"""
function main()
    target = resample_log(CEM_I_TARGET, N_RESIDUALS_COUPLED)
    target_fine = resample_log(CEM_I_TARGET, N_RESIDUALS_SURROGATE)
    holdout = resample_log(CEM_I_HOLDOUT, N_RESIDUALS_COUPLED)
    θ0 = prior_vector()

    rule("the data")
    println("target : ", CEM_I_TARGET)
    println("holdout: ", CEM_I_HOLDOUT)
    @printf(
        "source : %s\n         doi %s, %s\n",
        CEM_I_TARGET.meta.source, CEM_I_TARGET.meta.doi, CEM_I_TARGET.meta.license
    )
    @printf(
        "both records open at ~1 h, with %.1f J/g already gone by then (the\n",
        CEM_I_TARGET.meta.q_before_start
    )
    println("depositors' estimate, not a measurement), so model and record are both")
    println("referred to the first sample rather than to mixing.")

    rule("the published parameters, untouched")
    t0 = time()
    Q_prior = forward_Q(θ0, target; mode = :coupled)
    t_coupled = time() - t0
    @printf("one coupled solve on %d instants: %.1f s\n", length(target.t), t_coupled)
    @printf(
        "Q(%.0f h): %.1f J/g computed, %.1f J/g measured  (%+.1f %%)\n",
        target.t[end] / 3600, Q_prior[end], target.Q[end],
        100 * (Q_prior[end] / target.Q[end] - 1)
    )
    @printf("RMSE over the curve: %.2f J/g\n", calorimetry_loss(θ0, target; mode = :coupled))
    if !isempty(target.Qref)
        println("for comparison, the depositors' own fitted affinity model, on this very")
        @printf(
            "record: Q(end) %.1f J/g (%+.1f %%), RMSE %.2f J/g\n",
            target.Qref[end], 100 * (target.Qref[end] / target.Q[end] - 1),
            sqrt(mean(abs2, target.Qref .- target.Q))
        )
    end
    println("The total heat is already right to a couple of percent, with parameters")
    println("fitted in 1984 to other cements. What a calibration has to earn is the")
    println("timing, and that is where the residual lives.")

    rule("why the assemblage cannot be imposed")
    t0 = time()
    Q_surr = forward_Q(θ0, target_fine; mode = :surrogate)
    t_surr = time() - t0
    @printf("one surrogate solve: %.3f s — %.0f× cheaper\n", t_surr, t_coupled / t_surr)
    println("   instant      surrogate      coupled     measured")
    for (lab, h) in (("6 h", 6), ("24 h", 24), ("end", 0))
        i = h == 0 ? lastindex(target.t) :
            something(findfirst(>=(h * 3600 - 1), target.t), lastindex(target.t))
        j = h == 0 ? lastindex(target_fine.t) :
            something(findfirst(>=(h * 3600 - 1), target_fine.t), lastindex(target_fine.t))
        @printf("   %-8s %11.1f %12.1f %12.1f\n", lab, Q_surr[j], Q_prior[i], target.Q[i])
    end
    println("Imposing Jennite plus portlandite, and sending the aluminate to hydrogarnet")
    println("because the sulfate cannot stretch, misses the level by a margin no rate")
    println("constant is responsible for. Fit its SHAPE and the parameters saturate their")
    println("bounds instead, which says the imposed assemblage produces a curve the")
    println("Parrott-Killoh family cannot make. So the surrogate is a demonstration, not")
    println("a preconditioner, and everything below runs on the coupled model.")
    rS = calibrate(
        target_fine; mode = :surrogate, θ0, spec = CALIB_SPEC_FULL,
        maxiters = 800, shape = true,
    )
    report_fit(rS, target_fine; label = "surrogate, shape only, six parameters")

    rule("what the measurement can constrain")
    @printf(
        "%d coupled solves for the six-parameter sensitivity — the expensive diagnostic\n",
        2 * length(CALIB_SPEC_FULL)
    )
    idF = local_identifiability(
        prior_vector(CALIB_SPEC_FULL), target; mode = :coupled, spec = CALIB_SPEC_FULL,
    )
    report_identifiability(idF, CALIB_SPEC_FULL)
    println("\n   parameter   ‖∂Q/∂log θ‖ (J/g)   weight in the leading 3 directions")
    for (i, p) in enumerate(CALIB_SPEC_FULL)
        @printf("   %-10s %14.4g %25.2f\n", p.name, norm(idF.J[:, i]), sum(abs2, idF.V[i, 1:3]))
    end
    println("\nfitted: ", join(string.(getfield.(CALIB_SPEC, :name)), ", "))
    println(
        "held at published values: ", join(
            string.(
                setdiff(
                    getfield.(CALIB_SPEC_FULL, :name), getfield.(CALIB_SPEC, :name)
                )
            ), ", "
        )
    )

    rule("the calibration")
    @printf(
        "budget: %d iterations × ≈%.0f s per solve\n", COUPLED_ITERS, t_coupled
    )
    r = calibrate(target; mode = :coupled, θ0, maxiters = COUPLED_ITERS)
    report_fit(r, target; label = "coupled")
    println("θ̂ = ", repr(r.θ), "   # paste into CALIBRATED_THETA")

    rule("identifiability at the optimum")
    id = local_identifiability(r.θ, target; mode = :coupled)
    report_identifiability(id)

    rule("the holdout, predicted and not fitted")
    @printf(
        "%s: w/b %.2f (target %.2f), Blaine %.0f (target %.0f) m²/kg, %.0f h\n",
        holdout.meta.cement, holdout.meta.wb, target.meta.wb,
        holdout.meta.blaine, target.meta.blaine, holdout.t[end] / 3600
    )
    for (label, θ) in (("published", θ0), ("calibrated", r.θ))
        @printf(
            "   %-11s RMSE on the holdout: %8.2f J/g\n",
            label, calorimetry_loss(θ, holdout; mode = :coupled)
        )
    end
    println("\nThe only machinery carrying the fit across the two mixes is")
    println("`powers_alpha_max(w/b)` on the hydration cap and `blaine_factor` on the rate.")

    return (;
        target, target_fine, holdout, prior = θ0, fit = r, surrogate_fit = rS,
        identifiability = id, identifiability_full = idF,
        Q_prior, Q_surrogate = Q_surr,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
