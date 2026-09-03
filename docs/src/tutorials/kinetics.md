# [Chemical Kinetics](@id sec-kinetics)

This tutorial covers the kinetics module of ChemistryLab.jl, which implements
mineral dissolution and precipitation kinetics coupled to fast aqueous speciation,
following the methodology of [Leal2017](@cite).

## Background

The kinetics algorithm solves:

```math
\frac{d n_k}{dt} = r_k(t), \quad k \in \text{kinetic minerals}
```

where ``r_k`` [mol/s] is the net rate of reaction ``k`` (positive = dissolution,
negative stoichiometric coefficient for the mineral so `dn/dt < 0`).
At each evaluation of this ODE right-hand side, the aqueous speciation is optionally
re-equilibrated, providing accurate activity coefficients for the saturation ratio
``\Omega = \text{IAP}/K``.

## [Which route: two ways to advance in time](@id sec-kinetics-routes)

There are two, they answer different questions, and choosing wrongly is the most
common mistake. Both take the same [`KineticReaction`](@ref) objects and the same
rate laws.

| | implicit step — [`kinetic_step`](@ref) | ODE integration — [`integrate`](@ref) |
|:--|:--|:--|
| what advances | one Gibbs minimization per step, extents included | `dn/dt = ν' r(n)`, handed to SciML |
| kinetics–equilibrium coupling | **strong**: the aqueous phase is at equilibrium *at the end of the step* | the aqueous phase is re-speciated between evaluations |
| the hydrate assemblage | **thermodynamic** — whatever minimizes `G` | **imposed** by the stoichiometry you wrote |
| time stepping | `Δt` is yours | adaptive, with a choice of stiff solvers |
| overshooting saturation | impossible, at any `Δt` | controlled by the integrator's tolerances |
| several reactions per mineral | yes | yes |

**Use the implicit step** when the products should be decided by thermodynamics —
a mineral dissolving into a solution, carbonation, an assemblage you do not want
to prescribe. This is the algorithm of [Leal2017](@cite), the one Reaktoro
implements, and it is unconditionally stable: with a rate law `r = k(1 − Ω)` on
calcite, a step of `10⁶ s` at `k = 10⁻⁵ mol/s` would dissolve 10 mol explicitly —
a thousand times the calcite present — and the implicit step lands at `Ω =
0.999989`, approaching saturation from below and never crossing it.

**Use the ODE route** when the reactions themselves are the model — a
stoichiometric hydration scheme in the form of [Lavergne2018](@cite), where
`C₃S + H₂O → C-S-H + CH` is written out and the products are prescribed, not
proposed. Every product then follows the coefficients you gave, several pathways
may share a mineral, and the integrator chooses the step.

### Prescribed products AND strong coupling: `coupling = :species`

The third combination, which neither Reaktoro nor the ODE route offers: a
solid-to-solid scheme whose **products are imposed by the stoichiometry**, solved
inside the same Gibbs minimization as the aqueous equilibrium.

`coupling = :reactions` (the default) puts one constraint per reaction, so the
minimization still rearranges the solids and the assemblage stays thermodynamic.
`coupling = :species` puts one per kinetic species, `nᵢ = nᵢ(0) + Σⱼ νᵢⱼ Δξⱼ`, so
every product follows the coefficients you wrote while the aqueous phase still
minimizes `G` under element conservation.

```julia
kss = KineticStepSolver(cs, DiluteSolutionModel(), reactions; coupling = :species)
```

Measured on two C₃A pathways — ettringite and monosulfoaluminate, one thousand
seconds, constant rates — against what the stoichiometry demands:

| | `:reactions` | `:species` |
|:--|--:|--:|
| extents `Δξ` | `Δt·M·r`, exactly | `Δt·r`, exactly |
| C₃A, gypsum, AFt, AFm | rearranged by `G`; AFm at zero | the stoichiometry, to 8–9 digits |
| element balance | `2×10⁻¹¹` mol | `2×10⁻¹⁴` mol |
| certificate | proved optimal | proved optimal |

!!! note "How `:species` is solved, and why it is not a constraint"
    The pinned species are **eliminated**, not constrained. Their amounts are an
    explicit affine function of the extents, `nᵢ = nᵢ(0) + Σⱼ νᵢⱼ Δξⱼ`, so they
    need not be unknowns at all: they are removed from the system, their element
    content is subtracted from the budget, and a plain equilibrium is solved over
    what remains, with a Newton on the `nr` extents around it.

    Holding them by a linear row instead — which is the obvious implementation —
    was tried and stalls. The species' stationarity row stays in the system,
    satisfied by a multiplier that must reach the mineral's own chemical
    potential, of order 10²–10³ in `RT` units, and the Newton sits at a fixed
    point its line search cannot leave: an element balance of `6.1×10⁻⁷` mol
    whatever `maxit`, `tol` or the number of active-set updates, against
    `2×10⁻¹⁴` after elimination.

    The Newton on the extents is an outer loop, deliberately. It is the price of
    imposing an assemblage, and `nr` is a handful where the composition is dozens;
    each of its evaluations is one exact, well-conditioned equilibrium solve.

Use `:reactions` when the assemblage should come from thermodynamics, the ODE
route when no aqueous coupling is needed, and `:species` when you want a
stoichiometric scheme with the pore solution still at equilibrium.

### The implicit step, in full

```julia
using ChemistryLab, DynamicQuantities, OptimaSolver, OrderedCollections

rxn  = Reaction(OrderedDict(cs["Cal"] => 1.0),
                OrderedDict(cs["Ca+2"] => 1.0, cs["CO3-2"] => 1.0);
                symbol = "calcite dissolution")
rate = KineticFunc((T, P, t, n, lna, n0) -> 1.0e-6, NamedTuple(), u"mol/s")
kr   = KineticReaction(cs, rxn, rate)

kss = KineticStepSolver(cs, DiluteSolutionModel(), [kr])

Δξ  = Ref(Float64[])
st1 = kinetic_step(kss, st0, 100.0u"s"; parameters = Δξ)
Δξ[]        # the reaction extents the step found, in mol
```

What the step solves is one problem, not two:

```math
\min_n G(n) \quad\text{s.t.}\quad
\begin{cases}
A n = b_0 \\
K^{\mathsf T} n - \Delta\xi = \xi_0 \\
\Delta\xi - \Delta t\, M\, r(n) = 0, \quad M = K^{\mathsf T}K \\
n \ge 0
\end{cases}
```

The extents are unknowns beside the amounts and the element potentials, and the
rate is evaluated at the **end-of-step** composition — backward Euler. Two things
follow. There is no frozen speciation in a right-hand side, so no lag between the
kinetic and the equilibrium species. And the reactivity constraint is *linear*, so
it joins the conservation block: the algebraic cost of kinetics is the number of
**reactions**, not of species.

`K` carries each reaction's stoichiometry restricted to the **non-aqueous**
participants, since an aqueous product re-speciates and pinning it would stop
pinning the mineral. On calcite that leaves `K = [−1]`, hence `M = 1` and
`Δξ = Δt·r`: measured, the calcite left after 1000 s at 1 µmol/s is
`n₀ − kΔt` to the last bit, with the element balance at `10⁻¹⁴`.

### Choosing the step length, and a case where the stiff ODE route is wrong

[`kinetic_step_adaptive`](@ref) takes one step of `Δt` and two of `Δt/2`, compares
the extents, and accepts the finer pair when the difference is within tolerance —
Richardson's estimate, which for a first-order method IS the error of the coarse
step. Three implicit solves per accepted step is the price.

Reaktoro has nothing equivalent: its kinetics adds a single initial step to the
equilibrium options, the step is the caller's, and backward Euler being first
order, one ten times too large is ten times less accurate with nothing to say so.

Measured on calcite dissolving under `r = k(1 − Ω)` with `k = 10⁻⁴ mol/s` over
`10⁵ s`, against the equilibrium the trajectory converges to:

| route | steps | result |
|:--|--:|:--|
| `Rodas5P`, and `OrdinaryDiffEq`'s default polyalgorithm | 6 | **extent −457 mol**, `retcode = Success` |
| `Tsit5`, explicit | 85 626 | correct, in 519 s |
| one `kinetic_step` of `10⁵ s` | 1 | wrong — and its certificate says so |
| `kinetic_step_adaptive` | **7** | correct to eight digits |

The first row is the one to know about, and the natural explanation is the wrong
one. It is **not** a missing Jacobian term: the ODE route's residual reads the
speciation frozen at the last accepted step, so `∂(du)/∂bₑ = 0` is exact for the
system being integrated. What goes wrong is that the frozen speciation makes the
right-hand side inconsistent with the state *within* a step — an implicit method
steps past the point where `Ω` crosses one, the rate changes sign, and the run
enters a branch it never leaves.

Three measurements settle it. Bounding `dtmax` to `10³ s`, which takes 103 steps
instead of 6, returns the **identical** wrong value to seven digits, so the error
is not in the time discretization. Removing the re-speciation returns a sane
answer. And routing the re-speciation through the certified route moves −457 to
−383, so the quality of the partition is not the cause either.

So the fix is not to freeze the speciation, which is exactly what the implicit
step does. The extension warns when a trajectory ends on amounts no chemistry can
produce — 457 mol of calcite from a budget of 55.6 — and that is what the ODE
route can offer: a rate law that reads the solution belongs on the implicit
route.

A step is accepted on the estimate **and** on the certificate, never on the
estimate alone. Richardson's difference measures the disagreement between two
resolutions, so it is blind to an error the two share: measured on calcite over
`10⁵ s`, the coarse step and both half-steps each dissolve the entire mineral,
their extents agree to `5×10⁻¹¹`, and the estimator reports `9×10⁻⁶` — a step it
should have refused, graded excellent. The certificate is not fooled, because a
composition that dissolved everything violates `Δξ − Δt·M·r(n) = 0` by the whole
extent. With both conditions the same march refuses `10⁵ s`, comes down to
`1.6×10³ s`, and grows back to `5×10⁴ s` over seven steps.

The tolerance is relative to the amount each reaction acts on, not to the extent.
Scaling by `Δξ` is the obvious thing to write and does not work: `Δξ ∝ Δt`, so
that tolerance vanishes with the step while the equilibrium solve's own noise does
not, the measured error behaves as `noise/(reltol·Δt)` and **grows** as the step
shrinks. Measured, the controller then halved to the floor without advancing, and
the final answer got worse as the tolerance was tightened. An ODE integrator's
`reltol` multiplies the solution for exactly this reason.

### Solid solutions in the step, and the two settings the certificate decides

A solid solution takes part like any other phase — the products of a kinetic
dissolution may be a C-S-H solution, decided by thermodynamics at the end of each
step. On C₃S dissolving at 20 µmol/s into a four-member CSHQ solution, the step
comes out certified at a stationarity of `2×10⁻¹³`, with `Δξ = Δt·r` exactly, all
four end-members present, and ten times the step giving ten times the C-S-H.

Two settings are decided by the certificate rather than by the caller, because
neither answer works on both kinds of problem:

  - `warm_start` (default `true`) equilibrates the starting **guess** when the
    system carries solid solutions, leaving the component totals untouched. A
    mixing phase is admitted by a tangent-plane test, and from a composition where
    the phase is absent that admission fails: a cold start left one end-member at
    `2.7×10⁻⁹` with a stationarity residual of 6.5. The failure belongs to the
    cold start and not to the kinetics — a plain equilibrium from the same guess
    fails identically.
  - `pin_minerals` (default `:auto`) says whether the kinetic minerals are held in
    the active set. Measured: on two C₃A pathways, **not** pinning is certified at
    `1.8×10⁻¹²` while pinning gives 7.9 and no certificate; on C₃S into a C-S-H
    solution, not pinning runs the step all the way to equilibrium — `C3S = 0`,
    `Δξ` twenty-five times too large, because the mineral is exhausted at
    equilibrium and the drop rule removes it, after which nothing enforces its
    reactivity row. `:auto` tries both and keeps the answer the certificate
    **proves**, which is exact rather than heuristic: the problem is convex, so
    its KKT conditions decide.

Pass `certificate = Ref{Any}(nothing)` to see the proof. It is taken on the
augmented problem, whose reactivity rows carry the multipliers that make a
kinetically held mineral's stationarity satisfiable — the same composition tested
against the *unconstrained* equilibrium reports a residual of 7 `RT`, because a
mineral held back by a rate law is supersaturated by construction. That is what
being held back means.

## Rate functions: KineticFunc and StateView

All rate functions are encapsulated as [`KineticFunc`](@ref) objects.  They share the
**positional six-argument calling convention**:

```julia
kf(T, P, t, n::StateView, lna::StateView, n_initial::StateView) -> Real [mol/s]
```

where

| Argument | Unit | Description |
| :-- | :-- | :-- |
| `T` | K | Temperature (plain `Real` or `ForwardDiff.Dual`) |
| `P` | Pa | Pressure |
| `t` | s | Current integration time |
| `n` | mol | Moles of all species — named access via [`StateView`](@ref) |
| `lna` | — | Log-activities — named access via `StateView` |
| `n_initial` | mol | Initial moles — named access via `StateView` |

### What a rate law may depend on

Everything in that table, by species name. A rate law is an ordinary Julia
function, so there is no expression language to learn and no restriction on what
it may read:

```julia
# 1. a constant rate
KineticFunc((T, P, t, n, lna, n0) -> 1.0e-6, NamedTuple(), u"mol/s")

# 2. depending on how much is left — first order in the mineral
KineticFunc((T, P, t, n, lna, n0) -> k * n["C3S"], NamedTuple(), u"mol/s")

# 3. depending on the SOLUTION — pH, or any activity
KineticFunc(
    (T, P, t, n, lna, n0) -> begin
        a_H = exp(lna["H+"])            # `lna` is the natural log of the activity
        k * a_H^0.5
    end, NamedTuple(), u"mol/s")

# 4. depending on time explicitly
KineticFunc((T, P, t, n, lna, n0) -> k * exp(-t / 86400), NamedTuple(), u"mol/s")

# 5. depending on the degree of reaction, since `n0` is there too
KineticFunc(
    (T, P, t, n, lna, n0) -> begin
        α = 1 - n["C3S"] / n0["C3S"]
        k * (1 - α)^(2 / 3)
    end, NamedTuple(), u"mol/s")
```

### A reaction that waits for another species to run out

This is the pattern [Lavergne2018](@cite) uses for the aluminate sequence, and it
needs nothing special: read the species you are waiting on and multiply by a gate.

In that model C₃A takes three successive routes, and which one runs depends on
what is left:

| stage | while | C₃A reacts with | forming |
|--:|:--|:--|:--|
| 1 | gypsum remains | gypsum | ettringite (AFt) |
| 2 | gypsum gone, ettringite remains | the ettringite already formed | monosulfoaluminate (AFm) |
| 3 | both gone | water alone | C₃AH₆ |

Three [`KineticReaction`](@ref) objects on the same mineral, each gated on the
reactant its stage consumes:

```julia
# Stage 1 — runs while there is gypsum
r_aft = KineticFunc(
    (T, P, t, n, lna, n0) -> begin
        Gp = n["Gp"]
        Gp <= 0 && return zero(Gp)
        k1 * n["C3A"] * _gate(Gp, n0["Gp"])
    end, NamedTuple(), u"mol/s")

# Stage 2 — takes over once the gypsum is exhausted, and consumes the AFt
r_afm = KineticFunc(
    (T, P, t, n, lna, n0) -> begin
        ett = n["ettringite"]
        ett <= 0 && return zero(ett)
        k2 * n["C3A"] * (1 - _gate(n["Gp"], n0["Gp"])) * _gate(ett, n0["Gp"])
    end, NamedTuple(), u"mol/s")

# Stage 3 — only once BOTH are gone
r_hydrogarnet = KineticFunc(
    (T, P, t, n, lna, n0) -> begin
        k3 * n["C3A"] *
            (1 - _gate(n["Gp"], n0["Gp"])) * (1 - _gate(n["ettringite"], n0["Gp"]))
    end, NamedTuple(), u"mol/s")

# A smooth gate: 1 while the reactant is plentiful, 0 once it is spent.
_gate(x, xref) = tanh(max(x, zero(x)) / (1e-3 * xref))
```

Note that the gates are written from the amounts, and that stage 2 consumes a
species stage 1 produced. Both work because the reaction extents are carried in
the state, so a non-kinetic amount moves during the step — see the two rules
below.

Two rules make gates like this behave:

  - **Return a zero of the right type**, `zero(CH)` rather than `0.0`. The rate is
    differentiated for the stiff solver's Jacobian and for parameter sensitivity,
    and a hard `Float64` zero breaks the dual-number chain.
  - **Prefer a smooth gate** where you can, as `_gate` above does. A `min` or a
    hard cut-off is a kink, and an adaptive integrator will shorten its step to
    crawl over it; `tanh` costs nothing and integrates cleanly. Keep the hard
    `<= 0 && return zero(x)` guard as well — it is what keeps a rate from being
    computed on a negative amount the integrator proposed and rejected.

A gate on a species the reactions consume only works because the reaction extents
are carried in the ODE state, which is what lets non-kinetic amounts move during
the step. Before that they were pinned at their initial values inside the
residual, so a gate never closed and the reaction ran past depletion — while the
kinetic mass balance stayed exact and the solver reported success. That failure
mode is pinned by a test; see
[Rate laws that depend on a consumed reactant](@ref kinetics-frozen-species).

!!! tip "Several reactions on one mineral"
    Write one [`KineticReaction`](@ref) per pathway, each with its own rate law,
    and pass them all. The ODE state holds one entry per unique mineral and the
    contributions accumulate, so C₃A reacting through both an ettringite and a
    monosulfoaluminate pathway is two reactions sharing a species — see
    [Multiple kinetic reactions per mineral](@ref). Gating one pathway on the
    exhaustion of the sulfate is how the switch between them is written.

[`StateView`](@ref) provides O(1) named access to a species vector via a pre-built
dictionary (`sv["C3S"]` — no per-step dict allocation):

```julia
index = Dict("C3S" => 1, "C2S" => 2)
data  = [2.5, 1.8]          # moles
sv    = StateView(data, index)

sv["C3S"]          # → 2.5
haskey(sv, "C3S")  # → true
```

The `index` dict is built **once** at `KineticsProblem` construction; the same dict is
shared by all `StateView` wrappers inside the ODE hot path.

## Setting up rate constants

Rate constants are built as [`NumericFunc`](@ref) objects using the Arrhenius factory:

```julia
using ChemistryLab, ForwardDiff

# Arrhenius rate constant for calcite acid dissolution ([PalandriKharaka2004](@cite))
k_acid = arrhenius_rate_constant(5.012e-1, 14400.0)   # k₀ [mol/(m²s)], Ea [J/mol]

k_acid(; T = 298.15)    # → 0.5012 mol/(m²s)
k_acid(; T = 310.0)     # → higher value at elevated T

ForwardDiff.derivative(T -> k_acid(; T = T), 298.15)  # AD-compatible
```

## Cement clinker hydration: [ParrotKilloh1984](@cite) model

[`parrot_killoh_avrami`](@ref) is the factory for the [ParrotKilloh1984](@cite)
kinetic model. It returns a [`KineticFunc`](@ref) that uses the `StateView`-based
calling convention and is fully AD-compatible.

```julia
using ChemistryLab, DynamicQuantities

# Predefined NamedTuple parameters for the four main clinker phases
#   PK84_PARAMS_C3S, PK84_PARAMS_C2S, PK84_PARAMS_C3A, PK84_PARAMS_C4AF
# Each has keys: k₁, n₁, k₂, k₃, n₃, Ea, T_ref (with units)

pk_C3S = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S")    # → KineticFunc

# Water availability limit of [Powers1948](@cite) for w/c = 0.40
WC    = 0.40
α_max = powers_alpha_max(WC)          # 0.952 at w/c = 0.40
pk_C3S_wc = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S"; α_max = α_max)

# Fineness correction — the rate scales as B / 385 m²/kg
pk_C3S_fine = parrot_killoh_avrami(
    PK84_PARAMS_C3S, "C3S"; α_max = α_max, blaine = 380u"m^2/kg",
)
```

!!! warning "`parrot_killoh` is deprecated — use `parrot_killoh_avrami`"
    ChemistryLab also exports [`parrot_killoh`](@ref), a smoothed variant with its
    own parameter sets `PK_PARAMS_*`. Its formulas are **not** those of
    [ParrotKilloh1984](@cite) and its parameters match no published set, so the
    attribution has been withdrawn. With `PK_PARAMS_*` the diffusion branch takes
    over within the first few percent of hydration, and C₃S, C₂S and C₃A then all
    land on α(7 d) = 0.2386 whatever their `K₁` — against the ≈ 0.61 the
    literature reports for a CEM I at w/c = 0.40. See
    [Two Parrot–Killoh variants](@ref pk-variants). Supplementary cementitious
    materials do not follow Parrot–Killoh at all: use [`waller`](@ref).

Rate evaluation:

```julia
index = Dict("C3S" => 1)
n_sv  = StateView([0.5], index)    # 0.5 mol C3S remaining
n0_sv = StateView([1.0], index)    # 1.0 mol initial
lna   = StateView([0.0], index)

r = pk_C3S(293.15, 1.0e5, 0.0, n_sv, lna, n0_sv)   # mol/s
```

AD smoke-test:

```julia
using ForwardDiff

drdT = ForwardDiff.derivative(T -> pk_C3S(T, 1.0e5, 0.0, n_sv, lna, n0_sv), 293.15)
@assert isfinite(drdT) && drdT > 0
```

## Transition-state theory rate model

For models based on solution chemistry (calcite, quartz, …),
[`transition_state`](@ref) builds a multi-mechanism TST [`KineticFunc`](@ref):

```julia
k_neutral = arrhenius_rate_constant(1.549e-6, 23500.0)
k_acid    = arrhenius_rate_constant(5.012e-1, 14400.0)

surface = BETSurfaceArea(90.0)   # 90 m²/kg

tst = transition_state(
    [
        RateMechanism(k_neutral, 1.0, 1.0),
        RateMechanism(k_acid,    1.0, 1.0, [RateModelCatalyst("H+", 1.0)]),
    ],
    cs,
    cs.dict_reactions["calcite dissolution"],
    surface,
)

# tst is a KineticFunc — same calling convention as parrot_killoh
kr = KineticReaction(cs, cs.dict_reactions["calcite dissolution"], tst)
```

!!! warning "`cs.dict_reactions` holds only the *kinetic* reactions"
    The lookup above works only if calcite was declared kinetic on the
    `ChemicalSystem`. For any other reaction, build the `Reaction` yourself and
    pass it directly — `dict_reactions` is not a catalog of the database.

The `transition_state` factory captures the ΔₐG⁰ callables from the reaction and
recomputes `ln K(T)` at every ODE step, so the saturation ratio Ω is correct even
when the temperature changes (semi-adiabatic calorimetry).

For a single mechanism without catalysts, use [`first_order_rate`](@ref) as a
convenience wrapper.

## Defining kinetic reactions

[`KineticReaction`](@ref) associates a [`Reaction`](@ref) with a `KineticFunc`:

```julia
pk = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S")

# Convenience: look up "C3S" in cs by symbol/formula, build minimal dissolution Reaction
kr = KineticReaction(cs, "C3S", pk)

# Explicit Reaction object (e.g. from cs.dict_reactions)
rxn = cs.dict_reactions["calcite dissolution"]
kr  = KineticReaction(cs, rxn, tst_calcite)

# Reaction-centric (kinetics stored in rxn.properties)
rxn[:rate] = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S")
kr = KineticReaction(cs, rxn)
```

## [Reaction-centric workflow](@id reaction-centric-workflow)

Kinetics can be attached directly to `Reaction` objects via their `properties` dict:

### Properties table

| Key | Type | Required | Description |
| --- | --- | :---: | --- |
| `:rate` | `KineticFunc` or `(T,P,t,n,lna,n₀)->Real` | ✓ | Net dissolution rate [mol/s] |
| `:ΔᵣH⁰` | `AbstractFunc` or `Number` | | Reaction enthalpy [J/mol] (thermodynamic convention: negative = exothermic). Computed automatically from species `ΔₐH⁰` for balanced reactions; set explicitly for custom species without `ΔₐH⁰`. |

When `:rate` is a plain callable (not a `KineticFunc`), it is automatically wrapped
in a `KineticFunc` with empty `refs`.

### Example: C₃S hydration (reaction-centric)

```julia
using ChemistryLab, OrdinaryDiffEq, DynamicQuantities

sp(name) = cs[name]

pk_C3S = parrot_killoh_avrami(PK84_PARAMS_C3S, "C3S"; α_max = powers_alpha_max(WC))

rxn_C3S = Reaction(
    OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 3.33),
    OrderedDict(sp("Jennite") => 0.167, sp("Portlandite") => 1.5);
    symbol = "C₃S hydration",
)
rxn_C3S[:rate] = pk_C3S   # heat computed from species ΔₐH⁰ (CEMDATA18)

kp = KineticsProblem(cs, [rxn_C3S, rxn_C2S, rxn_C3A, rxn_C4AF], state0, tspan)
```

## Unit-aware API

All constructors support DynamicQuantities `Quantity` values.  Plain `Real` → SI assumed.

| Constructor | Parameter | Default SI unit | Example |
| :-- | :-- | :-- | :-- |
| `arrhenius_rate_constant` | `k₀` | mol/(m²·s) | `0.5u"mmol/(m^2*s)"` |
| `arrhenius_rate_constant` | `Ea` | J/mol | `62.0u"kJ/mol"` |
| `parrot_killoh_avrami` params | `k₁, k₂, k₃` | s⁻¹ | `1.5u"1/d"` |
| `parrot_killoh_avrami` params | `Ea` | J/mol | `42.0u"kJ/mol"` |
| `parrot_killoh_avrami` params | `T_ref` | K | `293.15u"K"` |
| `FixedSurfaceArea` | `A` | m² | `500.0u"cm^2"` |
| `BETSurfaceArea` | `A_specific` | m²/kg | `0.09u"m^2/g"` |
| `Reaction` | `:ΔᵣH⁰` (custom species) | J/mol | `NumericFunc((T,) -> -36100.0, (:T,), u"J/mol")` |
| `KineticsProblem` | `tspan` | s | `(0.0u"s", 7.0u"d")` |
| `IsothermalCalorimeter` | `T` | K | `298.15u"K"` |
| `SemiAdiabaticCalorimeter` | `Cp` | J/K | `4.0u"kJ/K"` |
| `SemiAdiabaticCalorimeter` | `T_env` | K | `293.15u"K"` |
| `SemiAdiabaticCalorimeter` | `L` | W/K | `500.0u"mW/K"` |
| `SemiAdiabaticCalorimeter` | `T0` | K | `293.15u"K"` |

## Running a simulation

```julia
using OrdinaryDiffEq   # activates KineticsOrdinaryDiffEqExt

kp  = KineticsProblem(cs, [kr], state0, (0.0, 7200.0))
sol = integrate(kp)   # default Rodas5P

ks  = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1e-8, abstol = 1e-10)
sol = integrate(kp, ks)
```

## Isothermal calorimetry

```julia
cal = IsothermalCalorimeter(298.15u"K")   # T = 25 °C
kp = KineticsProblem(cs, reactions, state0, tspan; calorimeter = cal)
sol = integrate(kp, ks)

t, Q    = cumulative_heat(sol, cal)   # Q(t) [J]
t, qdot = heat_flow(sol, cal)         # q̇(t) [W]
```

## Calorimetry under partial equilibrium

The two calorimeters above take their heat from [`heat_rate`](@ref), which sums
`rᵢ(−ΔᵣH⁰ᵢ)` over the **kinetic** reactions. That is exact when those reactions
produce the hydrates. It is not, and cannot be, when a `equilibrium_solver` is
attached: the kinetic reactions then only dissolve the anhydrous phases into
ions, and the hydrates are precipitated by the Gibbs minimization, whose heat
that sum does not see. On an ordinary Portland cement, driving a semi-adiabatic
cell from it gave a temperature rise of 207 K against the few tens of kelvin a
Langavant test measures — so that combination now warns rather than returning the
number in silence.

Use [`heat_release`](@ref) instead. Enthalpy is a state function, so the heat
released between two states at the same temperature is simply their difference,
with reactants, ions and hydrates each counted once and no reaction stoichiometry
to write down — Eqs. (17)–(21) of [Lavergne2018](@cite):

```math
-\delta Q \;=\; \mathrm{d}H
\;=\; \Bigl(\sum_i n_i C^\circ_{p,i}(T)\Bigr)\mathrm{d}T
\;+\; \sum_i \Delta_f H_i(P,T)\,\mathrm{d}n_i .
```

```julia
kp  = KineticsProblem(cs, reactions, state0, tspan;
                      equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()))
sol = integrate(kp, ks)

t, Q, q̇ = heat_release(sol, kp; times = my_times)   # J and W
```

It reads the **certified** speciations of [`speciated_states`](@ref), not the
composition the integrator carries. The in-run minimization is warm-started and
uncertified, and one hydrate is worth hundreds of kilojoules: read that way the
curve came out at 12.7, 145, 1174, 936 and 631 J/g at 1 h, 6 h, 12 h, 1 d and
2 d — heat that rises and then falls, which no calorimeter has ever measured.

[`enthalpy`](@ref) and [`heat_capacity`](@ref) give the same sums for a single
state, and [`missing_enthalpy`](@ref) lists the species that carry no `ΔₐH⁰` and
are therefore absent from the balance — worth checking, because a heat curve
missing one hydrate is not visibly wrong.

## Semi-adiabatic calorimetry [Lavergne2018](@cite)

The semi-adiabatic calorimeter solves:

```math
\frac{dT}{dt} = \frac{\dot{q}(t) - \varphi(T(t) - T_{\rm env})}{C_p + \sum_i n_i C^\circ_{p,i}(T)}
```

The denominator uses the **variable total heat capacity** `Cp_total = Cp + Σᵢ nᵢ Cp°ᵢ(T)`,
where `Cp°ᵢ(T)` are the molar heat capacities from the thermodynamic database
([Lavergne2018](@cite)).

[`SemiAdiabaticCalorimeter`](@ref) bundles hardware parameters and initial temperature:

```julia
using ChemistryLab, DynamicQuantities

# Quadratic heat loss ([Lavergne2018](@cite): φ = a·ΔT + b·ΔT²)
cal = SemiAdiabaticCalorimeter(;
    Cp        = (1.0 * 800.0 + WC * 4186.0 + 1.0 * 900.0) * u"J/K",   # ≈ 3449 J/K
    T_env     = 293.15u"K",
    heat_loss = ΔT -> 0.3 * ΔT + 0.003 * ΔT^2,
    T0        = 293.15u"K",
)

# Shorthand for linear Newton cooling (L keyword)
cal_lin = SemiAdiabaticCalorimeter(;
    Cp = 4000.0u"J/K", T_env = 293.15u"K", L = 0.5u"W/K", T0 = 293.15u"K",
)

kp = KineticsProblem(cs, reactions, state0, tspan; calorimeter = cal)
sol = integrate(kp, ks)

t, T_profile = temperature_profile(sol, cal)   # T(t) [K]
t, qdot      = heat_flow(sol, cal)             # q̇(t) [W]
t, Q         = cumulative_heat(sol, cal)       # Q(t) [J]
```

## Full workflow: OPC clinker hydration with KineticsProblem

Complete chain — [`ChemicalSystem`](@ref) → [`ChemicalState`](@ref) →
[`KineticReaction`](@ref) → [`KineticsProblem`](@ref) → [`integrate`](@ref):

```julia
using ChemistryLab, OrdinaryDiffEq, DynamicQuantities, Printf

# ── 1. ChemicalSystem from CEMDATA18 ────────────────────────────────────────
DATA_FILE = datapath("cemdata18-thermofun.json")
substances = build_species(DATA_FILE)

input_species = split(
    "C3S C2S C3A C4AF " *
    "Portlandite Jennite ettringite monosulphate12 C3AH6 C3FH6 H2O@",
)
species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])
cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

# ── 2. Initial state: 1 kg OPC (CEM I 52.5 R), w/c = 0.40 ─────────────────
WC          = 0.40
COMPOSITION = (C3S=0.619, C2S=0.165, C3A=0.080, C4AF=0.087)

state0 = ChemicalState(cs)
for (name, frac) in pairs(COMPOSITION)
    set_quantity!(state0, string(name), frac * u"kg")
end
set_quantity!(state0, "H2O@", WC * u"kg")

# ── 3. [ParrotKilloh1984](@cite) rate functions with Powers α_max ───────────────────────
α_max   = powers_alpha_max(WC)
BLAINE  = 380.0u"m^2/kg"
pk_C3S  = parrot_killoh_avrami(PK84_PARAMS_C3S,  "C3S";  α_max, blaine = BLAINE)
pk_C2S  = parrot_killoh_avrami(PK84_PARAMS_C2S,  "C2S";  α_max, blaine = BLAINE)
pk_C3A  = parrot_killoh_avrami(PK84_PARAMS_C3A,  "C3A";  α_max, blaine = BLAINE)
pk_C4AF = parrot_killoh_avrami(PK84_PARAMS_C4AF, "C4AF"; α_max, blaine = BLAINE)

# ── 4. Kinetic reactions (reaction-centric) ─────────────────────────────────
# Reactions follow [LothenbachWinnefeld2006](@cite) — Jennite = Ca₉Si₆O₁₈(OH)₆·8H₂O
# Balanced hydration reactions — ΔᵣH⁰ computed from species ΔₐH⁰.
sp(name) = cs[name]

rxn_C3S = Reaction(
    OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 103/30),
    OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 4/3);
    symbol = "C₃S hydration",
)
rxn_C3S[:rate] = pk_C3S

rxn_C2S = Reaction(
    OrderedDict(sp("C2S") => 1.0, sp("H2O@") => 73/30),
    OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 1/3);
    symbol = "C₂S hydration",
)
rxn_C2S[:rate] = pk_C2S

rxn_C3A = Reaction(
    OrderedDict(sp("C3A") => 1.0, sp("H2O@") => 6.0),
    OrderedDict(sp("C3AH6") => 1.0);
    symbol = "C₃A hydration",
)
rxn_C3A[:rate] = pk_C3A

rxn_C4AF = Reaction(
    OrderedDict(sp("C4AF") => 1.0, sp("Portlandite") => 2.0, sp("H2O@") => 10.0),
    OrderedDict(sp("C3AH6") => 1.0, sp("C3FH6") => 1.0);
    symbol = "C₄AF hydration",
)
rxn_C4AF[:rate] = pk_C4AF

# ── 5. Problem + semi-adiabatic calorimeter ─────────────────────────────────
cal = SemiAdiabaticCalorimeter(;
    Cp        = (1.0 * 800.0 + WC * 4186.0 + 1.0 * 900.0) * u"J/K",
    T_env     = 293.15u"K",
    heat_loss = ΔT -> 0.30 * ΔT + 0.003 * ΔT^2,
    T0        = 293.15u"K",
)

kp = KineticsProblem(
    cs, [rxn_C3S, rxn_C2S, rxn_C3A, rxn_C4AF], state0, (0.0, 7.0 * 86400.0);
    calorimeter = cal,
    equilibrium_solver = nothing,
)

# ── 6. Integrate and post-process ───────────────────────────────────────────
ks  = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1e-6, abstol = 1e-9)
sol = integrate(kp, ks)

_, T_vec = temperature_profile(sol, cal)
_, Q_vec = cumulative_heat(sol, cal)

n0_kin = [sol.prob.p.n_initial_full[i] for i in kp.idx_kinetic]
n_kin  = [[u[i] for u in sol.u] for i in eachindex(n0_kin)]

function phase_alpha(cs, kp, n0_kin, n_kin, name)
    sp_idx = findfirst(s -> ChemistryLab.symbol(s) == name, cs.species)
    pos    = findfirst(==(sp_idx), kp.idx_kinetic)
    isnothing(pos) && return fill(NaN, length(sol.t))
    return 1.0 .- n_kin[pos] ./ n0_kin[pos]
end

α_C3S = phase_alpha(cs, kp, n0_kin, n_kin, "C3S")

@printf "ΔT_max = %.2f °C   Q_7d = %.1f kJ/kg   α(C3S) = %.4f\n" \
    maximum(T_vec .- 273.15) - 20.0  Q_vec[end]/1000  α_C3S[end]
```

!!! tip "Choosing `equilibrium_solver`"
    Setting `equilibrium_solver = nothing` skips the Gibbs minimization, which is
    appropriate for the [ParrotKilloh1984](@cite) model: its rate closure ignores
    the `lna` argument, so re-speciation cannot change `α(t)` or the calorimetry.
    It does change what the *products* are — with `nothing`, the hydrate
    assemblage is whatever the hand-written reaction stoichiometry says, not what
    thermodynamics gives. For `transition_state` models, whose rate depends on the
    saturation ratio `Ω = IAP/K`, an `EquilibriumSolver` is required for the rates
    themselves to mean anything.

    The solver may be passed on the [`KineticsProblem`](@ref) or on the
    [`KineticsSolver`](@ref); both work. If it is set on both, the one on the
    problem is used and the conflict is reported.

## The equilibrium–kinetics coupling

The coupling is the partitioned formulation of [Leal2017](@cite), the one
Reaktoro implements. Species are split into a **kinetic partition** — the
minerals carrying a rate law — and an **equilibrium partition**, everything
else: the aqueous phase and any mineral free to precipitate or dissolve
instantaneously.

The ODE state is `(bₑ, nₖ)`: the element amounts held by the equilibrium
partition, and the moles of the kinetic minerals. It advances as

```math
\frac{\mathrm{d} n_k}{\mathrm{d} t} = \nu_k^{\mathsf T} r,
\qquad
\frac{\mathrm{d} b_e}{\mathrm{d} t} = A_e \, \nu_e^{\mathsf T} r ,
```

and the composition of the equilibrium partition is recovered at each step by

```math
n_e = \varphi(b_e) \;=\; \arg\min_{n} \; G(n)
\quad \text{s.t.} \quad A_e\, n = b_e , \; n \ge 0 .
```

Three points are worth stating, because each is a way the coupling can be got
wrong and look plausible:

**The minimization runs over the equilibrium partition only.** Posing it on the
whole system would equilibrate the kinetic minerals instantaneously, which is
what a kinetic description exists to prevent. `KineticsProblem` therefore builds
a sub-system restricted to the partition, sharing the parent's primary species
so that `bₑ`, `dbₑ/dt` and the solve all live in one conservation basis.

**`bₑ` is integrated, not `nₑ`.** Along the way an individual species may want
to go negative — the generated dissolution reactions are written in `H⁺`, and a
cement paste contains no acid — and it is the minimizer, not the caller, that
redistributes the elements over a feasible set. Element amounts are what is
conserved; species amounts are what is solved for.

**The rate sign is fixed by the stoichiometry.** Each kinetic reaction is
normalized so its controlling mineral carries `ν = −1`: a positive rate is a
dissolution. A reaction generated from the nullspace comes out with an arbitrary
orientation, and taken as-is the ODE grows the clinker instead of consuming it.

!!! note "How the two are coupled in time — operator splitting"
    The equilibrium is **not** solved inside the ODE right-hand side. The ODE
    advances the kinetic minerals with the speciation held frozen, and
    `respeciate!` re-equilibrates the equilibrium partition once per accepted
    step, wired as a `DiscreteCallback`.

    This is deliberate. A stiff solver evaluates the right-hand side many times
    per step and differentiates it to build its Jacobian; a Gibbs minimization in
    there makes the cost per step unpredictable and leaves the Jacobian
    describing a different model from the one being integrated. Splitting the two
    is first-order accurate in the step size — the standard arrangement in
    reactive transport — and keeps residual and Jacobian consistent.

    The element amounts `bₑ` carried by the ODE state are handed to the solver
    as the constraint of the sub-problem, not derived from a starting
    composition. The composition passed alongside is a starting guess only.

!!! warning "A failed re-speciation is reported, not hidden"
    If the equilibrium solve fails, that step keeps its frozen composition, the
    first failure is reported with its exception, and `integrate` states how many
    steps were affected. A run in which re-speciation never succeeded must not
    look like a healthy one.

!!! tip "Calorimetry and ΔᵣH⁰"
    The calorimeter computes the heat generation rate as `q̇ = Σ rᵢ × (−ΔᵣH⁰ᵢ)`,
    where `ΔᵣH⁰ᵢ(T)` is the reaction enthalpy (thermodynamic convention:
    negative = exothermic). It is built automatically from species `ΔₐH⁰`
    properties via `complete_thermo_functions!` — **reactions must be
    mass-balanced** for this computation to be correct. For **custom species**
    that lack a `ΔₐH⁰` entry (GGBS, MK, …), set `:ΔᵣH⁰` directly on the
    reaction: `rxn[:ΔᵣH⁰] = NumericFunc((T,) -> -36_100.0, (:T,), u"J/mol")`.

## Multiple kinetic reactions per mineral

Following [Leal2017](@cite), **reactions** — not species — carry kinetics.
A single mineral can appear in multiple `KineticReaction` objects
(e.g. C₃A via an early ettringite pathway and a late monosulphate pathway).
The ODE state contains **one entry per unique mineral**; contributions accumulate
in `du[j]`.

```julia
pk_c3a = parrot_killoh_avrami(PK84_PARAMS_C3A, "C3A"; α_max)

kr_C3A_ett  = KineticReaction(cs, rxn_C3A_ettringite,   pk_c3a)
kr_C3A_mono = KineticReaction(cs, rxn_C3A_monosulphate, pk_c3a)

kp = KineticsProblem(
    cs, [kr_C3S, kr_C2S, kr_C3A_ett, kr_C3A_mono, kr_C4AF], state0, tspan,
)
# length(build_u0(kp)) == 4  (one entry per unique mineral: C3S, C2S, C3A, C4AF)
```

!!! note "Species classification [Leal2017](@cite)"
    `KineticsProblem` automatically applies the [Leal2017](@cite) classification:
    - **Kinetic species** (tracked in ODE state `u`): `AS_CRYSTAL` species with non-zero
      stoichiometry in any kinetic reaction.
    - **Equilibrium species**: aqueous species re-equilibrated by `equilibrium_solver`.
    - **Inert species**: all others.

## Parameter sensitivity and optimization

The rate laws themselves are ForwardDiff-clean, and tested to be: no `Float64`
casts in the evaluation path, branch guards through `_primal`, and
[`arrhenius_rate_constant`](@ref) differentiable in `k₀`, `Ea` and `T_ref`.

```julia
using ForwardDiff

pk(Ea) = parrot_killoh_avrami(merge(PK84_PARAMS_C3S, (Ea = Ea,)), "C3S")
idx = Dict("C3S" => 1)
f(Ea) = pk(Ea)(300.0, 1.0e5, 3600.0, StateView([0.9], idx), StateView([0.0], idx),
    StateView([1.0], idx))
ForwardDiff.derivative(f, 42_000.0)
```

!!! warning "Differentiating *through* `integrate` with respect to parameters does not work"
    The ODE right-hand side is AD-clean in the state `u` and in the time `t` —
    deliberately, because `Rodas5P` evaluates it with a dual `t` for the time
    gradient. It is **not** AD-clean in the parameters, and a dual number baked
    into a rate closure cannot reach the integrator:

      - `build_u0` returns a `Vector{Float64}`, so `∂/∂n₀` and `∂/∂T₀` die at the
        initial condition;
      - `build_kinetics_params` casts `T`, `P`, the initial amounts, the
        stoichiometric matrices and the calorimeter's `Cp` and `T_env` to
        `Float64`;
      - the outer constructors of [`FixedSurfaceArea`](@ref) and
        [`BETSurfaceArea`](@ref) cast to `Float64` although the structs are
        declared `{T<:Real}`, so `∂/∂A` is blocked;
      - `_strip_heat_per_mol` casts an explicit `heat_per_mol` to `Float64`, so
        `∂/∂ΔᵣH` is blocked;
      - `system_enthalpy` accumulates into a `0.0`;
      - and under partial equilibrium `respeciate!` writes into a `Float64` buffer
        by construction, so the coupled branch is `Float64`-only end to end.

    An earlier version of this section claimed the opposite and showed a
    derivative through `integrate`. There was no test for it, and the lines above
    are why. Until that changes, a parameter study needs a derivative-free or
    finite-difference outer loop —
    [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/) with
    `AutoFiniteDiff()` and `NelderMead()`, as
    [the calibration example](@ref ex-hydration-calibration) does — or finite
    differences by hand.

    What *does* differentiate is the equilibrium map on its own:
    `EquilibriumSolver` implements the implicit function theorem on the KKT system
    for `ForwardDiff.Dual` inputs, so `∂n*/∂θ` through a Gibbs minimization is
    exact.

## [Two Parrot–Killoh variants](@id pk-variants)

ChemistryLab ships **two** implementations of the Parrot & Killoh clinker
hydration model. They are not interchangeable, and their parameter sets are not
transferable between them.

| | [`parrot_killoh`](@ref) (deprecated) | [`parrot_killoh_avrami`](@ref) |
|:--|:--|:--|
| Nucleation–growth | `(K₁/N₁)(1-ξ)^N₁ / (1 + B·ξ^N₃)` | `(k₁/n₁)(1-ξ)(-ln(1-ξ))^(1-n₁)` (Avrami) |
| Second mechanism | `K₂(1-ξ)^N₂` | `k₂(1-ξ)^(2/3) / (1-(1-ξ)^(1/3))` (Jander) |
| Third mechanism | `3K₃(1-ξ)^(2/3) / (N₃(1-(1-ξ)^(1/3)))` | `k₃(1-ξ)^n₃` (power law) |
| Combination | `min(max(r_NG, r_I), r_D)` | `min` of all three |
| Parameters | `PK_PARAMS_C3S` …, provenance unestablished | [`PK84_PARAMS_C3S`](@ref) … |
| Status | deprecated, warns once | **use this one** |

`parrot_killoh_avrami` is the **canonical 1984 formulation**, with the parameters
of [ParrotKilloh1984](@cite) as reported by [LothenbachWinnefeld2006](@cite) and
used by [Lavergne2018](@cite). Use it when you want the α(t) curves of the cement
literature. A useful signature of a correct transcription: with those parameters
C₂S has no nucleation–growth stage and C₃S no diffusion-controlled stage.

`parrot_killoh` is kept only so that existing scripts keep running. Its
nucleation–growth term carries no Avrami logarithm, `K₃` sits where the canonical
form has `k₂`, `N₁ = 3.3` is the canonical `n₃`, and `k₃ = 1.1` has no counterpart
at all — so the two are different models, not two parameterizations of one. The
primary source (*British Ceramic Proceedings* **35**, 41–53, 1984) has no DOI and
could not be consulted, so the smoothed variant is not attributed to it. With
`PK_PARAMS_*` the diffusion branch takes over at α ≈ 0.003 (C₂S), 0.013 (C₃S) and
0.057 (C₃A); those three then follow the closed form
`α(t) = α_max·[1 − (1 − √(2·K₃·t / N₃))³]`, giving **α(7 d) = 0.2386 for all
three**, while C₄AF is limited by its own nucleation branch and reaches only
0.193. Converting the four clinker demos to the canonical variant moved a CEM I
paste at w/c = 0.40 from a mean degree of 0.234 to 0.628, an adiabatic
temperature rise from 2.0 to 14.2 °C, and the heat released at seven days from
115 to 308 kJ/kg of cement.

```julia
pk = parrot_killoh_avrami(
    PK84_PARAMS_C3S, "C3S";
    α_max    = powers_alpha_max(0.5),      # Powers water availability limit
    blaine   = 380u"m^2/kg",               # rate ∝ B/385 m²/kg
    humidity = 0.95,                       # β_h reduction, 0 below 80 % RH
)
```

The three corrections are exported separately — [`powers_alpha_max`](@ref),
[`blaine_factor`](@ref) and [`humidity_factor`](@ref) — and multiply the rate.
`humidity` also accepts a callable `t -> h(t)` for a drying history.

Supplementary cementitious materials do not follow Parrot & Killoh at all. Their
pozzolanic or latent-hydraulic reaction follows [`waller`](@ref), a sigmoid in
log-time, with [`WALLER_PARAMS_FLY_ASH`](@ref), [`WALLER_PARAMS_SILICA_FUME`](@ref)
or [`WALLER_PARAMS_SLAG`](@ref).

## [Rate laws that depend on a consumed reactant](@id kinetics-frozen-species)

The ODE state carries the **extents of reaction** `ξ` alongside the kinetic
moles, integrating `dξ/dt = r`. The residual therefore reconstructs *every*
species from `n = n(0) + νᵀ ξ` before evaluating the rates, so a rate closure may
read any amount — including species that are not themselves kinetic.

That is what makes reaction sequencing expressible. Sulfate available for
ettringite, portlandite available for a pozzolanic reaction, a product inhibiting
its own formation: all are written directly.

```julia
# C3A reacts with gypsum while sulfate lasts, then by the sulfate-free route.
# "Gp" is consumed but is not a kinetic species — reading it here is correct.
ε = 1.0e-3
sulfated  = (T, P, t, n, lna, n0) -> pk(T, P, t, n, lna, n0) * n["Gp"] / (n["Gp"] + ε)
unsulfated = (T, P, t, n, lna, n0) -> pk(T, P, t, n, lna, n0) * ε / (n["Gp"] + ε)
```

Keep such a gate **smooth**. A hard `n["Gp"] > 0 ? r : 0` is a discontinuity in
the residual, which a stiff solver will either step over or grind against; the
`x/(x+ε)` form above costs nothing and stays differentiable.

!!! note "With an equilibrium solver, the non-kinetic partition is piecewise constant"
    When re-speciation is active the equilibrium partition is owned by the
    equilibrium solve, which runs once per accepted step by operator splitting.
    Those amounts are therefore refreshed between steps but held constant *within*
    a step, and are not reconstructed from `ξ`. A gate on such a species still
    works; it simply lags by at most one step.

## [Post-processing a kinetics run](@id kinetics-postprocessing)

The ODE state vector carries only the **kinetic** species. Every other amount is
held in a single buffer that the integrator mutates in place, so after a run it
reflects the last accepted step and nothing else. Recovering the composition at an
arbitrary time therefore means replaying the stoichiometry, which is what
[`reaction_extents`](@ref) and [`state_at`](@ref) do.

```julia
sol = integrate(kp, ks)

α  = degrees_of_hydration(sol, kp)          # Dict(symbol => α(t))
ᾱ  = mean_degree_of_hydration(sol, kp)      # mass-weighted binder average

ξ  = reaction_extents(sol, kp)              # extent of each reaction, mol
st = state_at(sol, kp, 28 * 86400.0)        # full ChemicalState at 28 days
```

[`state_at`](@ref) rebuilds **every** species from `n = n₀ + νᵀξ`, so the result
satisfies `A(n - n₀) = (Aνᵀ)ξ` to machine precision — whatever the quadrature
error. How well the elements themselves balance is a property of the reactions
you wrote, not of the reconstruction; check it with `maximum(abs, A * ν')`.

!!! warning "Re-speciation is not replayed"
    When the run used an equilibrium solver, the aqueous partition was
    re-speciated at every accepted step, and that redistribution cannot be
    recovered from the stoichiometry. `state_at` then returns the purely kinetic
    reconstruction; call [`equilibrate`](@ref) on it to speciate.

### From moles to volume fractions

[`volume_fractions`](@ref) turns a state into the input a mean-field
homogenization scheme consumes, using the standard molar volumes `V⁰` that
CEMDATA18 supplies for every cement phase.

```julia
groups = [
    "anhydrous" => ["C3S", "C2S", "C3A", "C4AF"],
    "C-S-H"     => "Jennite",
    "CH"        => "Portlandite",
    "AFt"       => "ettringite",
    "water"     => "H2O@",
]
f = volume_fractions(st, groups; reference = state0)
```

Passing `reference` selects the **sealed-curing** convention of
[Lavergne2018](@cite): fractions are referred to the initial volume, held fixed,
and the deficit left by the reactions appears as a `"void"` phase. That void is
the chemical shrinkage — hydration products occupy less space than the reactants
they consume — and it is a genuine phase of the microstructure, not a rounding
error. Without `reference`, fractions are referred to the current volume and the
void does not exist.

```julia
sum(values(f))                              # 1.0, void included
f["void"]                                   # Le Chatelier contraction
chemical_shrinkage(st, state0)              # the same thing, as a volume
```

!!! note "Species without a molar volume are invisible"
    A species carrying no `V⁰` contributes nothing to [`volume`](@ref),
    [`porosity`](@ref) or [`volume_fractions`](@ref) — a custom species added by
    hand, for instance. Call [`missing_molar_volumes`](@ref) before trusting a
    volume balance:

    ```julia
    isempty(missing_molar_volumes(st)) || @warn "incomplete volume balance"
    ```

Individual fractions may be slightly **negative** for aqueous solutes, whose
standard partial molar volumes are negative (electrostriction contracts the
solvent around an ion). They are kept so that `volume_fractions` and
[`volume`](@ref) always agree; grouping folds them back into the liquid.
