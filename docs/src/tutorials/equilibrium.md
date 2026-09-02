# [Chemical Equilibrium](@id sec-equilibrium)

ChemistryLab computes thermodynamic equilibrium by minimizing the Gibbs free energy of the system subject to element-conservation constraints. The workflow always follows the same four steps:

1. Build a [`ChemicalSystem`](@ref) (species + stoichiometric matrix).
2. Create an initial [`ChemicalState`](@ref) (temperature, pressure, initial amounts).
3. Call [`equilibrate`](@ref) (or use [`EquilibriumSolver`](@ref) explicitly).
4. Inspect the resulting [`ChemicalState`](@ref).

---

## Minimal workflow

The convenience function [`equilibrate`](@ref) handles everything with sensible defaults.
The example below computes the equilibrium state of calcite (CaCO₃) dissolving in mildly acidic water — a standard geochemical benchmark.

```@example eq_setup
using Optimization, OptimizationIpopt
using ChemistryLab
using DynamicQuantities

substances = build_species(datapath("slop98-inorganic-thermofun.json"))

# Select the carbonate-system species, calcite and its dissolution product Ca²⁺
dict = Dict(symbol(s) => s for s in substances)
species = [dict[sym] for sym in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")]

cs = ChemicalSystem(species, ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
```

```@example eq_setup
state = ChemicalState(cs)

# 1 mmol calcite dissolved in 1 L of acidic water (initial pH ≈ 4)
set_quantity!(state, "Cal",  1e-3u"mol")
set_quantity!(state, "H2O@", 1.0u"kg")

V = volume(state)
set_quantity!(state, "H+",  1e-4u"mol/L" * V.liquid)   # pH = 4
set_quantity!(state, "OH-", 1e-10u"mol/L" * V.liquid)  # charge seed

state_eq = equilibrate(state)
```

!!! tip "Quick shortcut"
    Calling `equilibrate(state)` with no extra arguments uses sensible defaults and is usually sufficient for aqueous geochemical problems.

---

## Choosing a solver

ChemistryLab provides two solver extensions. Load whichever fits your workflow:

| Extension packages | Default solver | When to use |
|:--|:--|:--|
| `Optimization`, `OptimizationIpopt` | `IpoptOptimizer` | general-purpose, robust |
| `OptimaSolver` | `OptimaOptimizer` | preferred when available |

When both are loaded, `OptimaSolver` provides the default single back end. On the
cement equilibria this package targets both reach an element balance of the order
of `1e-14`, and OptimaSolver is 3 to 26 times faster.

### `equilibrate(state)` does not pick one of them — it proves the answer

Neither back end is reliable on its own, so the one-argument `equilibrate` solves
by **every** route available and keeps the answer
[`optimality_certificate`](@ref) proves optimal. Measured, with the element
balance judged row by row against each row's own budget:

| case | interior point | dual Newton | certified route |
|:--|--:|--:|--:|
| calcite in water, 10–40 °C | **3.0e-2** | 3.0e-12 | certified |
| calcite + 1 mmol CO₂ | 1.0e-3 | 6.3e-13 | certified |
| calcite + 50 mmol CO₂ | 1.3e-16 | 8.5e-14 | certified |
| pure water | 1.3e-16 | 1.3e-16 | certified |
| CEM I paste, w/c 0.45 and 0.60 | 1.5e-14 | 2.0e-14 | certified |
| CEM I paste, w/c 0.30 | 7.2e-16 | 3.3e-10, **not** certified | certified |

Ten of ten certified through both routes; seven through the interior point alone,
nine through the dual Newton alone. The interior point's 3 % on the first row is
not a tolerance: it is a **charge balance wrong in the second digit**, and the
reason it cannot be improved is that the fraction-to-boundary rule caps its step
at 15 % of a correction the next iteration re-poses — traced over twenty-three
iterations with the residual frozen at 3.0e-6 and `‖dn‖` decaying at `1 − α`.

```julia
state_eq = equilibrate(state)                   # certified (the default)
state_eq = equilibrate(state; certify = false)  # one back end, as before
eq, cert = equilibrate_certified(state)         # when the proof itself is wanted
```

Pass a solver explicitly — `equilibrate(state, OptimaOptimizer())` — to use that
one back end and nothing else.

## Proving that an answer is the answer

An interior-point method minimizes `G` by walking the interior of the feasible
set, and on a cement equilibrium it stops on `MaxIters` — at any tolerance.
Whether the point it returns is the minimum is then an open question, and the
package can now settle it rather than assume it.

### Why the question has an answer

Write the Gibbs energy in `RT` units as `G(n) = Σᵢ nᵢ μᵢ(n)`. Its ideal part

```math
\varphi(n) = \sum_i n_i \ln\frac{n_i}{N}, \qquad N = \sum_j n_j
```

has Hessian ``\operatorname{diag}(1/n_i) - \tfrac1N \mathbf{1}\mathbf{1}^\top``,
and for any ``v``

```math
v^\top \nabla^2\varphi\, v = \sum_i \frac{v_i^2}{n_i} - \frac{1}{N}\Bigl(\sum_i v_i\Bigr)^2 \;\ge\; 0
```

by Cauchy–Schwarz applied to ``v_i = (v_i/\sqrt{n_i})\sqrt{n_i}``. A pure phase
has unit activity, so it contributes a term **linear** in its amount. Hence `G`
is convex, the feasible set ``\{An=b,\ n\ge 0\}`` is a polyhedron, and — the
constraints being affine, so that the linearity constraint qualification holds
everywhere — the KKT conditions are **necessary and sufficient**.

Two consequences follow. The minimizer is unique, so a solver returning different
answers from different starting points is not finding local minima but stopping
short of stationarity. And optimality can be *checked*: a composition satisfying
the KKT conditions is proved globally optimal.

### The certificate

[`optimality_certificate`](@ref) checks the three conditions, on any composition
and whatever produced it. Writing ``u = -A^\top y`` for the element potentials:

| condition | on which species | meaning |
|:--|:--|:--|
| ``\mu_i + (A^\top y)_i = 0`` | interior (`n > floor`) | stationarity |
| ``An = b`` | — | conservation of matter |
| ``u_i \le g_i`` | a **pure** phase at its bound | that phase undersaturated |
| ``\ln \sum_i \exp(u_i - g_i - \ln\gamma_i) \le 0`` | a **mixing** phase held entirely absent | that solution cannot form |

The last row is Michelsen's tangent-plane measure, and it is a separate test
because a mixing phase needs one: its members are never exactly zero while it
exists, so they are neither interior nor at a bound, and a solid solution left out
of the assemblage used to pass the certificate **unexamined**. The trial
composition is refined against the phase's own activity model, so the test is not
the ideal approximation.

Two subtleties decide whether the check is meaningful.

A species **at its bound** obeys the inequality, not the equality. Imposing the
equality on an amount held at `1e-16` whose mass-action value is `e⁻³⁰⁰`
misstates its log-activity by 263 `RT` units, and the check then reports a
residual of 74 for a composition solved to `5e-12`.

A species carrying a **vanished component** is absent by the *constraint*, not by
thermodynamics, and its saturation index is meaningless — the element potential
of a component nobody supplies is determined by nothing. The test for that is not
`bₖ ≈ 0` but `bₖ ≈ 0` **with the non-zero entries of row `k` sharing a sign**:
only then does ``\sum_i A_{ki} n_i = 0`` with ``n \ge 0`` force each term to
vanish. The `H⁺` row carries `+1` for `H⁺` and `−1` for `OH⁻`, so its zero total
is the ordinary state of pure water; treating it as degenerate kills the entire
acid–base system and returns pH 7.000 with the calcite undissolved.

### The certifying solver

[`DualEquilibriumSolver`](@ref) solves the KKT system directly, in element
potentials. From ``\mu_i + (A^\top y)_i = 0`` an aqueous species obeys the
mass-action law ``a_i = \exp(u_i - g_i)``, and a pure phase is present exactly
when ``u_i = g_i``, absent when undersaturated — the classical phase-stability
criterion.

Two levels. The inner one inverts the **solutes'** mass-action laws at fixed
potentials and fixed solvent amount; the outer is a Newton on ``1 + m + |P|``
unknowns — the solvent, the `m` element potentials, and the amounts of the
active phases. Parameterizing the solutes by ``\ln n`` makes their positivity
automatic, which is what removes the fraction-to-boundary limit that caps the
interior-point step at every iteration.

The solvent is deliberately **not** inverted through its own mass-action law: its
activity is a mole fraction, so ``\ln a_w \le 0`` always, and an arbitrary `y` can
demand more, for which no finite composition exists. It belongs to the outer
system, where the balance determines it.

```julia
des  = DualEquilibriumSolver(cs, HKFActivityModel())
ipm  = equilibrate(state, OptimaOptimizer())      # into the neighbourhood
dual = solve(des, ipm; b = b)                     # to the KKT conditions
cert = optimality_certificate(des, dual; b = b)
cert.optimal    # true: a proof, for a convex problem
```

!!! note "What it buys, measured"
    On calcite in pure water the certified pH is **9.90** against an
    interior-point 6.96 — not an imprecision but a wrong answer, and one nothing
    in that solver's output reveals. On the Reaktoro reference (calcite, CO₂ and
    water) both routes now agree with Reaktoro on every species: above `10⁻⁵` mol
    to `10⁻³` relative, the trace ions to 5 %, the worst being `CaOH⁺` at ×1.032
    on 1.6 nmol. That reference used to carry a `@test_broken` for `CaOH⁺` at
    ×2.47; what closed it was the convergence test moving to the true KKT error at
    `μ = 0` (`OptimaSolver` 0.4.1), and `test/equilibrium_reference.jl` is now 26
    plain assertions.

    [`speciated_states`](@ref) certifies every instant it replays and names any it
    cannot. On a full ordinary Portland cement over 28 days, all forty replayed
    instants are certified, with element balances between `1e-11` and `1e-13` mol.

## [Constraints other than fixed T and P](@id sec-equilibrium-constraints)

A closed system at given temperature and pressure is one case among several. A
vessel that exchanges no heat has its temperature *determined* by the reaction; a
titration holds the pH and lets the amount of acid follow. Both are equilibrium
problems, and both are stated by saying what is held and what is unknown.

!!! compat "Needs OptimaSolver 0.5"
    The blocks in this section are shown rather than executed: the documentation
    environment resolves `OptimaSolver` from the registry, and the parameter block
    these constraints ride on arrived in 0.5.0. The figures quoted below are the
    ones `test/equilibrium_constraints.jl` asserts.

```julia
using ChemistryLab, DynamicQuantities, OptimaSolver

# The same calcite system, but adiabatic: no heat leaves, so T is an unknown
st = ChemicalState(cs)
set_quantity!(st, "Cal",  1e-2u"mol")
set_quantity!(st, "H2O@", 1.0u"kg")
set_quantity!(st, "H+",   1e-2u"mol")
set_quantity!(st, "OH-",  1e-10u"mol")

des   = DualEquilibriumSolver(cs, DiluteSolutionModel())
eq_ad = solve(des, st; constraint = Adiabatic())

temperature(eq_ad)      # 298.1996 K — the reaction warms it by 0.0496 K
pH(eq_ad)               # 6.4323, against 6.4329 at fixed temperature
```

| constraint | held | unknown | vehicle |
|:--|:--|:--|:--|
| [`FixedTP`](@ref) (default) | `T`, `P` | — | — |
| [`Adiabatic`](@ref) | the enthalpy of the initial state | `T` | parameter |
| [`FixedEnthalpy`](@ref) | a prescribed `H` | `T` | parameter |
| [`FixedVolume`](@ref) | a prescribed `V` | `P` | parameter |
| [`SealedVolume`](@ref) | the volume of the initial state | `P` | parameter |
| [`FixedpH`](@ref) | `−log₁₀ a(H⁺)` | the titrant amount | column |
| [`FixedActivity`](@ref) | `a` of any species | the titrant amount | column |

### The two vehicles

They are not interchangeable, and the distinction is the same one Reaktoro draws.

A **prescribed property** — an enthalpy, a volume — adds one *parameter* and one
equation to the solver's own square system. Nothing loops around the equilibrium
solve: the temperature is an unknown beside the amounts and the element
potentials, and the enthalpy balance is one more row of the same Newton system.
So an adiabatic solve costs one equation, not a solve per trial temperature.

A **prescribed chemical potential** — a pH, an activity — adds one *column* to the
conservation matrix instead: an unknown amount of a titrant the system may draw
on. The linear rows become `A n − A[:, titrant] q = b`, so the system is open to
that one substance and closed to everything else, and `q` comes back as part of
the answer — it is the reagent consumed, which is what a titration measures.

```julia
q = Ref(Float64[])
eq_ph = solve(des, st; constraint = FixedpH(7.0), parameters = q)
q[][1]      # 3.4673e-3 mol of H+ had to be added to hold pH 7
```

On 10 mmol of calcite in a kilogram of water, whose free pH is 9.90, holding the
pH takes more acid the lower the target, and at pH 6 the calcite is gone
entirely:

| pH held | H⁺ added (mol) | calcite left (mol) |
|--:|--:|--:|
| 6 | 1.689e-2 | 0 |
| 7 | 3.467e-3 | 7.063e-3 |
| 8 | 8.648e-4 | 9.149e-3 |
| 9 | 2.504e-4 | 9.728e-3 |

### What is checked, and what is refused

The adiabatic answer is validated two ways. Against physics: `H⁺ + OH⁻ → H₂O`
comes out at **−55.85 kJ/mol** at every amount tested, against the accepted
−55.8, with nothing fitted to it — the enthalpies come from the database and the
temperature is an unknown of the system. And against itself: solving at a *fixed*
temperature equal to the one the adiabatic solve found returns the same
composition to `1e-12`.

A volume constraint on a **condensed** system is **refused**, with the lever it
measured named in the error. The molar volumes of water and of the minerals in the
shipped databases are exactly pressure-independent — `V⁰(1 bar) = V⁰(100 bar)` to
the last bit for `H2O@` and `Cal` — and only a few aqueous ions vary, `OH⁻` by 8 %
over 100 bar. The relative lever `(∂V/∂P)·P/V` is then about `1e-6`, meaning some
**9 600 bar to change the volume by one percent**. That is the physics, not a
solver limitation: the volume of an incompressible condensed system is fixed by
its composition. Declare a gas phase to give the pressure something to do. To
report the volume change of a sealed specimen at fixed pressure — what a hydrating
binder actually needs — use [`porosity`](@ref) with a reference state instead.

!!! note "Constraints need the certifying solver"
    The prescribed property or potential is an unknown of the dual Newton's own
    system, so a constraint other than `FixedTP` requires `OptimaSolver` and a
    system with an aqueous phase and `H2O@`. Asking for one otherwise raises,
    rather than being silently ignored. `equilibrate(state; constraint = ...)`
    routes through the certified path, so the answer still comes with its proof.

## Differentiating an equilibrium

A `ChemicalState` carries whatever number type its amounts have, so a
composition built from `ForwardDiff.Dual` values propagates through the
speciation, and `pH`, `pOH`, `porosity` and `saturation` come back as duals too.

Crossing the **solve** works as well, and without asking any solver to iterate on
dual numbers — Ipopt is a C library and never could. The equilibrium is solved
once at the primal values, and the sensitivities come from the optimality
conditions, the implicit-function-theorem route:

```math
\begin{bmatrix} H & A^{\mathsf T} \\ A & 0 \end{bmatrix}
\begin{bmatrix} \dot n \\ \dot y \end{bmatrix}
=
\begin{bmatrix} -\partial_\theta \nabla G \\ \dot b \end{bmatrix},
\qquad H = \nabla^2 G(n^\star),
```

restricted to the species actually present. One factorization serves every
partial derivative, and the result is exact — no step size to choose.

```julia
using ChemistryLab, DynamicQuantities, ForwardDiff, OptimaSolver

f(x) = begin
    n = Any[fill(0.0u"mol", length(cs.species))...]
    n[i_h2o] = 55.5u"mol";  n[i_cal] = 0.05u"mol";  n[i_co2] = x * u"mol"
    eq = equilibrate(ChemicalState(cs, n), OptimaOptimizer())
    ustrip(us"mol", eq.n[i_ca])
end

ForwardDiff.derivative(f, 0.01)     # → 0.15193
```

!!! note "Why the complementarity conditions cannot be skipped"
    The stationarity conditions are `∇G − Aᵀy − z = 0`, `A n = b`, `nᵢzᵢ = 0`,
    with `z ≥ 0` the stability multipliers. The last block partitions the
    species, and dropping it does not degrade the answer gently — it destroys
    it. On calcite + CO₂ in water with a gas phase declared, the unreduced
    system puts the whole response into the **absent** gas species
    (`n = 5×10⁻¹¹` mol), returning a sensitivity that satisfies the element
    balance to `4×10⁻¹⁶` and means nothing.

    No back-end returns `z`, so the active set is recovered internally: a
    species negligible on the scale of the system that nonetheless takes a
    leading share of the response is pinned, and the system re-solved.

!!! note "Verified against Reaktoro"
    On calcite + CO₂ + water, `∂n/∂(CO₂)` from this route agrees with the
    package's own finite differences to `9×10⁻⁵` — the finite-difference
    truncation error — and with [Reaktoro](https://reaktoro.org) 2.13 reading
    **the same Cemdata18 file**, over the same eleven species, under the same
    (ideal) activity model:

    | species | ChemistryLab (AD) | Reaktoro (FD) | rel. diff |
    |:--|--:|--:|--:|
    | H₂O | −0.181336 | −0.181397 | 3.4×10⁻⁴ |
    | Ca²⁺ | +0.151920 | +0.151987 | 4.4×10⁻⁴ |
    | HCO₃⁻ | +0.333308 | +0.333427 | 3.6×10⁻⁴ |
    | CO₂(aq) | +0.818655 | +0.818600 | 6.7×10⁻⁵ |
    | Ca(HCO₃)⁺ | +0.029333 | +0.029338 | 1.6×10⁻⁴ |
    | calcite | −0.181251 | −0.181324 | 4.0×10⁻⁴ |

    The equilibrium amounts agree to the same order (`Ca²⁺` 3.5281×10⁻³ against
    3.52902×10⁻³). Reaktoro's own spread across `h ∈ {10⁻³, 10⁻⁴, 10⁻⁵}` is
    `7.2×10⁻⁴`, so the residual difference sits below the oracle's truncation
    error on every species.

    The absent gas species gets exactly zero from the active-set treatment,
    against `2×10⁻⁹` by finite differences.

!!! warning "A cross-code comparison has three knobs, not one"
    Database, species list and activity model all have to match, and each is
    worth tens of percent here. Reading Cemdata18 in both codes but leaving
    Reaktoro on its HKF activity model against this package's default
    `DiluteSolutionModel()` moves `∂Ca²⁺/∂(CO₂)` from `+0.1520` to `+0.2179` —
    a 35 % gap that says nothing about either code.

    The species list matters just as much. Dropping the aqueous calcium
    complexes leaves free `Ca²⁺` as the only aqueous home for calcium, so
    `∂Ca²⁺/∂(CO₂)` and `∂calcite/∂(CO₂)` mirror each other exactly. Restoring
    `CaOH⁺`, `Ca(CO₃)@` and `Ca(HCO₃)⁺` breaks that mirror — the difference is
    what `Ca(HCO₃)⁺` takes up — and **both codes break it the same way**. The
    element balance closes to `2×10⁻¹⁶` either way.

### Explicit solver (always works)

Pass the solver as the **second positional argument**:

```julia
using Optimization, OptimizationIpopt
state_eq = equilibrate(state, IpoptOptimizer())

using OptimaSolver
state_eq = equilibrate(state, OptimaOptimizer())
```

### Default shortcut

When only one extension is loaded:

```julia
using Optimization, OptimizationIpopt
state_eq = equilibrate(state)   # → IpoptOptimizer

using OptimaSolver
state_eq = equilibrate(state)   # → OptimaOptimizer
```

With both loaded, `OptimaSolver` always wins:

```julia
using Optimization, OptimizationIpopt
using OptimaSolver
state_eq = equilibrate(state)   # → OptimaOptimizer (priority)
```

---

## Inspecting the equilibrium state

The returned [`ChemicalState`](@ref) carries all derived thermodynamic quantities:

```@example eq_setup
println("pH      = ", pH(state_eq))
println("pOH     = ", pOH(state_eq))
println("porosity   = ", porosity(state_eq))
println("saturation = ", saturation(state_eq))
```

Phase volumes and mole amounts are accessible via named tuples:

```@example eq_setup
v = volume(state_eq)
println("V liquid = ", v.liquid)
println("V solid  = ", v.solid)
println("V total  = ", v.total)

m = moles(state_eq)
println("n liquid = ", m.liquid)
println("n solid  = ", m.solid)
```

Individual species amounts (in mol):

```@example eq_setup
cs_eq = state_eq.system
for (i, sp) in enumerate(cs_eq.species)
    n_i = state_eq.n[i]
    println(rpad(symbol(sp), 20), ustrip(n_i), " mol")
end
```

---

## Scaling and normalization

It is often useful to express a composition relative to a reference amount — per mole, per kilogram, or per cubic meter of system. Two mechanisms are provided.

### Scalar multiplication

A [`ChemicalState`](@ref) can be multiplied or divided by a real number. All molar amounts are scaled proportionally; temperature, pressure, and the chemical system are unchanged. The operation is **non-mutating** — a new state is returned:

```julia
state2  = state_eq * 2.0    # double all amounts
state_m = state_eq / 1000   # millimolar scale
```

### `rescale!` — rescale to a target total

[`rescale!`](@ref) scales all molar amounts **in-place** so that the total of the matching physical quantity equals `target`:

| `target` dimension | Quantity brought to `target` |
|:-------------------|:-----------------------------|
| mol                | `moles(state).total`         |
| kg (mass)          | `mass(state).total`          |
| m³ (volume)        | `volume(state).total`        |

All derived quantities (pH, porosity, volume, …) are recomputed automatically after scaling.

```@example eq_setup
# Express the equilibrium composition per kilogram of total system
state_pkg = copy(state_eq)
rescale!(state_pkg, 1.0u"kg")

println("Ca²⁺ = ", moles(state_pkg, "Ca+2"), "  mol/kg")
println("pH   = ", pH(state_pkg))   # intensive quantities are invariant
```

!!! note "Intensive quantities"
    pH, porosity, and saturation are **intensive** — they are invariant under homothety and remain unchanged after `rescale!` or scalar multiplication.

---

## Controlling the solver

### Variable space: `:linear` vs `:log`

[`equilibrate`](@ref) accepts a `variable_space` keyword that selects the optimization variable space:

| `variable_space`        | Variables | Recommended when |
|:-----------------|:----------|:----------------|
| `Val(:linear)`   | mole amounts `nᵢ ≥ 0` | most systems, default |
| `Val(:log)`      | `log nᵢ` | systems spanning many orders of magnitude |

```julia
state_eq_log = equilibrate(state; variable_space=Val(:log))
```

!!! warning "Convergence"
    Solving a system of equations in chemistry can be a difficult undertaking. The orders of magnitude can vary greatly, and convergence is not guaranteed.

---

### Tolerances

Tighter tolerances are passed directly as keyword arguments and forwarded to the underlying Ipopt solver:

```julia
state_eq_tight = equilibrate(state; abstol=1e-12, reltol=1e-12)
```

---

## Using `EquilibriumSolver` explicitly

For batch calculations where many different initial states share the same system and activity model, construct an [`EquilibriumSolver`](@ref) once and reuse it:

```julia
using Optimization, OptimizationIpopt

opt = IpoptOptimizer(
    acceptable_tol        = 1e-12,
    dual_inf_tol          = 1e-12,
    acceptable_iter       = 1000,
    constr_viol_tol       = 1e-12,
    warm_start_init_point = "no",
)

solver = EquilibriumSolver(
    cs,
    DiluteSolutionModel(),
    opt;
    variable_space = Val(:linear),
    abstol  = 1e-10,
    reltol  = 1e-10,
)
```

Once built, `solver` is called with any compatible `ChemicalState`:

```@example eq_setup
using Optimization, OptimizationIpopt #hide

opt = IpoptOptimizer( #hide
    acceptable_tol        = 1e-12, #hide
    dual_inf_tol          = 1e-12, #hide
    acceptable_iter       = 1000, #hide
    constr_viol_tol       = 1e-12, #hide
    warm_start_init_point = "no", #hide
) #hide

solver = EquilibriumSolver( #hide
    cs, #hide
    DiluteSolutionModel(), #hide
    opt; #hide
    variable_space = Val(:linear), #hide
    abstol  = 1e-10, #hide
    reltol  = 1e-10, #hide
) #hide
state_eq2 = solve(solver, state)
```

!!! note "Performance"
    The potential function `μ(n, p)` is compiled once during `EquilibriumSolver` construction. Repeated calls to `solve(solver, ...)` with different states reuse it, avoiding redundant compilation overhead.

---

## Temperature dependence (10–30 °C)

Calcite solubility varies with temperature. Using the `solver` built above, we sweep from 10 to 30 °C and track pH, dissolved calcium and remaining solid calcite:

```@example eq_setup
using Plots


temperatures = 10:30   # °C

pH_vals   = Float64[]
nCa_vals  = Float64[]  # mmol
nCal_vals = Float64[]  # mmol

i_Ca  = findfirst(sp -> symbol(sp) == "Ca+2", cs.species)
i_Cal = findfirst(sp -> symbol(sp) == "Cal",  cs.species)

# Start from the charged state built above, not from `ChemicalState(cs)`:
# a fresh state holds no matter at all, so every element balance would be zero
# and the sweep would return the same trivial solution at all 21 temperatures.
for θ in temperatures
    s_T = deepcopy(state)
    set_temperature!(s_T, (273.15 + θ) * u"K")
    s_eq = solve(solver, s_T)
    push!(pH_vals,   pH(s_eq))
    push!(nCa_vals,  ustrip(s_eq.n[i_Ca]) * 1e3)
    push!(nCal_vals, ustrip(s_eq.n[i_Cal]) * 1e3)
end
```

The figures can then be drawn.

```@example eq_setup
p1 = plot(collect(temperatures), pH_vals,
    xlabel = "T (°C)", ylabel = "pH", label = "pH",
    marker = :circle, linewidth = 2, title = "pH")
p2 = plot(collect(temperatures), nCa_vals,
    xlabel = "T (°C)", ylabel = "n (mmol)", label = "Ca²⁺",
    marker = :circle, linewidth = 2, title = "Dissolved species")
plot!(p2, collect(temperatures), nCal_vals,
    label = "Cal", marker = :square, linewidth = 2)
plot(p1, p2, layout = (1, 2), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm, size = (900, 400))
```

!!! note "Retrograde Kₛₚ, and yet more dissolved calcium"
    Calcite is **retrograde soluble**: its solubility product falls as temperature
    rises. The sweep above reproduces that — the ionic product ``[\text{Ca}^{2+}]
    [\text{CO}_3^{2-}]`` it settles on goes from ``10^{-8.411}`` at 10 °C to
    ``10^{-8.517}`` at 30 °C, matching the ``K_{sp}`` the database itself gives to
    within 0.003 log units, and ``-8.48`` at 25 °C is the accepted value for
    calcite.

    The dissolved calcium nevertheless **increases**, from 0.146 to 0.158 mmol/L.
    That is not a contradiction: this system is **closed**, with 1 mmol of total
    carbonate and no CO₂ reservoir. As temperature rises the pH drops from 9.83
    to 9.43, which shifts carbonate to bicarbonate and takes free CO₃²⁻ from 26.5
    down to 19.2 µmol/L — more than enough to offset the smaller ``K_{sp}``. The
    familiar statement that retrograde solubility means *less* dissolved calcium
    holds for a system buffered at a fixed CO₂ partial pressure, not for a sealed
    one. Fix the partial pressure instead (add a gas phase) and the sign flips.

---

## Activity models

All activity models inherit from [`AbstractActivityModel`](@ref). Three built-in
models are provided, covering ideal behavior through to the extended Debye-Hückel
level used by standard geochemical codes.

### Choosing a model

| Model | Formula | Valid range | Parameters needed |
| :-- | :-- | :-- | :-- |
| [`DiluteSolutionModel`](@ref) | Raoult / Henry | I ≪ 1 mol/kg | none |
| [`HKFActivityModel`](@ref) | B-dot extended Debye-Hückel | I ≲ 1 mol/kg | `A`, `B`, `Ḃ` (defaults at 25 °C) |
| [`DaviesActivityModel`](@ref) | Davies equation | I ≲ 0.5 mol/kg | `A`, `b` (defaults at 25 °C) |

---

### `DiluteSolutionModel` (ideal dilute solution)

| Phase | Law | Expression |
| :-- | :-- | :-- |
| Solvent (H₂O) | Raoult | `ln a = ln xₛ` |
| Aqueous solutes | Henry | `ln a = ln(cᵢ / c°)`, `c° = 1 mol/L` |
| Crystals | Pure solid | `ln a = 0` |
| Gas | Ideal mixture | `ln a = ln xᵢ` |

```julia
state_eq = equilibrate(state)   # DiluteSolutionModel is the default
```

---

### `HKFActivityModel` (extended Debye-Hückel B-dot)

Implements the extended Debye-Hückel model of Helgeson (1969) and
Helgeson, Kirkham & Flowers (1981), identical to the model used by
PHREEQC and EQ3/6.

**Ion activity coefficient:**
```
log₁₀ γᵢ = −A zᵢ² √I / (1 + B åᵢ √I)  +  Ḃ I
```

**Neutral aqueous species (salting-out):**
```
log₁₀ γᵢ = Kₙ I
```

**Water activity** is computed from the osmotic coefficient via Gibbs-Duhem
(not Raoult), which is accurate up to `I ≈ 1 mol/kg`.

**Ionic radius lookup** (priority order):
1. `sp[:å]` — explicit value set in species properties.
2. [`REJ_HKF`](@ref) — Helgeson et al. (1981) Table 3 (27 common ions).
3. [`REJ_CHARGE_DEFAULT`](@ref) — fallback by formal charge.
4. `model.å_default` (default: 3.72 Å).

**Usage:**

```julia
# Fixed A, B at 25 °C / 1 bar (fast — suitable for isothermal calculations)
state_eq = equilibrate(state; model=HKFActivityModel())

# Temperature-dependent A and B (recomputed from T, P at each equilibrium solve)
state_eq = equilibrate(state; model=HKFActivityModel(temperature_dependent=true))

# Custom parameters
model = HKFActivityModel(A=0.52, B=0.33, Ḃ=0.04)
```

The A and B parameters depend on the water dielectric constant and density
and can be computed explicitly via [`hkf_debye_huckel_params`](@ref):

```julia
ab = hkf_debye_huckel_params(298.15, 1e5)   # → (A=0.5114, B=0.3288)
```

!!! note "Valid range"
    The B-dot model is reliable for `I ≲ 1 mol/kg`. For higher ionic
    strengths (brines, evaporites), use the Pitzer model (planned future extension).

---

### `DaviesActivityModel` (Davies equation)

Simpler alternative with no species-specific ionic radii. Suitable when
ionic radii data are unavailable or for rapid screening calculations.

**Ion activity coefficient:**
```
log₁₀ γᵢ = −A zᵢ² (√I / (1 + √I)  −  b I)
```

**Water activity** uses the Raoult (mole fraction) approximation.

```julia
state_eq = equilibrate(state; model=DaviesActivityModel())

# Temperature-dependent A
state_eq = equilibrate(state; model=DaviesActivityModel(temperature_dependent=true))
```

---

### Custom activity models

To implement a custom activity model, define a new subtype and extend `activity_model`:

```julia
struct MyModel <: AbstractActivityModel
    # model parameters
end

function ChemistryLab.activity_model(cs::ChemicalSystem, ::MyModel)
    # Precompute species indices and constants here (called once)
    idx_solvent = only(cs.idx_solvent)
    # ...

    # Return a closure lna(n, p) -> Vector compatible with ForwardDiff
    function lna(n::AbstractVector, p)
        # p contains at minimum: p.ΔₐG⁰overT, p.T, p.P, p.ϵ
        # n is dimensionless mole vector, same indexing as cs.species
        out = zeros(eltype(n), length(n))
        # ... fill log-activities ...
        return out
    end
    return lna
end
```

Pass your model to `equilibrate` or `EquilibriumSolver`:

```julia
state_eq = equilibrate(state; model=MyModel(...))
```

---

## Solid solutions

Pure crystalline species have activity `ln a = 0`. **Solid solutions** are mineral phases
with variable composition (e.g. C-S-H, AFm, hydrogarnet), where the activity of each
end-member depends on its mole fraction within the phase.

### Defining end-members and phases

End-member species must carry `aggregate_state = AS_CRYSTAL`.
[`SolidSolutionPhase`](@ref) automatically requalifies any end-member whose class is not
already `SC_SSENDMEMBER`, so database species with `SC_COMPONENT` can be passed directly.

**Workflow A — pass database species directly:**

```julia
using ChemistryLab

substances = build_species(datapath("cemdata18-thermofun.json"))
dict = Dict(symbol(s) => s for s in substances)

# SolidSolutionPhase requalifies SC_COMPONENT → SC_SSENDMEMBER automatically
ss_afm = SolidSolutionPhase("AFm", [dict["Ms"], dict["Mc"]])
```

**Workflow B — automated via [`build_solid_solutions`](@ref) and a TOML file:**

```julia
# Load all phases defined in the TOML
ss_phases = build_solid_solutions(datapath("solid_solutions.toml"), dict)
```

See the [Databases](@ref sec-databases) tutorial for the TOML format and the
pre-built `data/solid_solutions.toml` file shipped with ChemistryLab.

Then pass `solid_solutions` as a keyword to `ChemicalSystem`:

```julia
cs = ChemicalSystem(
    [H2O_sp, dict["Ms"], dict["Mc"], ...],
    ["H2O@", "Al+3", ...];           # primaries
    solid_solutions = [ss_afm],      # or solid_solutions = ss_phases
)
```

### Activity models for solid solutions

| Model | Formula | Notes |
| :-- | :-- | :-- |
| [`IdealSolidSolutionModel`](@ref) | `ln aᵢ = ln xᵢ` | Default, any number of end-members |
| [`RedlichKisterModel`](@ref) | `ln aᵢ = ln xᵢ + ln γᵢ` (Margules) | Binary only (2 end-members), parameters in J/mol |

The solid-solution activity is computed **inside** the aqueous activity closure — no
separate activity model is needed. The existing `equilibrate(state)` call handles
solid solutions automatically.

### Ideal solid solution

```julia
ss = SolidSolutionPhase("AFm", [em_ms, em_mc])   # IdealSolidSolutionModel() by default
cs = ChemicalSystem([...]; solid_solutions=[ss])
state_eq = equilibrate(state)
```

### Non-ideal binary: Redlich-Kister

```julia
# Interaction parameters for monosulfoaluminate-monocarboaluminate (example values)
rk = RedlichKisterModel(a0 = 3000.0, a1 = 500.0)          # a2 defaults to 0.0
# or 3-parameter:  RedlichKisterModel(a0 = 3000.0, a1 = 500.0, a2 = 50.0)
ss = SolidSolutionPhase("AFm", [em_ms, em_mc]; model=rk)
```

Activity coefficients (Guggenheim / ThermoCalc convention):

```math
\begin{aligned}
\ln \gamma_1 &= \frac{x_2^2}{RT}\bigl[a_0 + a_1(3x_1 - x_2) + a_2(x_1 - x_2)(5x_1 - x_2)\bigr] \\[4pt]
\ln \gamma_2 &= \frac{x_1^2}{RT}\bigl[a_0 - a_1(3x_2 - x_1) + a_2(x_2 - x_1)(5x_2 - x_1)\bigr]
\end{aligned}
```

!!! note "Valid range"
    `RedlichKisterModel` requires exactly 2 end-members. For ternary or
    higher-order solid solutions, use the ideal model (`IdealSolidSolutionModel`).

!!! note "Integration with aqueous models"
    Solid-solution activities are computed independently of the aqueous activity model.
    You can combine `HKFActivityModel()` for the aqueous phase with any solid-solution
    model — the same `equilibrate` call handles both.
