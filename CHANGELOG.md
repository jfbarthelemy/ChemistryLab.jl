# Changelog

## v0.5.0 — the bridge from a kinetics run to a microstructure

Everything here is additive except the Rosenbrock fix below; no existing API changes.

### From moles to volume fractions

The pieces to turn a hydration run into the input of a mean-field homogenization
scheme were all present — `V⁰` from CEMDATA18, `volume(state, species)`,
`porosity` — but nothing joined them, and the sealed-volume balance was missing
entirely.

- **`volume_fractions(state)`** and **`volume_fractions(state, groups)`**, the
  latter aggregating species into the phase families a homogenization scheme
  consumes (`"C-S-H"`, `"AFt"`, `"anhydrous"`, …). A species listed in two groups
  is an error; species in no group are collected rather than dropped.
- **`chemical_shrinkage(state, state₀)`** and the `"void"` phase. Passing
  `reference` to `volume_fractions` selects the sealed-curing convention: the
  fractions are referred to the initial volume, held fixed, and the deficit left
  by the reactions becomes an explicit gas-filled phase. Without it the fractions
  sum to less than one and the microstructure is wrong.
- **`missing_molar_volumes(state)`** — a species with no `V⁰` is silently absent
  from every volume computation in the package; this makes that visible.

Aqueous solutes have negative standard partial molar volumes, so an individual
fraction can be negative. Those contributions are kept, which is what makes
`volume_fractions` and `volume` agree exactly.

### Post-processing a kinetics solution

The ODE state carries only the kinetic species, and `n_full` is a buffer mutated
in place — after a run it holds the last accepted step and nothing else. There
was no way to recover the composition at an arbitrary time.

- **`reaction_extents(sol, kp)`** — the extent of each reaction, by quadrature on
  the dense output. It integrates on the union of the requested grid and the
  solver's own accepted steps, so a coarse output grid cannot step over the early
  transient.
- **`extent_residual`** — the mass-balance residual of that quadrature, so its
  accuracy is measurable rather than assumed.
- **`state_at(sol, kp, t)`** — the full `ChemicalState` at any instant, every
  species rebuilt from `n = n₀ + νᵀξ` so the result conserves exactly what the
  reactions conserve.
- **`degrees_of_hydration`** and **`mean_degree_of_hydration`**, replacing the
  `phase_alpha` closure copy-pasted into three shipped scripts.

### Parrot & Killoh, canonical formulation

The shipped `parrot_killoh` implements a smoothed variant whose parameters are
not those of the 1984 paper — its diffusion branch uses the same `K₃` and `N₃`
for all four clinker phases. It is unchanged, and now documented as one of two
variants.

- **`parrot_killoh_avrami`** with **`PK84_PARAMS_C3S/C2S/C3A/C4AF`** — the
  canonical Avrami / Jander / power-law form, `α̇ = min(α̇₁, α̇₂, α̇₃)`. With these
  parameters C₂S has no nucleation–growth stage and C₃S no diffusion-controlled
  stage, which is a sharp check on a transcription.
- **`waller`** with **`WALLER_PARAMS_FLY_ASH/SILICA_FUME/SLAG`** — supplementary
  cementitious materials do not follow Parrot & Killoh; `blended_cement_kinetics.jl`
  had to invent PK parameters for slag and metakaolin for want of this.
- **`blaine_factor`**, **`humidity_factor`**, **`powers_alpha_max`** — the three
  rate corrections, previously either absent or retyped inline in every script.

The Avrami branch vanishes at `α = 0`, so `α ≡ 0` solves the ODE and hydration
never starts. Parrot & Killoh's own discrete scheme escapes this by integrating
over the first time step; a continuous solver cannot, so the argument is floored
at `PK_AVRAMI_SEED`.

### Fixed

- **A time-dependent rate law broke every Rosenbrock solver.** The ODE right-hand
  side typed its rate vector from `eltype(u)` alone, but Rosenbrock methods
  (`Rodas5P`, `Rodas4`, `Rosenbrock23`, …) need a *time* gradient, which they take
  by calling the residual with a dual `t` and a plain `u`. Any rate depending on
  `t` then failed with "First call to automatic differentiation for time gradient
  failed". `parrot_killoh` ignores `t`, so nothing exposed it until `waller`. The
  type is now promoted with `typeof(t)`.
- **`reaction_extents` integrated from `times[1]`, not from the start of the run.**
  Asking for a single late instant silently returned a zero extent. The quadrature
  now always starts at `kp.tspan[1]`.

### Bibliography

Every entry of `docs/src/refs.bib` was checked against Crossref (six DOIs) or,
for the four entries that have no DOI, against the publisher or an authoritative
catalog record. Three defects were corrected:

- The DOI of `Lavergne2018` pointed at `10.1016/j.cemconres.2017.11.007`, which
  Crossref resolves to a **different** paper (Machner et al., *CCR* **105**, 1–17).
  Corrected to `10.1016/j.cemconres.2017.10.018`.
- `Lothenbach2015` held the Cemdata18 paper, published in **2019**, while
  `data/solid_solutions.toml` uses the tag `Lothenbach2015` for the genuinely
  different Lothenbach & Nonat (2015) paper and `Lothenbach2019` for Cemdata18.
  The entry is now keyed `Lothenbach2019`, and the real
  Lothenbach & Nonat (2015), *CCR* **78**, 57–70, is added.
- The author of `ParrotKilloh1984` is **Parrott**, with two t's. The citation key
  is unchanged, being referenced throughout the sources and documentation.

## v0.4.0 — Equilibrium actually coupled to kinetics, and dual numbers everywhere

### Breaking

- **`ChemicalState` gained a fourth type parameter**, `R <: Real`, carrying the
  number type of the dimensionless diagnostics: the signature is now
  `ChemicalState{C, S, Q, R}`. Code that spelled the type out in full has to be
  updated; code that used `ChemicalState` without parameters is unaffected.
- `pH`, `pOH`, `porosity` and `saturation` no longer return `Float64`
  unconditionally. They return whatever number type the composition carries, so
  that `∂pH/∂composition` and `∂porosity/∂(w/c)` are obtainable.

### Fixed — the equilibrium was never actually solved during a kinetics run

Four independent defects, each sufficient on its own to disable re-speciation.
They combined to make the whole coupling dead code, and the errors were hidden.

- **`equilibrate(state, ::EquilibriumSolver)` throws.** That is the call
  `build_kinetics_ode` made. `equilibrate` expects a SciML *algorithm* and wraps
  it into an `EquilibriumSolver`; handing it one that is already built wraps it
  twice, and the solve has no applicable method. The entry point for a prebuilt
  solver is `solve(esolver, state)`, which is what is called now.
- **The failure was swallowed by a bare `catch`**, so a run in which
  re-speciation never once succeeded looked exactly like a healthy one.
  Failures are now counted, the first is reported with its exception, and the
  total is reported at the end of `integrate`.
- **`KineticsSolver(; equilibrium_solver = …)` was never read.** `integrate`
  built its parameters from the `KineticsProblem` alone, so the documented way
  of passing a solver did nothing. It is honored now; a solver set on the
  problem still wins, and the conflict is reported.
- **The integrated element amounts `bₑ` were discarded.** The ODE integrates
  `dbₑ/dt = Aₑ νₑᵀ r`, but the re-speciation rebuilt its input from the
  persistent buffer instead, so the constraint of the equilibrium sub-problem
  (Leal et al. 2017, Eq. 54) was not the quantity being integrated. The previous
  speciation is now projected onto `bₑ` through `pinv(Aₑ)` before the solve.

### Fixed — three defects that made the `kinetic_species` API produce nonsense

Found while implementing the partitioned coupling, all on the 3-argument
`KineticsProblem(cs, state, tspan)` path:

- **The controlling mineral was mis-identified.** `_find_mineral_idx` takes the
  first crystalline reactant and, failing that, the first reactant of any kind.
  On a reaction generated from the nullspace — where the kinetic mineral may sit
  on either side — that is the solvent, so the ODE integrated **water** as the
  dissolving phase. A species the system declares kinetic now wins over any
  guess.
- **The stoichiometry had an arbitrary orientation**, so the kinetic mineral
  could carry `ν = +1` and the ODE grew the clinker: α(C₃S) came out at −4.2
  over seven days. Reactions are now normalized to `ν = −1` on their controlling
  mineral, which also delivers the `1/|νₖ|` scaling the `ChemicalSystem`
  docstring already promised.
- **The conservation basis was inconsistent.** `Aₑ` was cut from the canonical
  element matrix `CSM.A` while the solve is posed on `SM.A`, the matrix with
  respect to the primaries; and the equilibrium sub-system re-derived its own
  primaries instead of inheriting the parent's. Either mismatch makes every
  equilibrium solve fail — 89 failures out of 89 steps in one intermediate
  state. `Aₑ` now comes from the sub-system itself, which shares the parent's
  primaries.

With the three fixed, seven days of C₃S/C₂S hydration at `w/c = 0.4` gives
α(C₃S) = 0.277, portlandite and C-S-H in comparable amounts, an alkaline pore
solution at millimolar calcium, and `‖Aₑnₑ − bₑ‖∞ = 8.7e-8`.

### Changed — the coupling is now the partitioned formulation of Leal et al.

`respeciate!` solves `nₑ = φ(bₑ)` over the **equilibrium partition only**, on a
sub-system built for the purpose, with the integrated element amounts `bₑ`
handed to the solver as the constraint rather than derived from a composition.
Running the minimization over the whole system, as the previous code did, would
equilibrate the kinetic minerals instantaneously.

**Validated against Reaktoro.** Calcite dissolving at a constant rate makes the
kinetic half analytic, so every difference of rate-law convention drops out and
only `nₑ(t) = φ(bₑ(t))` is under test. The two codes agree to 0.3 % at `t = 0`
and 4.3 % at 3600 s, the largest deviation always on `CO₂(aq)` at `2.3e-8` mol.
Assertions in `test/coupling_reference.jl`, oracle in
`test/reference/reaktoro_coupling.py`.

### Changed — operator splitting instead of a solve inside the right-hand side

The equilibrium sub-problem is no longer solved inside the ODE right-hand side.
It is solved once per accepted step, by the new `respeciate!`, wired as a
`DiscreteCallback`. The right-hand side is evaluated many times per step and
differentiated to build the Jacobian; solving an optimization problem in there
made the cost unpredictable and the Jacobian inconsistent with the model being
integrated. With the solve outside, residual and Jacobian see the same
speciation.

### Fixed — the conservation matrix was rebuilt by finite differences

`OptimaSolverExt` did not pass `A` and `b` through the problem parameters, so
`OptimaSolver` fell back to reconstructing the constraint Jacobian by finite
differences with a `1e-7` step. The element balance was capped at ~2e-6 mol and
did **not** respond to the requested tolerance. `A` is known exactly; it is
handed over now. Measured on calcite/water: feasibility `2.17e-06 → 4.02e-14`.

This also settles the choice of default back-end. Once the matrix is exact,
OptimaSolver reaches machine-level feasibility *and* is 3 to 26 times faster
than Ipopt on the cement equilibria, so it remains the high-priority default.

### New — differentiating *through* the equilibrium solve

`ForwardDiff` now crosses an `equilibrate`. No solver is asked to iterate on
dual numbers — Ipopt is a C library and never could. The equilibrium is solved
once at the primal values and the sensitivities come from the optimality
conditions of [Leal2017](@cite), the problem Reaktoro solves: one saddle-point
factorization serves every partial derivative, and the answer is exact rather
than a finite difference.

The complementarity block is what partitions the species, and skipping it does
not degrade the answer gently. On calcite + CO₂ in water with a gas phase
declared, the unreduced system puts the whole response into the **absent** gas
species (`n = 5e-11` mol) — a sensitivity that satisfies the element balance to
`4e-16` and means nothing. No back-end returns the stability multipliers `z`, so
the active set is recovered internally: a species negligible on the scale of the
system that nonetheless takes a leading share of the response is pinned and the
system re-solved.

Verified against the package's own finite differences (`9e-5`, the truncation
error) and against **Reaktoro 2.13 reading the same Cemdata18 file**, over the
same eleven species and under the same activity model: the two agree to
`4.4e-4` relative or better on every species, below Reaktoro's own `7.2e-4`
spread across finite-difference step sizes. The equilibrium amounts agree to
the same order. The absent gas species gets exactly zero from the active-set
treatment, against `2e-9` by finite differences.

A cross-code comparison has three knobs — database, species list, activity
model — and each is worth tens of percent. Leaving Reaktoro on HKF against this
package's default `DiluteSolutionModel()` moves `∂Ca²⁺/∂(CO₂)` from `+0.1520`
to `+0.2179`; dropping the aqueous calcium complexes makes `∂Ca²⁺/∂(CO₂)` and
`∂calcite/∂(CO₂)` mirror each other, and restoring them breaks that mirror in
*both* codes alike. The element balance closes to `2e-16` throughout.

### Fixed — trace species now agree with Reaktoro

On calcite + CO₂, every species except `CaOH⁺` agrees with Reaktoro to 5 % or
better, and `pKw` comes out at 13.979 instead of 13.40. Before, `CaOH⁺` was
20.8× high, `OH⁻` 3.19× and `H⁺` 1.24×.

The cause was in the back-end's convergence test, which excluded variables
judged "at their bound" using a threshold scaled by the *largest* amount in the
problem. With a solvent at 55 mol that threshold was `5.5e-7`, so every trace
ion below it was declared to sit on a bound of `1e-16` and its stationarity was
never enforced. `OptimaSolver` 0.2.5 judges it against the variable's own bound.

`CaOH⁺` remains 2.5× high at `4e-9` mol. Ipopt lands on the same value (×2.46),
so it is attributable to neither back-end, and at that amount it is chemically
inconsequential. The assertion is `@test_broken`.

### Changed — the analytic gradient is handed to the back-end

`∂G/∂nᵢ = μᵢ` exactly, the remaining term vanishing by Gibbs–Duhem. The
extension now supplies it instead of leaving the back-end to differentiate
`dot(n, μ(n))`, where that cancellation is only exact in theory: evaluated, the
solvent term alone is `55 × 0.018 ≈ 1` against a `μ` of order 10. This does not
resolve the defect above, but a Gibbs energy is stationary exactly where its
gradient is, and steering on a gradient wrong by ten percent is not defensible
whatever else is true.

### Added — documentation

- **Coupling kinetics and equilibrium** — the partitioned formulation derived,
  including why the state is `(bₑ, nₖ)` rather than `(nₑ, nₖ)`, and what to
  check when a run finishes.
- **The hydrating paste, end to end** — a worked example computing the hydrate
  assemblage instead of imposing it, with measured output.
- **Validation against Reaktoro** — what agrees, what does not, and the three
  knobs a cross-code comparison has to match.

### Fixed — the solver's return code was ignored

Neither extension looked at the optimizer's `retcode`; a non-converged iterate
was written into the state as though it were the equilibrium. It is now checked
and warned about. The warning, rather than an error, is deliberate: the flag is
unreliable in both directions on these problems — pure water comes back flagged
`MaxIters` while giving `[H⁺]/[OH⁻] = 1.000003`. Set
`ChemistryLab.STRICT_CONVERGENCE[] = true` to raise instead.

### Fixed — dual numbers could not enter a `ChemicalState`

The constructor took its element type from the *temperature*
(`Q = typeof(T_q)`), which silently cast the composition to `Float64`: nothing
downstream of a `ChemicalState` was differentiable. The type is now the
promotion of `T`, `P` **and the amounts**. The dimensionless diagnostics follow
(see Breaking above); only their *display* strips the dual part.


## v0.3.1 — Maintenance

- `[compat]` upper bounds raised: `OrdinaryDiffEq` to `"6, 7"`, `TimerOutputs`
  to `"0.5, 1"`.
- CI badge restored; Runic badge.
- Installation instructions updated for registration in Julia's General
  registry; documented the optional optimization backend required to solve
  equilibria (`Optimization`+`OptimizationIpopt`, or `OptimaSolver`).
- Confirmed each GitHub Release keeps archiving automatically to Zenodo's
  existing concept DOI `10.5281/zenodo.17756074` via the native
  GitHub↔Zenodo integration (no workflow or token needed).
- MPCM-Registry deprecated in favor of the General registry.

## v0.3.0 — Chemical kinetics module

### New features

**`src/kinetics/` — mineral dissolution/precipitation kinetics**

- `KineticReaction` — couples a reaction to a rate law; supports
  `transition_state`, `first_order_rate`, and the empirical
  Parrot–Killoh (1984) model for cement clinker hydration
- `KineticsProblem` / `KineticsSolver` — ODE problem formulation
  following the SciML `(u, p, t)` convention; integrated via
  `KineticsOrdinaryDiffEqExt` (weakdep, activated by `using OrdinaryDiffEq`)
- `SemiAdiabaticCalorimeter` / `IsothermalCalorimeter` — coupled
  heat-balance ODE; `temperature_profile`, `heat_flow`,
  `cumulative_heat` post-processing helpers
- `StateView{T,I}` — O(1) named access to species data vectors in the
  ODE hot path (no per-step allocation)
- `RateMechanism`, `RateModelCatalyst`, `BETSurfaceArea`,
  `FixedSurfaceArea` — building blocks for custom rate closures
- `parrot_killoh(params, mineral_name)` factory with built-in
  Schindler & Follliard (2005) Arrhenius correction and default
  parameters for C₃S, C₂S, C₃A, C₄AF
- ForwardDiff-compatible throughout; `KineticFunc` and `transition_state`
  closures accept `Dual` numbers

### Infrastructure
- Relicensed to LGPL-2.1-or-later
- Registered in MPCM-Registry (OptimaSolver resolved via registry)
- GitHub Actions workflows: CI, Documentation, Register, CompatHelper, Format, TagBot
- Multi-version documentation deployment (`docs/deploy_docs.jl`)

## v0.2.3 — Activity models & solid solutions

- Extended Debye–Hückel and Davies activity models
- Redlich–Kister solid solution phases
- `SolidSolutionPhase`, `build_solid_solutions` from TOML
- `with_class` for end-member requalification

## v0.2.0 — Equilibrium solver integration

- `EquilibriumProblem` / `EquilibriumSolver` / `equilibrate` API
- `OptimizationIpoptExt` and `OptimaSolverExt` extensions
- Implicit-differentiation sensitivity via OptimaSolver
- HKF aqueous solute thermodynamic model (`NumericFunc`)
- ThermoFun JSON and PHREEQC `.dat` database parsers
