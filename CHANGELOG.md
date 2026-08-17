# Changelog

## v0.7.1 — the element balance of a re-speciation, measured and enforced

Bug fixes only. No API was removed or renamed, and no documented behavior
changed, so a downstream `[compat] ChemistryLab = "0.7"` accepts this release
unchanged.

### The reported element-balance residual could not see a violated element

`_respeciate_solve!` scaled the whole residual `|Aₑnₑ − bₑ|` by
`maximum(abs, bₑ)` — in a cement paste, by the water budget. Water is 34 mol
against 0.27 mol of sulfate, so a violation of 0.465 mol of sulfate, i.e. 174 %
of the sulfate present, was reported as `1.4e-2` and read as a converged solve.

The residual is now taken row by row and scaled by each element's own budget,
or by the matter flowing through that row when the budget is zero — as it is for
the charge row, where dividing by `bᵢ` alone turned a rounding error into a
reported residual of 2·10³.

This is a diagnostic change, but not a cosmetic one: every convergence study
run against the old number was measuring the wrong quantity.

### A non-converged speciation was handed on as the next warm start

`eq_warm` was set after any solve that did not throw, non-converged ones
included, so a single bad point seeded every step after it and the error was
locked in for the rest of the run. A speciation is now handed on only when its
per-element residual is within `EQ_RESIDUAL_TOL`; otherwise the step falls back
to the cold stoichiometric reconstruction.

### Re-speciation could start outside its own feasible set

The Gibbs minimization is posed with `Aₑn = bₑ` as a hard equality, and the
warm start is the equilibrium of the *previous* `bₑ`. Once an element has been
spent — the sulfate of an ordinary Portland cement, after the gypsum is gone —
that guess demands more of it than now exists, and the interior-point solve
starts infeasible.

Two guards now run before the solve, both of which leave a feasible guess
untouched:

- `_budget_clip!` caps each species at what the totals can supply,
  `nⱼ ≤ minᵢ bᵢ/Aᵢⱼ`, over the rows it consumes. Rows with a negative total —
  H⁺ in a hydrating cement, met by the hydroxides — bound nothing.
- `_restore_feasibility!` then lands the guess in `{Aₑn = bₑ, n ≥ 0}` by
  alternating projection, ending on the positivity clamp so no small negative
  amount survives for the barrier to reject.

On an isolated ordinary-Portland-cement equilibrium — alite, belite, aluminate
and gypsum dissolved, with ettringite, monosulphate, katoite, portlandite and
C-S-H free to precipitate — this takes the achieved element balance from 174 %
to `2e-15`, machine precision, and returns the textbook assemblage: AFm 0.2323,
ettringite 0.0116, C-S-H 1.8556, portlandite 2.3021, at pH 12.58. Both budgets
close on the last digit — sulfate `0.2323 + 3x0.0116 = 0.2672` and aluminum
`2x(0.2323 + 0.0116) = 0.4879` — against the 0.2672 and 0.4879 available.

### The cold start was supersaturated in every phase at once

The cold-start guess added the stoichiometric reconstruction `Ne^T xi`, which
places every dissolved element in solution with no hydrates at all. For an
aqueous-only system that is harmless; for a cement it is close to the worst
possible start, supersaturated in every phase simultaneously, with an H+ entry
so negative that it is clamped to the floor and loses the acidity the hydroxides
must balance. The guess is now the composition the specimen was cast with, which
`_restore_feasibility!` then carries onto the current `be`.

### What this unblocks

A full ordinary Portland cement now runs end to end: alite, belite, aluminate,
ferrite and gypsum, over 28 days, 202 accepted steps, `retcode = Success`. The
pore solution holds at pH 12.58, and the aluminate sequence comes out of the
thermodynamics with no sequencing rule written anywhere -- ettringite forms
early, peaks at 6 hours and converts to monosulphate once the sulfate is spent,
the AFm settling at 0.26719 mol, exactly the sulfate budget at one SO4 per
formula.

Note that `p.n_full` is a scratch buffer rewritten at every right-hand-side
evaluation, Jacobian differences and rejected steps included; it is not the
accepted composition at `t_end`. Recover a speciation from a solution by taking
`be` from the ODE state and re-equilibrating, as `state_at` already documents.

## v0.7.0 — the real porosity of a setting binder, and a warm-started coupling

### Breaking

- Being a minor bump below 1.0, a downstream `[compat] ChemistryLab = "0.6"`
  will not accept 0.7 and must be widened. No API was removed and no documented
  behavior changed; the additions below are new methods.

### The porosity of a cement was not available

`porosity(state)` returns `(V_liquid + V_gas) / V_total`, and both ends of that
ratio are wrong for a hydrating binder: the denominator is the *current* volume,
which shrinks as the reactions proceed, while a sealed specimen keeps the volume
it was cast with; and the numerator has no gas term, so the empty porosity left
by the chemical shrinkage is structurally invisible. The errors compound — on a
w/c = 0.5 paste at 28 days the method returns 0.327 where the porosity referred
to the specimen is 0.375, the volume having shrunk 7.2 %.

The method is unchanged, and correct for a fixed-volume aqueous system, but now
carries that warning. Two new methods give the right calculation:

- **`porosity(state, reference)`** → `(; liquid, void, total)`, referred to the
  fresh material and counting the Le Chatelier contraction as empty porosity.
- **`saturation(state, reference)`** — a sealed paste desaturates as it hydrates
  though no water ever leaves it.

`porosity(state, ref).void` is exactly the `"void"` entry of `volume_fractions`
called with the same reference, so a transport or micromechanical model fed by
either sees the same thing.

### Fixed — the equilibrium sub-solve was restarted cold at every step

`respeciate!` rebuilt its starting guess from the reaction extents each time,
placing every dissolved element in solution with zero hydrates — a wildly
supersaturated composition. Harmless for an aqueous-only system, nearly the worst
possible start for a cement. It now warm-starts from the speciation left by the
previous accepted step, which is an equilibrium for a nearby `bₑ`.

## v0.6.0 — the coupled path made trustworthy

The v0.5.0 release opened the door to coupled dissolution/precipitation modeling.
Building a cement model on it exposed five defects, four of which were silent.

### Breaking

- **`EquilibriumSolver` gained a `model` field**, so its signature is now
  `EquilibriumSolver{F, S, V, M}`. Code constructing the raw struct positionally
  must be updated; the documented `(cs, model, solver)` constructor is unchanged.
  A new accessor `activity_model(::EquilibriumSolver)` returns it.
- Being a minor bump below 1.0, a downstream `[compat] ChemistryLab = "0.5"`
  will not accept 0.6 and must be widened.

### Fixed — the activity model was silently discarded in a coupled run

A coupled run must rebuild the equilibrium solver for the equilibrium
*sub-system*, because the compiled potential does not carry over. Having no
record of the model it was built from, it rebuilt with `kp.activity_model` —
`DiluteSolutionModel()` by default. Asking for `HKFActivityModel` on a cement
pore solution at I ≈ 0.5 mol/kg therefore gave an infinitely dilute solve, with
no warning. The solver now remembers its model, the rebuild uses it, and a
disagreement with the problem's own model is reported.

### Fixed — the equilibrium sub-solve started on its own bound

`respeciate!` floored its starting guess at `p.ϵ` (1e-30), which
`EquilibriumProblem` raises to exactly the `1e-16` lower bound. An interior-point
method started on its bound stalls: the package's own Reaktoro reference case
reported six non-converged solves out of eight steps for this reason alone. The
guess is now floored strictly inside the box.

**Do not "fix" this class of stall by loosening the optimizer tolerance.**
Measured against Reaktoro 2.13 on that same case, the worst species error grows
from 4.3 % at the default `tol = 1e-10` to 13 % at `1e-9`, 38 % at `1e-8` and
252 % at `1e-7`. A green retcode bought that way costs a factor of sixty in
accuracy.

### Fixed — a bad speciation was invisible

A non-success retcode was a `@warn` at `maxlog = 1` whose result was used anyway,
and it never reached the failure count reported by `integrate`, which only saw
solves that *threw*. Over thousands of steps that is one warning for any number
of bad speciations.

- `NONCONVERGED` counts them, and `integrate` reports the total.
- More usefully, `integrate` also reports the worst `|Aₑn − bₑ|∞` over the run.
  That is the criterion with physical meaning here: a solve can stop short of the
  optimizer's tolerance and still satisfy the element balance to machine
  precision, which is exactly what the remaining stalls do (9e-12 on the
  reference case).

### Fixed — `integrate(kp; reltol = …)` threw a `MethodError`

The documented shortcut forwarded its keywords to a concrete method that accepted
none. Precedence is now explicit: call site beats the solver's own settings,
which beat the defaults.

### Fixed — two solid solutions were dead entries

`data/solid_solutions.toml` declared AFm on `Ms`/`Mc` and Hydrotalcite on
`Ht_OH`/`Ht_CO3`. None of those four symbols exists in either shipped database,
so `build_solid_solutions` skipped both with a warning. AFm being the only
Redlich-Kister entry, the non-ideal mixing path had no live case at all. They are
now `monosulphate12`/`monocarbonate` and, for hydrotalcite, the OH/CO₃ couple at
matching Mg:Al = 2 (`hydrotalcite`/`Mg2AlC0.5OH`). All five phases build.

### Documentation

- `saturation_ratio` stated `ln Ω = Σνᵢlnaᵢ + ln K`, contradicting both the next
  line and the code. The sign is corrected.
- `examples/cement_carbonation.md` claimed `cemdata18-thermofun.json` does not
  contain calcite. It does; the merged database is needed for the phase volumes.
- `man/kinetics.md` used `cs.dict_reactions` as though it were a catalog of the
  database, when it holds only the declared kinetic reactions.

## v0.5.0 — the bridge from a kinetics run to a microstructure

### Breaking

- **The ODE state vector gained the extents of reaction.** Its layout is now
  `[bₑ, nₖ, ξ, (calorimeter)]`, where `ξ` integrates `dξ/dt = r`, one entry per
  kinetic reaction. Code indexing `sol.u` positionally for anything other than
  `bₑ` or `nₖ` must be updated; the calorimeter's slot is addressed from the end
  of the vector and is unaffected, as are `temperature_profile`,
  `cumulative_heat` and `heat_flow`.
- Being a minor bump below 1.0, a downstream `[compat] ChemistryLab = "0.4"`
  will not accept 0.5 and must be widened regardless.

### A rate law may now depend on any species

This is what the state change buys, and it is the point of the release.

Previously only the kinetic species evolved inside the ODE residual: every other
amount was pinned to its initial value for the whole run, because the buffer
holding them was refreshed only by `respeciate!`, which runs only when an
equilibrium solver is present. A rate closure gating on a consumed reactant
therefore never saw it move, and the failure was silent — the kinetic mass
balance stayed exact and the solver reported success while the gated species went
negative. In the cement model that motivated this release, 0.44 mol of ettringite
formed out of 0.27 mol of gypsum.

The residual now reconstructs every species from `n = n(0) + νᵀ ξ` before
evaluating the rates, so reaction sequencing — sulfate available for ettringite,
portlandite available for a pozzolanic reaction, product inhibition — is written
directly. With an equilibrium solver the equilibrium partition is still owned by
the equilibrium solve and refreshed once per accepted step.

As a side effect, [`reaction_extents`](@ref) and [`state_at`](@ref) read the
extents from the solution instead of re-integrating the rates by quadrature: both
are now exact to the solver's tolerance, and the `nsub` keyword is gone.

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

`n_full` is a buffer mutated in place — after a run it holds the last accepted
step and nothing else — so there was no way to recover the composition at an
arbitrary time.

- **`reaction_extents(sol, kp)`** — the extent of each reaction, read from the
  ODE state and therefore exact at any instant.
- **`extent_residual`** — the drift between `nₖ` and `νₖᵀ ξ`, which are redundant
  by construction, so the integrator's own consistency is measurable.
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
- **Non-kinetic amounts were frozen inside the residual** — see the section above,
  which is the substance of this release.

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
