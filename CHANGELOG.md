# Changelog

## v0.14.0 — an equilibrium that comes with a proof, and a kinetic step that is one problem

### Breaking changes

- `equilibrate(state)` now solves through **every** loaded back end and returns
  the answer `optimality_certificate` proves optimal. It is slower and, on some
  systems, gives a different — correct — answer. `equilibrate(state; certify =
  false)` restores the single-back-end behavior.
- `OptimaSolver = "0.5"` is required. The parameter block the new constraints and
  the kinetic step ride on arrived there; against 0.5's predecessor the extension
  fails with a `MethodError` on `DualNewtonProblem`.
- The `equilibrate` keyword `model` is now named explicitly in the one-argument
  form, and `certify` and `constraint` are new keywords. Code passing positional
  arguments is unaffected.

### An interior-point answer whose charge balance was wrong in the second digit

On 1 mmol of calcite dissolving in 1 kg of water, the interior-point path returned
`‖An − b‖∞ = 3×10⁻⁶`, which looks like round-off. It is not: the rows of a
conservation matrix do not share a scale, and against its own budget that residual
is **3 % on the charge row** and 0.3 % on the carbonate, while the water row at
55.5 mol sits at `1.8×10⁻⁸`. The pH came out 9.5210 against a correct 9.5274.

The measure is fixed in `OptimaSolver` 0.5.0. The answer is fixed here, and the
iteration is not the place: traced over twenty-three iterations, `α` is pinned at
its ceiling of 0.15 and `‖dn‖` decays at `1 − α` with the residual frozen — the
fraction-to-boundary rule lets through 15 % of a correction the next iteration
re-poses. Thirty thousand iterations change nothing, nor does a tolerance of
`10⁻⁶`.

What does is that the problem is convex, so its KKT conditions are sufficient and
`optimality_certificate` **decides** rather than ranks. `equilibrate` therefore
solves by every route and keeps a proved answer:

| case | interior point | dual Newton | certified |
|:--|--:|--:|--:|
| calcite in water, 10–40 °C | 3.0e-2 | 3.0e-12 | proved |
| calcite + 1 mmol CO₂ | 1.0e-3 | 6.3e-13 | proved |
| calcite + 50 mmol CO₂ | 1.3e-16 | 8.5e-14 | proved |
| pure water | 1.3e-16 | 1.3e-16 | proved |
| CEM I, w/c 0.45 and 0.60 | 1.5e-14 | 2.0e-14 | proved |
| CEM I, w/c 0.30 | 7.2e-16 | 3.3e-10, **not** proved | proved |

Ten of ten, against seven through the interior point alone and nine through the
dual Newton alone.

One bug found in the process, and it mattered: `solve_certified` recomputed the
component totals from **each** starting point. A start that violates the balance
therefore posed a different problem, and the dual solve certified the answer to
the shifted one — two "certified" compositions 0.2 % apart on dissolved calcium,
which a convex problem with one minimum cannot have. `b` is now fixed once, from
the state as given.

### Equilibrium under something other than fixed (T, P)

`Adiabatic`, `FixedEnthalpy`, `FixedVolume`, `SealedVolume`, `FixedpH` and
`FixedActivity`, passed as `equilibrate(state; constraint = ...)`.

The prescribed quantity is an unknown of the solver's own square system, not an
outer loop: an adiabatic solve costs one equation, not a solve per trial
temperature. Two vehicles, and they are not interchangeable — a prescribed
**property** adds a parameter and its residual, a prescribed **chemical
potential** adds a *column* to the conservation matrix, an unknown amount of a
titrant the system may draw on, and that amount comes back as part of the answer.

Validated against a published number: adiabatic `H⁺ + OH⁻ → H₂O` gives
**−55.85 kJ/mol** at every amount tested, against the accepted −55.8, with
nothing fitted to it. `ΔT = 6.62 K` for 0.5 mol in a kilogram of water, against
27.9 kJ over 4.18 kJ/K. Enthalpy conserved to `10⁻¹²` relative, and solving at a
*fixed* temperature equal to the one found returns the same composition to
`10⁻¹²`.

A volume constraint on a condensed system is **refused**, with the lever it
measured named in the error. The molar volumes of water and of the minerals in
these databases are exactly pressure-independent, so the relative lever
`(∂V/∂P)·P/V` is about `10⁻⁶` — some 9 600 bar to change the volume by one
percent. That is the physics: the volume of an incompressible condensed system is
fixed by its composition. Declare a gas phase, or use `porosity(state, reference)`
for a sealed specimen at fixed pressure.

### A kinetic step as Leal poses it: one problem, fully implicit

`KineticStepSolver` and `kinetic_step`. The reaction extents are unknowns of the
same Gibbs minimization as the amounts and the element potentials, and the rate is
evaluated at the **end-of-step** composition. No frozen speciation in a right-hand
side, and no lag between the kinetic and the equilibrium species. The reactivity
constraint being linear, it joins the conservation block, so the algebraic cost of
kinetics is the number of **reactions**, not of species.

Measured against an analytic answer: at a constant 1 µmol/s, `Δξ = Δt·r` to
`10⁻¹⁰` relative and the calcite left is `n₀ − kΔt` to the last bit, with the
element balance at `10⁻¹⁴`.

Measured on the feedback path nothing exercised before: with `r = k(1 − Ω)` at
`k = 10⁻⁵ mol/s`, an explicit step of `10⁶ s` would dissolve 10 mol — a thousand
times the calcite present. The implicit step lands at `Ω = 0.999989`, approaching
saturation from below and never crossing it, at any step size.

A **solid solution** takes part: C₃S dissolving into a four-member CSHQ solution
comes out certified at a stationarity of `2×10⁻¹³`, with all four end-members
present and ten times the step giving ten times the C-S-H. Two settings are
decided by the certificate rather than guessed, because neither answer works on
both kinds of problem — `warm_start` for the admission of a mixing phase from a
cold start, `pin_minerals` for whether the kinetic minerals are held in the active
set. Both are documented with the measurements that decide them.

The certificate of a kinetic step is taken on the **augmented** problem. A mineral
held back by a rate law is supersaturated by construction, so tested against the
unconstrained equilibrium it reports a residual of 7 `RT` and a correct answer is
called wrong.

`coupling = :species`, which would impose the assemblage instead of proposing it,
**raises**. It was implemented and does not converge: pinning a pure phase by a
linear row leaves its stationarity row in the system as well, made trivially
satisfiable by that row's own multiplier, and the Jacobian degenerates — the
extents came out 4.3 times short. For a model whose assemblage IS the
stoichiometry, the ODE route already does exactly that, and the tutorial now opens
with a table saying which route answers which question.

### Choosing the step length, and a case the stiff ODE route gets wrong

`kinetic_step_adaptive` chooses `Δt` from a Richardson estimate on the extents —
one step against two half-steps, which for a first-order method is the error of
the coarse one. Reaktoro has no equivalent: the step is the caller's and nothing
reports what taking it cost.

Measured on calcite under `r = k(1 − Ω)`, `k = 10⁻⁴ mol/s`, over `10⁵ s`:

| route | steps | result |
|:--|--:|:--|
| `Rodas5P`, and OrdinaryDiffEq's default polyalgorithm | 6 | **extent −457 mol**, `retcode = Success` |
| `Tsit5`, explicit | 85 626 | correct, 519 s |
| one `kinetic_step` of `10⁵ s` | 1 | wrong, and its certificate says so |
| `kinetic_step_adaptive` | **7** | correct to eight digits |

The cause is not a missing Jacobian term, which is the natural guess. The ODE
route's residual reads the speciation **frozen** at the last accepted step, so
`∂(du)/∂bₑ = 0` is exact for the system being integrated. What goes wrong is that
the frozen speciation makes the right-hand side inconsistent with the state within
a step: an implicit method steps past the point where `Ω` crosses one, the rate
changes sign, and the run enters a branch it never leaves. Three measurements
settle it — bounding `dtmax` to `10³ s`, 103 steps against 6, returns the
identical wrong value to seven digits; removing the re-speciation returns a sane
answer; and routing the re-speciation through the certified route moves −457 to
−383. So the fix is not to freeze, which is what the implicit step does.

The `OrdinaryDiffEq` extension now warns when a trajectory ends on amounts no
chemistry can produce, naming the species and the budget it exceeded. That is
what the ODE route can offer; a rate law that reads the solution belongs on the
implicit one.

### Re-speciation escalates to the certified route when it fails the balance

`respeciate!` went through the interior point alone, and where that is wrong it is
wrong by percent. It now escalates: when the plain solve's absolute element-balance
residual exceeds the retry tolerance, the partition is re-solved through every
back end and the answer the KKT certificate proves optimal is kept. The cost of a
normal run is unchanged, because the escalation fires only on the solves that need
it.

### `ode_solver = :auto`

Hands the choice of integrator to `OrdinaryDiffEq`'s default polyalgorithm, which
detects stiffness at run time. It is offered rather than made the default: on
these problems it lands on the same stiff method as `Rodas5P` and returns the same
answers, so switching would change nothing while removing a reproducible choice.

A step is accepted on the estimate **and** on the certificate. Richardson's
difference is blind to an error the two resolutions share: on calcite over
`10⁵ s` the coarse step and both half-steps each dissolve the entire mineral,
their extents agree to `5×10⁻¹¹`, and the estimator grades that step excellent.
The certificate is not fooled, since such a composition violates
`Δξ − Δt·M·r(n) = 0` by the whole extent. This also closed a hole in
`kkt_certificate`, which checked stationarity, the linear rows and phase
stability but not the **nonlinear** parameter residual, so a step violating its
own rate equation could be proved optimal (`OptimaSolver` 0.5.0).

One further design point, because the obvious choice fails. The tolerance
is relative to the amount each reaction acts on, **not** to the extent. Scaling by
`Δξ` makes the tolerance vanish with the step while the equilibrium solve's noise
does not, so the measured error grows as the step shrinks: the controller halved
to the floor without advancing, and the answer got worse as the tolerance was
tightened — 2.5e-5 at `reltol = 1e-3` against 9.2e-5 at `1e-5`.

### `coupling = :species`: prescribed products inside the strong coupling

The third combination, which neither Reaktoro nor the ODE route offers. One
constraint per kinetic species, `nᵢ = nᵢ(0) + Σⱼ νᵢⱼ Δξⱼ`, so a solid-to-solid
scheme in the form of [Lavergne2018] imposes its products while the aqueous phase
still minimizes `G` under element conservation. Measured on two C₃A pathways, the
extents come out at `Δt·r` exactly and every species follows the stoichiometry to
eight or nine digits, ettringite and monosulfoaluminate included, where
`:reactions` leaves the AFm at zero because the minimization is free to rearrange
the solids.

It was refused in an earlier draft of this release, and finding out why turned up
a genuine defect in `OptimaSolver`: a pinning row for a product that starts absent
has one positive entry and a zero budget, exactly the shape `degenerate_components`
reads as "this component is absent from the system". It pinned that row's
multiplier and declared the species dead — stationarity residual 458, extents 4.3
times short. `conservation_rows` now restricts that criterion to rows that
conserve something.

The pinned species are **eliminated**, not constrained, and that is what makes it
work. Holding them by a linear row — the obvious implementation — stalls: the
species' stationarity row stays in the system, satisfied by a multiplier that must
reach the mineral's own chemical potential, of order 10²–10³ in `RT` units, and
the Newton sits at a fixed point its line search cannot leave. Measured, an
element balance of `6.1×10⁻⁷` mol whatever `maxit`, `tol` or the number of
active-set updates. Their amounts being an explicit affine function of the
extents, they are removed from the system instead, their element content is
subtracted from the budget, and a plain equilibrium is solved over what remains
with a Newton on the `nr` extents around it: the balance goes to `2×10⁻¹⁴` and the
step is certified.

### A single step far beyond the relaxation time can find the other root

For a rate law that vanishes at equilibrium the implicit step has two solutions,
and the second is the composition with the mineral wholly dissolved: it satisfies
the element balance and the reactivity row while violating `Δξ − Δt·M·r(n) = 0`
by the entire extent. The certificate refuses it, and `kinetic_step` reports the
step uncertified rather than passing it off.

Which root the Newton finds is a property of the build, not of the chemistry.
Measured on calcite under `r = k(1 − Ω)` with `k = 10⁻⁵ mol/s`: steps of `10⁴ s`
and `10⁶ s` converge to the right root on Julia 1.12 and to the other one on
1.13.0-rc4 — which is how this was found, a suite green on one and five failures
on the other. `kinetic_step_adaptive` is build-independent, because it refuses an
uncertified step and halves until one certifies; it reaches the equilibrium
values to eight digits in seven steps on either.

The test now asserts the contract rather than the outcome of a Newton on a
particular build: a certified step never crosses saturation, an uncertified one is
flagged and its rate-equation residual is of order one, steps well inside the
relaxation time certify on any build, and the adaptive march reaches the right
answer. Asserting the physics unconditionally on a step of `10⁶ s` was asserting
which root a given LLVM finds.

### `pin_minerals = :auto` ranks the two formulations by margin

It used to return the first one whose certificate passed. That makes the choice a
**threshold** decision, and a threshold decision between two answers that straddle
the tolerance is settled by the last bits: the same commit then behaves differently
on a different machine. Found exactly that way — the suite green locally and five
failures in CI on the identical commit, in the test asserting that a saturating
rate law never crosses equilibrium.

Both formulations are now computed and ranked by a single margin —
stationarity, element balance, worst supersaturation among absent phases, and the
residual of the rate equation. Ranked, the two are not close: on the case that
exposed it the free branch has run to equilibrium and violates its own rate
equation by the whole extent while the pinned one satisfies it to `10⁻¹³`, fifteen
orders of magnitude apart. The adaptive march resolves `:auto` the same way.

### ForwardDiff through the certified route

A composition carrying dual numbers takes the implicit-function route: the
certified solve runs in real arithmetic and the derivative is attached at the
answer. This had to be built rather than inherited — the certified path converts
the component totals to `Float64` to fix them once, which silently dropped every
partial. `∂pH/∂n(CO₂)` now agrees with a central difference to `4.3e-9` relative,
and the primal value is the certified one.

Pushing duals through the search would be wrong in any case: an active set has a
discrete component, so the map `b ↦ n*(b)` is smooth only piecewise and the
derivative belongs at the solution with the active set frozen.

### Correctness of published values

- **The "maleic acid" titration is malonate.** The SLOP98 entries are
  `MALONIC-ACID,AQ` / `H-MALONATE,AQ` / `MALONATE,AQ`, the three-carbon diacid.
  Their `ΔₐG⁰` give pKa **2.851 / 5.696** against a tabulated malonic 2.83 / 5.69,
  and nowhere near maleic's 1.92 / 6.27. The database was right and the name was
  wrong. The script and the page are renamed, and both pKa are now derived rather
  than hard-coded — the old figures were plotted as dashed lines that crossed the
  curve at no half-equivalence, so the defect was visible in the published figure.
- **`parrot_killoh` is deprecated and no longer attributed to Parrott & Killoh.**
  Its nucleation term carries no Avrami logarithm, `K₃` sits where the canonical
  form has `k₂`, `N₁ = 3.3` is the canonical `n₃`, and `k₃ = 1.1` has no
  counterpart: two different models, not two parameterizations. With `PK_PARAMS_*`
  the diffusion branch takes over at α ≈ 0.003 (C₂S), 0.013 (C₃S) and 0.057
  (C₃A), and those three then land on **α(7 d) = 0.2386 whatever their `K₁`**,
  while C₄AF is limited by its own nucleation branch at 0.193. A CEM I at
  w/c = 0.40 is reported near 0.61. All demos and doc pages now use
  `parrot_killoh_avrami` with `PK84_PARAMS_*`: the CEM I paste moves from
  ᾱ(7 d) = 0.234 to **0.628**, ΔT from 2.0 to **14.2 °C**, and the heat released
  from 115 to **308 kJ/kg**. The slag and metakaolin laws move to `waller`.
- **The calorimeter's `Cp` was double-counted.** The denominator is
  `Cp + Σᵢ nᵢ Cp°ᵢ(T)` with the second term recomputed at every step, so adding
  the sample to `Cp` counts it twice and understates ΔT by a factor of about 1.75.
  Three scripts and one doc page corrected; the docstring, which contradicted
  itself, now carries the warning.
- **A documented temperature sweep ran on an empty state.** `ChemicalState(cs)`
  holds no matter, so every element balance was zero and all 21 points returned
  the same trivial answer — flat lines, with prose describing a curve that was not
  there. Repaired, it also shows the note was backwards: `Kₛₚ` is retrograde
  (`10^-8.411` at 10 °C to `10^-8.517` at 30 °C, matching the database to 0.003
  log units, and `−8.48` at 25 °C is the accepted value), and yet dissolved
  calcium **rises**, because the pH falls 0.4 and shifts carbonate to bicarbonate.
  The familiar statement holds for a system buffered at fixed pCO₂, not a sealed
  one.
- **The w/c page was out of navigation, and not by accident.** Its heavy blocks
  were never executed and the scan contradicts its analysis on three points:
  ettringite never forms (the sulfate all goes to `monosulphate12` at
  equilibrium), no clinker survives at any w/c including 0.30, and the porosity has
  neither minimum nor inflection — so no Powers threshold, which is kinetic and
  not thermodynamic. Rewritten with executed blocks, the sealed-curing porosity
  convention (`porosity(eq, fresh)`, 5.3 points and a 7.4 % volume shrinkage apart
  from the one-argument form at w/c = 0.50), its assumptions written out, and put
  back in the navigation.
- Stale documentation removed: the `@test_broken` on `CaOH⁺` at ×2.47 has been
  gone since the convergence test moved to the true KKT error at `μ = 0`, and the
  `pKw = 13.979` that came with it was computed from the pre-fix `OH⁻` — it is
  13.9994 against Reaktoro's 14.0001.

### New tests

`test/published_values.jl` pins the malonate pKa and the retrograde calcite sweep
against the database's own `Kₛₚ`; `test/certified_equilibrium.jl` the certified
route, including that `b` must be fixed once; `test/equilibrium_constraints.jl` the
adiabatic enthalpy of neutralization, the refused volume constraint and the
implicit titrant; `test/kinetics/test_implicit_step.jl` the analytic kinetic step,
the impossibility of crossing saturation, several reactions on one mineral, and a
solid solution inside the step.


## v0.13.0 — a data file is named, not located

### Breaking changes

No name was removed or renamed and no signature changed. Below 1.0 the resolver
treats a minor bump as breaking regardless, so a downstream
`[compat] ChemistryLab = "0.12"` must be widened to `"0.13"`. In this repository
group that is `MeanFieldHomogenization.jl`'s `docs/Project.toml`, which pins
`ChemistryLab = "0.12"` and will otherwise refuse the new version.

One behavior did change, and it can only turn a former failure into a success:
the database readers now fall back to the bundled `data/` directory when a path
does not resolve against the working directory. A call that already worked keeps
resolving to exactly the same file — see below.

### `datapath`, because a script was only runnable from one directory

`build_species("data/slop98-inorganic-thermofun.json")` resolves its argument
against the working directory. That is the repository root when the script is
launched from a shell sitting there, and almost never otherwise: running the
same file from an editor whose REPL started elsewhere failed with

    SystemError: opening file "data/slop98-inorganic-thermofun.json"

Half the scripts in `scripts/` were written that way and half were not; the test
suite had independently converged on `joinpath(pkgdir(ChemistryLab), "data", …)`
in eighteen places, and the documentation pages carried
`"../../../data/…"`, a depth that silently depended on where the page sat in
`docs/src/`.

The exported [`datapath`](@ref) names a bundled file without reference to any
working directory:

```julia
substances = build_species(datapath("cemdata18-thermofun.json"))
ss_phases  = build_solid_solutions(datapath("solid_solutions.toml"), dict)
readdir(datapath())                     # what ships with the package
```

Every script, documentation block and test now uses it, so the code a reader
copies out of a page is code that runs.

### The readers resolve a bundled name, so older code keeps working

`read_thermofun_database`, `build_solid_solutions`, `extract_primary_species` and
`merge_json` (for its two inputs, never its output) resolve their path argument
in this order: as given, then under the bundled `data/` directory, then relative
to the package root. The working directory comes **first**, which is what makes
this safe: a call that resolves today resolves to the very same file tomorrow, a
local database still takes precedence over a bundled one of the same name, and
only a name that *is* one of the shipped files can reach the fallback — so a
mistyped path to a file of your own still fails, now with an error listing what
is available.

The banner these functions print shows the path relative to the package
(`data/cemdata18-thermofun.json`) rather than the absolute one. Documenter
captures that banner into the built page, and an absolute path would have baked
a build machine's directories into the documentation.

### Also

- `scripts/README.md`, which did not exist, records how to run a script and the
  one rule that keeps the collection working: a script consumed by the
  documentation or the test suite (`ionic_hydration.jl`,
  `hydration_calibration.jl`) never calls `Pkg.activate`, because the active
  project is global process state and switching it mid-build breaks every
  `@example` block that follows, on every page.
- The `Usage:` headers of the scripts said `julia --project` while the scripts
  activated `scripts/`; they now agree, and mention that running the file
  directly from an editor works.
- `using Revise` is gone from the three `tutorial_*.jl` scripts — it belongs in
  a `startup.jl`, and it resolved only through the machine's global environment.
  `SymPy`, which `tutorial_symbolic_reactions.jl` imports, is now declared in
  `scripts/Project.toml`, so that script runs in the environment it activates.
- The `merge_json` snippet in `docs/src/man/advanced.md` named a `.dat` file that
  does not exist (`cemdata18.dat`, where the shipped one is
  `CEMDATA18-31-03-2022-phaseVol.dat`) and passed a bundled path as the output
  argument. Both corrected.

### A red test suite, and the version pin that hid it

`test/kinetics/test_calibration.jl` errored with
`UndefVarError: NelderMead not defined in Main`, taking the whole suite down
(496 passed, 1 errored). The cause is upstream: **OptimizationOptimJL 0.4.19
stopped re-exporting Optim's algorithm names**, and
`scripts/hydration_calibration.jl` reached `NelderMead` through that re-export.

What made it look like a local mystery is that the script kept working: a
`scripts/Manifest.toml` is gitignored, and the one on the author's machine still
pinned 0.4.18, where the re-export was there. The test environment resolves fresh
from the registry, gets 0.4.19, and fails — so the suite was red on continuous
integration while every local run of the script was fine.

The script now does `import OptimizationOptimJL: NelderMead`. The binding exists
in that module either way, exported or not, so the explicit import works on both
versions and no longer depends on which one an environment happens to resolve.

## v0.12.0 — calibrated against a measured heat curve, and an isothermal heat that is a heat

### Breaking changes

**The ODE state layout gained one slot for `IsothermalCalorimeter`.** Code that
reads `sol.u` directly, or carries its own offset arithmetic into the state
vector, has to be revisited. `cumulative_heat(sol, ::IsothermalCalorimeter)` and
`heat_flow(sol, ::IsothermalCalorimeter)` also change what they return — from a
reaction extent in moles to a heat in joules. See below.

No name was removed or renamed and no signature changed. Below 1.0 the resolver
treats a minor bump as breaking regardless, so a downstream
`[compat] ChemistryLab = "0.11"` must be widened to `"0.12"`.

### `cumulative_heat` was returning a number of moles

`IsothermalCalorimeter`'s docstring promised that "cumulative heat
`Q(t) = ∫₀ᵗ q̇(τ) dτ` [J]" was "tracked as an extra ODE state". It was not.
`build_u0` appended a trailing slot only for `SemiAdiabaticCalorimeter`, and the
right-hand side wrote `du[end]` under the same condition, so with an isothermal
device `u[end]` was the **last reaction extent ξ_M, in moles**. `cumulative_heat`
returned it as a heat and `heat_flow` finite-differenced it. The function's other
branch, `Q(t) = H(0) − H(t)` from `p.saved_H`, was unreachable: `saved_H` and
`saved_t` are read in exactly two places in the package and written in none.

The isothermal calorimeter now contributes a real state integrating `dQ/dt = q̇`,
so `n_extra_states`, `extend_u0` and `extend_ode!` — until now reachable only from
the test suite — describe what the integrator actually does. The dead branch is
gone. Every index keyed on the state's length was audited; the semi-adiabatic path
already went through `n_extra_states(cal)` and is unaffected, and `p.has_T`
remains false for an isothermal device, so nothing collides on `u[end]`.

A new test compares both routes to the heat on the stoichiometric formulation,
where both are valid, and requires them to agree — the assertion that would have
caught this. The state-layout docstrings were also wrong in a second way, omitting
the reaction-extent block entirely; both now describe all four segments.

One caveat is stated where it belongs rather than left to be discovered. `q̇` here
is `heat_rate`, the heat of the **kinetic** reactions. That is the heat of
hydration when those reactions produce the hydrates. Under partial equilibrium
they only dissolve the anhydrous phases into ions, the hydrates are precipitated
by the Gibbs minimization, and this sum cannot see their heat — worth hundreds of
joules per gram on an ordinary Portland cement. `integrate` now warns on that
combination, as it already did for the semi-adiabatic device, and points at
`heat_release`.

### Calibrating hydration kinetics against measured calorimetry

`scripts/ionic_hydration.jl` has always run a complete CEM I forward with the
*published* Parrott & Killoh parameters, and its `IONIC_CALIBRATION` dictionary
has always been an explicit, deliberately unused hook. Nothing in the repository
closed that loop: there was no experimental dataset, no loss function and no
inverse problem anywhere in `src/`, `scripts/` or `test/`.

`scripts/hydration_calibration.jl` and
[its documentation page](https://micropochemomechanics.github.io/ChemistryLab.jl/stable/examples/hydration_calibration/)
close it, on real measurements, fitting **kinetic parameters only** — the CEMDATA18
thermodynamics are measured quantities, and the Blaine fineness, w/b ratio and
temperature are reported by the experiment.

**The data.** Two CEM I records now ship under `data/experimental/`, subsets of the
CC-BY-4.0 Zenodo deposit of Šmilauer & Reiterman
([10.5281/zenodo.15212785](https://doi.org/10.5281/zenodo.15212785)): a CEM I
52.5 R at w/b 0.50 and Blaine 397 m²/kg over 262 h, fitted, and a CEM I 52.5 R at
w/b 0.45 and Blaine 415 m²/kg over 617 h, held out. They are the only source found
that is at once openly licensed for redistribution, real CEM I measured by
isothermal calorimetry, and documented well enough to *fix* the non-kinetic inputs
— every file states the fineness, the w/b ratio, the temperature and the
normalizing binder mass.

`regenerate.jl` in that directory rebuilds both files **byte for byte** from the
deposit, so nothing about them is unverifiable, and `--check` is an exact
comparison. The subset is the rows nearest to 500 log-spaced times, by nearest
index rather than by interpolation: every number in the vendored files is one the
calorimeter reported. A fourth column carries the depositors' own fitted curve,
read from their tabulated data rather than reimplemented from a functional form
recalled from memory, so a fit can be judged against somebody else's on the same
record. `data/experimental/LICENSE` states the terms — CC-BY-4.0, separate from
the package's LGPL — and the companion article, published CC BY-NC-ND 3.0, is
cited with nothing reproduced from it.

**The published parameters are already close.** With no adjustment at all the
coupled model gives Q(262 h) = 381.9 J/g against 376.0 J/g measured, +1.6 %. The
depositors' four-parameter affinity model, fitted to this very record, gives
388.9 J/g. So parameters fitted in 1984 to other cements, carried through
CEMDATA18 thermodynamics and a Gibbs minimization, land closer at the endpoint
than a purpose-fitted empirical curve. What a calibration can earn here is the
*timing*.

**Three parameters, not six, and the number is measured.** Six candidates were
written down, one per mechanism with a distinct signature on the curve. The
singular values of `∂Q/∂log θ` on the coupled model come out
`[422, 101, 60, 6.3, 1.4, 0.20]` — a factor of nine between the third and the
fourth — so the measurement determines three *combinations*. The correlation
matrix says which: `k₁_C3S` with `n₁_C3S` at −0.985, and `k₃_C3S` with `n₃_C3S`
at −0.977, the latter pair being essentially the whole leading singular direction.
One per pair is all that can be fitted, and the rate constant is kept in each; the
third is `k₁_C3A`. Belite's `k₃_C2S` carries a weight of 0.001 in the three
directions the data see, which is a quantitative way of saying that eleven days of
heat cannot see belite.

Activation energies are not fitted at all, and the reason is not caution: all these
records are at 20 °C, an activation energy is a temperature sensitivity, and a fit
that reported one would be reporting a number the data cannot contain.

**Two failures worth having in the record.** A cheap stand-in for the coupled
model — same rate laws, hydrate assemblage written down by hand — was tried and
does not work, in two distinct ways. First, `C₃A + 3 Gp + 26 H₂O → ettringite` is
stoichiometrically impossible in this mix: 11 % C₃A demands about 1.1 mol of gypsum
per kilogram of binder and 4.6 % gypsum supplies 0.27, so a rate law reading only
C₃A drives the gypsum **negative** and the released heat to 1609 J/g against a
measured 376. Second, with the aluminate sent to hydrogarnet instead, the level is
right to some 9 % but a fit on the *normalized* curve saturates whatever parameter
bounds it is given, whichever subset is fitted — the imposed assemblage produces a
curve shape the Parrott–Killoh family cannot make. Both are now documented sections
of the page and assertions in the test suite, and together they are the
quantitative argument for the model `ionic_hydration.jl` already implements.

**A docstring claim that does not survive measurement.** `parrot_killoh_avrami`
states, after Parrott & Killoh, that "C₃S has no diffusion-controlled stage". The
sensitivity of the released heat to `k₂` for alite was expected to be exactly zero
and is not: the Jander term `k₂(1-ξ)^{2/3}/(1-(1-ξ)^{1/3})` falls as ξ grows, so
past a high degree of hydration it does become the minimum of the three branches,
and `k₂` moves the heat by up to about an eighth of what `k₁` does — all of it
after the first day. The published statement describes the 1984 fit's intent; this
implementation of it has a diffusion-limited tail. Pinned by a test so the claim
cannot drift back.

**Parameter-space automatic differentiation does not work, and the manual said it
did.** `docs/src/man/kinetics.md` claimed "the entire chain is
ForwardDiff-compatible" and showed a derivative through `integrate`. The ODE
right-hand side is AD-clean in the state and in time — deliberately, because
`Rodas5P` needs the time gradient — but not in the parameters: `build_u0` returns
a `Vector{Float64}`; `build_kinetics_params` casts the temperature, the initial
amounts, the stoichiometry and the calorimeter constants to `Float64`; the
surface-area constructors cast to `Float64` despite being declared `{T<:Real}`;
`system_enthalpy` accumulates into a `0.0`; and `respeciate!` is `Float64`-only by
construction. There was no test through `integrate`. The claim is corrected, the
blocking lines are named, and the calibration uses `AutoFiniteDiff()` with
`NelderMead` — still the SciML interface, `Optimization.jl`, and honest about why.

**Reusable, not general.** The helpers live in the script, not in `src/`, on
purpose: the interface should be settled by use before it is exported. The seam for
somebody else's data is `read_calorimetry`, which reads a documented CSV and takes
the experimental conditions from its comment header, and `CALIB_SPEC`, a list of
(phase, rate-law field, prior, box) rows.

### The dormant period is now on by default in `scripts/ionic_hydration.jl`

`run_ionic_hydration` and `ionic_reactions` take `induction` and
`induction_phases`, and the default is **on**: `τ = 5 h`, `m = 2.5`, applied to
the two silicates. A CEM I has a dormant period; Parrot–Killoh does not.

The numbers are round on purpose. The calibration returns 5.6 h and 3.56 on one
record and then shows `τ` correlated with `k₁_C3S` at 0.994, so the data barely
determine either separately. Carrying 20 228 s into an uncalibrated example would
lend a precision the measurement does not support. `induction = nothing` recovers
the previous behavior exactly.

What it changes is the early heat and little else. Measured over 28 days on the
example's own formulation, `Q(6 h)` falls from 74.1 to 40.2 J/g — the point of the
exercise — while the 28-day pore solution (pH 12.754) and porosity (0.3635 against
0.3634) do not move, and the integrator takes 216 accepted steps instead of 200.

Re-measuring both configurations caught an unrelated drift.
`IONIC_DEFAULT_SYSTEM`'s docstring claimed 202 steps and pH 12.58, where the model
*without* any dormant period now gives 200 and 12.754 — that figure had gone stale
for some earlier reason, and the docstring now says so rather than letting the new
default take the blame.

### Downstream: MeanFieldHomogenization.jl

Two things follow for MFH, and the second is not optional.

1. **`docs/Project.toml` pins `ChemistryLab = "0.11"` and must be widened to
   `"0.12"`.** Below 1.0 the resolver treats a minor bump as breaking whatever the
   API did, so the MFH documentation will not resolve against this release until
   that bound moves. MFH's own `Project.toml` does not depend on ChemistryLab at
   all — the coupling lives only in its docs environment.
2. **MFH keeps its own copy of the ionic setup**, in
   `scripts/common/ionic_hydration.jl`, with its own `IONIC_CALIBRATION` and its
   own `parrot_killoh_avrami` call. The dormant period has to be added there, and
   in `scripts/common/stoichiometric_hydration.jl` for the stoichiometric route — on the
   clinker silicates only, since the silica fume there already goes through
   `waller`, whose sigmoid carries its own onset delay.

   It will move a *reported* number, not just an input.
   `scripts/common/paste_micromechanics.jl` returns `E = 0` until the hydrate foam
   percolates — the setting transition as a genuine zero of the self-consistent
   fixed point — and `45_ionic_hydration_micromechanics.jl` reports a setting time
   in hours from the first non-zero modulus. Without a dormant period the model
   accumulates hydrates from `t = 0`, so percolation, and the setting time with it,
   arrives too early.

### `heat_release`'s `states` keyword never worked

`heat_release` accepts `states` so a caller who already holds the certified replay
does not pay for it twice, and its docstring says so: *"it is the expensive part,
and a caller that needs the compositions anyway should not pay for it twice"*. The
line that implemented it was

```julia
states = something(states, speciated_states(sol, kp; times = times))
```

and `something` is an ordinary function, so Julia evaluated **both** arguments and
the replay ran whether or not the states were supplied. Measured on the ordinary
Portland cement of `scripts/ionic_hydration.jl` over 28 days, 60 instants:

| call | before | after |
|--- |--- |--- |
| `heat_release(...; states = sts)` | 121 s | **0.017 s** |
| `heat_release(...)`, replay included | 118 s | 118 s |

`docs/src/examples/ionic_hydration.md` computes its replay once and hands it to
both `ionic_phase_history` and `heat_release` for exactly this reason; it has been
doing the replay twice. A timing assertion in `test/kinetics/test_postprocessing.jl`
now guards it, deliberately — a redundant computation whose result is discarded
cannot be detected any other way.

Note what this does **not** speed up: a calibration loop calling `heat_release`
without `states` already paid one replay and still does. The replay is the price of
correctness, and the same docstring records why — reading the running composition
instead gave 1174 J/g where the certified answer was a few hundred.

Found by benchmarking, not by reading. The hypothesis under test was that the
enthalpy sum dominated; it costs 0.029 s. An optimization written for it was
reverted, because fifty lines and a fallback branch to save 26 ms is not an
improvement, and the measurement that refuted it is what exposed the real defect.

### Proper names out of the identifiers

Author names belong in the bibliography, not in constants. In the two example
scripts and the page that includes one of them:

| was | is |
|--- |--- |
| `LAVERGNE_MIX_C100`, `_VESSEL_CP`, `_LOSS_A`, `_LOSS_B` | `CALORIMETRY_MIX_C100`, `_VESSEL_CP`, `_LOSS_A`, `_LOSS_B` |
| `lavergne_semiadiabatic` | `semiadiabatic_cell` |
| `PK_L2018_C3S` and siblings | `PK_SMOOTHED_C3S` and siblings |

Nothing in `src/` is affected and no exported name changes. `Lavergne2018` and
every citation in prose stay: that is where the attribution belongs.

The names of *models* stay too, and the distinction is deliberate.
`parrot_killoh`, `parrot_killoh_avrami` and `waller` are how the literature
designates those rate laws, in the same way it says Arrhenius, Avrami, Jander,
Powers, Langavant or Blaine — they are technical terms, not credits, and they are
part of the public API. What has been removed is the name attached to a *source of
parameter values*, which is a citation wearing an identifier's clothes.

### Also

  - `ionic_reactions` and `run_ionic_hydration` take a `pk_params` keyword that
    overrides the published Parrott & Killoh sets, so the calibration varies the
    rate laws of the existing model instead of restating it. The default sets are
    now a single `IONIC_PK84` constant rather than a dictionary literal repeated
    per call.
  - `scripts/opc_semiadiabatic_calorimetry.jl` cited Lavergne et al. (2018) with
    the DOI `10.1016/j.cemconres.2017.11.007`, which Crossref resolves to a
    different paper (Machner et al., *CCR* **105**, 1–17). Corrected to
    `10.1016/j.cemconres.2017.10.018`. The same defect was found and fixed in
    `docs/src/refs.bib` in v0.4.0; the script was missed then.
  - Nine bibliography entries added, every DOI resolved against Crossref and the
    Šmilauer entries also against the publisher's own landing page. The dataset and
    the article have different author lists — two creators and three respectively —
    and each is credited for its own artifact.

## v0.11.0 — a proved answer, whichever route found it

### Breaking changes

One new exported name, `solve_certified`. Nothing was removed or renamed and no
existing signature changed, but below 1.0 the resolver treats a minor bump as
breaking regardless, so a downstream `[compat] ChemistryLab = "0.10"` must be
widened to `"0.11"`.

### `solve_certified`: several routes, one proof

`solve_certified(des, starts; b)` solves from each starting composition in turn and
returns the first answer `optimality_certificate` **proves** optimal, or — if none
is proved — the one with the smallest KKT error together with its certificate, so
the caller always sees what it is getting.

Offering several routes is rigorous here rather than opportunistic, and the reason
is the certificate. For a convex problem it does not rank answers, it decides them:
the proof is the same proof whichever start produced the point, and nothing about it
depends on having predicted the winner. What would be a fudge is choosing a route by
taste and reporting its output unproved.

It is also necessary, because no back end dominates. Measured on an LC³ equilibrium
solved COLD at four degrees of reaction, with the dual Newton started from each
interior-point back end in turn:

| degree of reaction | from `OptimaOptimizer` | from `IpoptOptimizer` |
|---|---|---|
| 0.05 | certified, 2.7e-12 | certified, 9.1e-13 |
| 0.25 | **not certified**, 7.2 | certified, 4.2e-12 |
| 0.50 | **not certified** | certified, 4.6e-12 |
| 1.00 | certified, 1.8e-11 | **not certified**, 3.5e-3 |

Each solves a case the other misses. With both offered, every degree of reaction is
certified from a cold start — where before, an intermediate one could only be
reached by continuation from a nearby solution. The starts are supplied by the
caller, so this adds no dependency: whichever back ends are loaded are the ones
available.

### Two assertions corrected because they were wrong

Neither was in the way; both stated something untrue.

The negative control on water autoprotolysis asserted that forcing the Schur
complement back on gives `[H⁺]/[OH⁻] > 2`, pinning a direction of error. Its point
is that the null-space step is what carries the correct ratio, so the assertion is
now that the answer is NOT one — it used to come out at 3.78 and now at 8.7e-5, both
far from unity, and pinning the sign made the test fail for a change that did not
touch what it is about.

The `solve_certified` testset was written asserting `9.5 < pH < 10.3` for a system
that carries 0.01 mol of CO₂. That is the pH of calcite in PURE water, from a
different test; with the CO₂ the answer is 6.43, which is Reaktoro's own value for
this system (`H⁺ = 3.69e-7` in the reference table). The certificate was right and
the assertion was not.

### `speciated_states` walks up to the first instant it is asked for

Every replayed speciation warm-starts from the previous one — and the FIRST one
requested has no previous one, so it started from the cast composition, which
carries no active set at all. The chain is now walked up to it through a few
earlier times whose compositions are discarded and whose only purpose is to hand
the guess an active set.

The consequence was not cosmetic. On the reference OPC replayed at forty instants,
the first of them came back with **56 interior species where the answer has 25** —
every candidate hydrate present, four of them at 1e-5 to 1e-6 mol, the signature of
an interior-point iterate that never reached a vertex — and neither the certifying
Newton nor its continuation recovers from that, because both inherit the start.
That instant is now certified, with stationarity 9.1e-13 and element balance
6.1e-15 mol: the whole replay is proved optimal where before one instant in forty
fell back silently to an uncertified composition. The 28-day heat of that paste
moves from 420.2 to 420.3 J/g, which is the size of the error the fallback was
hiding.

Two other routes to the same instant were implemented and measured, and neither is
in the code: retrying it from each of the thirty-nine certified neighbors, nearest
first, left it unproved with every one of them; and the continuation between
certified instants was rewritten to step forward adaptively rather than bisect —
correct in itself, since bisection lowers the upper end on failure and so abandons
the target for good, but it does not reach this instant either. What was missing was
the start, and only the run-up supplies it.

### The full Portland cement, through its pore solution

New example page and script, `scripts/ionic_hydration.jl`: a complete CEM I —
alite, belite, aluminate, ferrite, gypsum, limestone filler — dissolving into ions,
with the whole hydrate assemblage decided by Gibbs minimization at every accepted
step. The aluminate cascade that a stoichiometric model has to encode by hand comes
out as a result: run with 3.5 % limestone the ettringite survives to 28 days, run
without it the ettringite is depleted to zero, and no line of code distinguishes
the two cases.

The page carries the calorimetry with it, isothermal and semi-adiabatic
(NF EN 196-9), with every number of the cell written out: 420 J/g against 405 J/g
at 28 days, a rise of 19 K against an adiabatic 75 K. The heat is taken from the
enthalpy of the certified compositions, not from the kinetic reactions — those only
dissolve the clinker, so the sum `Σ rᵢ(−Δ_r H⁰ᵢ)` cannot see the precipitation heat
at all, and driving the cell from it gave a rise of 207 K.

## v0.10.0 — the heat of hydration, and solid solutions the solver can hold

### Breaking changes

Five new exported names — `enthalpy`, `heat_capacity`, `missing_enthalpy`,
`system_enthalpy` and `heat_release`. Nothing was removed or renamed and no
existing signature changed, but below 1.0 the resolver treats a minor bump as
breaking regardless, so a downstream `[compat] ChemistryLab = "0.9"` must be
widened to `"0.10"`.

`[compat] OptimaSolver` moves to `"0.4"`, which is required: the solid-solution
work below depends on `SolutionPhase(...; mole_fraction = true)`, introduced
there.

### Calorimetry that a coupled model can actually use

`enthalpy(state)` and `heat_capacity(state)` sum `nᵢ ΔₐH⁰ᵢ(T,P)` and
`nᵢ Cp⁰ᵢ(T,P)` over the composition, and `heat_release(sol, kp; times)` returns
the cumulative heat and the heat rate along a trajectory. This is Eqs. (17)–(21)
of Lavergne et al. (2018): enthalpy is a state function, so its drop between two
states at the same temperature is the heat given off, with reactants, ions and
hydrates each counted once and no reaction stoichiometry to write down.

That last point is what makes it necessary. `heat_rate` sums `rᵢ(−ΔᵣH⁰ᵢ)` over
the KINETIC reactions, which is right when those reactions produce the hydrates —
and wrong under partial equilibrium, where they only dissolve the anhydrous
phases into ions and the hydrates are precipitated by the Gibbs minimization.
Driving a semi-adiabatic cell from it put an ordinary Portland cement at a
temperature rise of 207 K. That combination now warns instead of returning the
number in silence.

`heat_release` reads the **certified** speciations of `speciated_states`, not the
composition the integrator carries. The in-run minimization is warm-started and
uncertified, and a single hydrate is worth hundreds of kilojoules: read that way
the curve came out at 12.7, 145, 1174, 936 and 631 J/g at 1 h, 6 h, 12 h, 1 d and
2 d — heat that rises and then falls. Certified, the same cement gives 185, 288,
347 and 420 J/g at 1, 3, 7 and 28 days, and the curve is monotone.

`missing_enthalpy(state)` lists the species that carry no `ΔₐH⁰` and are
therefore absent from the balance, because a heat curve missing one hydrate is
not visibly wrong.

### Every species now matches Reaktoro, `CaOH+` included

`equilibrium_reference.jl` recorded `CaOH+` as `@test_broken`: out by a factor 2.5
at 4e-9 mol, and left alone because Ipopt landed on the same point, which made it
look like a property of the problem rather than of one solver. It was neither. The
interior-point stage was reporting points it had reached without ever satisfying
the KKT conditions at `μ = 0` — see `OptimaSolver` v0.4.0 — and the trace species
are exactly where that shows. The assertion is now a plain `@test`.

### `s[:Cp⁰]` returned zero on a species that had one

The thermodynamic functions are built on demand. `getproperty` knew that;
`getindex` did not, and returned its not-found value `0` — an `Int64` — for any
of `:Cp⁰`, `:ΔₐH⁰`, `:S⁰`, `:ΔₐG⁰`, `:V⁰` on a species whose functions had not
been forced yet. Callers found out only when they tried to evaluate it. Returning
0 is right for a missing ATOM, which is what that fallback is for; it was never
right for a property the species can produce.

### A solid solution now has a reference worth the name

`DualEquilibriumSolver` took the first end-member of every solid solution as the
phase reference. The reference is the member whose stationarity the outer system
carries rather than inverting, so it has to be the one the phase is mostly made
of — for the aqueous phase that is the solvent, and not by convention. Picking a
minor end-member sends the outer unknown `ln x_ref` towards −∞ and the Jacobian
with it. It is now chosen by magnitude from the caller's own composition, and
solid solutions are declared to `OptimaSolver` as mole-fraction phases, which is
what lets them be solved at all.

## v0.9.0 — equilibria that are proved, not hoped for

### Breaking changes

Two new exported names, `DualEquilibriumSolver` and `optimality_certificate`.
Nothing was removed or renamed and no existing signature changed, but below 1.0
the resolver treats a minor bump as breaking regardless, so a downstream
`[compat] ChemistryLab = "0.8"` must be widened.

### The interior-point solver could not say whether its answer was the answer

`EquilibriumSolver` minimizes `G` by an interior-point method and, on a cement,
it stops on `MaxIters` — at any tolerance. That was known. What was not known is
how much it costs, because nothing in the package could check.

The Gibbs problem is **convex**: an ideal mixing entropy, which is convex, plus
terms that are LINEAR in the amounts of the pure phases, whose potential does not
depend on how much of them there is, over the polyhedron `{A n = b, n ≥ 0}`. The
minimizer is therefore unique and every stationary point is global — so a solver
returning different answers from different starting points is not finding local
minima, it is stopping short of stationarity, and the KKT conditions are
*sufficient*: they can be checked, and checking them is a proof.

**`optimality_certificate`** does that check, on any composition and whatever
produced it. It reports the stationarity of the interior species, the component
balance, and the worst supersaturation among absent phases. A species below the
floor is at its bound, where the condition is the inequality `μ + Aᵀy ≥ 0` and
not the equality — imposing the equality on an amount held at `1e-16` whose
mass-action value is `e⁻³⁰⁰` misstates its log-activity by 263 RT units.

### `DualEquilibriumSolver`

Newton on the KKT system in element-potential space — the Brinkley–Karpov
formulation of the geochemical Gibbs-minimization codes. Writing `u = −Aᵀy`, an
aqueous species obeys `aᵢ = exp(uᵢ − gᵢ)` and a pure phase is present exactly
when `uᵢ = gᵢ` and absent when undersaturated, which is the classical
phase-stability criterion.

**The algorithm is not in this package.** It is
`OptimaSolver.dual_newton_solve`, where it belongs: the active set, the
degeneracy test on the rows of `A`, the two-level Newton and the certificate are
statements about a convex program, not about chemistry. What this package
supplies is the four things that *are* chemistry — the conservation matrix and
the reference potentials `Δ_aG⁰/RT`; the activity model as the callback `h` with
`∇f = g + h(n)`; which species are strictly positive at any solution and which
may vanish; and which species is the solvent, whose activity is a mole fraction
and therefore bounded above, so that its stationarity cannot be inverted. This
release requires OptimaSolver 0.3.

Two levels. The inner one inverts the **solutes'** mass-action laws at fixed
potentials and fixed solvent amount; the outer is a Newton on `1 + m + |P|`
unknowns — the solvent, the `m` element potentials, the amounts of the active
phases — some fourteen numbers for a cement, against forty-seven species in the
interior-point route. Parameterizing the solutes by `ln n` makes their positivity
automatic, so the fraction-to-boundary rule that capped the interior-point step
at *every* iteration has nothing left to act on.

Six things had to be right, and each was found the hard way:

  - **the solvent cannot be inverted through its own mass-action law.** Its
    activity is a mole fraction, so `ln a_w ≤ 0` always, and an arbitrary `y` can
    demand `ln a_w > 0`, for which no finite composition exists. It belongs to the
    outer system, where the balance determines it;
  - **two phases declared saturated over-determine `y`.** Their stationarity rows
    are then jointly infeasible. With portlandite wrongly admitted at 1e-6 mol
    from the starting guess, the residual sat at 29.5 and fell by 1e-3 an
    iteration; dropping it *during* the Newton rather than after, the same solve
    reaches 5e-12 in eleven steps;
  - **an amount must not be clamped to `ϵ` on the way out**, for the reason given
    above;
  - **the active set must be updated during the Newton, not after it**, and one
    phase admitted at a time. Two phases both declared saturated over-determine
    `y` and their stationarity rows are jointly infeasible; admitting a batch
    feeds a cycle in which a phase is admitted, driven negative, dropped, and
    readmitted. On a cement without limestone that produced a solve converged to
    2e-12 — of the wrong subproblem, with an absent phase supersaturated by 10.9.
    Visited active sets are recorded, which makes the loop finite in practice as
    well as in theory;
  - **a component whose total has vanished must be removed from the unknowns.**
    Its element potential is determined by nothing, the Jacobian is singular in
    that direction, and this is the ordinary state early in a hydration run, where
    the iron, sulfur and aluminum totals are at machine zero. The test is not
    `bₖ ≈ 0` but `bₖ ≈ 0` **with the non-zero entries of row `k` sharing a sign** —
    only then does `Σ Aₖᵢnᵢ = 0` with `n ≥ 0` force each term to vanish. The `H⁺`
    row carries `+1` for `H⁺` and `−1` for `OH⁻`, so its zero total is the
    ordinary state of pure water; treating it as degenerate kills the entire
    acid–base system and returns pH 7.000 with the calcite undissolved;
  - **the replay must be seeded and, if need be, continued.** `speciated_states`
    tries a ladder of early instants until one certifies, then walks forward; when
    an instant still resists it bisects the interval in `bₑ` from the last
    certified one, which is a homotopy in the component totals and terminates.
    Without the ladder the first requested instant of a cement could not be
    proved — neither the interior-point answer nor the cast composition puts the
    Newton close enough.

### Measured

On the package's own Reaktoro reference — calcite, CO₂ and water, ideal
activities — the certified answer matches **every species to 1 %**, including
`CaOH+`, which `equilibrium_reference.jl` records as `@test_broken` because the
interior-point answer is 147 % high. On calcite in pure water the certified pH is
**9.90**; the interior-point answer is 6.96, which is not an imprecision but a
wrong answer, and nothing in that solver's output reveals it.

`speciated_states` now certifies each instant it replays, and **names any it
cannot**: falling back silently would leave the caller unable to tell a certified
trajectory from an uncertified one, and those are not the same object. On a full
ordinary Portland cement over 28 days, with and without limestone, **all forty
replayed instants of both mixes are certified**, with element balances between
1e-11 and 1e-13 mol — where the six-hour instant, the AFt peak at which the
assemblage switches, stood at 6.9e-2 mol before this work.

### What is not certifiable yet

**Solid solutions.** `DualEquilibriumSolver` refuses a system that declares one,
and does so explicitly rather than returning a number. A pure phase has unit
activity, hence `gᵢ = uᵢ` and a bound-constrained variable with an active set. An
end-member of a solid solution has `μᵢ = gᵢ + ln aᵢ(n)`, whose activity goes to
`−∞` as its mole fraction goes to zero: it is never exactly absent while the
phase exists, so the active set belongs at the level of the PHASE and not of the
species. That is a different algorithm. Treating an end-member as a pure phase
would not fail loudly — it would return a composition that looks reasonable and
is not the minimum — which is why the refusal is explicit and
[`speciated_states`](@ref) reports such instants as uncertified.

This matters for what comes next: slag and calcined-clay binders form C-A-S-H,
hydrotalcite and AFm solid solutions, and certifying those needs the phase-level
active set.

**The in-run speciation.** `respeciate!` still uses the interior-point solver, so
the compositions *inside* the integration are not certified; only the replay is.
Measured, that makes no difference here, because `bₑ` is integrated from the
rates alone and a Parrot–Killoh or Waller law reads only its own degree of
reaction. A rate law reading log-activities — a saturation ratio, or a
pH-dependent dissolution law — would feed the speciation back into the
trajectory, and would need this closed first.

### Still true

The interior-point solver is unchanged and still does not report convergence on a
cement. It is what gets the certifying Newton into the neighborhood — on 74 of
80 measured cement equilibria the certified answer is reached from its guess —
and that division of labor is deliberate.

## v0.8.2 — a replay that conserves matter, solid solutions that reach the solve

Bug fixes. No API changed and nothing was removed, so a downstream
`[compat] ChemistryLab = "0.8"` accepts this release unchanged.

The `docs/src/examples/coupled_hydration.md` page is rewritten on the real API:
it described four calls that never existed and had no executed block, so nothing
had ever checked it. Every block now runs when the documentation is built.

The `OptimaSolver` bound is relaxed from `"0.2.7"` to `"0.2"`. The patch-level
lower bound guarded against 0.2.6, where `OptimaOptimizer` discarded the
caller's `u0`; a fresh resolution always takes the newest patch, so the bound
only ever mattered when reviving an old manifest.

### The relative residual saturated on near-empty rows

`_row_residual` scaled each row by the larger of its own element total and the
matter flowing through it. For a row whose total is a millionth of the largest —
an element barely present, or the charge row early on — that divisor is
essentially zero, and a rounding error came out as an alarm. A full ordinary
Portland cement balanced to **1.4e-10 mol** at 28 days was reported at `3.2e-2`,
and the near-empty rows of the first steps saturated the measure at `1.0`, which
read as a 100 % violation and was nothing of the sort. Every row is now floored
at a millionth of the system scale, and the same run reports `3.9e-6`.

### A replayed speciation could stop on a broken balance

Two guards were added to [`speciated_states`](@ref), both judged on conservation
of matter rather than on a retcode:

- **A retry from a guess carrying no active set.** The warm start is the previous
  speciation, and it brings the previous *active set* with it. Where the
  assemblage switches — an ordinary Portland cement around six hours, when the
  ettringite peaks and the sulfate starts moving into AFm — that set is the wrong
  one and an interior-point method started inside it does not cross over. When
  the balance exceeds `1e-6 mol` the instant is solved again from the cast
  composition, and whichever result closes better is kept.

- **A feasibility restoration that actually converges.** Alternating projection
  converges linearly, at a rate set by the angle between the box and the affine
  set, and on a cement that angle can be small: at that same six-hour instant,
  where the iron row carries 0.013 mol across thirteen species, 200 sweeps left
  a residual of 8.4e-1 and 20 000 were needed to reach 6.7e-9. A replay runs a
  handful of times and now asks for the budget it needs.

  Inside the ODE right-hand side the budget stays where it was: raising it there
  made the run *worse* — the worst in-run balance went from 1.1 mol to 8.5 at
  2000 sweeps and to 41 at 100 000 — which is the same unpredictability the
  back-end shows elsewhere and is not understood.

Measured on a full OPC replay, with OptimaSolver 0.3: the six-hour instant —
the AFt peak, where the assemblage switches — goes from **6.9e-2 mol to
7.1e-15**, and every other instant lands between 1.8e-15 and 3.6e-9. Over the
201 accepted steps of the coupled run the worst balance is 6.0e-6 mol and the
median 2.7e-9. On the 40-instant grid the chapters use, the porosity and the
elastic modulus are both monotone across every interval and the pore solution
holds at pH 12.52–12.59.

### Solid solutions were dropped from the equilibrium partition

`_equilibrium_subsystem` rebuilds a `ChemicalSystem` for the partition the
equilibrium is actually solved on, and it did not carry `solid_solutions` over.
The loss was silent and total: a coupled run could declare CSHQ, AFm or
Hydrogarnet and the solve would treat their end-members as separate pure phases,
the mixing entropy never entering the Gibbs energy. Measured on alite and belite
with the four CSHQ end-members, the run was bit-identical with and without the
declaration. With the fix the two differ — the element balance goes from 1.5e-12
to 3.6e-15 and the pore solution from pH 12.48 to 12.37.

A solution survives into the partition only if **all** its end-members are in it;
one split across the kinetic and equilibrium sides is not a phase the equilibrium
can mix, and is dropped rather than passed truncated.

This makes the code path work. Whether a given solid solution is *stable* in a
given system is a separate, thermodynamic question: with these Cemdata18
end-members the small silicate system above still precipitates no C-S-H.

### The warning could not tell the trajectory from a probe

The right-hand side is evaluated far more often than the solution advances —
Jacobian finite differences and rejected steps included — and a poor speciation
on a *perturbed* `bₑ` never enters the result. Reporting the worst over all
evaluations alarmed about something that does not affect the answer.

`integrate` now reports the worst on the **accepted steps** first, and the
all-evaluations figure second. It also states when the distinction matters:
`bₑ` is integrated from the rates alone, so a rate law reading only its own
degree of reaction (Parrot–Killoh, Waller) gives a trajectory independent of the
speciation, while a law reading log-activities feeds it back in.

### The warning now gives moles

`integrate` reported only the relative figure, which cannot be judged without
knowing what it was divided by. It now leads with the absolute worst in moles —
`1e-10 mol` is machine precision whatever the system, `1e-2 mol` against a
0.3 mol sulfate budget is not — and says where such cases come from: the first
steps, where the paste has barely reacted. On the OPC above the honest reading
is 1.05 mol at the worst early step and machine precision from three days on.

## v0.8.0 — a public replay of the equilibrium partition

### Breaking changes

- **New exported name: `speciated_states`.** Below 1.0 Julia's resolver treats a
  minor bump as breaking whatever the API did, so a downstream
  `[compat] ChemistryLab = "0.7"` will **not** accept 0.8 and must be widened.
  Nothing was removed or renamed, and no existing behavior changed, so widening
  the bound is the only adjustment required.

### The composition at a given time was not recoverable

`state_at` returns the purely kinetic reconstruction `n(0) + νᵀξ` and says so:
the redistribution performed by the equilibrium solve is not recoverable from
the stoichiometry. For a cement that reconstruction is meaningless — every
dissolved element in solution, not one hydrate. The composition left in the
solver's buffers is no better: it is rewritten at every right-hand-side
evaluation, Jacobian differences and rejected steps included, so it is not the
accepted composition at any time. Reading it is what made a full OPC look as
though it held 0.244 mol of ettringite against 0.267 mol of sulfate.

**`speciated_states(sol, kp; times)`** does the recovery properly: the kinetic
species from the ODE state, and the equilibrium partition re-solved from the
element totals the run carried, walking the instants in order.

Two things it must do, both found the hard way:

- **cap and restore the guess.** The warm start is the equilibrium of the
  *previous* `bₑ`, so once an element has been spent it demands more of it than
  now exists and the interior-point solve begins outside its own feasible set.
- **use a back-end instance the integration has not touched.** Replaying through
  `p.eq_solver.solver` — the object that drove thousands of solves during the
  run — returned pH 14.2 with 0.31 mol of ettringite and no AFm, where a clean
  instance of the same type and settings gives pH 12.58 with the sulfate
  entirely in AFm, on the same trajectory, the same `bₑ` and the same guess.
  **The back-end carries state across solves**; this is a defect upstream, and
  until it is fixed a replay must not inherit it.

The first requested instant is the loose one, having only the cast composition
to start from; every instant after it lands at machine precision. Ask for a
first point within the first hours.

### The back-end defect is fixed upstream

The state the replay had to avoid was `OptimaOptimizer._cache`: with
`warm_start = true` the algorithm object started every solve from its previous
solution, discarding an explicit `u0` — including the guess this package builds
in `respeciate!`, so the caps and projections above were partly defeated during
the run itself. **OptimaSolver 0.2.7** consults the cache only when the caller
supplies no interior point, and `[compat]` now requires it. `speciated_states`
still builds a clean back-end instance, which costs nothing and keeps the replay
correct against an older back-end.

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
