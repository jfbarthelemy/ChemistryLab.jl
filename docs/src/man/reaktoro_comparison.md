# Validation against Reaktoro

[Reaktoro](https://reaktoro.org) solves the same problem this package does, from
the same literature ([Leal2017](@cite)), and is mature and widely used. It makes
a good oracle — provided the comparison is set up so that a disagreement means
something.

This page reports what agrees, what does not, and by how much. The scripts that
produce the reference values live in `test/reference/` and the assertions in
`test/equilibrium_reference.jl`, so every number below is reproducible and is
checked on every CI run.

## Setting up a comparison that means something

A cross-code comparison has **three** knobs, and each is worth tens of percent:

| knob | must match | what happens otherwise |
|:--|:--|:--|
| thermodynamic database | same file | different standard-state data, different answer |
| species list | same species | see below — this one is subtle |
| activity model | same convention | `∂Ca²⁺/∂(CO₂)` moves from `+0.152` to `+0.218` |

Reaktoro reads ThermoFun databases, so it can be pointed at **the very file this
package ships**, `data/cemdata18-thermofun.json`, with identical species names.
This package's default is `DiluteSolutionModel()`, matched on the Reaktoro side
by `ActivityModelIdealAqueous()`.

!!! warning "The species list is not a detail"
    Drop the aqueous calcium complexes and free `Ca²⁺` becomes the only aqueous
    home for calcium, so `∂Ca²⁺/∂(CO₂)` and `∂calcite/∂(CO₂)` mirror each other
    exactly. Restore `CaOH⁺`, `Ca(CO₃)@` and `Ca(HCO₃)⁺` and that mirror breaks
    — the difference is what `Ca(HCO₃)⁺` takes up. **Both codes break it the
    same way**, and the element balance closes to `2×10⁻¹⁶` either way. Comparing
    a system that mirrors against one that does not is comparing two different
    chemistries, not validating one.

## The test case

Calcite in water with dissolved CO₂ — small enough to reason about, and it
exercises a solid in equilibrium with an aqueous phase spanning ten orders of
magnitude in amount:

```
H2O@ 55.5 mol   ·   Cal 0.05 mol   ·   CO2@ 0.01 mol   ·   25 °C, 1 bar
species: H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 CaOH+ Ca(CO3)@ Ca(HCO3)+ Cal
```

## What agrees

**Major species** — everything present above about `10⁻⁵` mol — agrees to
`10⁻³` relative or better.

**Sensitivities.** `∂n/∂(CO₂)` obtained by differentiating *through* the solve
(see [Differentiating an equilibrium](@ref)) matches Reaktoro's finite
differences on every species:

| species | ChemistryLab (AD) | Reaktoro (FD) | rel. diff |
|:--|--:|--:|--:|
| H₂O | −0.181336 | −0.181397 | 3.4×10⁻⁴ |
| Ca²⁺ | +0.151920 | +0.151987 | 4.4×10⁻⁴ |
| HCO₃⁻ | +0.333308 | +0.333427 | 3.6×10⁻⁴ |
| CO₂(aq) | +0.818655 | +0.818600 | 6.7×10⁻⁵ |
| Ca(HCO₃)⁺ | +0.029333 | +0.029338 | 1.6×10⁻⁴ |
| calcite | −0.181251 | −0.181324 | 4.0×10⁻⁴ |

Reaktoro's own spread across `h ∈ {10⁻³, 10⁻⁴, 10⁻⁵}` is `7.2×10⁻⁴`, so the
residual difference sits **below the oracle's own truncation error**. This is as
close as a comparison against a finite-difference reference can get.

**Element balance.** Exact, and exactly where a finite difference is not:
calcium has no source in this system, so the sensitivities of all Ca-bearing
species must cancel. The implicit-function route gives `2×10⁻¹⁶`; a species
absent from the solution gets identically zero, against `2×10⁻⁹` by finite
differences.

## What does not agree

!!! danger "Trace species are wrong"
    Species below about `10⁻⁵` mol do **not** agree, and the discrepancy is not
    small:

    | species | ChemistryLab | Reaktoro | ratio |
    |:--|--:|--:|--:|
    | H⁺ | 4.590×10⁻⁷ | 3.693×10⁻⁷ | 1.24 |
    | OH⁻ | 8.623×10⁻⁸ | 2.706×10⁻⁸ | **3.19** |
    | CO₃²⁻ | 1.033×10⁻⁶ | 9.379×10⁻⁷ | 1.10 |
    | Ca(CO₃)@ | 5.645×10⁻⁶ | 5.548×10⁻⁶ | 1.02 |
    | CaOH⁺ | 3.295×10⁻⁸ | 1.586×10⁻⁹ | **20.8** |

    On this system the water autoprotolysis therefore still comes out at
    `pKw = 13.40` instead of `14.00`. **This matters for cement**: a pore
    solution sits at pH 13, and it is exactly these species that set it. Treat
    any pH, pOH or trace-ion figure from a system containing a solid as
    indicative until this is fixed.

    An **aqueous-only** system is a different matter — see below.

### Aqueous-only systems are fixed

Pure water used to come out at `[H⁺]/[OH⁻] = 3.78`, an electroneutrality
violation that no thermodynamic database can excuse. It now gives **1.0**, with
`pKw = 13.9994`.

The cause was the Newton step, not the chemistry. The back-end formed the Schur
complement `S = A H⁻¹ Aᵀ`, which requires `H` invertible — and a pure phase has
unit activity, hence `∂²G/∂nᵢ² = 0` exactly. `OptimaSolver` 0.2.5 computes the
step by the **nullspace** method instead: with `Z` a basis of `null(A)`,

```math
(Z^{\mathsf T} H Z)\, \delta z = -Z^{\mathsf T}(e_x + H\, \delta n_p),
```

in which `H` appears only as a product and is never inverted. This is what the
C++ Optima does by default; its `Rangespace` counterpart — the Schur complement
— is documented there as suitable for invertible diagonal Hessians only.

Mixed solid/aqueous systems are **unchanged** by this, which is why the table
above still stands for calcite + CO₂.

### What remains: conditioning, not thermodynamics

Three observations locate the residual cause:

1. **The standard-state data is right.** Pure water on its own gives
   `pKw = 13.99` and `[H⁺] = [OH⁻]` to six digits.
2. **The element balance is satisfied** to `4×10⁻¹⁶` in every case above. What
   fails is the *stationarity* — the solver returns feasible points that are
   not the Gibbs minimum.
3. **A second, independent solver hits the same wall.** Ipopt does markedly
   better but still does not reach Reaktoro:

   | | CaOH⁺ | OH⁻ | H⁺ | worst |
   |:--|--:|--:|--:|--:|
   | OptimaSolver (linear) | ×20.8 | ×3.19 | ×1.24 | 19.8 |
   | Ipopt (linear) | ×2.46 | ×1.31 | ×1.03 | 1.46 |

   That a mature interior-point code is also off rules out a defect specific to
   one back-end.

The difficulty is scale, and it can be measured rather than supposed. By
Gibbs–Duhem the Hessian of the Gibbs energy is the Jacobian of the chemical
potentials, ``\nabla^2 G = \partial\mu/\partial n``. Evaluated at Reaktoro's
solution for this system:

```
cond = 1.1e22        eigmin = 5.7e-14
H2O@   ∂²G/∂n² = 5.6e-6    (against 1/n = 0.018 — a factor 3200)
Cal    ∂²G/∂n² = 0         exactly, a pure phase has unit activity
CaOH+  ∂²G/∂n² = 6.3e8
```

A condition number of `10²²` is beyond what double precision (`10¹⁶`) carries.
That is the whole story: no interior-point method can be relied on to solve this
in the natural variables, which is why a second, mature solver is also off.

### What has been ruled out

Getting to the nullspace step took eight attempts, seven of which were disproved
by measurement. They are recorded so the next attempt on the *remaining* defect
does not repeat them:

| attempt | result |
|:--|:--|
| more iterations (300 → 200 000) | identical answer; the iteration is stalled, not slow |
| relative finite-difference step for the Hessian | worse (×20 → ×922) |
| exact Hessian `∂μ/∂n` by AD | fixes water, breaks the mixed case (×3750) |
| regularizing the pure-phase zero curvature | no effect at any magnitude tried |
| capping the inverse curvature spread | no effect over five orders of magnitude |
| supplying the analytic gradient `∂G/∂nᵢ = μᵢ` | correct in itself, changes nothing here |
| Ipopt instead of the default back-end | 8× closer, still short of Reaktoro |
| log-space parameterization | does not work in either back-end |
| **nullspace Newton step** | **fixes water, mixed case unchanged** |

The lesson that finally paid: a Hessian affects the *rate* of Newton's method,
not its limit, so a merely slow solver would still land on the right answer. It
did not — which placed the fault in the linear algebra of the step rather than
in the curvature, and that is where it was.

!!! tip "If you need trace species now"
    Use Ipopt, which is roughly eight times closer on this system:

    ```julia
    using Optimization, OptimizationIpopt
    state_eq = equilibrate(state, IpoptOptimizer())
    ```

    It is 3 to 26 times slower than the default on cement equilibria, which is
    why the default is what it is — but on this evidence the default is the
    faster wrong answer for anything below `10⁻⁵` mol.

## The kinetics/equilibrium coupling

**Not yet validated against Reaktoro.** The coupling is checked only against
itself: the partition constraint holds to `‖Aₑnₑ − bₑ‖∞ = 8.7×10⁻⁸` over a
seven-day hydration run with no failed re-speciation, and the resulting
assemblage is physically sensible. That is internal consistency, not agreement
with an independent code.

`test/reference/reaktoro_coupling.py` sets up the intended comparison. It
dissolves calcite at a *constant* rate so the kinetic half is analytic —

```math
n_{\text{Cal}}(t) = n_{\text{Cal}}(0) - k t,
\qquad
b_e(t) = b_e(0) + k t \, (\text{one CaCO}_3)
```

— which removes every difference of rate-law convention between the codes and
leaves only the question that matters: does the equilibrium partition follow
``n_e(t) = \varphi(b_e(t))``? Since the oracle then needs no kinetics integrator
at all, just one equilibrium solve per sample time, it is not blocked on
anything except being run.

## Reproducing

```sh
conda create -n reaktoro-env -c conda-forge reaktoro thermofun
conda run -n reaktoro-env python test/reference/reaktoro_calcite_co2.py
```

Paste the output into `REAKTORO` in `test/equilibrium_reference.jl`. The
assertions that currently fail are marked `@test_broken`, so a fix announces
itself as an unexpected pass rather than sitting unnoticed.
