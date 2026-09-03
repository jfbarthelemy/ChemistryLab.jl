# Coupling kinetics and equilibrium

Some reactions are fast enough to be treated as instantaneous, others are not.
A hydrating cement paste has both: aqueous speciation — protonation,
complexation, water autoprotolysis — reaches equilibrium in microseconds, while
alite dissolves over days. Integrating everything with rate laws would demand
kinetic parameters nobody measures for the fast reactions and would force the
integrator down to their timescale. Equilibrating everything would dissolve the
clinker instantly.

The way out is [Leal2017](@cite)'s partition, which is also what Reaktoro
implements. This page derives it; [The hydrating paste, end to end](@ref sec-coupled-hydration) puts
it to work.

## The partition

Split the species into two sets:

- the **kinetic partition** ``n_k`` — phases whose transformation is rate-limited
  (clinker, most minerals);
- the **equilibrium partition** ``n_e`` — everything assumed instantaneously
  equilibrated (the aqueous species, and any phase you accept as equilibrated).

Each kinetic reaction ``j`` has a stoichiometric row that touches both, so the
full stoichiometric matrix splits column-wise as ``\nu = [\nu_e \;\; \nu_k]``.

## Why the state is ``(b_e, n_k)`` and not ``(n_e, n_k)``

The obvious choice — integrate every species — does not work, and the reason is
worth stating because it explains the whole construction.

Take ``\dot n_e = \nu_e^{\mathsf T} r``: this advances the equilibrium species
along the kinetic reactions *without re-equilibrating them*. After one step the
aqueous phase is no longer at equilibrium, so the assumption that justified the
partition has been abandoned.

Re-equilibrating after each step does not rescue it either, because the natural
constraint for the equilibrium problem is not a composition — it is a set of
**element totals**. Along the way an individual species genuinely wants to go
negative: the dissolution reactions generated from the nullspace are written in
``\mathrm{H^+}``, and a cement paste contains no acid. What has physical meaning
is not "how much ``\mathrm{H^+}`` was consumed" but "how much hydrogen is in the
equilibrium partition", and it is the minimizer, not the caller, that
redistributes those elements over a feasible composition.

So the state carries the element amounts of the equilibrium partition,

```math
b_e = A_e \, n_e ,
```

where ``A_e`` is the conservation matrix restricted to that partition. These are
conserved by every fast reaction by construction, and only the kinetic reactions
move them.

!!! warning "`Aₑ` must be the matrix the solve is posed on"
    The formula matrix over the *canonical elements* and the matrix over the
    system's *primary species* are different matrices. Building `Aₑ` from one
    while the minimization is posed on the other makes every equilibrium solve
    infeasible — in one intermediate version of this package, 89 failures out of
    89 steps. `Aₑ` is taken from the equilibrium sub-system, which inherits the
    parent system's primaries.

## The system

With that state, the coupled problem is ([Leal2017](@cite), Eqs. 54–65):

```math
\frac{\mathrm{d} n_k}{\mathrm{d} t} = \nu_k^{\mathsf T} r(n, T, t),
\qquad
\frac{\mathrm{d} b_e}{\mathrm{d} t} = A_e \, \nu_e^{\mathsf T} r(n, T, t),
```

closed by the equilibrium map

```math
n_e = \varphi(b_e) \;=\; \arg\min_{n} \; G(n)
\quad \text{subject to} \quad A_e n = b_e, \;\; n \ge 0 .
```

Two ODEs and one constrained minimization. The rates ``r`` depend on the full
composition — a dissolution rate needs the saturation index, hence the aqueous
activities — so ``\varphi`` feeds back into the right-hand side.

## Solving it: operator splitting

``\varphi`` is not evaluated inside the ODE right-hand side. It is applied once
per **accepted** step, as a `DiscreteCallback`.

This is a deliberate choice, and both reasons matter:

- a stiff integrator evaluates the right-hand side many times per step and
  differentiates it to build a Jacobian. An optimization solve in there makes
  the cost unpredictable, and the Jacobian would be taken through an
  active-set-dependent map;
- with the solve outside, the residual and its Jacobian see the *same*
  speciation, so the Jacobian is consistent with the model actually being
  integrated.

The initial state is equilibrated before the first step, so the trajectory
starts on the constraint manifold rather than drifting onto it.

```
    u = (bₑ, nₖ)  ──▶  ODE step  ──▶  accepted?  ──▶  respeciate!  ──▶  u
                          ▲                              │
                          └───── frozen nₑ ◀─────────────┘
```

## What to check when it runs

Two diagnostics tell you whether the coupling is doing its job:

| check | meaning | healthy value |
|:--|:--|:--|
| `‖Aₑ nₑ − bₑ‖∞` | the speciation carries the integrated element totals | `~10⁻⁸` mol |
| re-speciation failure count | steps that kept a frozen composition | 0 |

The second is reported automatically at the end of `integrate`. It exists
because an earlier version of this package swallowed solve failures in a bare
`catch`: a run in which re-speciation **never once succeeded** looked exactly
like a healthy one. Silence is not success — if the count is nonzero, the
trajectory is not the one the equations describe.

!!! note "Validation status"
    Both halves are checked against Reaktoro: the **equilibrium** species by
    species, and the **coupling** along a constant-rate dissolution trajectory,
    where the two agree to 4.3 % or better — see
    [Validation against Reaktoro](@ref).

## See also

- [The hydrating paste, end to end](@ref sec-coupled-hydration) — a worked example
- [Validation against Reaktoro](@ref) — what agrees and what does not
- [Chemical Kinetics](@ref sec-kinetics) — rate laws and the `KineticFunc` interface
