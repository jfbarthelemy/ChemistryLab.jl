# The hydrating paste, end to end

A worked application of [Coupling kinetics and equilibrium](@ref): alite and
belite dissolve according to [ParrotKilloh1984](@cite), and the hydrate
assemblage that forms is **computed** by Gibbs minimization rather than imposed
from a stoichiometric recipe.

Read the theory page first if you have not — this one assumes the partition
``(b_e, n_k)`` and the map ``n_e = \varphi(b_e)``.

## What is different from a stoichiometric model

The usual approach writes reactions like

```
C₃S + 5.3 H → C₁.₇SH₄ + 1.3 CH
```

and integrates the extent. It is serviceable, but the products and their
proportions are decided in advance. Nothing tells you the pore solution pH, and
nothing adapts if you change the water/cement ratio, add limestone, or leach the
paste.

Here the clinker dissolution rate is prescribed — that is genuinely kinetic —
and everything downstream follows from thermodynamics: which hydrates appear, in
what amounts, and what the pore solution looks like.

Every block below is executed when this page is built, so the numbers are
whatever the code actually produced.

## 1. Species and system

Take the two silicate clinker phases, the hydrates Cemdata18 offers for them,
and the aqueous species `speciation` pulls in from the primaries.

```@example coupled
using ChemistryLab, DynamicQuantities, OptimaSolver, OrdinaryDiffEq, Printf

data = datapath("cemdata18-thermofun.json")
substances = build_species(data)

anhydrous = ["C3S", "C2S"]
hydrates = ["Portlandite", "Jennite"]

sp = speciation(substances, vcat(anhydrous, hydrates); aggregate_state = [AS_AQUEOUS])
cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)

@printf "%d species: %d aqueous, %d crystalline\n" length(cs.species) count(
    s -> aggregate_state(s) == AS_AQUEOUS, cs.species
) count(s -> aggregate_state(s) == AS_CRYSTAL, cs.species)
```

## 2. The initial state

A paste at w/c = 0.40: one kilogram of clinker, 65 % alite and 11 % belite by
mass, and 400 g of water.

```@example coupled
state0 = ChemicalState(cs)
set_quantity!(state0, "C3S", 0.65u"kg")
set_quantity!(state0, "C2S", 0.11u"kg")
set_quantity!(state0, "H2O@", 0.40u"kg")

for s in anhydrous
    @printf "%-5s %8.4f mol\n" s ustrip(us"mol", moles(state0, s))
end
```

## 3. Dissolution reactions

The clinker does not react *into hydrates* here — it **dissolves into ions**,
and ChemistryLab balances each reaction from the phase and the primaries. The
acid-driven form is what drives the pore solution alkaline.

```@example coupled
primaries = ["Ca+2", "SiO2@", "H2O@", "H+"]
reactions = ChemistryLab.Reaction[]
for (phase, pk) in (("C3S", PK84_PARAMS_C3S), ("C2S", PK84_PARAMS_C2S))
    rxn = Reaction([cs[phase]], [cs[p] for p in primaries]; symbol = "$phase dissolution")
    rxn[:rate] = parrot_killoh_avrami(pk, phase; α_max = 1.0, blaine = 380.0u"m^2/kg")
    push!(reactions, rxn)
    println("  ", rxn)
end
```

The negative H⁺ coefficient is the whole point: dissolving one mole of alite
consumes six moles of protons, which is what a pore solution at pH 12.5 and
above records.

## 4. The coupled problem

The reaction list is a **positional** argument, and `equilibrium_solver` belongs
on the problem. Without it the run is a pure kinetics integration and the
aqueous phase never re-speciates; with it, `respeciate!` solves ``\varphi(b_e)``
once per accepted step.

The activity model matters. A cement pore solution sits at an ionic strength of
0.1–0.7 mol/kg, where the dilute model is not defensible, so pass
`HKFActivityModel` — to the problem *and* to the solver, which must agree.

```@example coupled
model = HKFActivityModel()
kp = KineticsProblem(
    cs, reactions, state0, (0.0, 7 * 86400.0);
    activity_model = model,
    equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()),
)
sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-7, abstol = 1.0e-10))
@printf "%d accepted steps, retcode = %s\n" length(sol.t) sol.retcode
```

## 5. Reading the result

The composition at a given time is **not** `state_at`: that returns the purely
kinetic reconstruction ``n(0) + \nu^{\mathsf T}\xi``, and the redistribution
performed by the equilibrium solve cannot be recovered from the stoichiometry.
Use [`speciated_states`](@ref), which re-solves ``\varphi(b_e)`` from the element
totals the run carried, walking the instants in order.

```@example coupled
times = [1.0, 3.0, 7.0] .* 86400
states = speciated_states(sol, kp; times = times)

@printf "%6s %8s %10s %12s %10s\n" "t [d]" "pH" "Jennite" "Portlandite" "H2O"
for (t, st) in zip(times, states)
    @printf "%6.1f %8.2f %10.4f %12.4f %10.3f\n" t / 86400 something(pH(st), NaN) (
        ustrip(us"mol", moles(st, "Jennite"))
    ) ustrip(us"mol", moles(st, "Portlandite")) ustrip(us"mol", moles(st, "H2O@"))
end
```

Check the partition constraint — conservation of matter is one of the three
conditions a certificate rests on, and the one with an immediate meaning:

```@example coupled
p = sol.prob.p
be = collect(sol(times[end])[1:(p.n_be)])
ne = Float64[ustrip(us"mol", moles(states[end], symbol(cs.species[j]))) for j in kp.idx_equilibrium]
@printf "‖Aₑnₑ − bₑ‖∞ = %.2e mol\n" maximum(abs, p.Ae * ne .- be)
```

Degrees of hydration come from the kinetic amounts against their initial values:

```@example coupled
α = degrees_of_hydration(sol, kp; times = times)
for s in anhydrous
    @printf "α(%s) = %.3f\n" s α[s][end]
end
```

## 6. What the thermodynamics decided

Nothing above states that portlandite forms, or in what proportion to the
C-S-H. The dissolution reactions produce only Ca²⁺, SiO₂ and H⁺; the assemblage
is the minimizer of the Gibbs energy under the element totals the kinetics
delivered.

```@example coupled
st = states[end]
solids = [(symbol(s), ustrip(us"mol", moles(st, symbol(s)))) for s in cs.species
    if aggregate_state(s) == AS_CRYSTAL]
filter!(x -> last(x) > 1.0e-6, solids)
sort!(solids; by = last, rev = true)
for (name, amount) in solids
    @printf "  %-14s %8.4f mol\n" name amount
end
```

The pore solution comes with it:

```@example coupled
aq = [(symbol(s), ustrip(us"mol", moles(st, symbol(s)))) for s in cs.species
    if aggregate_state(s) == AS_AQUEOUS && symbol(s) != "H2O@"]
sort!(aq; by = last, rev = true)
for (name, amount) in aq[1:min(5, end)]
    @printf "  %-14s %.3e mol\n" name amount
end
```

## Caveats worth carrying

  - **The C-S-H is `Jennite`**, the Cemdata18 end-member normalized per silicon.
    A real C-S-H is a solid solution of varying Ca/Si; using the CSHQ solid
    solution inside a coupled run is not exercised by this package's tests.
  - **The interior-point optimizer rarely reports convergence** on a cement
    equilibrium, and its return code is not the thing to read. `integrate` reports
    the worst element balance instead, separating the accepted steps — the
    trajectory — from the Jacobian probes that never enter it. The compositions
    above do not rest on that solver: [`speciated_states`](@ref) passes each
    instant to [`DualEquilibriumSolver`](@ref), which solves the KKT system and
    **certifies** the result. The problem is convex, so stationarity of the
    interior species, the component balance, and undersaturation of every absent
    phase together prove global optimality.
  - **A Parrot–Killoh rate reads only its own degree of reaction**, so the
    trajectory here does not depend on the speciation at all; the speciation is
    what you read out of it. A rate law reading log-activities would feed the
    speciation back into the trajectory, and would need the balance above to be
    tight at every step, not just at the instants you ask for.

## 7. The certificate

Nothing above asks you to take the composition on trust. `G` is convex — an ideal
mixing entropy plus terms linear in the amounts of the pure phases, over a
polyhedron — so its minimizer is unique and the KKT conditions are sufficient.
[`optimality_certificate`](@ref) checks them.

```@example coupled
sub = ChemistryLab._equilibrium_subsystem(kp.system, kp.idx_equilibrium)
des = DualEquilibriumSolver(sub, model)
snames = symbol.(sub.species)
be = collect(sol(times[end])[1:(sol.prob.p.n_be)])
sub_state = ChemicalState(
    sub, [ustrip(us"mol", moles(states[end], s)) for s in snames] .* u"mol";
    T = sol.prob.p.T_q[], P = sol.prob.p.P_q[],
)
cert = optimality_certificate(des, sub_state; b = be)
@printf "optimal          %s\n" cert.optimal
@printf "stationarity     %.2e  (RT, over %d interior species)\n" cert.stationarity cert.n_interior
@printf "component balance %.2e mol\n" cert.balance
@printf "worst supersaturation %+.2e  (negative: every absent phase undersaturated)\n" cert.worst_supersaturation
```
