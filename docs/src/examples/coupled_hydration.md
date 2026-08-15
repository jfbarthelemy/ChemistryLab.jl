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

## 1. Species and system

Take the clinker phases, the hydrates Cemdata18 offers, and enough aqueous
species to carry the chemistry:

```julia
using ChemistryLab, DynamicQuantities, OptimaSolver, OrdinaryDiffEq

data = joinpath(pkgdir(ChemistryLab), "data", "cemdata18-thermofun.json")
dict = Dict(symbol(s) => s for s in build_species(data))

kinetic  = ["C3S", "C2S"]                       # dissolve over days
hydrates = ["Portlandite", "Jennite", "Tobermorite-14"]
aqueous  = ["H2O@", "H+", "OH-", "Ca+2", "CaOH+", "SiO2@", "HSiO3-", "Ca(HSiO3)+"]

cs = ChemicalSystem(
    [dict[s] for s in vcat(aqueous, hydrates, kinetic)],
    ["H2O@", "SiO2@", "Ca+2", "H+"],            # primaries
)
```

The primaries matter: the equilibrium sub-system inherits them, so `Aₑ` and the
minimization are posed on the same basis (see the warning on the theory page).

## 2. Initial state — a paste at w/c = 0.40

```julia
n = Any[fill(0.0u"mol", length(cs.species))...]
idx(s) = findfirst(sp -> symbol(sp) == s, cs.species)

n[idx("H2O@")] = 22.20u"mol"     # 400 g water
n[idx("C3S")]  = 3.29u"mol"      # 750 g alite
n[idx("C2S")]  = 0.87u"mol"      # 150 g belite

state = ChemicalState(cs, n; T = 298.15u"K", P = 1u"bar")
```

## 3. Rate laws

`parrot_killoh` returns a `KineticFunc` with the standard `(T, P, t, n, lna,
n₀)` signature, so it drops straight into a `KineticReaction`:

```julia
α_max = 1.0 - 3.33 * (0.40 - 0.40)          # Powers limit at this w/c

reactions = [
    kinetic_species(cs, "C3S", parrot_killoh(PK_PARAMS_C3S, "C3S"; α_max = α_max)),
    kinetic_species(cs, "C2S", parrot_killoh(PK_PARAMS_C2S, "C2S"; α_max = α_max)),
]
```

!!! note "`kinetic_species` normalizes the stoichiometry for you"
    The reaction it builds from the nullspace can come out with either sign on
    the controlling mineral. It is normalized to `ν = −1` there, so the ODE
    dissolves the clinker rather than growing it, and the `1/|νₖ|` scaling the
    rate law expects is applied.

## 4. Declare the partition and solve

Everything not declared kinetic is equilibrium. Attaching an equilibrium solver
is what turns a plain kinetics run into the coupled problem:

```julia
kp = KineticsProblem(
    cs, reactions, state, (0.0u"s", 7.0u"d");
    equilibrium_solver = EquilibriumSolver(cs, DiluteSolutionModel(), OptimaOptimizer()),
)

ks = KineticsSolver(; ode_solver = Rodas5P(),               # stiff: rates span decades
                    reltol = 1e-6, abstol = 1e-10)

sol = integrate(kp, ks)
```

The reaction list is a **positional** argument, and `equilibrium_solver` belongs
on the problem (a `KineticsSolver` accepts it too, and the problem wins if both
are given).

Without `equilibrium_solver` the run is a pure kinetics integration and the
aqueous phase never re-speciates. With it, `respeciate!` solves ``\varphi(b_e)``
once per accepted step.

## 5. Reading the result

`integrate` reports failed re-speciations at the end — **check that it says
nothing**, and check the partition constraint:

```julia
using LinearAlgebra
be, ne = sol.u[end][1:kp.n_be], final_speciation(sol, kp)
@show maximum(abs, kp.Ae * ne - be)          # → 8.7e-8 mol
```

Degrees of hydration come from the kinetic amounts against their initial values:

```julia
α(s) = 1 - amount(sol, s, 7.0u"d") / amount(sol, s, 0.0u"s")
α("C3S"), α("C2S")                            # → 0.277, 0.279
```

### Measured output, seven days at w/c = 0.40

| quantity | value |
|:--|--:|
| α(C₃S), α(C₂S) | 0.277, 0.279 |
| Jennite | 1.019 mol |
| Portlandite | 1.086 mol |
| H₂O remaining | 18.97 mol (from 22.20) |
| Ca²⁺ | 3.5 mmol |
| OH⁻ | 8.5 mmol |
| `‖Aₑnₑ − bₑ‖∞` | 8.7×10⁻⁸ mol |
| failed re-speciations | 0 of 89 steps |

C-S-H and portlandite in comparable amounts, water consumed, an alkaline pore
solution at millimolar calcium. None of it was imposed: change the water content
or add a species and the assemblage rearranges itself.

!!! note "How far to trust the pore-solution figures"
    The speciation now agrees with Reaktoro to 5 % or better on every species
    except `CaOH⁺`, which is 2.5× high at `4×10⁻⁹` mol, and `pKw` comes out at
    13.979. The `OH⁻` figure above and any pH derived from this run are usable;
    see [Validation against Reaktoro](@ref) for what was checked and how.

## 6. Differentiating the whole thing

Because the sensitivities come from the optimality conditions rather than from
solver iterations, `ForwardDiff` crosses the equilibrium solve. So the
composition — and quantities derived from it — can be differentiated with
respect to a formulation parameter:

```julia
using ForwardDiff

function portlandite_at_7d(wc)
    n = Any[fill(0.0u"mol", length(cs.species))...]
    n[idx("H2O@")] = (wc / 0.40 * 22.20)u"mol"
    n[idx("C3S")]  = 3.29u"mol"
    n[idx("C2S")]  = 0.87u"mol"
    sol = integrate(KineticsProblem(cs, reactions, ChemicalState(cs, n),
                                    (0.0u"s", 7.0u"d")), ks)
    return amount(sol, "Portlandite", 7.0u"d")
end

ForwardDiff.derivative(portlandite_at_7d, 0.40)
```

This is what makes the package usable inside an optimization or an
identification loop, rather than only as a forward simulator.

## See also

- [Coupling kinetics and equilibrium](@ref) — the theory behind this example
- [Validation against Reaktoro](@ref) — what is verified and what is not
- [Chemical Kinetics](@ref) — the rate-law catalog
