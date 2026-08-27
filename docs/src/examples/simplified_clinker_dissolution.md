# [Simplified Clinker Dissolution](@id sec-clinker-dissolution)

This example illustrates ChemistryLab's equilibrium workflow applied to cement chemistry: computing the hydration state of a simplified Portland clinker at a given water-to-cement ratio. ChemistryLab computes thermodynamic equilibrium by minimizing the Gibbs free energy of the system subject to element-conservation constraints. The workflow always follows the same four steps:

1. Build a [`ChemicalSystem`](@ref) (species + stoichiometric matrix).
2. Create an initial [`ChemicalState`](@ref) (temperature, pressure, initial amounts).
3. Call [`equilibrate`](@ref) (or use [`EquilibriumSolver`](@ref) explicitly).
4. Inspect the resulting [`ChemicalState`](@ref).

---

## Minimal workflow

The convenience function [`equilibrate`](@ref) handles everything with sensible defaults.
The example below reproduces a simplified clinker dissolution calculation.

```@setup eq_setup
using Optimization, OptimizationIpopt
using ChemistryLab
using DynamicQuantities

substances = build_species(datapath("cemdata18-thermofun.json"))
input_species = split("C3S C2S C3A C4AF Gp Anh Portlandite Jennite H2O@ ettringite monosulphate12 C3AH6 C3FH6 C4FH13")
species = speciation(substances, input_species; aggregate_state=[AS_AQUEOUS])

cs = ChemicalSystem(species, CEMDATA_PRIMARIES)
```

```julia
using ChemistryLab
using DynamicQuantities

substances = build_species(datapath("cemdata18-thermofun.json"))
input_species = split("C3S C2S C3A C4AF Gp Anh Portlandite Jennite H2O@ ettringite monosulphate12 C3AH6 C3FH6 C4FH13")
species = speciation(substances, input_species; aggregate_state=[AS_AQUEOUS])

cs = ChemicalSystem(species, CEMDATA_PRIMARIES)
```

```@example eq_setup
state = ChemicalState(cs)

# Clinker + gypsum composition (mass fractions, total = 1)
compo = ["C3S" => 0.678, "C2S" => 0.166, "C3A" => 0.04, "C4AF" => 0.072, "Gp" => 0.028]
c     = sum(last.(compo))
wc    = 0.4          # water-to-cement ratio
w     = wc * c
mtot  = c + w

for x in compo
    set_quantity!(state, x.first, x.second / mtot * u"kg")
end
set_quantity!(state, "H2O@", w / mtot * u"kg")

# pH-neutral seed for H⁺ and OH⁻
V = volume(state)
set_quantity!(state, "H+",  1e-7u"mol/L" * V.liquid)
set_quantity!(state, "OH-", 1e-7u"mol/L" * V.liquid)

state_eq = equilibrate(state)
```

!!! tip "Quick shortcut"
    If your `ChemicalState` was constructed from a `ChemicalSystem` built with `CEMDATA_PRIMARIES`, calling `equilibrate(state)` with no other arguments is usually sufficient for cement-chemistry problems.

---

## Inspecting the equilibrium state

The returned [`ChemicalState`](@ref) carries all derived thermodynamic quantities:

```julia
println("pH      = ", pH(state_eq))
println("pOH     = ", pOH(state_eq))
println("porosity   = ", porosity(state_eq))
println("saturation = ", saturation(state_eq))
```

Phase volumes and mole amounts are accessible via named tuples:

```julia
v = volume(state_eq)
println("V liquid = ", v.liquid)
println("V solid  = ", v.solid)
println("V total  = ", v.total)

m = moles(state_eq)
println("n liquid = ", m.liquid)
println("n solid  = ", m.solid)
```

Individual species amounts (in mol):

```julia
cs_eq = state_eq.system
for (i, sp) in enumerate(cs_eq.species)
    n_i = state_eq.n[i]
    println(rpad(symbol(sp), 20), ustrip(n_i), " mol")
end
```

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

```julia
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

## Activity models

All activity models inherit from [`AbstractActivityModel`](@ref). The only built-in model is [`DiluteSolutionModel`](@ref), which implements:

| Phase | Law | Expression |
|:------|:----|:----------|
| Solvent (H₂O) | Raoult | `ln a = ln xₛ` |
| Aqueous solutes | Henry | `ln a = ln(cᵢ / c°)`, `c° = 1 mol/L` |
| Crystals | Pure solid | `ln a = 0` |
| Gas | Ideal mixture | `ln a = ln xᵢ` |

To implement a custom activity model, define a new subtype and extend `activity_model`:

```julia
struct MyModel <: AbstractActivityModel
    # model parameters
end

function ChemistryLab.activity_model(::ChemicalSystem, ::MyModel)
    # return a closure lna(n, p) -> Vector{Float64}
    return (_, _) -> begin
        # compute and return log-activities
    end
end
```

Pass your model to `equilibrate` or `EquilibriumSolver`:

```julia
state_eq = equilibrate(state; model=MyModel(...))
```
