# [ChemicalSystem and ChemicalState](@id sec-system-state)

These two types are the bridge between species definitions and equilibrium calculations:

- **[`ChemicalSystem`](@ref)** — immutable container that groups species, stoichiometric matrices, and index maps. Built once; shared across many states.
- **[`ChemicalState`](@ref)** — mutable snapshot of a system at a given temperature, pressure, and molar composition. Built from a `ChemicalSystem`; mutated in place.

---

## ChemicalSystem type

### What it holds

```
ChemicalSystem
├── species           → AbstractVector{<:AbstractSpecies}    (ordered list)
├── CSM               → CanonicalStoichMatrix                (elements × species)
├── SM                → StoichMatrix                         (primaries × species)
├── reactions         → Vector{<:AbstractReaction}           (derived or provided)
├── idx_aqueous       → Vector{Int}  (aqueous species)
├── idx_crystal       → Vector{Int}  (crystalline species)
├── idx_gas           → Vector{Int}  (gas-phase species)
├── idx_solvent       → Vector{Int}  (solvent, i.e. H₂O@)
├── idx_solutes       → Vector{Int}  (SC_AQSOLUTE)
├── idx_components    → Vector{Int}  (SC_COMPONENT)
└── solid_solutions   → Nothing | Vector{SolidSolutionPhase}
```

All derived fields (stoichiometric matrices, index maps) are **computed once at construction** and remain consistent for the lifetime of the object.

### Minimal construction

```@example cs_minimal
using ChemistryLab

# Three aqueous species — no primaries specified → all species are primaries
H2O  = Species("H2O";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
Hp   = Species("H+";   aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
OHm  = Species("OH-";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)

cs = ChemicalSystem([H2O, Hp, OHm])
length(cs)
```

```@example cs_minimal
cs.idx_solvent      # index of the solvent
```

```@example cs_minimal
cs.idx_solutes      # indices of solutes
```

```@example cs_minimal
pprint(cs.CSM; label = :symbol)
```

### Specifying primary species

Primary species (independent components) determine which species are "dependent" and how the stoichiometric matrix `SM` is built. Balanced reactions are then extracted from the null space of `SM`.

```@example cs_primaries
using ChemistryLab

H2O  = Species("H2O";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
Hp   = Species("H+";   aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
OHm  = Species("OH-";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
CO2  = Species("CO2";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
HCO3 = Species("HCO3-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)

# Primaries: H₂O, H⁺, HCO₃⁻  (a chemist's choice for the carbonate system)
primaries = [H2O, Hp, HCO3]

cs = ChemicalSystem([H2O, Hp, OHm, CO2, HCO3], primaries)
pprint(cs.SM; label = :symbol)
```

```@example cs_primaries
# Independent reactions derived from the null space
rxns = reactions(cs.SM)
for r in rxns
    println(r.equation)
end
```

### Filtered views

`ChemicalSystem` is an `AbstractVector{<:AbstractSpecies}`. Filtered views return sub-vectors without copying data:

```@example cs_views
using ChemistryLab

H2O   = Species("H2O";   aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
Hp    = Species("H+";    aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
Cal   = Species("Cal";   aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
CO2g  = Species("CO2";   aggregate_state = AS_GAS,     class = SC_GASFLUID)

cs = ChemicalSystem([H2O, Hp, Cal, CO2g])

println("aqueous:  ", symbol.(aqueous(cs)))
println("crystal:  ", symbol.(crystal(cs)))
println("gas:      ", symbol.(gas(cs)))
println("solvent:  ", symbol(solvent(cs)))
println("solutes:  ", symbol.(solutes(cs)))
```

### Species lookup

```@example cs_lookup
using ChemistryLab

H2O  = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
Hp   = Species("H+";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
cs   = ChemicalSystem([H2O, Hp])

cs["H2O"]       # lookup by symbol string
```

```@example cs_lookup
cs[2]            # lookup by index
```

```@example cs_lookup
# Find the index of a species
findfirst(s -> symbol(s) == "H+", cs.species)
```

### Getting reactions

When the system was built from a database (with thermodynamic data attached to species), each dependent species has a corresponding dissolution/formation reaction:

```julia
# From a database workflow
cs = ChemicalSystem(species, primaries)
r_cal = get_reaction(cs, "Cal")   # reaction for calcite
r_cal.logK⁰(T = 298.15)
```

---

## ChemicalState type

### Construction

A `ChemicalState` is created from a `ChemicalSystem`. All mole amounts start at zero; temperature defaults to 298.15 K and pressure to 10⁵ Pa.

```@example cst_basic
using ChemistryLab
using DynamicQuantities

H2O  = Species("H2O";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
Hp   = Species("H+";   aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
OHm  = Species("OH-";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
cs   = ChemicalSystem([H2O, Hp, OHm], [H2O, Hp])

state = ChemicalState(cs)
ustrip(state.T[])      # temperature in K
```

```@example cst_basic
ustrip(state.P[])      # pressure in Pa
```

Override temperature and pressure at construction:

```@example cst_basic
state2 = ChemicalState(cs; T = 350.0u"K", P = 2e5u"Pa")
ustrip(state2.T[])
```

### Setting molar amounts

`set_quantity!` accepts any compatible unit — moles, kilograms, grams, or concentration (mol/L) multiplied by a volume:

```@example cst_setq
using ChemistryLab
using DynamicQuantities

H2O  = Species("H2O";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
Hp   = Species("H+";   aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
OHm  = Species("OH-";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
Cal  = Species("Cal";  aggregate_state = AS_CRYSTAL,  class = SC_COMPONENT)
cs   = ChemicalSystem([H2O, Hp, OHm, Cal], [H2O, Hp])

state = ChemicalState(cs)

# 1 litre of water (by mass)
set_quantity!(state, "H2O", 1.0u"kg")

# pH 4 water: 10⁻⁴ mol/L × volume of liquid phase
V = volume(state)
set_quantity!(state, "H+",  1e-4u"mol/L" * V.liquid)
set_quantity!(state, "OH-", 1e-10u"mol/L" * V.liquid)

# 1 mmol of calcite
set_quantity!(state, "Cal", 1e-3u"mol")

# Inspect moles
ustrip.(state.n)
```

### Changing temperature and pressure

```@example cst_setq
set_temperature!(state, 350.0u"K")
set_pressure!(state,    2e5u"Pa")
ustrip(state.T[])
```

### Derived quantities

All derived quantities are **recomputed automatically** after any `set_*!` call:

```@example cst_derived
using ChemistryLab
using DynamicQuantities

H2O  = Species("H2O";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
Hp   = Species("H+";   aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
OHm  = Species("OH-";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
cs   = ChemicalSystem([H2O, Hp, OHm], [H2O, Hp])

state = ChemicalState(cs)
set_quantity!(state, "H2O", 1.0u"kg")
# H⁺ and OH⁻ are auto-seeded at neutral pH — no manual seeding needed

pH(state)
```

```@example cst_derived
pOH(state)
```

!!! info "Automatic H⁺/OH⁻ seeding"
    When water (the aqueous solvent, `SC_AQSOLVENT`) is added to a `ChemicalState`
    that contains H⁺ and OH⁻ species, ChemistryLab automatically seeds them at
    the neutral pH concentration ``c = 10^{-pK_w/2}`` based on the water
    autoprotolysis constant ``K_w(T, P)``. This ensures a chemically consistent
    initial state without manual intervention.

    The auto-seeding only triggers when **both** H⁺ and OH⁻ are at zero — if you
    set them explicitly (e.g. to impose a specific initial pH), your values are
    preserved.

```@example cst_derived
v = volume(state)
println("V liquid = ", v.liquid)
println("V solid  = ", v.solid)
println("V total  = ", v.total)
```

```@example cst_derived
m = moles(state)
println("n liquid = ", m.liquid)
println("n solid  = ", m.solid)
```

```@example cst_derived
porosity(state)
```

### Copying a state

`copy` creates a new `ChemicalState` that **shares** the underlying `ChemicalSystem` (no duplication) but has its own independent mole vector, temperature, and pressure:

```@example cst_copy
using ChemistryLab
using DynamicQuantities

H2O = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
Hp  = Species("H+";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
cs  = ChemicalSystem([H2O, Hp])

s1 = ChemicalState(cs)
set_quantity!(s1, "H2O", 1.0u"kg")

s2 = copy(s1)
set_temperature!(s2, 350.0u"K")   # only s2 is modified

ustrip(s1.T[]), ustrip(s2.T[])
```

---

## Full workflow example

Building a minimal carbonate system from scratch — no database required:

```@example full_example
using ChemistryLab
using DynamicQuantities

# 1. Declare species
H2O  = Species("H2O";   aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
Hp   = Species("H+";    aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
OHm  = Species("OH-";   aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
CO2  = Species("CO2";   aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
HCO3 = Species("HCO3-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
CO3  = Species("CO3-2";  aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)

# 2. Build the system (H₂O, H⁺, CO₃²⁻ as primaries)
cs = ChemicalSystem([H2O, Hp, OHm, CO2, HCO3, CO3], [H2O, Hp, CO3])

println("Species: ", join(symbol.(cs.species), ", "))
println("Aqueous: ", join(symbol.(aqueous(cs)), ", "))
```

```@example full_example
# 3. Build an initial state
state = ChemicalState(cs)
set_quantity!(state, "H2O",   1.0u"kg")
set_quantity!(state, "H+",    1e-7u"mol/L" * volume(state).liquid)
set_quantity!(state, "OH-",   1e-7u"mol/L" * volume(state).liquid)
set_quantity!(state, "CO2",   1e-3u"mol")

println("pH       = ", pH(state))
println("n liquid = ", moles(state).liquid)
```

!!! tip "Next step: equilibrium"
    Once you have a `ChemicalSystem` and a `ChemicalState`, pass the state to [`equilibrate`](@ref) to find the thermodynamic equilibrium. See the [Equilibrium](@ref sec-equilibrium) page for details.
