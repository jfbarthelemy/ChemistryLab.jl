# Chemical Reactions

In ChemistryLab it is possible to build chemical reactions and manipulate them. A reaction is constructed as a structure, "a composite data type that allows you to store multiple values in a single object". The `struct` is organized as follows:

```julia
struct Reaction{SR<:AbstractSpecies, TR<:Number, SP<:AbstractSpecies, TP<:Number}
    equation::String
    colored::String
    reactants::OrderedDict{SR, TR}
    products::OrderedDict{SP, TP}
    equal_sign::Char
    properties::OrderedDict{Symbol,PropertyType}
end
```

## Parsing reactions

`Reaction` is a composite type `struct` which can be build from:

- a string containing [species](./databases.md#species)

```julia
equation = "13H⁺ + NO₃⁻ + CO₃²⁻ + 10e⁻ = 6H₂O@ + HCN@"
reac, prod, equal_sign = parse_equation(equation)
```

- a string containing [cement species](./databases.md#species)

```julia
eqC3S = "C₃S + 5.3H = 1.3CH + C₁.₇SH₄"
rC3S = CemReaction(eqC3S)
```

- an operation on species

```@example
using ChemistryLab
C3S = CemSpecies("C3S", symbol="C₃S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
H = CemSpecies("H", symbol="H₂O@", aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)
CH = CemSpecies("CH", symbol="C₃S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
CSH = CemSpecies("C1.7SH4", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
r = C3S + 5.3H ↔ 1.3CH + CSH
typeof(r)
```

- a balance calculation

```@example
using ChemistryLab
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
CSH = CemSpecies("C1.7SH4")
r = Reaction([C3S, H, CH, CSH]; equal_sign='→')
pprint(r)
```

- a balance calculation with symbolic numbers

```@example
using ChemistryLab
using Symbolics
@variables a b g
CSH = CemSpecies(Dict(:C => a, :S => one(Num), :H => g))
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
r = Reaction([CSH, C3S, H, CH]; equal_sign='→')
```

```@example
using ChemistryLab
using Symbolics
using PrettyTables
@variables a b g
CSH = CemSpecies(Dict(:C => a, :S => one(Num), :H => g))
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
r = map(simplify, Reaction([C3S, H], [CH, CSH]; equal_sign='→'))
SM = StoichMatrix([C3S], [CSH, H, CH])
pprint(SM)
```

!!! note "Collection of arrow symbols"
    In ChemistryLab, there are collections of arrow symbols used in chemical reaction notation, such as:
    - `>`, `→`, `↣`, `↦`, `⇾`, `⟶`, `⟼`, `⥟`, `⥟`, `⇀`, `⇁`, `⇒`, `⟾` for reaction directionality from reactants to products;
    - `<`, `←`, `↢`, `↤`, `⇽`, `⟵`, `⟻`, `⥚`, `⥞`, `↼`, `↽`, `⇐`, `⟽` for reaction directionality from products to reactants;
    - `↔`, `⟷`, `⇄`, `⇆`, `⇌`, `⇋`, `⇔`, `⟺` for reversible reactions and equilibrium states;
    - `=`, `≔`, `⩴`, `≕` to separate reactants from products in balanced equations.

## Thermodynamic properties of reactions

When the species involved in a reaction carry thermodynamic data (loaded from a database), the reaction automatically exposes temperature-dependent thermodynamic functions. These are computed lazily on first access and stored in the reaction's `properties` dict:

| Property | Description |
|:---------|:------------|
| `r.ΔᵣCp⁰` | Heat capacity of reaction (J mol⁻¹ K⁻¹) |
| `r.ΔᵣH⁰`  | Enthalpy of reaction (J mol⁻¹) |
| `r.ΔᵣS⁰`  | Entropy of reaction (J mol⁻¹ K⁻¹) |
| `r.ΔᵣG⁰`  | Gibbs free energy of reaction (J mol⁻¹) |
| `r.logK⁰` | Decimal logarithm of the equilibrium constant |

Each property is a `SymbolicFunc` callable with a keyword argument `T` (temperature in K):

```julia
using ChemistryLab

# Load reactions from a database-built stoichiometric matrix
all_species = build_species(datapath("cemdata18-merged.json"))
species = speciation(all_species, split("Cal H2O@");
              aggregate_state=[AS_AQUEOUS], exclude_species=split("H2@ O2@ CH4@"))
dict_species = Dict(symbol(s) => s for s in species)
candidate_primaries = [dict_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_species, s)]
cs = ChemicalSystem(species, candidate_primaries)
list_reactions = reactions(cs.SM)
dict_reactions = Dict(r.symbol => r for r in list_reactions)

r_cal = dict_reactions["Cal"]   # calcite dissolution reaction

# Evaluate at 25 °C
r_cal.logK⁰(T = 298.15)          # log₁₀ K at 25 °C
r_cal.ΔᵣG⁰(T = 298.15)           # ΔᵣG° at 25 °C  (J/mol)
```

Properties can also be set manually on a reaction:

```julia
r = Reaction("H2 + O2 = H2O")
r[:ΔᵣH⁰] = -241800.0   # J/mol
```

## Reactions from a stoichiometric matrix

The most common source of reactions in a database workflow is `reactions(SM)`, which derives all independent reactions from a `StoichMatrix`:

```julia
list_reactions = reactions(cs.SM)
dict_reactions = Dict(r.symbol => r for r in list_reactions)
```

Each reaction symbol matches the dependent species it describes. Temperature-dependent thermodynamic functions are computed automatically when the constituent species carry thermodynamic data.

---

## Algebraic reaction balancing

The null space of the stoichiometric matrix automatically yields **balanced chemical reactions** — no manual coefficient guessing needed. This is the algebraic counterpart to the chemical insight that mass must be conserved.

### Principle

Given a set of species, `ChemicalSystem` builds the stoichiometric matrix **A** (rows = elements, columns = species). Any vector **ν** in the null space of **A** satisfies **Aν = 0**, i.e. it represents a balanced reaction. `reactions(cs.SM)` extracts one independent balanced reaction per dependent species.

### Example: combustion of methane

```@example combustion
using ChemistryLab

# Declare gas-phase species
CH4 = Species("CH4"; aggregate_state = AS_GAS, class = SC_GASFLUID)
O2  = Species("O2";  aggregate_state = AS_GAS, class = SC_GASFLUID)
CO2 = Species("CO2"; aggregate_state = AS_GAS, class = SC_GASFLUID)
H2O = Species("H2O"; aggregate_state = AS_GAS, class = SC_GASFLUID)

# CH4, O2, H2O are the primary (independent) components; CO2 is the derived species
primaries = [CH4, O2, H2O]
cs = ChemicalSystem([CH4, O2, CO2, H2O], primaries)

rxns = reactions(cs.SM)
pprint(rxns[1])
```

!!! note "Reading the matrix"
    The stoichiometric matrix for this system (elements as rows, species as columns) is:

    |       | C | H | O |
    |:------|:--|:--|:--|
    | CH₄   | 1 | 4 | 0 |
    | O₂    | 0 | 0 | 2 |
    | CO₂   | 1 | 0 | 2 |
    | H₂O   | 0 | 2 | 1 |

    With CH₄, O₂, H₂O as primaries, CO₂ is the dependent species.
    `reactions(cs.SM)` expresses CO₂ as: CO₂ = CH₄ + 2 O₂ − 2 H₂O,
    yielding the balanced equation CH₄ + 2 O₂ → CO₂ + 2 H₂O.

### Example: combustion of alkanes (CₙH₂ₙ₊₂)

```@example combustion_alkanes
using ChemistryLab

C2H6 = Species("C2H6"; aggregate_state = AS_GAS, class = SC_GASFLUID)
C3H8 = Species("C3H8"; aggregate_state = AS_GAS, class = SC_GASFLUID)
O2   = Species("O2";   aggregate_state = AS_GAS, class = SC_GASFLUID)
CO2  = Species("CO2";  aggregate_state = AS_GAS, class = SC_GASFLUID)
H2O  = Species("H2O";  aggregate_state = AS_GAS, class = SC_GASFLUID)

primaries = [CO2, H2O, O2]

# Ethane: C₂H₆ + 7/2 O₂ → 2 CO₂ + 3 H₂O
# reactions(SM) returns the reverse here; negate to get the forward combustion direction
cs_ethane = ChemicalSystem([C2H6, O2, CO2, H2O], primaries)
pprint(-reactions(cs_ethane.SM)[1])
```

```@example combustion_alkanes
# Propane: C₃H₈ + 5 O₂ → 3 CO₂ + 4 H₂O
cs_propane = ChemicalSystem([C3H8, O2, CO2, H2O], primaries)
pprint(-reactions(cs_propane.SM)[1])
```

The general formula for saturated alkanes CₙH₂ₙ₊₂ is:
`CₙH₂ₙ₊₂ + (3n+1)/2 O₂ → n CO₂ + (n+1) H₂O`

### Example: combustion of alkenes (CₙH₂ₙ)

```@example combustion_alkenes
using ChemistryLab

C2H4 = Species("C2H4"; aggregate_state = AS_GAS, class = SC_GASFLUID)
C3H6 = Species("C3H6"; aggregate_state = AS_GAS, class = SC_GASFLUID)
O2   = Species("O2";   aggregate_state = AS_GAS, class = SC_GASFLUID)
CO2  = Species("CO2";  aggregate_state = AS_GAS, class = SC_GASFLUID)
H2O  = Species("H2O";  aggregate_state = AS_GAS, class = SC_GASFLUID)

primaries = [CO2, H2O, O2]

# Ethylene: C₂H₄ + 3 O₂ → 2 CO₂ + 2 H₂O
# reactions(SM) returns the reverse here; negate to get the forward combustion direction
cs_ethylene = ChemicalSystem([C2H4, O2, CO2, H2O], primaries)
pprint(-reactions(cs_ethylene.SM)[1])
```

```@example combustion_alkenes
# Propylene: C₃H₆ + 9/2 O₂ → 3 CO₂ + 3 H₂O
cs_propylene = ChemicalSystem([C3H6, O2, CO2, H2O], primaries)
pprint(-reactions(cs_propylene.SM)[1])
```

Alkenes CₙH₂ₙ follow: `CₙH₂ₙ + 3n/2 O₂ → n CO₂ + n H₂O`.

### Example: multiple simultaneous reactions

When more than one dependent species is present, `reactions(cs.SM)` returns one balanced equation per dependent species:

```@example multi_rxns
using ChemistryLab

CH4  = Species("CH4";  aggregate_state = AS_GAS, class = SC_GASFLUID)
C2H6 = Species("C2H6"; aggregate_state = AS_GAS, class = SC_GASFLUID)
C2H4 = Species("C2H4"; aggregate_state = AS_GAS, class = SC_GASFLUID)
O2   = Species("O2";   aggregate_state = AS_GAS, class = SC_GASFLUID)
CO2  = Species("CO2";  aggregate_state = AS_GAS, class = SC_GASFLUID)
H2O  = Species("H2O";  aggregate_state = AS_GAS, class = SC_GASFLUID)

primaries = [CO2, H2O, O2]
cs = ChemicalSystem([CH4, C2H6, C2H4, O2, CO2, H2O], primaries)

rxns = reactions(cs.SM)
for r in rxns
    println((-r).equation)
end
```

!!! tip "When to use algebraic balancing"
    Algebraic balancing via `reactions(cs.SM)` is most useful when:
    - The system has many species and manual balancing would be error-prone.
    - You load species from a database and want all independent reactions automatically.
    - You want to verify that a set of species is stoichiometrically consistent.
