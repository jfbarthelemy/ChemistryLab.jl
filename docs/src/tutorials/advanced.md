# Advanced Topics

This page covers advanced usage patterns and techniques for working with ChemistryLab: complex formula transformations, cement species, reaction algebra, and programmatic database operations.

## Advanced Formula Manipulation

### Working with Quantities (units)

If your project uses `DynamicQuantities`, you can attach units to stoichiometric coefficients and apply transformations that preserve dimensions:

```julia
using ChemistryLab, DynamicQuantities
using OrderedCollections

# Create formulas with quantities
comp_with_units = OrderedDict(:H => 2.0u"g/mol", :O => 1.0u"g/mol")
f = Formula(comp_with_units)

# Apply a function across coefficients (units preserved where possible)
f_scaled = apply(x -> x * 2, f)
```

### Arithmetic with fractional stoichiometry

ChemistryLab preserves rational coefficients when parsing fractional formulas (use `//` for rational notation) and when doing arithmetic.

As an example, Rankinite (Ca₃Si₂O₇) can be written as a half-unit formula where each coefficient is a ratio:

```@example advanced
using ChemistryLab #hide
f = Formula("Ca3//2Si1O7//2")   # half the formula unit of Rankinite (Ca₃Si₂O₇)
composition(f)  # OrderedDict(:Ca => 3//2, :Si => 1, :O => 7//2)
```

```@example advanced
f2 = f * 2   # full Rankinite formula
composition(f2)  # rational coefficients become integer-valued Rationals; plain integers stay Int
```

### Converting between notations

Switch between Phreeqc (plain text), Unicode (pretty), and colored output:

```@example advanced
f = Formula("Fe+3")

expr(f)       # "Fe+3"
phreeqc(f)    # "Fe+3"
unicode(f)    # "Fe³⁺"
colored(f)    # colored terminal string (ANSI)
```

```@example advanced
# Convert between notations programmatically
phreeqc_to_unicode("SO4-2")   # "SO₄²⁻"
unicode_to_phreeqc("SO₄²⁻")  # "SO4-2"
```

---

## Advanced Species Operations

### CemSpecies and oxide-to-atom decomposition

`CemSpecies` represents a species in cement nomenclature (as oxides, e.g., C3S = 3CaO·SiO₂) and automatically converts to/from atomic composition:

```@example advanced
# Create a cement species from oxide formula
c3s = CemSpecies("C3S")
oxides(c3s)   # OrderedDict(:C => 3, :S => 1)  — cement components (C=CaO, S=SiO₂)
atoms(c3s)    # C3S = Ca₃SiO₅: OrderedDict(:Ca => 3, :Si => 1, :O => 5)
```

```@example advanced
# Convert Species → CemSpecies
# C3S = Ca₃SiO₅  (3 CaO + 1 SiO₂)
s = Species("Ca3SiO5")
cem_s = CemSpecies(s)  # automatically decomposes to oxides
```

### Type conversion

Convert numeric types in a `Species` using `convert` on its `Formula`:

```julia
s1 = Species("H2O")              # Species{Int64} by default
s2 = Species(convert(Float64, formula(s1)))  # Species{Float64}
```

### Custom properties and thermodynamic data

Attach arbitrary properties to species (molar mass is auto-calculated; add more):

```julia
water = Species("H2O")
water[:Cp] = 75.3       # J/(mol·K)
water[:ΔfH0] = -285.8   # kJ/mol

# Access via bracket or dot notation
water[:Cp]           # 75.3
water.Cp             # 75.3 (same)
haskey(water, :Cp)   # true
```

---

## Advanced Reaction Operations

### Reaction algebra (arithmetic)

Combine and transform reactions using operator overloading:

```@example advanced
r1 = Reaction("H2O = H+ + OH-")
r2 = Reaction("H2O = H2 + 1//2 O2")

r_sum  = r1 + r2     # combine all species and coefficients
r_diff = r1 - r2     # subtract r2 from r1 (reverse coefficients of r2)

# Scalar multiplication — note: scalar must be on the left
r_scaled = 2 * r1    # double all coefficients
r_half   = (1/2) * r1   # halve all coefficients
```

### Reaction simplification

Eliminate species appearing on both sides (cancel them). For example, combining two reactions and simplifying:

```@example advanced
r_water = Reaction("H+ + OH- = H2O")
r_elec  = Reaction("H2O = H2 + 1//2 O2")

r_combined = r_water + r_elec    # H⁺ + OH⁻ + H₂O = H₂O + H₂ + ½O₂
r_simple   = simplify_reaction(r_combined)  # H₂O cancels → H⁺ + OH⁻ = H₂ + ½O₂
```

### Building reactions from species lists

Construct a balanced reaction by passing a list of species; the stoichiometry is solved automatically:

```@example advanced
h2o = Species("H2O")
h2  = Species("H2")
o2  = Species("O2")

# ChemistryLab solves for the integer stoichiometry
r = Reaction([h2o, h2, o2]; equal_sign='→')
```

---

## Advanced Stoichiometric Matrix Operations

### Mass-based stoichiometric matrices

By default, stoichiometric matrices use atom counts. `mass_matrix` returns a version with coefficients scaled by molar masses:

```@example advanced
species  = [Species("H2O"), Species("H2"), Species("O2")]
SM       = StoichMatrix(species)        # atom count coefficients
SM_mass  = mass_matrix(SM)             # mass-weighted coefficients
pprint(SM_mass)
```

### Extracting independent and dependent species

`StoichMatrix` automatically identifies independent components (primaries) and expresses all species as combinations of them:

```julia
species = [Species("Ca+2"), Species("OH-"), Species("CaOH+")]
SM = StoichMatrix(species)

SM.primaries   # basis components (independent)
SM.species     # all species (columns of the matrix)
SM.A           # stoichiometric matrix (primaries × species)
SM.N           # nullspace matrix
```

### Canonical form and redox handling

`CanonicalStoichMatrix` reorders species automatically and handles charged species (charge balance via `:Zz`):

```julia
species = [Species("Fe+2"), Species("Fe+3"), Species("H+")]
CSM = CanonicalStoichMatrix(species)
# CSM includes :Zz column for charge balance
```

---

## Advanced Database Operations

### Standard workflow: `build_species` + `speciation`

The recommended entry point is `build_species(filename)`, which reads a ThermoFun JSON database and returns a `Vector{Species}` with compiled thermodynamic functions. Use `speciation` to filter to a relevant chemical sub-space:

```julia
using ChemistryLab

# Load all species from the database
all_species = build_species(datapath("cemdata18-merged.json"))

# Filter to species compatible with a seed set (Ca–C–H–O system)
species = speciation(all_species, split("Cal H2O@ CO2");
                     aggregate_state=[AS_AQUEOUS],
                     exclude_species=split("H2@ O2@ CH4@"))
```

To load only a specific subset of species by symbol, pass a list as the second argument to `build_species`:

```julia
# Load only the listed symbols (faster for large databases)
selected = build_species(datapath("cemdata18-merged.json"), split("Cal H2O@ H+ OH- Ca+2 CO3-2"))
```

### Low-level access: `read_thermofun_database`

For direct DataFrame access (e.g. to inspect raw fields or apply custom filters), use `read_thermofun_database`. Note that `aggregate_state` is stored as a single-entry `Dict` — use `only(values(...))` to extract its string value:

```julia
df_elements, df_substances, df_reactions = read_thermofun_database(datapath("cemdata18-thermofun.json"))

# Filter aqueous substances
aqueous = filter(row -> only(values(row.aggregate_state)) == "AS_AQUEOUS", df_substances)

# Filter by charge
charged_species = filter(row -> row.charge != 0, df_substances)

# Convert the filtered DataFrame to Species objects
species_list = build_species(aqueous)
```

### Merging databases

`merge_json` combines a ThermoFun JSON file with a Phreeqc `.dat` file (phase definitions) into a single merged JSON:

```julia
merge_json(
    datapath("cemdata18-thermofun.json"),            # bundled ThermoFun database
    datapath("CEMDATA18-31-03-2022-phaseVol.dat"),   # bundled Phreeqc phase definitions
    "cemdata18-merged.json",                         # output: an ordinary path of your choosing
)
```

The two inputs are resolved against the bundled `data/` directory, so the call
works from any working directory. The third argument is a file to be *written*,
and is never resolved that way: give it the path where you want the merged
database. A pre-merged `cemdata18-merged.json` already ships with the package.

The merged file can then be loaded with `build_species` as usual.

### Cemdata .dat parsing and extraction

Extract primary species from a Phreeqc / Cemdata `.dat` file:

```julia
# Extract primary aqueous species (SOLUTION_MASTER_SPECIES section)
df_primary = extract_primary_species("path/to/file.dat")
```

---

## Tips for complex workflows

1. **Batch operations on Species/Formulas**: Use `apply(func, formula)` or `apply(func, species)` to transform stoichiometric values or properties in bulk.

2. **Debugging stoichiometric matrices**: Call `pprint(SM)` or `pprint(SM; label=:name)` to print a formatted, colored matrix table to the REPL.

3. **Type stability**: Prefer homogeneous numeric types (`Species{Float64}` or `Species{Rational}`) across collections for performance.

4. **Combining databases**: Use `read_thermofun_database` to load multiple sources, then concatenate DataFrames or merge reactions as needed.

5. **Custom properties as metadata**: Attach source information, uncertainty, or application-specific data as custom properties — they are preserved during conversions and copies.

## Next steps

- Explore the reference documentation for full API signatures and options.
- See `docs/src/examples/` for complete worked examples (cement hydration calculations, equilibrium problems, etc.).
- If you encounter edge cases or need custom workflows, open an issue or extend the examples here.
