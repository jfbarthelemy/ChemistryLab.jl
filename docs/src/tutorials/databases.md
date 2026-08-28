# [Database Interoperability](@id sec-databases)

So far, we have looked at the possibility of creating and manipulating any species, whether they exist or not. If we wanted to create a H₂O⁺⁴ molecule, it would not be a problem. However, you will admit that it is a little strange...

This is why ChemistryLab relies on existing databases, in particular [Cemdata18](https://www.empa.ch/web/s308/thermodynamic-data) and [PSI-Nagra-12-07](https://www.psi.ch/en/les/thermodynamic-databases). Cemdata18 is a chemical thermodynamic database for hydrated Portland cements and alkali-activated materials. PSI-Nagra is a Chemical Thermodynamic Database. The formalism adopted for these databases is that of [Thermofun](https://thermohub.org/thermofun/thermofun/) which is a universal open-source client that delivers thermodynamic properties of substances and reactions at the temperature and pressure of interest. The information is stored in json files.

## [Naming a bundled database](@id sec-datapath)

Several databases ship with ChemistryLab, in its `data/` directory. Name them
with [`datapath`](@ref), which returns an absolute path and therefore does not
depend on the working directory:

```@example datapath
using ChemistryLab

path = datapath("cemdata18-thermofun.json")
isfile(path)
```

The value is an absolute path, so it is not shown here: it depends on where the
package is installed, and printing it would put this machine's directories into
the page. What it points at does not depend on the machine:

```@example datapath
relpath(path, pkgdir(ChemistryLab))
```

```@example datapath
readdir(datapath())
```

This matters in practice: a script written as `build_species("data/cemdata18-thermofun.json")`
only runs when the working directory happens to be the repository root, which is
not the case in an editor whose REPL started elsewhere, nor in a documentation
build. Written with `datapath`, the same line runs from anywhere.

!!! note "Resolution of a bundled name"
    The readers are tolerant, so older code keeps working: a path is tried
    against the working directory **first**, and only then against the bundled
    `data/` directory. A call that already resolves therefore keeps resolving to
    exactly the same file, and a local database always takes precedence over a
    bundled one of the same name. The fallback can only succeed for a name that
    *is* one of the bundled files, so a mistyped path to a file of your own
    still fails, with an error listing what is available.

## Loading species from a database

The simplest way to load species from a ThermoFun-compatible JSON file is `build_species`, which reads the file and directly returns a `Vector{Species}` with compiled thermodynamic functions:

```julia
using ChemistryLab
all_species = build_species(datapath("cemdata18-merged.json"))
```

Each species already carries its molar mass and temperature-dependent thermodynamic functions (Cp⁰, ΔₐH⁰, ΔₐS⁰, ΔₐG⁰, logK⁰) as `SymbolicFunc`s and `NumericFunc`s.

!!! note "Low-level access"
    If you need the raw DataFrames (e.g. to inspect metadata or filter on database columns), the lower-level function `read_thermofun_database` is still available and returns three DataFrames `(df_elements, df_substances, df_reactions)`:
    ```julia
    df_elements, df_substances, df_reactions = read_thermofun_database(datapath("cemdata18-merged.json"))
    ```
    `build_species(df_substances)` can then be called on the filtered DataFrame.

## Filtering species with `speciation`

In practice, only a small subset of the database is relevant to a given problem. `speciation` filters a species list to those whose atomic composition is a subset of the atoms found in a set of *seed* species:

```@example database
using ChemistryLab #hide
all_species = build_species(datapath("cemdata18-merged.json")) #hide
# Keep only species that can form from the calcite / water system
species_calcite = speciation(all_species, split("Cal H2O@");
                             aggregate_state=[AS_AQUEOUS],
                             exclude_species=split("H2@ O2@ CH4@"))
dict_species_calcite = Dict(symbol(s) => s for s in species_calcite)
```

The `aggregate_state` keyword restricts results to aqueous species (the seed species themselves — `Cal` and `H2O@` — are always kept through `include_species` internally). For example, the properties of Ca(HCO₃)⁺ can then be read as:

```@example database
dict_species_calcite["Ca(HCO3)+"]
```

### `speciation` signatures

`speciation` accepts seed arguments in three forms:

| Seed argument | Description |
|:--------------|:------------|
| `Vector{Symbol}` | Explicit list of atom symbols |
| `Vector{<:AbstractSpecies}` | Species objects — their union of atoms defines the space |
| `Vector{<:AbstractString}` | Species symbol strings — looked up in `species_list` |

Common keyword arguments:

| Keyword | Default | Description |
|:--------|:--------|:------------|
| `aggregate_state` | all states | restrict to `[AS_AQUEOUS]`, `[AS_CRYSTAL]`, etc. |
| `class` | all classes | restrict to `[SC_AQSOLUTE]`, etc. |
| `exclude_species` | `[]` | species (or symbols) to always exclude |
| `include_species` | `[]` | species to always include regardless of composition |

## Primary species extraction

It is also possible to retrieve primary species from the Cemdata18 database. Primary species are a minimal subset such that every other species can be expressed as their linear combination.

```julia
df_primaries = extract_primary_species(datapath("CEMDATA18-31-03-2022-phaseVol.dat"))
show(df_primaries, allcols=true, allrows=true)
```

---

## Building solid solution phases from a TOML file

Solid solution phases (e.g. C-S-H gel, AFm) group several end-member species with a
mixing model. Because database species carry `class = SC_COMPONENT` rather than
`class = SC_SSENDMEMBER`, they cannot be passed directly to [`SolidSolutionPhase`](@ref).

[`build_solid_solutions`](@ref) automates the full pipeline:

1. Reads phase definitions from a TOML file.
2. Looks up each end-member in a pre-built species dictionary.
3. Requalifies end-members to `SC_SSENDMEMBER` via [`with_class`](@ref).
4. Constructs and returns a `Vector{SolidSolutionPhase}`.

### TOML format

Each `[[solid_solution]]` entry specifies a phase name, the list of end-member
symbols (as they appear in the database), and the mixing model:

```toml
# Ideal solid solution (any number of end-members)
[[solid_solution]]
name        = "CSHQ"
end_members = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD"]
model       = "ideal"
source      = "Lothenbach2015"

# Binary Redlich-Kister (exactly 2 end-members, parameters in J/mol)
[[solid_solution]]
name        = "AFm"
end_members = ["Ms", "Mc"]
model       = "redlich_kister"
a0          = 3000.0
a1          = 500.0
a2          = 0.0
source      = "Lothenbach2019"
```

The file `data/solid_solutions.toml` shipped with ChemistryLab.jl contains
pre-calibrated entries for the main cemdata18 solid solutions (CSHQ, AFm, Hydrogarnet,
Ettringite_ss, Hydrotalcite).

### Usage

```julia
using ChemistryLab, DynamicQuantities

substances = build_species(datapath("cemdata18-thermofun.json"))
dict       = Dict(symbol(s) => s for s in substances)

# Build all solid solution phases defined in the TOML
ss_phases = build_solid_solutions(datapath("solid_solutions.toml"), dict)

# Use them directly in ChemicalSystem
cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions = ss_phases)
```

Phases whose end-members are not found in `dict` are skipped with a warning
(pass `skip_missing = false` to raise an error instead).

!!! note "Manual requalification"
    If you only need one or two end-members, you can requalify them individually
    with [`with_class`](@ref) instead of going through the TOML file:

```julia
em = with_class(dict["CSHQ-TobD"], SC_SSENDMEMBER)
```

---
