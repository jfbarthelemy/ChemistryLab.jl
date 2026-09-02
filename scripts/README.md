# ChemistryLab.jl demo and tutorial scripts

Runnable companions to the documentation. Each one is self-contained: it
activates this directory's environment and names every data file through
`datapath`, so it runs from any working directory.

## Running one

```bash
julia --project=scripts scripts/titration_acetic_acid.jl
```

From a REPL, `include("scripts/titration_acetic_acid.jl")` works as well. In
VS Code, opening the file and running it is enough — nothing needs to be
configured, because no path in these scripts depends on the working directory.

The first run of a fresh clone needs the environment instantiated once:

```bash
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
```

## Two rules that keep these scripts working

**Name data files with `datapath`.** The databases live in the package's `data/`
directory, and `datapath("cemdata18-thermofun.json")` returns an absolute path to
one of them. A relative `"data/cemdata18-thermofun.json"` only resolves when the
working directory happens to be the repository root — which is not the case in an
editor whose REPL started elsewhere, nor inside a documentation build. The
readers do fall back to the bundled directory, so older code still runs, but new
code should say `datapath` and mean it.

**A script consumed by the documentation or the test suite never calls
`Pkg.activate`.** The active project is global process state: activating an
environment halfway through a documentation build changes it for every
`@example` block that runs afterwards, on every page. Two scripts are in that
position today, and neither activates anything:

| script | consumed by |
|---|---|
| `ionic_hydration.jl` | `docs/src/examples/ionic_hydration.md` |
| `hydration_calibration.jl` | `docs/src/examples/hydration_calibration.md`, `test/kinetics/test_calibration.jl` |

Both are run with `julia --project=scripts` instead. Every other script in this
directory activates `scripts/` itself.

## The scripts

### Aqueous chemistry and equilibrium

| script | what it shows |
|---|---|
| `aqueous_equilibration_basic.jl` | smallest complete equilibration: species, system, state, solve |
| `titration_acetic_acid.jl` | pH curve of a weak monoacid by NaOH, against Henderson–Hasselbalch |
| `titration_malonic_acid.jl` | the same for a diacid, with both pKa derived from the database |
| `clinker_hydrate_equilibration.jl` | equilibrium assemblage of a clinker with its hydrates |
| `slag_clinker_equilibration.jl` | the same for a slag-blended system |

### Kinetics and calorimetry

| script | what it shows |
|---|---|
| `cement_clinker_kinetics.jl` | Parrot–Killoh dissolution of the four clinker phases |
| `blended_cement_kinetics.jl` | clinker plus supplementary cementitious materials |
| `opc_semiadiabatic_calorimetry.jl` | aggregated solid → solid reactions under a semi-adiabatic calorimeter |
| `ionic_hydration.jl` | hydration driven by the pore solution: no sequencing rule, the assemblage is a result |
| `hydration_calibration.jl` | fitting rate parameters against measured calorimetry |
| `kinetic_species_api_demo.jl` | the `KineticReaction` / rate-model API surface |

### Structures, databases, symbolics

| script | what it shows |
|---|---|
| `tutorial_chemical_structures.jl` | formulas, species, reactions across five databases |
| `tutorial_database_workflows.jl` | loading and inspecting a ThermoFun database |
| `tutorial_symbolic_reactions.jl` | symbolic manipulation of reactions (needs SymPy) |
