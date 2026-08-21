# Experimental data

Measured data, vendored so that the test suite and the documentation build
without network access. Unlike everything else under `data/`, these files are
**not** thermodynamic databases and **not** covered by the package's license —
see [`LICENSE`](LICENSE) in this directory.

| file | cement | w/b | Blaine | *T* | duration | final *Q* | role |
|--- |--- |--- |--- |--- |--- |--- |--- |
| `smilauer2025-122-cemI-52.5R-cizkovice.csv` | CEM I 52.5 R, Čížkovice | 0.50 | 397 m²/kg | 20 °C | 262 h | 376 J/g | calibration target |
| `smilauer2025-116-cemI-52.5R-ladce.csv` | CEM I 52.5 R, Ladce | 0.45 | 415 m²/kg | 20 °C | 617 h | 355 J/g | holdout |

Both come from the CC-BY-4.0 Zenodo deposit
[10.5281/zenodo.15212785](https://doi.org/10.5281/zenodo.15212785) of Šmilauer
and Reiterman, and are used by
[`scripts/hydration_calibration.jl`](../../scripts/hydration_calibration.jl).

## Format

A `#`-commented provenance block, then a plain comma-separated table with the
header `time_h,heat_flow_W_per_g,heat_J_per_g`. `read_calorimetry` in the script
reads it, and reads the metadata lines it needs out of the comment block.

Any file in this shape can be substituted: that is the seam for calibrating
against your own measurements. The four comment lines the reader looks for are
`blaine`, `wb`, `temperature` and `Released heat up to`.

## Why these two records

Of the 65 cements in the deposit, 14 are plain CEM I. Record `122` was chosen as
the calibration target because its w/b of 0.50 and Blaine of 397 m²/kg sit within
a few percent of the CEM I 52.5 N formulation that
[`scripts/ionic_hydration.jl`](../../scripts/ionic_hydration.jl) already runs
after Lavergne et al. (2018) — so the one input the deposit does not report, the
clinker phase composition, is the least of an extrapolation available.

Record `116` is the holdout: same nominal strength class, but a different w/b
**and** a different Blaine. Predicting it with parameters fitted on `122`, without
refitting, is a direct test of the two corrections that claim to carry a
calibration across mixes — `powers_alpha_max(w/b)` and `blaine_factor(blaine)`.

## Caveats, all of them

These matter for interpreting a fit and are stated in the documentation page too.

- **The records do not start at *t* = 0.** Both are truncated at the calorimeter's
  thermal equilibration — 1.02 h for `122`, 1.21 h for `116`. The header line
  `Released heat up to 45 minutes (J/g of binder)` gives the heat already gone
  before the record begins (12.0 J/g, "estimated", for both). A model compared
  against these curves must be offset by that amount, not zeroed at the first
  sample.
- **No phase composition, no oxide analysis.** The deposit reports the cement
  type, the Blaine fineness, the w/b ratio and the temperature — not the clinker
  mineralogy. The composition therefore has to be assumed, and only the products
  of a fitted rate multiplier and its assumed phase fraction are identifiable.
- **One temperature.** Everything was measured at 20 °C, so activation energies
  are not identifiable from these data and must be held at published values.
- **Filename and internal header disagree for the Ladce records.** The file named
  `116-CEM I 52.5 R Ladce-415.csv` in the deposit carries
  `Cement name: CEM I 42.5R Ladce-415` internally; the same is true of `115`.
  Which of the two is right is not resolvable from the deposit, which is why
  `116` is used as a plausibility check rather than as a validation.
- **The depositors' own reference fit is quoted twice, inconsistently.** The
  provenance block reports the affinity-model parameters from the deposit's
  `Fit-affinity.zip`, which give `Ea 38300 J/mol` for every cement. The header of
  the source data file for `122` carries a different line, with `AEn 40000
  J/mol`. The value from the fit archive is the one kept, since it is the one that
  goes with the fitted `B1`, `B2` and `eta`.

## Regenerating

```
julia data/experimental/regenerate.jl            # download and rebuild
julia data/experimental/regenerate.jl --check    # rebuild and compare only
```

Needs network access and `unzip`; neither is needed to *use* the files. The
subset is deterministic — nearest row to each of 500 log-spaced times — so
`--check` is an exact byte comparison, and it is how you tell that nothing here
was hand-edited.
