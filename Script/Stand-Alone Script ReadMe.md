# Stand-Alone Script — README

A standalone R port of the Finomena app's BAM pipeline. Runs end-to-end from
raw Zebrabox CSVs to a corrected master results table, with no Python /
PyQt app needed.

The script bundles the raw-data aggregation step that the app's
`data_loader.py` normally performs, then feeds the aggregated long-format
data into the same modeling and correction logic used by
`App/R/scripts/TweedieAR1 BAM.R`.

---

## What it does

1. Reads three input files (`design.csv`, `time_frame.csv`, plus one
   `plate_layout.csv` per plate subdirectory).
2. For each plate, reads its Zebrabox CSV(s), drops the marker rows
   (`SOUND`, `BACK_LIGHT`, `TOP_LIGHT`, `Not identify`), bins to 1-second
   pixel-difference sums per well, and attaches that plate's well-to-condition
   map.
3. Concatenates all plates, assigns phase / group from the `time_frame`, and
   splits each Condition on `"+"` into one column per design variable.
4. Fits a Tweedie BAM with AR(1) autocorrelation (two passes: pilot fit to
   estimate rho, refit with AR(1)).
5. Runs `emmeans` contrasts for the per-variable effect families and the
   Full_Interaction family, then applies the chosen multiple-testing
   correction.
6. Writes `aggregated_data.csv`, `master_results.csv`, and
   `summary_statistics_by_group.csv` to `OUTPUT_DIR`.

---

## Quick start

1. Open `Stand-Alone Script.R` and edit the four top-of-file constants if
   needed (defaults assume the directory layout described below):

   ```r
   INPUT_DIR  <- "Stand-Alone Inputs"
   OUTPUT_DIR <- "Stand-Alone Output"
   ALL_PAIRS  <- FALSE
   CORRECTION_SPEC <- list(
     strategy   = "flat",
     error_rate = "FDR",
     threshold  = 0.05,
     scope      = "per_family"
   )
   ```

2. Make sure the required R packages are installed:

   ```r
   install.packages(c("tidyverse", "data.table", "mgcv", "emmeans",
                      "jsonlite", "gratia", "MASS"))
   ```

3. Run the script from R / RStudio (`source("Stand-Alone Script.R")`) or
   from the command line:

   ```
   Rscript "Stand-Alone Script.R"
   ```

---

## Directory layout

```
<INPUT_DIR>/
├── design.csv                  <-- declares design variables + reference levels
├── time_frame.csv              <-- declares phases and phase-groups
├── Plate 1/                    <-- a "Plate" group (same drug map across replicates)
│   ├── Replicate 1/
│   │   ├── plate_layout.csv    <-- this replicate's well-to-condition map
│   │   └── *.xls               <-- raw Zebrabox files for this replicate
│   └── Replicate 2/
│       ├── plate_layout.csv
│       └── *.xls
├── Plate 2/
│   └── Replicate 1/
│       ├── plate_layout.csv
│       └── *.xls
└── ...
```

- A "plate" is **any** directory under `INPUT_DIR` (at any depth) that
  contains a `plate_layout.csv`. The discovery walk is recursive — flat or
  nested hierarchies both work.
- The plate's name (used as the `plate` factor level in the model) is its
  path **relative to `INPUT_DIR`**, with separators normalized to `_`. So
  `Plate 1/Replicate 1` becomes `"Plate 1_Replicate 1"` and
  `Plate 1/Replicate 2` becomes `"Plate 1_Replicate 2"` — distinct factor
  levels, even when their drug maps happen to be identical.
- Per-replicate layouts are always read independently. There is no
  inheritance from `Plate 1/plate_layout.csv` (if one exists) down to
  `Plate 1/Replicate 1/plate_layout.csv`; each layout file is parsed on its
  own. Different replicates of the same Plate can have **different** layouts.
- Each plate can have its own size and conditions — the script doesn't
  assume any uniformity across plates.
- Multiple raw data files (`.csv`, `.xls`, or `.xlsx`) sitting **directly**
  in a plate's directory are concatenated before binning. The Zebrabox often
  emits tab-separated text files with a `.xls` extension; these are handled
  natively (`fread` auto-detects the delimiter). Data files in
  sub-subdirectories of a plate dir are **not** picked up — keep raw files
  alongside the layout.

---

## Input files

### `design.csv`

Lists the experimental design variables and the reference level (the
"control") for each. Order matters: it sets the order in which variable
values are joined by `"+"` to form the canonical Condition string.

```csv
variable,reference
Genotype,WT
Drug,DMSO
```

The reference Condition string is automatically derived by joining the
reference levels with `"+"` — for the example above it would be `"WT+DMSO"`.

### `time_frame.csv`

Declares the recording's phases and phase-groups. There is **no `Trans`
column** — list every phase you want included in the analysis (including
transition phases if you want them modeled).

```csv
Group,Phase,Start,End
1,Dark,1,56
1,Trans 1,57,72
1,Light,73,115
1,Trans 2,116,145
1,Dark,146,180
2,600Hz,181,210
2,200Hz,211,240
```

- `Group` is the phase-group integer used by the model's `Group_factor`.
- `Phase` is a label (it doesn't enter the model directly, but is preserved
  on the aggregated data for inspection).
- `Start` / `End` are inclusive bounds in seconds, matching the integer
  `time_sec` produced by the aggregation step.
- Any second not covered by any window is dropped from the analysis.

### `plate_layout.csv` (one per plate subdirectory)

A grid mapping wells to Condition strings.

- First column: row letter (`A`, `B`, `C`, ...). The header cell for this
  column is literally `row` — this label is **not** required by the script
  (which reads the column by position), but giving it a name keeps Excel
  from dropping the column on save when you edit the file as a spreadsheet.
- First row (header): `row, 1, 2, 3, ...` — `"row"` then the well-column
  numbers.
- Each interior cell is a Condition string built by `"+"`-joining the
  per-variable values, in the same order declared in `design.csv`.
- Use `"empty"` (or leave the cell blank) for unused wells. Empty wells are
  dropped before modeling.

Example for an 8 × 12 plate with two design variables (Genotype, Drug):

```csv
row,1,2,3,4,5,6,7,8,9,10,11,12
A,WT+DMSO,WT+DMSO,MUT+Mut1,MUT+Mut1,MUT+Mut2,MUT+Mut2,MUT+Mut3,MUT+Mut3,MUT+Mut4,MUT+Mut4,MUT+DMSO,MUT+DMSO
B,WT+DMSO,WT+DMSO,MUT+Mut1,MUT+Mut1,MUT+Mut2,MUT+Mut2,MUT+Mut3,MUT+Mut3,MUT+Mut4,MUT+Mut4,MUT+DMSO,MUT+DMSO
...
```

**Editing tips:**
- If you edit `plate_layout.csv` in Excel, watch out for cells whose contents
  look numeric — Excel may convert text like `1E02` into the scientific
  number `1.00E+02`. If your conditions don't include numeric-looking
  strings you'll be fine; otherwise pre-format the affected cells as Text
  before pasting.
- A plain text editor (VS Code, Notepad, etc.) is the safest way to edit
  the layout — no auto-formatting surprises.

Each non-empty cell must split into exactly N parts on `"+"`, where N is
the number of rows in `design.csv`. The script will refuse to start if any
plate layout has the wrong arity.

---

## Configuration knobs (top of the script)

| Variable | Type | Purpose |
|----------|------|---------|
| `INPUT_DIR` | path | Root directory holding `design.csv`, `time_frame.csv`, and the plate subdirectories. Relative paths are resolved against the script's directory. |
| `OUTPUT_DIR` | path | Where to write `aggregated_data.csv`, `master_results.csv`, and `summary_statistics_by_group.csv`. Created if absent. |
| `ALL_PAIRS` | `TRUE` / `FALSE` | Controls the Full_Interaction contrast set. `TRUE` = all pairwise. `FALSE` = experimental conditions vs. the reference condition only (smaller family, tighter correction). |
| `CORRECTION_SPEC` | list | Multiple-testing correction config (see below). |
| `HELPERS_DIR` | path | Where to find `corrections.R` and `tree_metadata.R`. Default points at `../App/R/scripts/`. |

### `CORRECTION_SPEC` shape

Mirrors the `correction.json` sidecar used by the in-app BAM widget.

**Flat correction (BH for FDR, Holm for FWER) — the default:**

```r
CORRECTION_SPEC <- list(
  strategy   = "flat",
  error_rate = "FDR",       # or "FWER"
  threshold  = 0.05,
  scope      = "per_family"
)
```

**Tree correction (TreeBH for FDR, Holm-gatekeeping or graphical-MCP for FWER):**

```r
CORRECTION_SPEC <- list(
  strategy   = "tree",
  error_rate = "FDR",
  threshold  = 0.05,
  # If `tree` is omitted, the default rescue tree from tree_metadata.R is used.
  tree       = NULL
)
```

For custom tree topologies, see `App/R/scripts/tree_metadata.R` for the
`tree_spec` list format.

---

## Outputs

All written to `OUTPUT_DIR`:

- **`aggregated_data.csv`** — long-format frame after binning, condition
  assignment, and phase join. Useful for inspecting what actually entered
  the model.
- **`master_results.csv`** — one row per emmeans contrast, with raw and
  adjusted p-values, log2-scaled effect sizes, and (for tree corrections)
  tree gating metadata.
- **`summary_statistics_by_group.csv`** — counts of significant contrasts
  per `Test_Family × Group` plus mean absolute effect size.

---

## How this compares to the in-app pipeline

| Step | In-app | Stand-alone |
|------|--------|-------------|
| Raw data loading | `data_loader.py` (Python / PyQt thread) | inline in `Stand-Alone Script.R` |
| Plate layout | GUI grid widget, per-plate | `plate_layout.csv` per plate subdir |
| Experimental design | Conditions widget | `design.csv` |
| Time frame | Plate format widget | `time_frame.csv` |
| Contrast selection | `contrast_selection.json` sidecar | `ALL_PAIRS` TRUE/FALSE toggle |
| Correction | `correction.json` sidecar | `CORRECTION_SPEC` R list |
| Distribution family | Family-selection widget | hardcoded to Tweedie |
| Modeling | `App/R/scripts/TweedieAR1 BAM.R` | same algorithm, inlined |
| Multiple testing | `apply_correction()` from `App/R/scripts/corrections.R` | sources and calls the same function |

The standalone script is intentionally narrower than the app: no per-condition
color palette, no role labels (RC / EXP), no ablation diagnostics, no
posterior-equivalence test, no architecture or random-basis switching, no
family selection. It uses the production-winning configuration (HGAM,
Tweedie, sum-to-zero random basis, AR(1) reset at animal-and-phase
boundaries) directly.

---

## Troubleshooting

- **"Missing R packages"** — run the `install.packages()` snippet in the
  Quick start section.
- **"Could not find correction helper scripts"** — `HELPERS_DIR` is wrong
  for your filesystem. Point it at the `App/R/scripts/` directory in the
  Finomena repo (absolute path is fine).
- **"X cell(s) do not split into N '+'-separated parts"** — at least one
  cell in a `plate_layout.csv` doesn't match the variable arity from
  `design.csv`. The error message names the offending plate and shows a
  few example cells.
- **"Reference level X for variable Y not found in the data"** — the
  reference declared in `design.csv` never appears in any plate layout
  (likely a typo, or the plate layouts were edited without updating the
  design).
- **"Reference condition X not found among Conditions"** — same idea but
  for the combined reference: `paste(reference_values, collapse = "+")`
  is missing from the plate layouts.
- **Pass 1 fails / "falling back to fallback rho = 0.20"** — the pilot fit
  errored. Often a sign the data is too thin (very few animals, very few
  time points). The script proceeds with a reasonable default rho but the
  AR(1) correction will be less well-tuned.
