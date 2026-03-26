# Finomena: Zebrafish Behavioral Analysis Suite

A desktop application for processing, visualizing, and statistically analyzing
zebrafish locomotor data captured by the Zebrabox system. The app takes raw
Zebrabox output files, lets you define your experimental design, and then runs
two complementary analyses: a Generalized Additive Mixed Model (GAMM) for
confirmatory statistics, and CATCH22/24 time-series feature clustering for
unsupervised pattern discovery.

---

## Table of Contents

1. [Installation & Requirements](#installation--requirements)
2. [Zebrabox Input Data Format](#zebrabox-input-data-format)
3. [Naming Conditions — The `+` Convention](#naming-conditions--the--convention)
4. [Tab-by-Tab Walkthrough](#tab-by-tab-walkthrough)
   - [Tab 1: Experimental Design](#tab-1-experimental-design)
   - [Tab 2: Data Loading](#tab-2-data-loading)
   - [Tab 3: GAMM Analysis](#tab-3-gamm-analysis)
   - [Tab 4: Catch22 Clustering](#tab-4-catch22-clustering)
5. [Why GAMM? — Statistical Rationale](#why-gamm--statistical-rationale)
6. [P-Value Corrections — What They Mean and When to Use Each](#p-value-corrections--what-they-mean-and-when-to-use-each)
7. [CATCH22/24 Feature Reference](#catch2224-feature-reference)
8. [Output Files](#output-files)

---

## Installation & Requirements

### Python dependencies

```
pip install PySide6 pandas numpy scikit-learn seaborn matplotlib pymupdf pycatch22 pyqtdarktheme
```

### R (bundled automatically)

The packaged app includes a complete portable R installation — **no system R is
needed on the target machine**. The build scripts (`build_windows.bat`,
`build_mac.sh`, `build_linux.sh`) automatically detect your system R, copy it
into the distribution folder, and install the required R packages (`tidyverse`,
`data.table`, `mgcv`, `emmeans`) into the bundled library.

R is only needed on the **build machine** during packaging. The person running
`Finomena.exe` (or `Finomena.app`) does not need R installed.

### Running the app (development)

```
python App/app.py
```

### Building the self-contained package

```
# Windows
build_windows.bat

# macOS
bash build_mac.sh

# Linux
bash build_linux.sh
```

The output in `dist/Finomena/` (or `dist/Finomena.app` on Mac) is fully
self-contained and can be distributed as-is.

---

## Zebrabox Input Data Format

The app reads the raw tab-delimited `.csv` (or `.xls`/`.xlsx`) files exported
directly from the Zebrabox ViewPoint software. Each file must contain at least
these three columns (exact names, case-sensitive):

| Column     | Description |
|------------|-------------|
| `time`     | Timestamp in **microseconds** since recording start |
| `location` | Well identifier string, e.g. `"location_1"`, `"location_48"` |
| `data1`    | Raw pixel difference value for that time point and well |

Rows with `time <= 0` are automatically discarded.

**Directory structure:**
Each directory you add in the Data Loading tab represents one **plate replicate**.
The app will recursively search every subfolder for `.csv`/`.xls`/`.xlsx` files
and treat all files found within a single directory tree as belonging to the same
plate. You can have one file per well, one file for the whole plate, or any
combination — the app concatenates them.

---

## Naming Conditions — The `+` Convention

### How conditions are split

Every condition label you define in the app is split on the `+` character to
produce two sub-categories:

```
Condition label → Genotype + Drug
e.g.  "WT+DMSO"    →  Genotype = "WT",  Drug = "DMSO"
      "KO+compound" →  Genotype = "KO", Drug = "compound"
```

This split drives the entire statistical model: Genotype and Drug become the
two main effects, and their combination becomes the interaction term.

### Rules for condition names

- **Always include exactly one `+`** in every condition name.
- **The part before `+` is Genotype, the part after is Drug.**
  If your experiment only has one factor, use a placeholder:
  `WT+vehicle`, `mutant+vehicle`, etc.
- **Do not use `+` for any other purpose** in a condition name. A name like
  `"WT+high+dose"` will only use `"WT"` as Genotype and `"high+dose"` (or just
  `"high"`) as Drug, depending on how the split is applied, which will cause
  incorrect model setup.
- **Be consistent** — `"WT"` and `"wt"` are treated as different genotypes.
  Choose a casing convention and stick to it across all plates.

### Examples

| Condition label | Genotype | Drug |
|-----------------|----------|------|
| `WT+DMSO`       | WT       | DMSO |
| `KO+DMSO`       | KO       | DMSO |
| `WT+compound`   | WT       | compound |
| `KO+compound`   | KO       | compound |

The reference controls (set in the Reference Controls sub-tab) should match
these exact strings — e.g., Reference Genotype = `WT`, Reference Drug = `DMSO`.

---

## Tab-by-Tab Walkthrough

### Tab 1: Experimental Design

Complete this tab **before** loading data. The configuration here drives
everything downstream.

#### Sub-tab 1 — Experimental Setup (Phases)

Define the temporal structure of your experiment. Each row is one **phase**:

| Column   | Meaning |
|----------|---------|
| End (s)  | The time in seconds at which this phase ends |
| Phase    | A descriptive name for the phase (e.g., "Baseline", "Dark 1", "Light 1") |
| Group    | An **integer** that groups related phases into a single statistical model. All rows sharing the same Group number are analyzed together in one GAMM. |

**Critical:** The Group column controls how the R script splits the data for
modeling (see [Why GAMM?](#why-gamm--statistical-rationale)). All phases that
share the same Group integer are pooled into a single contiguous time window
and fitted with one GAMM. Group phases by **stimulus modality** — all
light/dark phases together, all color stimuli together, all vibration stimuli
together — so that each model captures a coherent behavioral context and the
smooth terms can fit the data without having to bridge across fundamentally
different stimulus types.

A realistic three-modality experiment looks like this:

| Group | Phase Name       | Length (s) | Start (s) | End (s) |
|-------|------------------|-----------|-----------|---------|
| 1     | Dark 1           | 300       | 0         | 300     |
| 1     | Light 1          | 300       | 300       | 600     |
| 1     | Dark 2           | 300       | 600       | 900     |
| 1     | Light 2          | 300       | 900       | 1200    |
| 2     | White            | 300       | 1200      | 1500    |
| 2     | Blue             | 300       | 1500      | 1800    |
| 2     | Green            | 300       | 1800      | 2100    |
| 3     | None             | 300       | 2100      | 2400    |
| 3     | Vibration 200 Hz | 300       | 2400      | 2700    |
| 3     | Vibration 600 Hz | 300       | 2700      | 3000    |

In this example:
- **Group 1** — all light/dark phases are modeled together. The GAMM for this
  group fits a nonlinear activity trajectory over the full 0–1200 s window,
  with condition differences estimated across the repeated light and dark
  transitions simultaneously.
- **Group 2** — white, blue, and green light stimuli are modeled together.
  These share the same stimulus modality (chromatic light) so their behavioral
  dynamics are comparable and a single model can capture the full 1200–2100 s
  window coherently.
- **Group 3** — silence baseline, 200 Hz vibration, and 600 Hz vibration are
  modeled together. Pooling the vibration conditions alongside the silent
  baseline allows the model to directly estimate how each frequency deviates
  from the no-stimulus response within the same time window.

The Length (s) column is computed automatically by the app (= End − Start of
the preceding phase). You can edit Length or Start directly in the table and
the other columns will update accordingly.

**Save/Load config:** Use the "💾 Save Config…" and "📂 Load Config…" buttons
at the top of the Experimental Design tab to save your entire experimental
design (phases, conditions, plate layout, reference controls) to a `.json` file
for reuse across sessions.

#### Sub-tab 2 — Experimental Conditions

Add one row per condition. Each condition needs:
- A **name** following the `Genotype+Drug` convention described above
- A **color** used to paint wells on the plate layout

Click "Add Condition" to add a row; click the color swatch to change the color.

#### Sub-tab 3 — Plate Format

1. Select your plate format (96-well, 48-well, 24-well, 6-well) from the
   dropdown.
2. Select a condition from the palette on the right.
3. Click individual wells or click-and-drag to paint them.
4. Right-click a well to clear it.
5. Unpainted wells are treated as empty and excluded from all analysis.

The plate layout you define here is applied to **every** plate replicate loaded
in the Data Loading tab. All replicate plates must use the same layout.

#### Sub-tab 4 — Reference Controls

Set the baseline levels used as the reference in the GAMM contrasts:

- **Reference Genotype:** The genotype all others are compared against (e.g., `WT`)
- **Reference Drug:** The drug/vehicle all others are compared against (e.g., `DMSO`)
- **Reference Condition:** The full condition used as baseline in catch22
  comparisons (e.g., `WT+DMSO`)

These values are passed directly to the R script and to the Catch22 tab.

---

### Tab 2: Data Loading

#### Sub-tab 1 — Load Data

1. Click **"Add Directory…"** to add one directory per plate replicate.
   The order in the list determines the plate number (top = Plate 1).
2. Use **"Move Up" / "Move Down"** to reorder plates if needed.
3. Click **"Load & Process Data"**.

The app will:
- Scan all files recursively in each directory
- Parse `time`, `location`, and `data1` columns
- Map each well's `location` number to its coordinate (A1–H12, etc.) and
  the condition you painted in the plate layout
- Convert microsecond timestamps to seconds
- Sum pixel-difference values per well per second
- Assign Phase and Group labels based on your phase table (temporal join)
- Split each Condition into Genotype + Drug

Progress is shown per-file so you can see the app is working.

#### Sub-tab 2 — Sanity Check

After loading, the processed data table is shown here. Verify:
- The `Condition` column has no unexpected blanks or `NaN` values
- `Phase` and `Group` columns are populated for all rows
- `plate` numbers match the number of directories you added
- Row count is reasonable (a 96-well plate for ~60 min at 1-second bins
  ≈ 3,600 rows per well × number of conditions)

Click the column headers to sort. Data does not proceed to analysis until
you confirm it looks correct here.

---

### Tab 3: GAMM Analysis

#### Section 1 — Statistical Correction Methods

Set both correction methods **before** running the analysis (see the
[P-Value Corrections](#p-value-corrections--what-they-mean-and-when-to-use-each)
section for full guidance).

**Global Multiple-Testing Correction** — applied across all Groups × contrasts:
- **Benjamini-Hochberg (BH)** — controls the False Discovery Rate
- **Holm** — controls the Family-Wise Error Rate (more conservative)

**Contrast Correction** — applied within each emmeans contrast call:
- **Dunnett's Test** — many-to-one, specifically designed for treatment vs.
  a single reference
- **Dunnett-Šidák** — closely related alternative using the Šidák inequality

#### Section 2 — Export Data for R

Click **"Export to CSV…"** to save the processed data in the format expected
by the R script. This CSV has these columns:

`time_sec`, `location`, `loc_coord`, `pixel_diff`, `Condition`, `Phase`,
`Group`, `animal_id`, `plate`

`animal_id` is constructed as `plate_location` (e.g., `1_location_5`) to
ensure each fish is uniquely identified across plates.

#### Section 3 — Run GAMM Analysis

Click **"Run R Analysis"**. You will be asked to choose an output directory
where the PDF figures and summary CSV will be saved. The R output log streams
live so you can monitor progress. A typical 8-group experiment with 96 wells
and ~60 minutes of data takes 10–30 minutes depending on your CPU.

Click **"Stop"** to terminate the R process early.

#### Section 4 — R Output Log

The raw stdout/stderr from Rscript is shown here. Key things to look for:
- The `Optimal Rho` value printed for each Group (should be between 0 and 1;
  values close to 1 mean strong autocorrelation, which is common in zebrafish
  locomotor data)
- Any `Warning:` messages from mgcv (usually about `k` basis dimension; if
  you see `k-index < 1`, the smooth may be under-specified — contact the
  developer)
- The final `CLEAN MASTER TABLE` showing all contrasts and corrected p-values

#### Section 5 — Results

PDF figures rendered inline after the run completes. Use the arrow buttons to
page through figures. Click **"Save All"** to copy the figures to a directory
of your choice.

---

### Tab 4: Catch22 Clustering

#### Section 1 — Data Status / Condition Comparison

After data is loaded, the condition dropdowns populate automatically. Set:
- **Condition A:** The "treatment" condition whose top-driving features you
  want to identify
- **Condition B:** The "reference" condition to compare against (defaults to
  your Reference Condition from the Experimental Design tab)

#### Section 2 — Run Analysis

**Include mean & variance features (CATCH24)** checkbox:
- **Checked (default):** Extracts 24 features per phase — the standard 22
  plus the signal mean and variance. Use this when overall activity level
  is expected to differ between conditions and you want the clustering to
  reflect that.
- **Unchecked:** Extracts only the 22 dynamics-based features. Use this when
  mean activity levels are similar between conditions and you want the
  clustering to focus purely on the temporal *pattern* of movement rather
  than the magnitude.

Click **"Run catch22 Analysis"**. The progress bar tracks feature extraction
per animal and per phase. A typical run on 96 wells × 8 phases takes 2–10
minutes.

#### Section 3 — Phase Clustermaps

One clustermap per Group is generated. Each heatmap shows:
- **Rows:** CATCH22/24 features (in their dendrogram-clustered order)
- **Columns:** Conditions (in their dendrogram-clustered order)
- **Color:** Robust Z-score (median-centered). Blue = below-median activity for
  that feature; red = above-median. The color intensity indicates how far from
  the center that condition falls for that feature.

Conditions that cluster together have similar temporal behavioral profiles for
that phase. Conditions far apart have distinct dynamics.

Use the **◀ / ▶** navigation buttons to step through phases. Click **"Save All"**
to export all phase clustermaps as PNG files.

#### Section 4 — Top Driving Features

For the two conditions you selected, this panel reports the 5 features with
the largest absolute difference in their robust Z-score centroids (using
the Chebyshev / L∞ distance, the same metric used in the clustermaps).
The report is shown for both directions (A vs B and B vs A) so you can see
which features most strongly characterize each condition relative to the other.

---

## Why GAMM? — Statistical Rationale

### Why not a t-test per phase group?

A t-test (or ANOVA) at each time point or phase assumes:
1. **Observations are independent** — but a single fish measured repeatedly
   across hundreds of seconds has strongly correlated consecutive measurements.
   Ignoring this inflates the effective sample size and produces falsely small
   p-values.
2. **The relationship with time is linear** — but zebrafish locomotor responses
   to light/dark transitions are highly nonlinear (rapid acceleration, gradual
   adaptation, oscillations).
3. **The residuals are normally distributed** — but pixel-difference data is
   non-negative, zero-inflated (fish rest), and right-skewed with occasional
   large bursts. A Gaussian model will fit poorly.

### Why a GAMM with Tweedie family and AR1?

**Generalized Additive Model (GAM):** The smooth terms `s(time_in_group)` and
`s(time_in_group, by = Geno_Drug)` capture the nonlinear trajectory of activity
over time without assuming any particular shape. The `mgcv` package selects
the smoothness automatically using penalized regression.

**Mixed model (the "M" in GAMM):** The random smooth `s(time_in_group, animal_id, bs="fs")`
treats each fish as its own random trajectory, accounting for between-animal
variability without overfitting. This is equivalent to a repeated-measures
design but far more flexible.

**Tweedie family:** The Tweedie distribution natively handles the mixture of
zeros (fish at rest) and positive continuous values (fish moving) in a single
model. It is parameterized by a power parameter `p` that is estimated from the
data, and it naturally accommodates the right-skewed, non-negative nature of
pixel-difference data. A Gaussian model would require log-transformation and
still could not handle the zeros properly.

**AR1 (first-order autoregressive) correlation:** The `rho` parameter in `bam()`
models the correlation between consecutive residuals within each fish. Without
this, the model treats each second of data as independent, which it is not —
a fish that was active 1 second ago is more likely to be active now. The app
estimates the optimal `rho` for each Group separately (Pass 1), then fits the
final model with that AR1 structure (Pass 2). This gives valid standard errors
and p-values.

### Why split the data into separate models by Group?

Each Group receives its own independently fitted model because the behavioral
responses zebrafish produce under different stimulus modalities are
fundamentally distinct — not just quantitatively different, but qualitatively
incomparable in ways that make a single joint model biologically inappropriate.

1. **Different stimulus types produce data that are not statistically comparable.**
   The model estimates a single set of parameters — including how
   spread-out the data are, how many zeros to expect, and the shape of the
   activity trajectory over time — that are assumed to apply consistently
   throughout. When you mix light/dark cycles with vibration trials in the
   same model, those parameters end up being a compromise that fits none of
   the conditions well. Instead of detecting differences between your
   experimental groups, the model ends up using most of its flexibility just
   to bridge the jump between one stimulus type and the next. Keeping each
   modality in its own model means the statistical estimates are grounded in
   data that are actually comparable to each other.

2. **The temporal dynamics within each modality are self-contained.** Within
   a Group, the fish are adapting, habituating, or sensitizing to stimuli of
   the same type over a coherent time window. The smooth terms in the GAMM
   capture that arc meaningfully. Pooling across modalities would force the
   model to fit an incoherent trajectory — one that crosses between entirely
   different behavioral states — and the smooth would have no interpretable
   biological meaning.

3. **The autocorrelation structure differs by modality.** The AR1 `rho`
   parameter is estimated fresh for each Group. Phases with high arousal (e.g.,
   active vibration response) have different within-animal serial correlation
   than quiet baseline or adapted light phases. A single global `rho` would
   misrepresent the dependence structure in every Group simultaneously.

The consequence of this design is that we are testing the same hypotheses
(drug effect, genotype effect, interaction) independently in each Group —
which is the right thing to do scientifically, but creates a multiple testing
problem that the two-stage p-value correction is designed to address.

---

## P-Value Corrections — What They Mean and When to Use Each

### The two-stage correction structure

Every p-value in the final output has been corrected **twice**, at two
different levels:

**Stage 1 — Contrast correction (within each Group):**
Dunnett's test (or Dunnett-Šidák) is applied inside each `emmeans` contrast
call. Within a single Group, you are comparing each treatment condition against
the reference simultaneously. For example, if you have 3 drug conditions, there
are 3 drug contrasts in that Group, all using the same DMSO reference. Dunnett's
exploits the shared correlation structure (all contrasts share the reference
group) to be more powerful than Bonferroni while still controlling the
within-Group false positive rate.

**Stage 2 — Global correction (across all Groups):**
After Stage 1, all contrast p-values from all Groups are pooled and a second
round of correction is applied across the entire table. This controls for the
fact that you are testing the same biological hypothesis in 8 separate models
(one per Group). Without this second correction, with 8 independent tests you
would expect roughly a 34% chance of at least one false positive at p < 0.05,
even if nothing were truly significant.

This sequential design is conservative by intent: it controls multiplicity at
both the within-phase and across-phase levels.

### Choosing the contrast correction

**Dunnett's Test (recommended)**
Designed specifically for comparing multiple treatment groups against a single
reference. It uses the exact joint distribution of all treatment-vs-reference
contrasts, which accounts for their positive correlation (they all share the
reference group mean). This makes Dunnett's more powerful than Bonferroni or
Holm for the many-to-one design while controlling the probability of even one
false positive across those contrasts below 5%.
*Use when:* You have a clear reference/vehicle control and all hypotheses are
of the form "is treatment X different from the control?" — which is the standard
pharmacology or genetic rescue design.

**Dunnett-Šidák**
Achieves the same goal as Dunnett's — all contrasts are against the reference
only, and the family-wise false positive rate is kept below 5% — but uses the
Šidák inequality instead of the exact multivariate t-distribution. In practice
the difference is negligible for most datasets (Dunnett-Šidák is very slightly
more conservative). Appropriate when you need a widely auditable correction
that is still anchored to the reference group.
*Use when:* You want a recognized alternative to Dunnett's for reporting or
review purposes.

### Choosing the global correction

**Benjamini-Hochberg (BH) — False Discovery Rate (default)**
Among all contrasts you call significant across all Groups, at most 5% are
expected to be false positives *on average*. The probability of making even
one false positive is typically higher than 5%, but BH is substantially more
powerful than FWER methods — it will detect more true effects. The trade-off
is that a small fraction of your significant hits may be false, and follow-up
validation is important.
*Use when:* You are exploring a large number of comparisons and expect real
effects to exist. Appropriate for initial discovery analyses and when
false positives can be followed up experimentally.

**Holm (step-down Bonferroni) — Family-Wise Error Rate**
The probability of making *even one* false positive across the entire table
is kept below 5%. This is substantially more conservative than BH — you will
miss more real effects, but nearly every result you call significant will be
a true positive.
*Use when:* A single false positive has serious consequences, such as claiming
a drug is efficacious when it is not, or when the analysis is part of a
pre-specified confirmatory study.

---

## CATCH22/24 Feature Reference

CATCH22 is a set of 22 time-series features selected from the hctsa feature
library to be maximally non-redundant and broadly informative across diverse
datasets. CATCH24 adds two simple features (mean and variance) to make 24.

In the app, features are extracted per animal per phase, then z-scored using
a RobustScaler (median and IQR-based) to reduce sensitivity to outliers.

The clustermaps use **Chebyshev distance** (maximum absolute difference across
all features) and **average linkage** hierarchical clustering, applied separately
to both rows (features) and columns (conditions).

### Feature dictionary

| Feature name | What it captures in zebrafish behavior |
|---|---|
| `DN_HistogramMode_5` | The most common activity level during the phase. High values mean the fish spent most of its time at a higher baseline of movement. |
| `DN_HistogramMode_10` | Same idea as above but with finer resolution — distinguishes between, e.g., a fish that rests at near-zero vs. one that idles with low-level movement. |
| `CO_f1ecac` | How quickly a fish's activity becomes unpredictable from one moment to the next. A short memory means rapid, erratic transitions between movement and rest; a long memory means sustained, smooth bouts. |
| `CO_FirstMin_ac` | The dominant rhythm of the fish's movement — the time interval at which activity tends to repeat. Captures regular swim-rest oscillations if they exist. |
| `CO_HistogramAMI_even_2_5` | Whether the fish's current activity level tells you something about what it was doing 2 seconds ago beyond what a simple correlation would suggest. Picks up on structured behavioral patterns (e.g., stereotyped swim bursts) that linear statistics would miss. |
| `CO_trev_1_num` | Whether the fish tends to speed up faster than it slows down (or vice versa). Asymmetry here can reflect startle-then-freeze patterns or gradual onset/abrupt offset of swimming bouts. |
| `CO_Embed2_Dist_tau_d_expfit_meandiff` | How complex or structured the overall pattern of movement is. Simple behavior (steady swimming or steady rest) scores low; complex, multi-state behavior (switching between distinct activity modes) scores high. |
| `IN_AutoMutualInfoStats_40_gaussian_fmmi` | The time scale of the fish's most structured behavioral rhythm, measured in a way that captures non-linear repeating patterns. Similar to `CO_FirstMin_ac` but sensitive to more subtle periodicities. |
| `MD_hrv_classic_pnn40` | How frequently the fish makes sudden large changes in activity from one second to the next — essentially, how "bursty" the movement is. A high value means frequent sharp transitions between rest and vigorous swimming. |
| `SB_BinaryStats_mean_longstretch1` | The average duration of active bouts (periods where the fish is moving above its typical level). Longer stretches indicate sustained swimming; shorter stretches indicate fragmented, start-stop locomotion. |
| `SB_BinaryStats_diff_longstretch0` | The longest continuous period where the fish's activity did not change — essentially the longest freeze or rest episode during the phase. |
| `SB_MotifThree_quantile_hh` | How predictable the transitions between high-activity states are. Low entropy means the fish follows a regular high-activity pattern; high entropy means the timing of active bouts is disorganized. |
| `SB_TransitionMatrix_3ac_sumdiagcov` | How likely the fish is to stay in the same behavioral state (active, moderate, or resting) from one moment to the next rather than switching. High values indicate stable, persistent behavioral states; low values indicate rapid, unpredictable switching. |
| `SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1` | Whether activity fluctuations are structured across long time scales. A high value means the fish has slow, coordinated waves of activity over tens of seconds — not just random noise. Often elevated in animals with altered arousal regulation. |
| `SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1` | An alternative measure of long-range temporal structure, complementary to the feature above. Together they indicate whether the fish's locomotor pattern has a "fractal" quality (self-similar at different time scales) vs. being purely random. |
| `SP_Summaries_welch_rect_area_5_1` | How much of the fish's movement energy is in slow oscillations (cycles longer than a few seconds) vs. rapid fluctuations. High values indicate slow, rhythmic patterns; low values indicate fast, erratic movement. |
| `SP_Summaries_welch_rect_centroid` | The dominant speed of activity fluctuations — whether the fish's movement varies mostly on a scale of seconds (fast centroid) or tens of seconds (slow centroid). A drug that induces tremor-like behavior would shift this toward faster frequencies. |
| `FC_LocalSimple_mean1_tauresrat` | How predictable the fish's next movement is from its recent activity. Low values mean the fish behaves consistently; high values mean its activity is erratic and hard to forecast even from the immediate past. |
| `FC_LocalSimple_mean3_stderr` | Similar to above but over a slightly longer window. Captures how variable the fish's short-term behavioral patterns are — high values indicate inconsistent or disorganized locomotion. |
| `DN_OutlierInclude_p_001_mdrmd` | The structure of the fish's most extreme bursts of activity. Characterizes how the highest-intensity movements are distributed — whether extreme bursts are isolated spikes or part of a broader high-activity state. |
| `DN_OutlierInclude_n_001_mdrmd` | The structure of the fish's deepest rest periods. Characterizes whether the most extreme inactivity occurs as brief pauses or as prolonged, deep freezing episodes. |
| `PD_PeriodicityWang_th0_01` | Whether the fish shows a dominant repeating cycle in its locomotion, and how long that cycle is. Detects regular behavioral rhythms such as periodic swim-rest oscillations, which can be altered by drugs or genetic background. |
| `mean` *(CATCH24 only)* | The average activity level across the entire phase. Directly reflects overall locomotor output — a drug that suppresses movement will lower this; a hyperactivity phenotype will raise it. |
| `variance` *(CATCH24 only)* | How much the fish's activity fluctuates around its average during the phase. High variance means the fish alternates between very active and very quiet periods; low variance means consistent activity (whether consistently active or consistently still). |

### When to turn off mean and variance (CATCH24 → CATCH22)

The mean and variance features measure *how much* a fish moves on average. If your
conditions differ primarily in overall activity level, these two features will 
dominate the clustermap and drive the clustering entirely, potentially obscuring 
differences in the *temporal dynamics* of movement — for example, whether fish 
show longer bouts of activity, more oscillatory patterns, or stronger adaptation 
to light transitions.

**Use CATCH24 (default) when:**
- You want an overall behavioral fingerprint that captures both level and dynamics
- Conditions may differ in mean activity and you want that reflected in the clustering

**Use CATCH22 (uncheck the box) when:**
- You are specifically investigating temporal pattern differences (periodicity,
  autocorrelation structure, transition dynamics) rather than magnitude

---

## Output Files

### GAMM Analysis

All files are written to the output directory you choose when running the analysis.

| File | Description |
|------|-------------|
| `forest_plot_drug_effect.pdf` | Forest plot of drug effect sizes (treatment vs. reference drug) for each Group × Genotype. Points are colored by corrected significance. |
| `forest_plot_genotype_effect.pdf` | Forest plot of genotype effect sizes (treatment vs. reference genotype) for each Group × Drug. |
| `forest_plot_interactions.pdf` | Forest plot of key Genotype × Drug interaction contrasts for each Group. |
| `heatmap_diverging_v1.pdf` | Heatmap of signed −log₁₀(corrected p-value) × sign(effect). Blue = decreased activity vs. reference; red = increased. Color intensity = significance. Effect size printed inside significant cells. |
| `rescue_assessment_context_faceted.pdf` | Line plots of effect sizes over Groups for key comparisons involving the reference condition. Useful for assessing whether a compound "rescues" the disease phenotype toward the wild-type. |
| `summary_statistics_by_group.csv` | Machine-readable summary: per Group, number of significant contrasts and mean absolute effect size, organized by test family (Drug, Genotype, Interaction). |

### Catch22 Clustering

Clustermaps and top-driving-feature reports are displayed inline. Use the
**"Save All"** button in the figure viewer to export all phase clustermaps as
PNG files to a directory of your choice.
