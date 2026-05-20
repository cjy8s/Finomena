# Finomena: Zebrafish Behavioral Time-Series Analysis Suite

A desktop application for processing, visualizing, and statistically analyzing
zebrafish locomotor data captured by the Zebrabox system. Finomena reads raw
Zebrabox output files, lets you define your experimental design, and runs a
Tweedie-family Generalized Additive Mixed Model with AR1 correlation correction
(BAM) — paired with optional Bayesian posterior equivalence testing — to
quantify how each condition's behavioral trajectory differs from the controls.

The app ships as a single self-contained installer for Windows and macOS. End
users do **not** need Python, R, or any other dependencies installed.

---

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Zebrabox Input Data Format](#zebrabox-input-data-format)
4. [Naming Conditions — The `+` Convention](#naming-conditions--the--convention)
5. [Tab-by-Tab Walkthrough](#tab-by-tab-walkthrough)
   - [Tab 1: Experimental Design](#tab-1-experimental-design)
   - [Tab 2: Data Loading](#tab-2-data-loading)
   - [Tab 3: Time Series BAM](#tab-3-time-series-bam)
6. [Why This Statistical Approach?](#why-this-statistical-approach)
7. [P-Value Corrections — What They Mean and When to Use Each](#p-value-corrections--what-they-mean-and-when-to-use-each)
8. [Posterior Equivalence Testing — Detailed Reference](#posterior-equivalence-testing--detailed-reference)
9. [Output Files](#output-files)
10. [Building from Source](#building-from-source)
11. [Troubleshooting](#troubleshooting)

---

## Installation

### For end users (no Python or R required)

| Platform | Download | What it does |
|---|---|---|
| Windows  | `Finomena_Setup_Windows.exe` from [GitHub Releases](https://github.com/cjy8s/Finomena/releases) | Wizard installer. Installs to `%LocalAppData%\Finomena`. **No admin rights required.** Adds a Start Menu shortcut and an uninstaller. |
| macOS    | `Finomena_macOS.dmg` from [GitHub Releases](https://github.com/cjy8s/Finomena/releases) | Drag `Finomena.app` to your Applications folder. |

Both installers bundle R 4.5.3 (downloaded from CRAN at build time) and every
required R package. Nothing else needs to be installed on the target machine.

### For developers / building from source

See the [Building from Source](#building-from-source) section near the end.

---

## Quick Start

If you just want to know the minimum sequence of clicks to get from raw data
to results:

1. **Tab 1 → Experimental Conditions:** add your variable names and levels;
   the app generates one condition per combination (e.g., `WT+DMSO`,
   `MUT+drugA`).
2. **Tab 1 → Metadata Assignment:** mark one condition as **Reference Control**
   and (optionally) one as **Positive Control**. Everything else stays
   **Experimental**.
3. **Tab 1 → Plate Format:** paint your plate wells with the conditions. Save
   the layout if you want to reuse it.
4. **Tab 1 → Experimental Setup:** define the temporal phases of your
   experiment and assign each row a Group integer (one model is fit per Group).
5. **Tab 1 → Save Config…** (recommended) — exports the entire design to JSON.
6. **Tab 2 → Load Data:** add one directory per plate replicate and click
   **Load & Process Data**.
7. **Tab 2 → Set Output Directory:** pick where the results should go.
8. **Tab 2 → Sanity Check:** verify the processed table looks right.
9. **Tab 3 → Family Selection (Optional):** click **Run** to test which
   distribution (Tweedie / Gamma / NegBinomial) fits best. If you skip this,
   Tweedie is used by default.
10. **Tab 3 → Contrast Selection:** the defaults already include every
    *Experimental vs Control* pair. Untick any pair you don't care about.
    Optionally turn on **Posterior Equivalence Testing** here.
11. **Tab 3 → BAM Analysis → Run R Analysis:** wait for the script to finish
    (typically 5–30 minutes depending on data size). Figures appear inline
    when complete.

---

## Zebrabox Input Data Format

The app reads the raw tab-delimited `.csv` (or `.xls`/`.xlsx`) files exported
directly from the Zebrabox ViewPoint software. Each file must contain at
least these three columns (exact names, case-sensitive):

| Column     | Description |
|------------|-------------|
| `time`     | Timestamp in **microseconds** since recording start |
| `location` | Well identifier string, e.g. `"location_1"`, `"location_48"` |
| `data1`    | Raw pixel-difference value for that time point and well |

Rows with `time <= 0` are automatically discarded.

**Directory structure.** Each directory you add in the Data Loading tab
represents one **plate replicate**. The app recursively searches every
subfolder for `.csv` / `.xls` / `.xlsx` files and treats all files found
within a single directory tree as belonging to the same plate. You can have
one file per well, one file per plate, or any combination — the app
concatenates them.

---

## Naming Conditions — The `+` Convention

Each condition is a combination of one or more **variable values** joined
with the `+` character. You define your variables once in the Experimental
Conditions sub-tab; the app generates every valid combination as a condition.

```
Two-variable design (Genotype × Drug):
    Variables: Genotype = {WT, MUT}, Drug = {DMSO, drugA, drugB}
    Conditions: WT+DMSO, WT+drugA, WT+drugB,
                MUT+DMSO, MUT+drugA, MUT+drugB

Three-variable design (Genotype × Drug × Pretreatment):
    Conditions: WT+DMSO+naive, WT+drugA+naive, …, MUT+drugB+stressed
```

The number of variables is set in the Experimental Conditions sub-tab and can
be 1, 2, 3, or more.

### Rules

- **Use exactly one `+` between adjacent variable values.** A condition with
  N variables has exactly N−1 `+` characters.
- **Do not use `+` inside a single variable value.** A variable value of
  `"high+dose"` will be split incorrectly. Use a different separator (e.g.,
  `high_dose`).
- **Be consistent with casing.** `WT` and `wt` are treated as different
  values. Pick a casing convention and stick with it across all plates.
- **Avoid leading/trailing whitespace.** The app strips it, but it's easy to
  accidentally introduce mismatches when typing.

The variable names you choose (e.g., `Genotype`, `Drug`) become the predictor
columns in the statistical model and the labels on every output figure.

---

## Tab-by-Tab Walkthrough

### Tab 1: Experimental Design

Complete this tab **before** loading data. The configuration here drives
every downstream step.

The four sub-tabs are arranged in the order you should fill them in.

#### 1. Experimental Conditions

Define your experiment's variables and the levels each variable can take.

- Click **Add Variable** for each independent variable in your design (e.g.,
  add a `Genotype` variable with levels `WT, MUT`, then add a `Drug` variable
  with levels `DMSO, drugA, drugB`).
- The app automatically generates every combination of variable levels as a
  condition (e.g., `WT+DMSO`, `WT+drugA`, …, `MUT+drugB`).
- Each generated condition has a **color**; click the color swatch to change it
  (used to paint wells on the plate layout).

#### 2. Metadata Assignment

Assign a **role** to each condition. Roles are used to (a) auto-derive the
reference levels passed to the BAM model and (b) drive the default
*Experimental vs Control* contrast selection.

| Role | Abbrev. | What it does |
|---|---|---|
| Reference Control | RC  | The baseline that every effect is measured against in the BAM. Exactly one condition must be assigned this role; the app auto-promotes the first row if you forget. |
| Positive Control  | PC  | A "known responder" condition (e.g., a validated reference compound). |
| Experimental      | EXP | The default for any condition you are testing. |
| Potential Rescue  | RES | Optional label for conditions you predict will rescue a mutant phenotype back to WT. |

The app uses these roles in two specific ways:

1. The Reference Control's variable values become the **reference levels** in
   the BAM (e.g., if `WT+DMSO` is the Reference Control, then `WT` becomes
   the reference level for `Genotype` and `DMSO` becomes the reference level
   for `Drug`).
2. In the **Contrast Selection** sub-tab (Tab 3), pairs are checked by
   default if and only if one side is Experimental and the other is a control
   (Reference or Positive). All other pairs default to unchecked.

#### 3. Plate Format

1. Pick a plate format (`6-well`, `12-well`, `24-well`, `48-well`, `96-well`,
   `384-well`).
2. Click a condition button on the left to make it active.
3. Paint wells on the grid by clicking individual wells, column headers, or
   row headers.
4. Click **Clear Plate Assignments** to reset.
5. Unpainted wells are treated as empty and excluded from analysis.

**Layout library.** You can save multiple named plate layouts (e.g.,
`Plate_1_design`, `Plate_2_with_swap`) using the layout dropdown at the top.
When you load multiple plate replicates in Tab 2, you assign a saved layout
to each replicate independently — useful if your replicates have different
plate maps.

#### 4. Experimental Setup

Define the temporal structure of your experiment. Each row is one **phase**:

| Column   | Meaning |
|----------|---------|
| End (s)  | The time in seconds at which this phase ends |
| Phase    | A descriptive name (`Baseline`, `Dark 1`, `Light 1`, …) |
| Group    | An **integer** that groups related phases into a single statistical model. All rows sharing the same Group number are modeled together. |

The Length column is computed automatically. You can edit Length or Start
directly and the others update.

**Group numbers matter.** All phases that share the same Group integer are
pooled into a single contiguous time window and fitted with one BAM. Group
phases by **stimulus modality** — all light/dark phases together, all color
stimuli together, all vibration stimuli together — so that each model
captures a coherent behavioral context. See
[Why split by Group](#why-split-by-group) below for the statistical
rationale.

A realistic three-modality experiment:

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

#### Save / Load Config

Use **Save Config…** and **Load Config…** at the top of the tab to persist
the entire experimental design — phases, conditions, roles, plate layouts,
and contrast selection — to a single JSON file. Reload it next session to
skip Tab 1 entirely.

---

### Tab 2: Data Loading

#### 1. Load Data

1. Click **Add Directory…** to add one directory per plate replicate. The
   order in the list determines the plate number (top = Plate 1).
2. Use **Move Up** / **Move Down** to reorder if needed.
3. For each directory, pick which saved plate layout applies (from the
   library you built in Tab 1).
4. Click **Load & Process Data**.

The app will:

- Recursively scan all files in each directory
- Parse `time`, `location`, and `data1` columns
- Map each well's `location` number to its plate coordinate (e.g., `A1`–`H12`)
  and the condition you painted in the corresponding plate layout
- Convert microsecond timestamps to seconds
- Sum pixel-difference values per well per second
- Assign Phase and Group labels via a temporal join against the phase table
- Split each condition string into its component variable values

Progress is reported per-file in the log.

#### 2. Set Output Directory

Pick the directory where all R results will be saved. Required before you
can run the BAM. The directory you pick is also passed to the Family
Selection sub-tab so its output lands alongside the BAM results.

#### 3. Sanity Check

After loading, the processed table appears here. Verify that:

- The `Condition` column has no unexpected blanks or `NaN` values.
- `Phase` and `Group` columns are populated for all rows.
- `plate` numbers match the number of directories you added.
- Row count is reasonable (e.g., a 96-well plate × 60 min × 1-second bins
  ≈ 3,600 rows per non-empty well).

Click column headers to sort. The data does **not** propagate to the
analysis tabs until you scroll through and confirm it looks correct (the
viewer emits the `data_passed_through` signal once you've engaged with it).

---

### Tab 3: Time Series BAM

This is the analysis tab proper. It has three sub-tabs that are intended to
be run in order, although Family Selection is optional.

#### 1. Family Selection (Optional)

Compares three candidate response distributions on a simplified version of
the model:

| Family       | When it fits best |
|---|---|
| **Tweedie**     | Pixel-difference data with genuine zero-activity periods (most common). |
| **Gamma (log link)** | Continuous, strictly positive data with minimal zeros (e.g., already-aggregated activity). A small additive shift is applied to ensure positivity if needed. |
| **Negative Binomial** | Discrete count-like data with overdispersion. |

Click **Run Family Selection**. The script fits each family with a fixed
ρ = 0.25 (for speed), reports AIC / BIC / dispersion / deviance explained /
zero-proportion match, and recommends a winner. The winner is automatically
passed to the BAM Analysis sub-tab.

If you skip this step, **Tweedie** is used by default — the right call for
nearly every Zebrabox dataset.

#### 2. Contrast Selection

Pick which **pairwise condition contrasts** to estimate in the BAM.

By default, the app pre-checks every pair where one side is **Experimental**
and the other is a **control** (Reference Control or Positive Control). All
other pairs (Experimental vs Experimental, control vs control, anything
involving Potential Rescue) default to **unchecked**. Tick or untick anything
you want to override the defaults.

**Why the default?** With many conditions, the pairwise comparison count is
combinatorial (C(N, 2)). Statistical power evaporates if you ask the
multiple-testing correction to handle pairs you never actually cared about.
The defaults reflect the typical research question — *how does each treatment
compare to the controls?* — and exclude pairs that are usually noise.

**Posterior Equivalence Testing (Bayesian).** Optional checkbox at the bottom
of this sub-tab. If enabled, the script will compute a posterior distribution
of the maximum log₂ fold change between each kept pair's trajectories, and
report Pr(M < δ) — a direct positive equivalence statement. See
[Posterior Equivalence Testing](#posterior-equivalence-testing--detailed-reference)
for the math and interpretation.

| Setting | Default | Meaning |
|---|---|---|
| δ (log₂ fold change) | 1.0 | Equivalence margin: trajectories never differ by more than 2× at any time point. |
| Posterior draws | 10,000 | Number of samples used to estimate the M distribution. ~0.5% Monte Carlo precision. Lower = faster, noisier. |

If a BAM has already been run, a red banner appears here ("BAM has been run.
Changing contrasts here will require re-running the BAM analysis.") to remind
you that any change invalidates the prior results.

#### 3. BAM Analysis

The main event.

##### 3a. Statistical Correction Methods

Set both correction methods **before** running. See
[P-Value Corrections](#p-value-corrections--what-they-mean-and-when-to-use-each)
for full guidance.

- **Global Multiple-Testing Correction** (across all Groups × contrasts):
  Benjamini-Hochberg (FDR, default) or Holm (FWER).
- **Contrast Correction** (within each emmeans contrast call): Dunnett or
  Dunnett-Šidák.

##### 3b. Family Display

The currently selected distribution family is shown here (italic gray text
above the Run button). It updates automatically when Family Selection runs;
otherwise it shows "Tweedie — default."

##### 3c. Run R Analysis

Click **Run R Analysis** to launch the R script as a subprocess. The output
log streams live in the panel below. Typical runtime for a 96-well plate ×
60 min × 8 phase groups is **5–30 minutes** depending on CPU.

If the contrast selection has changed since the last run, a red banner above
the Run button reminds you that prior figures are stale.

The **Stop** button kills the entire R process tree (including any worker
processes spawned by `parallel::`). Closing the app also kills any in-flight
R subprocess automatically.

##### 3d. R Output Log

Raw stdout/stderr from Rscript. Key lines to watch for:

- `Optimal Rho: <value>` per Group — should be in [0, 1]; values near 1 mean
  strong autocorrelation (typical for resting fish).
- `[Group N] Estimated grid: <rows> × <coefs> ≈ <MB>` — the size and memory
  cost of the emmeans reference grid for that Group. Mostly informational;
  the script auto-bumps `rg.limit` to fit.
- `Successfully fitted Group <N>!` — model fit completed.
- `Contrast selection sidecar found: <N> pair(s) requested.` — confirms the
  contrast list reached R.
- Any `Warning:` from mgcv (most are harmless; "k-index < 1" suggests a
  smooth might benefit from a higher basis dimension).
- `Analysis complete. All outputs saved to: <dir>` — done.

##### 3e. Results

PDF/PNG figures are rendered inline after the run completes. Use the arrow
buttons to navigate; **Save All** copies them to a directory of your choice.

---

## Why This Statistical Approach?

Behavioral time-series data from zebrafish are unlike most biology datasets
in three crucial ways:

1. Measurements are **densely repeated** within each animal (one row per
   second, hundreds of rows per fish per phase).
2. The responses are **non-linear in time** (e.g., light-on responses are not
   straight lines but curves with rapid acceleration and gradual decay).
3. There is **massive zero inflation** (resting fish contribute long runs of
   zero pixel difference).

Naïve statistical methods — t-tests per phase, repeated-measures ANOVA, even
non-parametric tests — make assumptions that all three of these properties
violate. Finomena uses a Generalized Additive Mixed Model fitted via
`mgcv::bam()` with AR1 correlation and a Tweedie response distribution,
paired with optional Bayesian posterior equivalence testing. The rest of
this section explains why.

### The problem with pairwise testing — even with repeated-measures correction

A common analysis is "average each fish's activity within each phase, run a
t-test (or Wilcoxon) per phase, Bonferroni-correct across phases." This
sounds reasonable but has several deep problems:

1. **Aggregating to a per-phase mean discards the time information that
   defines the phenotype.** A drug that causes a *delayed* response and one
   that causes a *blunted* response can produce identical per-phase means
   but reflect different mechanisms. Time-series data exists precisely so
   you can distinguish these — averaging it away erases the very signal you
   paid to collect.

2. **Repeated-measures ANOVA assumes sphericity** — that within-subject
   variances and covariances follow a specific structure. Time-series data
   essentially never satisfies sphericity (consecutive measurements are far
   more correlated than distant ones). Greenhouse-Geisser corrections are
   crude band-aids, not solutions.

3. **Pairwise testing at every time point produces enormous multiplicity
   costs.** With 300 seconds × 5 phases × 6 conditions, you have thousands
   of comparisons. Bonferroni-correcting across this hyperinflates false
   negatives — most "significant" differences disappear, and most real
   differences are missed.

4. **Linear models cannot capture nonlinear dynamics.** Zebrafish startle,
   accelerate, then exponentially decay back to baseline. No t-test or
   linear regression captures the *shape* of this trajectory, only point
   estimates of differences in level.

5. **Pixel-difference data is non-negative and zero-inflated.** A
   Gaussian-based test (t-test, ANOVA) assumes normally distributed residuals
   around the mean. Pixel-difference values cluster at zero (resting fish)
   with a long right tail (occasional bursts). A Gaussian test will
   systematically misestimate the variance and produce wrong p-values.

### Why a Generalized Additive Mixed Model

The BAM/GAMM approach addresses each of those problems directly:

- **Smooth terms** — `s(time_in_group)` and `s(time_in_group, by = Condition_Combo)`
  fit a flexible curve through the time axis for each condition, *without*
  assuming a specific shape. The model learns the trajectory from the data
  using penalized regression. Curves can have rapid rises, plateaus, and
  slow decays — all in one fit.

- **Per-animal random smooths** — `s(time_in_group, animal_id, bs = 'fs')`
  treats each fish as having its own unique trajectory drawn from a
  population distribution. This properly handles between-fish variability
  *without* the sphericity assumption that ANOVA requires. Each animal
  contributes information about both the population effect and its own
  deviation from it.

- **AR1 correlation** explicitly models the correlation between consecutive
  residuals within each fish. The model fits the autocorrelation parameter
  ρ from the data (separately for each Group, since arousal-state
  autocorrelation differs from baseline) and uses it to compute valid
  standard errors. Without this, the model would effectively double-count
  each fish's time series.

- **Tweedie family** is a compound Poisson-Gamma distribution that natively
  accommodates a mass at zero (resting fish) and a continuous positive tail
  (active fish), in a single response distribution. No log-transformation
  needed, no zero-handling tricks.

- **Effects are tested as smooth differences over time**, not as single
  scalar contrasts. The result for "drug X vs reference at Genotype WT in
  Phase Group 1" is a *trajectory* of effect sizes — the difference between
  the WT+drugX smooth and the WT+reference smooth as a function of time.
  emmeans aggregates this into a forest-plot contrast for reporting, but the
  underlying model carries the full time-resolved information used by the
  posterior equivalence test.

The net effect is a model that captures the structure of the experiment,
properly accounts for the dependency structure of the data, and yields
valid (not inflated, not deflated) standard errors and p-values.

### Why pair it with Bayesian posterior equivalence testing

A frequentist p-value can only **reject** the null hypothesis. It cannot
**support** it. Failing to reject (p > 0.05) is *not* evidence that two
conditions are the same — it's just absence of evidence for them being
different. This becomes a critical limitation when you want to make claims
like:

- "This drug rescues the mutant phenotype to wild-type levels"
- "This compound has no off-target effect on wild-type fish"
- "Treatment X is bioequivalent to treatment Y"

All of these are **equivalence claims** — they require demonstrating
sameness, not failing to demonstrate difference. Pairwise NHST cannot do
this. You can only say "we found no significant difference," which is a
weak and frequently-misinterpreted statement.

The optional **posterior equivalence test** addresses this directly. For
each kept pair of conditions, the model's posterior (every coefficient is
sample-able via mgcv's Bayesian interpretation of penalized smooths) is used
to draw thousands of plausible trajectory pairs. For each draw, the script
computes **M** — the maximum absolute log₂ fold change between the two
condition trajectories across the entire phase time window. The result is a
full posterior distribution of M, summarised as **Pr(M < δ)** — the posterior
probability that the trajectories never differ by more than a factor of δ
(default δ = 1.0 → factor of 2) at any time point.

This statistic is interpretable directly: "Pr(M < 1.0) = 0.92" means there
is a 92% posterior probability that the two trajectories never differ by
more than 2× anywhere in the phase window. That is a *positive* equivalence
statement, not a non-rejection of the null, and it uses **every time point**
in the data — exactly the strength a pairwise test cannot offer.

### Why split by Group

<a name="why-split-by-group"></a>Each Group receives its own independently
fitted model because the behavioral responses zebrafish produce under
different stimulus modalities are fundamentally distinct — not just
quantitatively different, but qualitatively incomparable in ways that make a
single joint model biologically inappropriate.

1. **Different stimulus types produce data that are not statistically
   comparable.** A single model estimates one set of parameters
   (dispersion, zero-inflation rate, smooth basis weights) that are assumed
   to apply consistently. Mixing light/dark cycles with vibration trials in
   one model produces a compromise that fits none of them well; the model
   spends most of its flexibility bridging the jumps between modalities
   instead of detecting differences within them. Splitting by Group means
   each model's estimates are grounded in data that are actually comparable.

2. **The temporal dynamics within each modality are self-contained.** Within
   a Group, the fish are adapting, habituating, or sensitizing to stimuli of
   the same type over a coherent time window. The smooth terms capture that
   arc meaningfully. Pooling across modalities forces the smooth to bridge
   entirely different behavioral states with no interpretable biology.

3. **The autocorrelation structure differs by modality.** The AR1 ρ
   parameter is estimated fresh for each Group. Phases with high arousal
   (e.g., vibration response) have different within-animal serial
   correlation than baseline. A single global ρ would misrepresent the
   dependence structure in every Group simultaneously.

The consequence is that the same biological hypotheses (drug effect,
genotype effect, interaction) are tested independently in each Group — which
is the correct thing to do scientifically, but creates a multiple-testing
problem that the two-stage p-value correction is designed to address (see
next section).

### When this matters (and when it doesn't)

This approach is overkill for a quick pilot with 4 conditions and 8 fish — a
basic t-test will tell you whether the assay is working. It becomes
essential when:

- You have **dose-response or large factorial designs** where Tweedie +
  smooth + AR1 properly controls each comparison's variance across the whole
  dataset.
- You are making **rescue or equivalence claims** — posterior equivalence is
  the only honest way to defend "the drug restores the wild-type
  phenotype."
- You have **complex stimulus protocols** (light, dark, vibration, all in
  one experiment) where Group-wise modeling captures the modality-specific
  dynamics.
- You will publish the results in a peer-reviewed venue where reviewers will
  ask why you used a t-test and whether you have considered autocorrelation.

If you only need a sanity check on whether a drug "does anything,"
practically any test will work. If you need to defend specific scientific
claims, the model in this app is the right tool.

---

## P-Value Corrections — What They Mean and When to Use Each

### The two-stage correction structure

Every p-value in the final BAM output has been corrected **twice**, at two
different levels:

**Stage 1 — Contrast correction (within each Group):** Dunnett's test (or
Dunnett-Šidák) is applied inside each `emmeans` contrast call. Within a
single Group, you are comparing each treatment condition against the
reference simultaneously. For example, with 3 drug conditions, there are 3
drug contrasts in that Group, all using the same DMSO reference. Dunnett's
exploits the shared correlation structure (all contrasts share the
reference) to be more powerful than Bonferroni while still controlling the
within-Group false-positive rate.

**Stage 2 — Global correction (across all Groups):** After Stage 1, all
contrast p-values from all Groups are pooled and a second round of
correction is applied across the entire table. This controls for the fact
that the same biological hypothesis is tested in every Group's separate
model. Without this second correction, with 8 independent tests you would
expect roughly a 34% chance of at least one false positive at p < 0.05 even
if nothing were truly significant.

This sequential design is conservative by intent: it controls multiplicity
at both the within-phase and across-phase levels.

### Choosing the contrast correction

**Dunnett's Test (recommended).** Designed specifically for comparing
multiple treatment groups against a single reference. It uses the exact
joint distribution of all treatment-vs-reference contrasts, accounting for
their positive correlation (they all share the reference group mean). This
makes Dunnett's more powerful than Bonferroni or Holm for the many-to-one
design while controlling the family-wise false-positive rate below 5%.
*Use when:* you have a clear reference/vehicle control and all hypotheses
are of the form "is treatment X different from the control?" — the standard
pharmacology or genetic-rescue design.

**Dunnett-Šidák.** Achieves the same goal as Dunnett's but uses the Šidák
inequality instead of the exact multivariate t-distribution. In practice the
difference is negligible (Dunnett-Šidák is very slightly more conservative).
*Use when:* you want a recognised alternative to Dunnett's for reporting
or review purposes.

### Choosing the global correction

**Benjamini-Hochberg (BH) — False Discovery Rate (default).** Among all
contrasts you call significant across all Groups, at most 5% are expected to
be false positives *on average*. BH is substantially more powerful than
FWER methods — it will detect more true effects — at the cost of allowing a
small fraction of significant hits to be false. Follow-up validation is
important.
*Use when:* exploring a large number of comparisons, expecting real effects,
and able to follow up false positives experimentally.

**Holm (step-down Bonferroni) — Family-Wise Error Rate.** The probability of
making *even one* false positive across the entire table is kept below 5%.
Substantially more conservative than BH — you will miss more real effects,
but nearly every result you call significant will be a true positive.
*Use when:* a single false positive has serious consequences (e.g., claiming
a drug is efficacious when it is not), or when the analysis is part of a
pre-specified confirmatory study.

---

## Posterior Equivalence Testing — Detailed Reference

This section assumes you've enabled posterior equivalence in the Contrast
Selection sub-tab. If you skip it, you can ignore this section entirely.

### The statistic: M

For each kept condition pair (A, B):

1. Draw `n_draws` samples (default 10,000) from the model's posterior.
2. For each sample, compute the predicted trajectory of A and the predicted
   trajectory of B over the full Group time window (excluding the per-animal
   random smooth — we want population-level inference).
3. Compute **M** = max over all time points t of `|log₂(μ_A(t) / μ_B(t))|`.
   This is the largest log-fold-change that occurs anywhere in the trajectory.
4. The posterior distribution of M summarises *how different* the two
   trajectories are at their most divergent point, marginalising over
   parameter uncertainty.

### The threshold: δ

You pick δ — a tolerance for "practically equivalent" expressed in log₂
units:

| δ | Means trajectories never differ by more than | Typical use |
|---|---|---|
| 0.5 | ~1.4× | Very strict. Use only for tightly characterised assays. |
| 1.0 (default) | 2× | Standard. "No more than a doubling/halving." |
| 1.5 | ~2.8× | Loose. Suitable for noisy assays or large dynamic range. |
| 2.0 | 4× | Very loose. Equivalence at this margin is weak evidence. |

### The summary: Pr(M < δ)

The script reports **Pr(M < δ)** for each pair — the posterior probability
that the two trajectories never differ by more than δ log₂ units anywhere in
the phase window. Rough interpretation guide:

| Pr(M < δ) | Interpretation |
|---|---|
| > 0.95 | Strong posterior support for equivalence |
| 0.80 – 0.95 | Moderate support |
| 0.50 – 0.80 | Inconclusive — data is consistent with equivalence but doesn't strongly support it |
| < 0.50 | Posterior favours non-equivalence; the conditions are likely meaningfully different |

A pair with Pr(M < 1.0) = 0.92 supports a defensible claim like:
*"There is a 92% posterior probability that drug X and DMSO never differ by
more than a factor of 2 in WT fish across the dark-light cycle."*

This is the kind of statement reviewers cannot dismiss as "absence of
evidence." It is **positive** evidence, with a quantified probability.

---

## Output Files

All files are written to the output directory you set in the Data Loading
tab.

### From Family Selection

| File | Description |
|---|---|
| `family_selection_results.csv` | Per-family AIC, BIC, dispersion, deviance explained, zero-proportion match — aggregated across phase groups (AIC/BIC summed, dispersion/dev_explained n-weighted). |
| `family_selection_by_group.csv` | Per-group × per-family breakdown of the same metrics. Useful for spotting groups where the family ranking disagrees with the overall winner. |
| `family_selection_winner.csv` | The recommended family + any additive shift applied for Gamma. |

### From the BAM run

| File | Description |
|---|---|
| `finomena_pre-processed_data.csv` | The full processed dataset that was sent to R. Contains `time_sec`, `location`, `loc_coord`, `pixel_diff`, `Condition`, `Phase`, `Group`, `animal_id`, `plate`. |
| `contrast_selection.json` | Sidecar file listing the kept pairs and posterior-equivalence settings. R reads this to know which contrasts to compute. |
| `forest_plot_<varname>_effect.pdf/png` | One forest plot per variable in your design (e.g., `forest_plot_drug_effect.pdf`, `forest_plot_genotype_effect.pdf`). Effect sizes vs the reference level, faceted by the other variable's levels. |
| `forest_plot_interactions.pdf/png` | Forest plot of the kept pairwise interactions across all Groups. |
| `heatmap_diverging_v1.pdf/png` | Heatmap of signed −log₁₀(corrected p) × sign(effect). Blue = decreased vs reference; red = increased. Effect sizes printed inside significant cells. |
| `rescue_assessment_context_faceted.pdf/png` | (2-variable designs only) Line plots of effect sizes over Groups for the comparisons that matter for rescue assessment. |
| `summary_statistics_by_group.csv` | Per Group: number of significant contrasts and mean absolute effect size, broken down by test family. Machine-readable. |

### From Posterior Equivalence (if enabled)

| File | Description |
|---|---|
| `posterior_equivalence_density.png` | Density plots of the posterior distribution of M for each kept pair, with δ marked as a vertical line. |
| `posterior_equivalence_heatmap.png` | Heatmap of Pr(M < δ) for every kept pair × Group combination. |
| `posterior_equivalence_summary.csv` | Machine-readable: per pair × Group, posterior mean of M, posterior SD, Pr(M < δ), and the δ used. |

---

## Building from Source

Most users do not need this section — use the installer above. This section
is for people who want to modify the app or build a fresh installer.

### Prerequisites

| Tool | Purpose | Where to get it |
|---|---|---|
| Python 3.10+ | App runtime | https://python.org or your distro |
| `pip install -r App/requirements.txt` | Python deps | `requirements.txt` lives in `App/` |
| Inno Setup 6 | Windows installer | https://jrsoftware.org/isinfo.php |
| `create-dmg` | macOS installer | `brew install create-dmg` |
| `pyinstaller` | Bundling | `pip install pyinstaller` (already in requirements) |

R itself is **not** needed on the build machine — the build scripts
download R 4.5.3 directly from CRAN.

### One-time dev environment setup

After cloning, run the platform-appropriate dev-setup script *once*. It
downloads R 4.5.3 from CRAN into `App/R/runtime/` (gitignored, ~250 MB) and
installs the required R packages into `App/R/library/` (which IS tracked in
git so the app's bundled library stays version-controlled).

```
# Windows
App\setup_dev.bat

# macOS
bash App/setup_dev.sh
```

After this, you can run the app in dev mode with system Python:

```
python App/app.py
```

The app will automatically prefer the bundled R at `App/R/runtime/` over any
system R.

### Building a distributable installer

```
cd App

# Windows  → Output\Finomena_Setup_Windows.exe
build_windows.bat

# macOS  → Finomena_macOS.dmg
bash build_mac.sh
```

Each build:

1. Runs PyInstaller to package the Python app into `dist/Finomena/` (or
   `dist/Finomena.app/`).
2. Downloads R 4.5.3 from CRAN (cached in `%TEMP%` / `/tmp` after the first
   run).
3. Extracts R into the dist as `R/runtime/`.
4. Runs `install_packages.R` against the bundled R to populate
   `R/library/`.
5. Wraps the result with Inno Setup (Windows) or `create-dmg` (macOS) into a
   single distributable installer file.

Upload the resulting installer to a [GitHub Release](https://github.com/cjy8s/Finomena/releases)
and the direct download URL from there is what end users use.

### Repo layout

```
App/
  app.py                      Entry point
  finomena.spec               PyInstaller spec
  finomena.iss                Inno Setup script (Windows)
  build_windows.bat           Windows build pipeline
  build_mac.sh                macOS build pipeline
  setup_dev.bat               Windows dev environment setup
  setup_dev.sh                macOS dev environment setup
  requirements.txt            Python dependencies
  Finomena/utils/             All Qt widgets
    paths.py                  Resource path resolution (dev vs frozen)
    bam_widget.py             Tab 3 sub-tab 3 — BAM Analysis
    family_selection_widget.py    Tab 3 sub-tab 1
    contrast_selection_widget.py  Tab 3 sub-tab 2
    data_loader.py            Tab 2
    experiment_design.py      Tab 1 (most sub-tabs)
    plate_format.py           Tab 1 sub-tab 3
    dataframe_viewer.py       Sanity-check viewer
    figure_viewer.py          Reusable image carousel
  R/
    install_packages.R        R-package bootstrap
    library/                  App's bundled R packages (tracked in git)
    runtime/                  Bundled R 4.5.3 (gitignored, populated by setup_dev)
    scripts/
      TweedieAR1 BAM.R        Main BAM + posterior-equivalence script
      FamilySelection.R       Family selection comparison script
```

---

## Troubleshooting

**"Could not locate Rscript" error on app launch.** The app couldn't find a
bundled or system R. Run `App\setup_dev.bat` (Windows) or
`bash App/setup_dev.sh` (macOS) to install the bundled R 4.5.3.

**Mangled characters in the R log (e.g., `Ã—` instead of `×`).** Should not
happen as of v1.0.0 — the subprocess decodes stdout as UTF-8. If you still
see them, your terminal/console is the problem, not the app.

**`The rows of your requested reference grid would be N` from emmeans.** The
app auto-bumps `rg.limit` to fit; this should not happen. If it does, the
diagnostic line right above ("Estimated grid: …") tells you which factors
caused the blow-up.

**BAM run takes hours.** Likely causes (in order of probability):
1. Many conditions × many time points (the smooth basis matrix grows fast).
2. AR1 estimation in Pass 1 is the bottleneck — use fewer phases per Group
   if it's intolerable.
3. Posterior equivalence testing adds ~1–2 minutes per kept pair × Group at
   the default 10,000 draws.

**App closes but R keeps running.** Should not happen as of v1.0.0 — the
app's `closeEvent` kills the subprocess tree. If you see orphaned
`Rscript.exe` processes in Task Manager after closing, please file an issue.

**Saved config loads but contrast selections look wrong.** If your saved
config has an empty `excluded_pairs` list, the app falls through to the
role-based defaults (Experimental vs Control) instead of "include
everything." This is intentional. Use **Select All** in the Contrast
Selection sub-tab if you really want every pair.

---

## License

[Add license text here]

## Contact

Questions, bug reports, feature requests: please open a GitHub issue at
https://github.com/cjy8s/Finomena/issues
