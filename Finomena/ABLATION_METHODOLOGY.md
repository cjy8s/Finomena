# Ablation Methodology — Choosing the Statistical Model for Finomena's BAM Pipeline

This document records how the current BAM model architecture (`Arch C — HGAM`, random-effect basis `sz`, by-variable `interaction(Condition_Combo, Group_factor)`, AR.start mode `phase`) was selected. It is intended to:

1. Make the methodological decision **auditable** — any reviewer can trace the chosen model back to the empirical evidence
2. Document the **trade-offs that were considered** so a future maintainer can re-evaluate if the data or scientific questions change
3. Anchor the **diagnostic tooling** (the "Architecture × Method Ablation" tab in the app) for future re-runs

The ablation was conducted using the `Architecture × Method Ablation` widget in the app. That tab is hidden by default after the analysis was completed but can be re-enabled via the toggle in the Experimental Design config bar.

---

## TL;DR

| Item | Value |
|---|---|
| Chosen architecture | **C — Hierarchical GAM (HGAM)** |
| Random-effect basis | `sz` (sum-to-zero constrained factor smooth over `time_in_group`) |
| By-variable for trajectory smooth | `interaction(Condition_Combo, Group_factor)` |
| AR.start reset mode | `phase` (resets at animal AND phase boundaries) |
| Held constant | `select=TRUE`, `discrete=TRUE`, `rho=0.2` (during ablation), `family=Tweedie` |
| Selection metric | 5-fold animal-stratified cross-validated holdout Tweedie deviance |
| CV deviance, winner | **1.956 × 10⁶** (C_sz_interaction_phase) |
| CV deviance, runner-up | 2.2 × 10⁶ (A_sz, per-phase reference architecture) |
| Margin | ~12% better holdout deviance for C over A |

---

## 1. The problem

The pre-ablation pipeline was producing **monstrously inflated raw p-values** in the BAM contrast outputs. This was a symptom that something was structurally wrong in the model — but inflated p-values themselves cannot be the *target* of method selection. Optimizing models to produce smaller p-values is p-hacking; the methodologically correct approach is to **find a structurally defensible model** and accept whatever p-values it honestly produces.

The question this ablation set out to answer is **not** "which model gives better p-values?" — it is **"which model is the correct statistical machinery for this experimental design?"**

---

## 2. The methodological framework

### Three candidate architectures

| Architecture | Description | Why it was a candidate |
|---|---|---|
| **A — Per-Phase BAMs** | One `bam()` fit per phase; everything (φ, ρ, smooths, random effects) is phase-local | The historical reference. Cleanly separated phases. Known to produce well-calibrated SEs within each phase. |
| **B — Fully Global** | One `bam()` fit across all phases, smooths on `time_sec`, phase enters as a parametric factor only | The "pool everything" baseline. Tests whether phase boundaries can be ignored if the model has enough flexibility. |
| **C — HGAM (Hierarchical GAM)** | One `bam()` fit, fixed-effect smooths on `time_in_group` (phase-local), random effects pool across phases as intercepts (Pedersen et al. 2019) | The methodologically defensible middle ground. Phase boundaries respected for fixed effects, animal identity pooled across phases for random effects. |

### The decision cascade

A three-stage cascade was used to combine cheap descriptive metrics, multi-criteria Pareto pruning, and out-of-sample cross-validation:

```
Stage 1 — fit every (arch × random-effect basis × by-variable × AR.start) variant
          and compute the full per-fit metric panel
   ↓
Stage 2 — within each architecture, apply adequacy filters then compute the
          Pareto front over the surviving variants
   ↓
Stage 3 — run 5-fold animal-stratified cross-validation on the Pareto union
          and rank by mean holdout deviance
```

This structure is borrowed from the hyperparameter-optimization literature (successive halving / Hyperband — Karnin et al. 2013, Li et al. 2017) crossed with classical multi-criteria decision analysis (Pareto efficiency — Pareto 1906; NSGA-II — Deb et al. 2002).

### Why this cascade, not a single-metric ranking

A single composite score (e.g. weighted AIC + dev.explained + SE) would require user-defined weights that are inherently subjective. The Pareto-then-CV cascade avoids subjective weighting:

- **Pareto** identifies the non-dominated set without weights — a candidate stays on the front only if no other variant is at-least-as-good on every metric AND strictly better on at least one
- **CV** then ranks the Pareto front by a single, principled, out-of-sample metric (holdout deviance) — there is no weight to choose

---

## 3. The metric panel (per Stage 1 fit)

Each Stage 1 variant fit produces a `run_metrics.json` with the following:

### Fit quality
- `AIC` (mgcv's EDF-aware AIC — Wood, Pya & Saefken 2016)
- `BIC`
- `fREML` score (`model$gcv.ubre`)
- `logLik`
- `dev_explained`
- `r_sq_adj`
- `tweedie_phi` (residual scale)
- `total_edf`

### Model adequacy
- `convergence_warnings` (mgcv warnings string)
- `k_index` (per-smooth k-check p-values from `mgcv::k.check`; pass rate aggregate)
- `ar1_residual_lag1` (lag-1 ACF of AR-corrected residuals — should be ≈0 if AR(1) captured autocorrelation)

### Identifiability + random-eats-fixed diagnostic
- `concurvity.max_overall` (worst pairwise concurvity, all smooths)
- `concurvity.max_random_vs_fixed` (worst concurvity specifically between random-effect smooths and fixed-effect smooths — the direct symptom of random absorbing fixed signal)
- `variance_components` (from `mgcv::gam.vcomp` — std devs per smooth)
- `edf.random_to_fixed` (sum of random-smooth EDF / sum of fixed-smooth EDF)
- `vcov_diagnostic.condition_num` (condition number of the contrast vcov matrix)

### Contrast variance
- `contrast_variance.<TestFamily>.mean_SE` / `median_SE` / `max_SE` — average estimated SE per Test_Family
- `mean_abs_est` — sanity check that estimates aren't all zero

### Cross-architecture summary
- For Arch A (per-phase) variants, per-group fits are aggregated into one record: AIC/BIC summed, dev.expl/φ averaged, concurvity/EDF taken max-across-phases

---

## 4. Stage 1 — adequacy filters

Variants that catastrophically fail one of the following are excluded from the Pareto front (their scores are still reported, but they cannot be selected). Thresholds are deliberately loose — only egregious failures are filtered:

| Filter | Threshold | Rationale |
|---|---|---|
| Convergence | Any mgcv warning containing "indefinite", "convergence", "boundary" | Hard model-fitting failure |
| `gam.check` k-index pass rate | < 30% of smooths pass | Most smooths basis-exhausted |
| Max pairwise concurvity | > 0.97 | Extreme identifiability failure |
| Max random-vs-fixed concurvity | > 0.95 | Direct random-eats-fixed failure |
| AR(1) residual lag-1 | > 0.7 | AR(1) catastrophically failed to capture autocorrelation |
| Tweedie φ relative to minimum across variants | > 100× | Dispersion catastrophe |

### Architecture B was eliminated at this stage

Every Arch B variant failed at least one adequacy filter. Most failed on **convergence warnings + k-index pass rate < 30%**. This is the expected, documented failure mode: Arch B's smooths over `time_sec` cannot capture the discontinuous protocol boundaries between phases, leaving large unmodeled residual structure that mgcv can't fit cleanly. This was a legitimate model rejection, not a tooling bug.

---

## 5. Stage 2 — Pareto front (within architecture)

After Stage 1's adequacy pruning, the Pareto front was computed within each surviving architecture using these metrics:

| Metric | Direction |
|---|---|
| AIC | minimize |
| BIC | minimize |
| dev_explained | maximize |
| r_sq_adj | maximize |
| max_random_vs_fixed (concurvity) | minimize |

`mean_contrast_SE` was deliberately **excluded** from the Pareto criteria — a candidate that minimizes SE could simply be overconfident, so SE belongs in reporting, not in selection.

### Result

- **Arch A**: one Pareto candidate — `A_sz`
- **Arch B**: no surviving candidates (all filtered at Stage 1)
- **Arch C**: one Pareto candidate — `C_sz_interaction_phase`

A Pareto front of one candidate per architecture is not a deficiency of the algorithm — it means within each architecture there is a single multi-criteria-dominant variant. This is efficient pruning, exactly what Pareto is supposed to do.

The random-eats-fixed diagnostic also did its job: `A_fs` (per-phase BAM with `bs='fs'` factor smooth over `time_in_group`) had an EDF random/fixed ratio of **15.5** and max random-vs-fixed concurvity of 1.0 — direct empirical confirmation that unconstrained factor smooths absorb fixed-effect signal in the per-phase scope. `A_fs` failed adequacy and was correctly excluded.

---

## 6. Stage 3 — Cross-validation

5-fold animal-stratified CV was run on the Pareto union. For each Pareto candidate:
- Animals were randomized into 5 folds with a fixed seed (deterministic across candidates)
- For each fold, the model was refit on 4/5 of animals and predicted on the held-out 1/5 at the **population level** (random-effect smooths excluded from the linear predictor via `exclude=`)
- Holdout Tweedie deviance was computed via the family's `dev.resids()`
- For Arch A, deviance was summed across all phase groups per fold (since A's "model" is the joint per-phase set)

Lower mean holdout deviance = better out-of-sample prediction.

### Result

| Variant | Architecture | Mean CV deviance | SE | n_folds |
|---|---|---|---|---|
| **C_sz_interaction_phase** | C (HGAM) | **1.956 × 10⁶** | 4.469 × 10⁴ | 5 |
| A_sz | A (Per-Phase) | 2.2 × 10⁶ | 8.648 × 10⁴ | 5 |

**C_sz_interaction_phase** is the chosen model — ~12% lower mean CV deviance with a smaller fold-to-fold SE (better stability).

---

## 7. The chosen model — full specification

```r
formula = pxl_diff ~ Var1 * Var2 * ... * Group_factor
                  + s(time_in_group, k = 30)
                  + s(time_in_group,
                      by = interaction(Condition_Combo, Group_factor),
                      k = 30)
                  + s(time_in_group, animal_id, bs = 'sz')
                  # (s(time_in_group, plate, bs='sz') if n_plates >= 2)

bam(formula  = formula,
    data     = gam_df,
    family   = chosen_family,           # Tweedie by default
    rho      = optimal_rho,              # Pass-1 estimated via itsadug start_value_rho-equivalent
    AR.start = gam_df$start_event,       # resets at animal AND phase boundaries
    select   = TRUE,
    method   = "fREML",
    discrete = TRUE,
    nthreads = usable_cores)
```

This is what runs in the main BAM Analysis tab when no ablation env vars are set. It matches the production formula in `App/R/scripts/TweedieAR1 BAM.R`.

---

## 8. Trade-offs and caveats

These are honest acknowledgments of what the chosen model gives up or relies on:

1. **Single global Tweedie φ.** Each phase doesn't get its own dispersion. If phases differ substantially in noise scale, this averaging is a misspecification. Per-phase φ would require Arch A (which lost on CV) or a separately-fit scale GAM (not currently implemented).

2. **Single global AR(1) ρ.** Estimated once (via itsadug's `start_value_rho` equivalent — lag-1 ACF of no-AR residuals) and plugged in as a constant. `bam()` does not jointly estimate ρ. A grid search over ρ is documented in `?bam` but not implemented.

3. **AR(1) on the residuals is a GEE approximation under Tweedie.** mgcv's `bam()` documents that AR(1) is technically supported only for Gaussian-identity models; for Tweedie + `discrete=TRUE` it becomes a working-residual approximation. mgcv source flags this as having "badly biased estimates for low count data" (Fletcher, *Biometrika*) — for continuous `pxl_diff` this is unlikely to bite, but it's a real limitation.

4. **The `interaction(Condition_Combo, Group_factor)` by-variable creates n_conditions × n_groups smooths**, each with its own smoothing parameter. This is more flexible than the Pedersen-canonical single-factor `by =` pattern. It does fit better empirically (won CV), but it's worth re-evaluating if the experimental design changes.

5. **AR.start resets at phase boundaries** — this is debatable. The principled argument *against* phase resets is that animals are continuous nervous systems across phases. The argument *for* is that protocol switches create real discontinuities. The ablation tested both modes; the `phase` mode won within the C_sz_interaction Pareto branch. Future experimental designs without sharp protocol boundaries may want to revisit.

6. **No simulation-based calibration of p-values.** CV measures predictive accuracy, not p-value calibration. If publication requires demonstrating Type-I error control, simulation-based calibration on permuted condition labels (within plate) would be the rigorous next step. This was out of scope for the present diagnostic.

---

## 9. Re-running the ablation

The diagnostic tooling stays in place. To re-run:

1. In the app, toggle on the diagnostic ablation tab from the Experimental Design config bar
2. Configure axes (architectures + random bases + by-variables + AR.start modes) in the "5. Architecture × Method Ablation" sub-tab
3. Run Stage 1 (~2.5 hr for the full 75-fit sweep, animal-stratified by group)
4. Inspect the Pareto front per architecture
5. Run Stage 3 CV on the union (~100 min)
6. Compare CV ranking to the historical winner

Outputs land under `<output_dir>/ablation/stage1/<run_id>/` and `<output_dir>/ablation/stage3_cv/<run_id>/`. Importing existing results is supported via the "Import existing results from output dir" button — no need to re-fit if you just want to re-view a prior run.

---

## 10. References

- **Pedersen, E. J., Miller, D. L., Simpson, G. L., & Ross, N.** (2019). Hierarchical generalized additive models in ecology: an introduction with mgcv. *PeerJ*, 7, e6876. https://peerj.com/articles/6876/
- **Wood, S. N., Pya, N., & Säfken, B.** (2016). Smoothing parameter and model selection for general smooth models. *JASA*, 111(516), 1548–1563. (basis for mgcv's EDF-aware AIC)
- **van Rij, J., Wieling, M., Baayen, R. H., & van Rijn, H.** (2017). itsadug: Interpreting time series and autocorrelated data using GAMMs. (R package providing `start_value_rho`, `acf_resid`)
- **Pareto, V.** (1906). *Manuale di economia politica*. (origin of Pareto efficiency)
- **Karnin, Z., Koren, T., & Somekh, O.** (2013). Almost optimal exploration in multi-armed bandits. *ICML*. (successive halving — basis for staged elimination)
