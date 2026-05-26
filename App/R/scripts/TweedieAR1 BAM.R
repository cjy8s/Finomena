# ── PACKAGE CHECK ─────────────────────────────────────────────────────────────
# Packages are managed by the Python app (installed into App/R/library/).
# R_LIBS is set by the app before launching this script, so R finds them there.
cat("R library path:", paste(.libPaths(), collapse = "\n               "), "\n")
required_pkgs <- c("tidyverse", "data.table", "mgcv", "parallel", "emmeans", "jsonlite",
                   "gratia", "MASS")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install them with:  Rscript App/R/install_packages.R"
  )
}

library(tidyverse)
library(data.table)
library(mgcv)
library(parallel)
library(emmeans)
library(jsonlite)
library(gratia)
# MASS::mvrnorm used directly without library() to avoid masking dplyr::select

# Null-coalescing operator (used when reading optional sidecar fields)
`%||%` <- function(a, b) if (is.null(a)) b else a

# Multiple-testing correction algorithms + tree metadata helpers
.script_dir <- tryCatch({
  args0 <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args0, value = TRUE)
  if (length(fa) > 0) dirname(sub("^--file=", "", fa[1])) else getwd()
}, error = function(e) getwd())
source(file.path(.script_dir, "corrections.R"))
source(file.path(.script_dir, "tree_metadata.R"))

# ── ARGUMENT PARSING ──────────────────────────────────────────────────────────
# Called from the app as:
#  Rscript "TweedieAR1 BAM.R" <input_csv> <output_dir> <var_names> <ref_values> <ref_condition> <global_corr> <contrast_adj> <roles> <family_name> <shift_val>
# var_names and ref_values are comma-separated strings, e.g. "Genotype,Drug" and "WT,DMSO"
# family_name: "Tweedie" | "Gamma" | "NegBinomial"  (default: "Tweedie")
# shift_val:   numeric constant added to pxl_diff before fitting (default: 0)
args <- commandArgs(trailingOnly = TRUE)
input_csv  <- args[1]
output_dir <- if (length(args) >= 2 && !is.na(args[2])) args[2] else dirname(input_csv)

if (is.na(input_csv) || !file.exists(input_csv)) {
  stop("ERROR: input_csv argument missing or file not found: ", input_csv)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# ── DYNAMIC VARIABLE PARSING ─────────────────────────────────────────────────
var_names_raw <- if (length(args) >= 3 && !is.na(args[3]) && nchar(args[3]) > 0) args[3] else "Genotype,Drug"
ref_values_raw <- if (length(args) >= 4 && !is.na(args[4]) && nchar(args[4]) > 0) args[4] else "WT,DMSO"
ref_condition  <- if (length(args) >= 5 && !is.na(args[5]) && nchar(args[5]) > 0) args[5] else "WT+DMSO"
# Legacy positional args — the Python BAM widget still passes "BH"/"none" here
# for back-compat with older R scripts, but the new correction architecture
# drives every step from correction.json (read further down) instead. These
# values aren't used by the analysis any more.
.legacy_global_correction <- if (length(args) >= 6 && !is.na(args[6])) args[6] else "BH"
.legacy_contrast_adjust   <- if (length(args) >= 7 && !is.na(args[7])) args[7] else "none"
roles_raw                 <- if (length(args) >= 8 && !is.na(args[8]) && nchar(args[8]) > 0) args[8] else ""
family_name       <- if (length(args) >= 9 && !is.na(args[9]) && nchar(args[9]) > 0) args[9] else "Tweedie"
shift_val         <- if (length(args) >= 10 && !is.na(args[10]) && nchar(args[10]) > 0) as.numeric(args[10]) else 0

# Build the family object from the name
chosen_family <- switch(family_name,
  Tweedie     = tw(),
  Gamma       = Gamma(link = "log"),
  NegBinomial = nb(),
  {
    warning("Unknown family '", family_name, "' — defaulting to Tweedie.")
    tw()
  }
)

var_names_raw_vec <- str_split(var_names_raw, ",")[[1]]
ref_values        <- str_split(ref_values_raw, ",")[[1]]
n_vars            <- length(var_names_raw_vec)

# Make variable names R-safe (replace spaces, special chars with dots)
# Keep original names for display (plot titles, labels)
var_names_display <- var_names_raw_vec
var_names         <- make.names(var_names_raw_vec, unique = TRUE)

# Build a named reference vector: c(Genotype = "WT", Drug = "DMSO", ...)
ref_map <- setNames(ref_values[1:n_vars], var_names)

# Build display-name map for plot labels
var_display_map <- setNames(var_names_display, var_names)

# Parse roles: "WT+DMSO=RC,KO+DMSO=EXP,..." → named vector c("WT+DMSO"="RC", ...)
role_map <- c()
if (nchar(roles_raw) > 0) {
  role_pairs <- str_split(roles_raw, ",")[[1]]
  for (pair in role_pairs) {
    parts <- str_split(pair, "=")[[1]]
    if (length(parts) == 2) {
      role_map[parts[1]] <- parts[2]
    }
  }
}

# Helper: append role abbreviation to a condition name for plot labels
label_with_role <- function(cond) {
  role <- role_map[cond]
  if (!is.na(role)) paste0(cond, "\n[", role, "]") else cond
}

# Smooth-basis dimension (k) used for every time-series smooth in the formula.
# Defined here (rather than just-in-time) so the run log can echo it back to
# the user alongside the other run parameters.
k_start <- 30

cat("Input CSV:        ", input_csv, "\n")
cat("Output dir:       ", output_dir, "\n")
cat("Variables:        ", paste(var_names_display, collapse = ", "), "\n")
# Only show the make.names()-cleaned form if it actually differs from the original
if (!identical(var_names_display, var_names)) {
  cat("Variables (R):    ", paste(var_names, collapse = ", "), "\n")
}
cat("Reference values: ", paste(ref_map, collapse = ", "), "\n")
cat("Ref Condition:    ", ref_condition, "\n")
# emmeans contrasts now always run with adjust = "none"; multiple-testing
# correction (flat BH/Holm or tree TreeBH/Holm-gatekeeping/graphicalMCP) is
# applied post-hoc against the raw p-values via the correction.json sidecar.
cat("Roles:            ", roles_raw, "\n")
cat("Family:           ", family_name, "\n")
cat("Shift applied:    ", shift_val, "\n")
cat("k (smooth basis): ", k_start, "\n")

# Find how many cores the computer has, and leave 1 free so the computer doesn't freeze
usable_cores <- max(1, detectCores() - 1)

# ── DATA LOADING ──────────────────────────────────────────────────────────────
# The Python app exports a CSV with columns:
#   time_sec, location, loc_coord, pixel_diff, Condition, Phase, Group, animal_id, plate
#   + individual variable columns (Genotype, Drug, etc.)
full_df <- read_csv(input_csv, show_col_types = FALSE) %>%
  rename(pxl_diff = pixel_diff) %>%
  mutate(
    Condition = as.character(Condition),
    Phase     = as.character(Phase),
    Group     = as.integer(Group),
    time_sec  = as.integer(time_sec),
    pxl_diff  = as.numeric(pxl_diff),
    plate     = as.character(plate),
    location  = as.character(location),
    loc_coord = as.character(loc_coord)
  ) %>%
  filter(!is.na(Condition), Condition != "")

# If individual variable columns are already in the CSV, use them;
# otherwise split Condition on "+"
# Check using the original (display) names since CSV columns come from Python
if (!all(var_names_display %in% colnames(full_df))) {
  full_df <- full_df %>%
    separate(Condition, into = var_names, sep = "\\+", remove = FALSE)
} else {
  # Rename columns from display names to R-safe names
  for (i in seq_along(var_names)) {
    if (var_names_display[i] != var_names[i] && var_names_display[i] %in% colnames(full_df)) {
      full_df <- full_df %>% rename(!!var_names[i] := !!var_names_display[i])
    }
  }
}

# ---------------------------------------------------------
# 1. FINAL DATA PREP FOR MGCV
# ---------------------------------------------------------
# Unified-model design:
#   - One bam() fit across all phase groups, not one per group.
#   - time_in_group resets at each phase boundary so per-(condition × phase)
#     smooths live entirely inside their own phase.
#   - time_sec is the global continuous-time variable used by the per-well
#     and per-plate random-effect smooths so they span the whole recording.
#   - Group_factor encodes phase identity as a categorical predictor for the
#     parametric expansion and the by-interaction smooth.
#   - start_event resets the AR(1) chain at BOTH well starts AND phase
#     boundaries within wells — otherwise mgcv would attempt to correlate
#     residuals across a discontinuous phase transition.
gam_df <- full_df %>%
  mutate(
    animal_id       = as.factor(paste(plate, location, sep = "_")),
    plate           = as.factor(plate),
    Condition_Combo = as.factor(Condition),
    Group_factor    = as.factor(Group)
  ) %>%
  group_by(Group) %>%
  mutate(time_in_group = time_sec - min(time_sec)) %>%
  ungroup() %>%
  arrange(animal_id, Group, time_in_group) %>%
  group_by(animal_id) %>%
  mutate(start_event = (row_number() == 1) | (Group != lag(Group, default = -1L))) %>%
  ungroup()

# Convert each design variable to factor and set reference level globally
for (v in var_names) {
  gam_df[[v]] <- as.factor(gam_df[[v]])
  gam_df[[v]] <- relevel(gam_df[[v]], ref = ref_map[[v]])
}
gam_df$Condition_Combo <- relevel(gam_df$Condition_Combo, ref = ref_condition)

unique_groups <- sort(unique(gam_df$Group))
n_groups      <- length(unique_groups)
cat("Phase groups detected:", n_groups, "  →  ", paste(unique_groups, collapse = ", "), "\n")

# Apply shift to pxl_diff if requested (required when using Gamma family)
if (shift_val != 0) {
  gam_df <- gam_df %>% mutate(pxl_diff = pxl_diff + shift_val)
  cat(sprintf("\nShift of %.4f applied to pxl_diff (family: %s).\n", shift_val, family_name))
}

# ---------------------------------------------------------
# 2. BUILD UNIFIED BAM FORMULA
# ---------------------------------------------------------
# Parametric:   Var1 * Var2 * ... * VarN * Group_factor
#               → all design-variable main effects + all interactions
#                 INCLUDING phase-group interactions; captures per-(condition,
#                 phase) intercept shifts and any cross-phase modulation.
# Smooths:
#   s(time_in_group, k = k_start)                                  — shared within-phase shape baseline
#   s(time_in_group, by = interaction(Condition_Combo, Group_factor), k = k_start)
#                                                                  — per-(condition × phase) trajectory deviation
#   s(time_sec, plate,     bs = 'fs', m = 1)                       — plate random smooth across all time
#   s(time_sec, animal_id, bs = 'fs', m = 1)                       — well random smooth across all time
#
# Transition spikes between phases are inferred from the early-time shape of
# the next phase's per-(condition, phase) smooth plus the parametric phase
# intercepts. The random-effect smooths over time_sec span the entire
# recording so per-well / per-plate variance is estimated once jointly.
n_plates <- length(unique(gam_df$plate))
if (n_plates >= 2) {
  plate_smooth_str <- " + s(time_sec, plate, bs = 'fs', m = 1)"
} else {
  plate_smooth_str <- ""
  cat("Only ", n_plates, " plate(s) detected — omitting the plate factor smooth.\n",
      sep = "")
}

parametric_str <- paste(c(var_names, "Group_factor"), collapse = " * ")
formula_str <- paste0(
  "pxl_diff ~ ", parametric_str,
  " + s(time_in_group, k = k_start)",
  " + s(time_in_group, by = interaction(Condition_Combo, Group_factor), k = k_start)",
  plate_smooth_str,
  " + s(time_sec, animal_id, bs = 'fs', m = 1)"
)
cat("\nBAM formula (unified across all phase groups):\n  ", formula_str, "\n\n")

bam_formula <- as.formula(formula_str)

# ---------------------------------------------------------
# B. PASS 1 — ESTIMATE AR(1) RHO FROM RESIDUALS
# ---------------------------------------------------------
# Fit a no-AR model, then compute lag-1 correlation of residuals within each
# (animal_id, Group) segment. Averaging across segments gives a data-driven
# rho estimate that Pass 2 plugs in. If anything goes wrong here, fall back
# to a fixed default so the run still completes.
cat("\n========================================\n")
cat("Pass 1: fitting no-AR model to estimate autocorrelation (rho)...\n")

fallback_rho <- 0.20
optimal_rho  <- tryCatch({
  model_no_ar <- bam(
    formula  = bam_formula,
    data     = gam_df,
    family   = chosen_family,
    select   = TRUE,
    method   = "fREML",
    discrete = TRUE,
    nthreads = usable_cores
  )
  gam_df$res <- resid(model_no_ar)
  rho_est <- gam_df %>%
    dplyr::group_by(animal_id, Group) %>%
    dplyr::summarise(
      animal_rho = cor(res, dplyr::lag(res), use = "pairwise.complete.obs"),
      .groups    = "drop"
    ) %>%
    dplyr::pull(animal_rho) %>%
    mean(na.rm = TRUE)
  # Clamp to a sane range — degenerate AR1 (|rho|>=1) blows mgcv up.
  rho_est <- max(min(rho_est, 0.95), -0.95)
  gam_df$res <- NULL
  rm(model_no_ar)
  rho_est
}, error = function(e) {
  warning("Pass 1 failed (", conditionMessage(e),
          ") — falling back to fixed rho = ", fallback_rho)
  fallback_rho
})

cat(sprintf("Pass 1 estimated rho: %.4f\n", optimal_rho))

# ---------------------------------------------------------
# C. PASS 2 — FIT THE MODEL WITH AR(1) CORRECTION
# ---------------------------------------------------------
cat(sprintf("\nPass 2: fitting unified model with AR1 correction (rho = %.4f)...\n",
            optimal_rho))

final_model <- bam(
  formula  = bam_formula,
  data     = gam_df,
  family   = chosen_family,
  rho      = optimal_rho,
  AR.start = gam_df$start_event,
  select   = TRUE,
  method   = "fREML",
  discrete = TRUE,
  nthreads = usable_cores
)

cat("Unified model fitted across all", n_groups, "phase groups.\n")

# ---------------------------------------------------------
# 3. GLOBAL MULTIPLE COMPARISONS: ACROSS ALL PHASE GROUPS
# ---------------------------------------------------------
# Dynamically generate per-variable contrasts and the full interaction contrast

# ── Optional contrast-selection sidecar ────────────────────────────────────
# The Python app may write contrast_selection.json into the output dir to
# narrow the Full_Interaction family to a user-chosen subset of pairs. The
# filter is applied BEFORE building the emmeans contrast set so the local
# (Sidak / Dunnett / Tukey) and global (FDR / Holm) corrections both operate
# on the trimmed family.
#   - file present, kept_pairs non-empty → build only those pairs
#   - file present, kept_pairs empty     → skip Full_Interaction entirely
#   - file absent                        → fall back to method = "pairwise"
contrast_sidecar <- file.path(output_dir, "contrast_selection.json")
if (file.exists(contrast_sidecar)) {
  contrast_spec <- tryCatch(
    fromJSON(contrast_sidecar, simplifyVector = FALSE),
    error = function(e) {
      warning("Could not parse contrast_selection.json: ", conditionMessage(e),
              " — falling back to full pairwise.")
      NULL
    }
  )
  if (!is.null(contrast_spec)) {
    cat("Contrast selection sidecar found: ",
        length(contrast_spec$kept_pairs), " pair(s) requested.\n", sep = "")
  }
} else {
  contrast_spec <- NULL
}

all_contrasts_list <- list()

# ── Reference-grid sizing ───────────────────────────────────────────────────
# Single unified model: the grid now expands across (var levels × Group_factor).
# Bump rg.limit to accommodate.
random_factor_names <- character(0)
for (sm in final_model$smooth) {
  if (inherits(sm, "fs.interaction")) {
    random_factor_names <- c(random_factor_names, sm$fterm)
  } else if (inherits(sm, "random.effect")) {
    random_factor_names <- c(random_factor_names, sm$term)
  }
}
random_factor_names <- unique(random_factor_names)

candidate_factor_vars <- intersect(
  c(var_names, "Condition_Combo", "Group_factor", "animal_id"),
  colnames(gam_df)
)
factor_sizes <- vapply(candidate_factor_vars, function(v) {
  if (is.factor(gam_df[[v]])) nlevels(gam_df[[v]]) else 1L
}, integer(1))
fixed_factor_vars  <- setdiff(candidate_factor_vars, random_factor_names)
fixed_factor_sizes <- factor_sizes[fixed_factor_vars]
grid_rows <- if (length(fixed_factor_sizes) > 0)
               prod(as.numeric(fixed_factor_sizes)) else 1
cat(sprintf("\nRef-grid fixed factors: %s\n",
            paste(sprintf("%s(%d)", fixed_factor_vars, fixed_factor_sizes),
                  collapse = " × ")))
cat(sprintf("Estimated emmeans grid: %d rows\n", as.integer(grid_rows)))
emm_options(rg.limit = max(50000, as.integer(grid_rows * 2)))

# ── Per-variable effects ────────────────────────────────────────────────────
# Group_factor joins the conditioning so the same trt-vs-ref contrast is
# produced once per (other-vars × Group) cell from the single unified model.
for (vi in seq_along(var_names)) {
  focal_var   <- var_names[vi]
  other_vars  <- var_names[-vi]
  ref_val     <- ref_map[[focal_var]]

  cond_vars <- c(other_vars, "Group_factor")
  spec_str  <- paste0("~ ", focal_var, " | ", paste(cond_vars, collapse = " + "))
  spec_formula <- as.formula(spec_str)

  emm_base <- emmeans(final_model, specs = spec_formula)
  emm_contrasts <- contrast(emm_base, method = "trt.vs.ctrl",
                            ref = ref_val, adjust = "none")

  df_contrasts <- as.data.frame(emm_contrasts) %>%
    mutate(
      Test_Family  = paste0(focal_var, "_Effect"),
      Tested_Level = as.character(contrast),
      Group        = as.integer(as.character(Group_factor))
    ) %>%
    select(-Group_factor)

  # Collapse remaining other-variable columns into a single Split_By string
  split_cols <- intersect(other_vars, colnames(df_contrasts))
  if (length(split_cols) > 0) {
    df_contrasts$Split_By <- apply(df_contrasts[, split_cols, drop = FALSE], 1,
                                    function(x) paste(x, collapse = " + "))
    df_contrasts <- df_contrasts %>% select(-all_of(split_cols))
  } else {
    df_contrasts$Split_By <- "All"
  }

  all_contrasts_list[[focal_var]] <- df_contrasts
}

# ── Full interaction: pairwise condition comparisons by Group ───────────────
# Use the DESIGN variables (Var1 * Var2 * ...) — NOT Condition_Combo — because
# emmeans detects Condition_Combo's levels as nested in the parametric
# Genotype*Drug*Group_factor expansion and strips it out of the reference
# grid. Asking for Condition_Combo by name then errors with
#   "No variable named Condition_Combo in the reference grid".
# Building the spec from the design variables sidesteps the nesting issue
# entirely; we reconstruct the "+"-joined condition string from the grid
# columns for matching against the sidecar's kept_pairs entries.
inter_spec <- as.formula(paste0(
  "~ ", paste(var_names, collapse = " * "), " | Group_factor"
))
emm_inter_base <- emmeans(final_model, specs = inter_spec)

# Pull a single by-level's sub-grid to learn the within-group row ordering;
# the weight vectors built below have length = nrow(that sub-grid) and are
# replayed by emmeans across every Group_factor level via by="Group_factor".
full_inter_grid <- emm_inter_base@grid
first_group     <- levels(gam_df$Group_factor)[1]
within_grid     <- full_inter_grid[
  as.character(full_inter_grid$Group_factor) == first_group,
  var_names, drop = FALSE
]
cond_strs <- apply(within_grid, 1,
                   function(r) paste(as.character(r), collapse = "+"))

emm_inter_contrasts <- NULL
if (!is.null(contrast_spec)) {
  kept_pairs <- contrast_spec$kept_pairs
  if (length(kept_pairs) == 0) {
    cat("Full_Interaction skipped (no pairs selected).\n")
  } else {
    custom_contrasts <- list()
    missing <- c()
    for (pair in kept_pairs) {
      lhs <- pair[[1]]; rhs <- pair[[2]]
      lhs_i <- which(cond_strs == lhs)
      rhs_i <- which(cond_strs == rhs)
      if (length(lhs_i) == 1L && length(rhs_i) == 1L) {
        v <- numeric(length(cond_strs))
        v[lhs_i] <-  1
        v[rhs_i] <- -1
        custom_contrasts[[paste(lhs, "-", rhs)]] <- v
      } else {
        missing <- c(missing, paste(lhs, "vs", rhs))
      }
    }
    if (length(missing) > 0) {
      warning("Couldn't locate pair(s) in the design-variable grid: ",
              paste(missing, collapse = "; "))
    }
    if (length(custom_contrasts) > 0) {
      emm_inter_contrasts <- contrast(emm_inter_base,
                                      method = custom_contrasts,
                                      by = "Group_factor",
                                      adjust = "none")
    }
  }
} else {
  # Legacy fallback (no contrast_selection.json sidecar): all-pairs.
  emm_inter_contrasts <- contrast(emm_inter_base,
                                  method = "pairwise",
                                  by = "Group_factor",
                                  adjust = "none")
}

if (!is.null(emm_inter_contrasts)) {
  df_inter <- as.data.frame(emm_inter_contrasts) %>%
    mutate(
      Test_Family  = "Full_Interaction",
      Tested_Level = as.character(contrast),
      Split_By     = "None",
      Group        = as.integer(as.character(Group_factor))
    ) %>%
    select(-Group_factor)
  all_contrasts_list[["interaction"]] <- df_inter
}

# ---------------------------------------------------------
# COMBINE + APPLY MULTIPLE-TESTING CORRECTION
# ---------------------------------------------------------
# The contrast loop above produced RAW p-values (emmeans `adjust = "none"`).
# Now apply the user's chosen correction method to those raw values per the
# correction.json sidecar:
#   - strategy   = "flat" | "tree"
#   - error_rate = "FDR"  | "FWER"
#   - contrast_set = "all_pairs" | "ref_only"   (affects graphicalMCP vs Holm)
#   - threshold  = q (FDR) or alpha (FWER), default 0.05
#   - tree       = list of level specs (only used when strategy=="tree")
#
# If no sidecar is present, default to flat BH per Test_Family at q = 0.05
# (matches the historical behaviour from before the new architecture).
master_results_df <- bind_rows(all_contrasts_list)

final_master_table <- master_results_df %>%
  mutate(
    raw_pvalue = p.value,
    estimate   = estimate / log(2),
    SE         = SE       / log(2),
    Group      = as.integer(Group)
  ) %>%
  select(Test_Family, Group, Split_By, Tested_Level, estimate, SE, raw_pvalue) %>%
  arrange(Test_Family, Group)

# Read correction.json if present
correction_sidecar <- file.path(output_dir, "correction.json")
correction_spec <- if (file.exists(correction_sidecar)) {
  tryCatch(
    fromJSON(correction_sidecar, simplifyVector = FALSE),
    error = function(e) {
      warning("Could not parse correction.json: ", conditionMessage(e),
              " — falling back to flat BH at q=0.05.")
      NULL
    }
  )
} else NULL

if (is.null(correction_spec)) {
  # Default: flat BH per-family at q=0.05
  correction_spec <- list(
    strategy   = "flat",
    error_rate = "FDR",
    threshold  = 0.05,
    scope      = "per_family"
  )
  cat("\nNo correction.json sidecar — using default: flat BH per family at q=0.05.\n")
} else {
  cat("\nCorrection.json found: strategy=", correction_spec$strategy %||% "flat",
      " error_rate=", correction_spec$error_rate %||% "FDR",
      " threshold=", correction_spec$threshold %||% 0.05, "\n", sep = "")
}

# If a tree strategy, annotate rows with tree metadata BEFORE correction.
if (identical(correction_spec$strategy, "tree")) {
  tree_spec <- correction_spec$tree
  if (is.null(tree_spec) || length(tree_spec) == 0) {
    cat("Tree strategy chosen but no tree spec provided — using default rescue tree.\n")
    tree_spec <- default_rescue_tree_spec(var_names, ref_map)
  }
  final_master_table <- assign_tree_metadata(final_master_table, tree_spec)
}

final_master_table <- apply_correction(final_master_table, correction_spec)

# Tidy up column order — tree columns may not have been added by flat methods
# but they're always present after apply_correction (NAs are fine).
final_master_table <- final_master_table %>%
  select(Test_Family, Group, Split_By, Tested_Level,
         estimate, SE,
         raw_pvalue, adjusted_pvalue, correction_method,
         tidyselect::any_of(c("tree_id", "tree_level",
                              "tree_parent_id", "tree_status")))

cat("\n--- CLEAN MASTER TABLE ---\n")
# Print via as.data.frame to avoid a partial-arg-matching crash when the
# table is a plain data.frame (mgcv select() loses the tibble class).
print(as.data.frame(final_master_table), row.names = FALSE)

# Persist the per-contrast results for downstream visualisation.
write_csv(final_master_table, "master_results.csv")
cat("Saved: master_results.csv\n")

# ---------------------------------------------------------
# SUMMARY TABLE FOR QUICK REFERENCE
# ---------------------------------------------------------

summary_table <- final_master_table %>%
  group_by(Test_Family, Group) %>%
  summarise(
    n_total = n(),
    n_sig = sum(adjusted_pvalue < 0.05, na.rm = TRUE),
    mean_effect = mean(abs(estimate)),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(n_sig, mean_effect),
    names_glue = "G{Group}_{.value}"
  )

write_csv(summary_table, "summary_statistics_by_group.csv")

# =============================================================================
# 4. POSTERIOR EQUIVALENCE TESTING (Bayesian)
# =============================================================================
# Optional. Runs when contrast_selection.json sets posterior_equivalence.enabled.
# For each kept pair (lhs, rhs) and each phase group, draws from the BAM
# posterior to compute
#     M = max_t |η_lhs(t) − η_rhs(t)|     (log_e fold change at each draw)
# and converts to log_2 by dividing by log(2). Reports the posterior of M and
# Pr(M < δ) — a direct Bayesian credibility statement of equivalence within ±δ.
# Implementation uses mgcv's lpmatrix posterior + MASS::mvrnorm coefficient
# draws. (gratia is loaded for users who want to do additional posterior work.)

# Helper: for one (model, full_df, lhs, rhs, group_id), return the posterior of M.
# Unified model produces per-(pair × group) draws by varying Condition_Combo
# AND Group_factor in the newdata frame.
.posterior_equivalence_pair <- function(model, full_df, lhs, rhs, var_names,
                                        group_id, delta_log2, n_draws,
                                        n_time = 200, seed = 42) {
  cond_levels  <- levels(full_df$Condition_Combo)
  group_levels <- levels(full_df$Group_factor)
  group_str    <- as.character(group_id)
  if (!(lhs %in% cond_levels) || !(rhs %in% cond_levels)) return(NULL)
  if (!(group_str %in% group_levels))                     return(NULL)

  parse_cond <- function(cond_str) {
    parts <- strsplit(cond_str, "+", fixed = TRUE)[[1]]
    if (length(parts) != length(var_names)) {
      stop("Cannot parse '", cond_str, "': expected ", length(var_names),
           " parts (one per variable), got ", length(parts))
    }
    setNames(as.list(parts), var_names)
  }
  lhs_vars <- parse_cond(lhs)
  rhs_vars <- parse_cond(rhs)

  # Use the time_in_group range for THIS phase group so the trajectory is
  # built only over the phase the comparison lives in.
  group_rows <- full_df[full_df$Group == group_id, ]
  if (nrow(group_rows) == 0) return(NULL)
  time_grid <- seq(0, max(group_rows$time_in_group), length.out = n_time)
  # time_sec is needed for the random-effect smooths even though those are
  # excluded from the linear predictor below — predict() still requires the
  # column to exist. Anchor it at the midpoint of the group's range.
  time_sec_anchor <- mean(range(group_rows$time_sec))

  build_nd <- function(cond_str, vars_list) {
    nd <- data.frame(
      time_in_group = time_grid,
      time_sec      = time_sec_anchor
    )
    for (v in var_names) {
      nd[[v]] <- factor(vars_list[[v]], levels = levels(full_df[[v]]))
    }
    nd$Condition_Combo <- factor(cond_str,  levels = cond_levels)
    nd$Group_factor    <- factor(group_str, levels = group_levels)
    # animal_id and plate are required by predict() but excluded from the
    # linear predictor (population-level trajectories).
    nd$animal_id <- factor(full_df$animal_id[1],
                           levels = levels(full_df$animal_id))
    nd$plate <- factor(full_df$plate[1],
                       levels = levels(full_df$plate))
    nd
  }
  nd_lhs <- build_nd(lhs, lhs_vars)
  nd_rhs <- build_nd(rhs, rhs_vars)

  # Linear-predictor design matrices, excluding the random-effect smooths so
  # we get population-level (not animal- or plate-specific) trajectories.
  # Random-effect smooths now live on time_sec (not time_in_group).
  excl <- c("s(time_sec,animal_id)", "s(time_sec,plate)")
  X_lhs <- predict(model, newdata = nd_lhs, type = "lpmatrix", exclude = excl)
  X_rhs <- predict(model, newdata = nd_rhs, type = "lpmatrix", exclude = excl)
  X_diff <- X_lhs - X_rhs                      # n_time × n_coef

  # Bayesian (large-sample Gaussian) posterior of β
  set.seed(seed)
  Cv <- vcov(model, unconditional = TRUE)
  beta_samples <- MASS::mvrnorm(n_draws, mu = coef(model), Sigma = Cv)  # n_draws × n_coef

  # Posterior draws of η_lhs(t) − η_rhs(t) on log_e scale → divide by log(2)
  diff_samples <- X_diff %*% t(beta_samples)   # n_time × n_draws
  M_log2 <- apply(abs(diff_samples), 2, max) / log(2)

  list(
    M_log2  = M_log2,
    median  = median(M_log2),
    lower95 = unname(quantile(M_log2, 0.025)),
    upper95 = unname(quantile(M_log2, 0.975)),
    pr_equiv = mean(M_log2 < delta_log2)
  )
}

if (!is.null(contrast_spec) &&
    !is.null(contrast_spec$posterior_equivalence) &&
    isTRUE(contrast_spec$posterior_equivalence$enabled) &&
    length(contrast_spec$kept_pairs) > 0) {

  delta_log2 <- as.numeric(contrast_spec$posterior_equivalence$delta %||% 1.0)
  n_draws    <- as.integer(contrast_spec$posterior_equivalence$n_draws %||% 10000)

  cat("\n========================================\n")
  cat("Posterior Equivalence Testing (Bayesian)\n")
  cat(sprintf("  delta = %.3f (log_2 fold change)\n", delta_log2))
  cat(sprintf("  posterior draws = %d\n", n_draws))
  cat("========================================\n")

  pe_summary_rows <- list()
  pe_draws_named  <- list()

  # Unified-model PE: iterate over (group, pair) directly against the single
  # final_model. Prediction grids carry both Condition_Combo and Group_factor.
  for (g in unique_groups) {
    g_str <- as.character(g)
    cat(sprintf("Group %s: posterior equivalence...\n", g_str))

    for (i in seq_along(contrast_spec$kept_pairs)) {
      pair <- contrast_spec$kept_pairs[[i]]
      lhs <- pair[[1]]; rhs <- pair[[2]]

      pe <- tryCatch(
        .posterior_equivalence_pair(final_model, gam_df, lhs, rhs, var_names,
                                    group_id = g, delta_log2 = delta_log2,
                                    n_draws = n_draws, seed = 42L + i),
        error = function(e) {
          warning(sprintf("Group %s, %s vs %s: %s", g_str, lhs, rhs,
                          conditionMessage(e)))
          NULL
        }
      )
      if (is.null(pe)) next

      pe_summary_rows[[length(pe_summary_rows) + 1]] <- data.frame(
        Group         = as.integer(g_str),
        lhs           = lhs,
        rhs           = rhs,
        pair          = paste(lhs, "vs", rhs),
        M_median_log2 = pe$median,
        M_lower95     = pe$lower95,
        M_upper95     = pe$upper95,
        Pr_equiv      = pe$pr_equiv,
        delta_log2    = delta_log2,
        stringsAsFactors = FALSE
      )
      pe_draws_named[[paste0("g", g_str, "__", lhs, "__", rhs)]] <- pe$M_log2
    }
  }

  if (length(pe_summary_rows) > 0) {
    pe_summary_df <- bind_rows(pe_summary_rows) %>%
      arrange(Group, pair) %>%
      as_tibble()
    cat("\n--- POSTERIOR EQUIVALENCE SUMMARY ---\n")
    # Print via as.data.frame to avoid a partial-arg-matching crash when the
    # table is a plain data.frame instead of a tibble.
    print(as.data.frame(pe_summary_df), row.names = FALSE)
    write_csv(pe_summary_df, "posterior_equivalence_summary.csv")

    # Long form of the raw posterior draws — persisted so the Python
    # Visualisations tab can render the PE density and heatmap interactively.
    draws_df <- bind_rows(lapply(seq_along(pe_draws_named), function(j) {
      key <- names(pe_draws_named)[j]
      m <- regmatches(key, regexec("^g(.+?)__(.+?)__(.+)$", key))[[1]]
      data.frame(
        Group  = as.integer(m[2]),
        pair   = paste(m[3], "vs", m[4]),
        M_log2 = pe_draws_named[[j]],
        stringsAsFactors = FALSE
      )
    }))
    write_csv(draws_df, "posterior_equivalence_draws.csv")
    cat("Saved: posterior_equivalence_draws.csv (",
        nrow(draws_df), " draws total)\n", sep = "")
  } else {
    cat("No posterior-equivalence results were produced (all pairs skipped).\n")
  }
}

cat("\nAnalysis complete. All outputs saved to:", output_dir, "\n")
