# =============================================================================
# Stand-Alone BAM Analysis Script
# -----------------------------------------------------------------------------
# A standalone R port of the Finomena app's BAM pipeline (App/R/scripts/
# TweedieAR1 BAM.R), bundled with the raw-data aggregation step so it can run
# end-to-end outside the app.
#
# Pipeline:
#   1. Read three input files (design.csv, time_frame.csv, per-plate
#      plate_layout.csv).
#   2. For each plate subdirectory, read the Zebrabox CSV(s), drop non-data
#      rows, bin to 1-second sums per well, attach Condition from the layout,
#      and tag plate + animal_id.
#   3. Join phase / group info via the time_frame.
#   4. Split each Condition string on "+" into the per-variable factor
#      columns declared in design.csv (same "+"-join convention as the app's
#      Experimental Conditions widget).
#   5. Fit a Tweedie BAM with AR(1) (two-pass: pilot fit -> rho -> refit).
#   6. Compute emmeans contrasts (per-variable effects + Full_Interaction),
#      apply the chosen multiple-testing correction, write master_results.csv
#      and summary_statistics_by_group.csv.
#
# See "Stand-Alone Script ReadMe.md" for input file shapes and usage.
# =============================================================================


# ── USER CONFIG ──────────────────────────────────────────────────────────────
# Paths can be absolute or relative to this script's location.

INPUT_DIR  <- "Stand-Alone Inputs"
OUTPUT_DIR <- "Stand-Alone Output"

# Contrast selection for the Full_Interaction family:
#   TRUE  - all pairwise Condition contrasts
#   FALSE - experimental-vs-reference-condition pairs only
ALL_PAIRS <- FALSE

# Multiple-testing correction (mirrors the correction.json sidecar shape used
# by the app). Strategy = "flat" or "tree"; error_rate = "FDR" or "FWER".
CORRECTION_SPEC <- list(
  strategy   = "flat",
  error_rate = "FDR",
  threshold  = 0.05,
  scope      = "per_family"
)

# Bayesian posterior equivalence testing (mirrors the app's
# contrast_selection.json posterior_equivalence block, but always uses
# experimental-vs-reference pairs in the standalone). For each kept pair ×
# phase group, draws coefficient samples from the BAM posterior, computes
# M = max_t |η_lhs(t) − η_rhs(t)| / log(2), and reports the posterior of M
# along with Pr(M < delta) — a direct Bayesian credibility statement of
# equivalence within ±delta log_2 fold change.
POSTERIOR_EQUIVALENCE <- list(
  enabled = TRUE,
  delta   = 1.0,
  n_draws = 10000
)

# Path (relative to this script) to the app's correction helper scripts. The
# standalone script sources them rather than re-implementing the algorithms.
HELPERS_DIR <- file.path("..", "App", "R", "scripts")


# ── ANCHOR ALL RELATIVE PATHS TO THIS SCRIPT'S LOCATION ─────────────────────
# Works whether the script is run via Rscript or sourced from RStudio.
.get_script_dir <- function() {
  args0 <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args0, value = TRUE)
  if (length(fa) > 0) return(normalizePath(dirname(sub("^--file=", "", fa[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(normalizePath(dirname(p)))
  }
  normalizePath(getwd())
}
SCRIPT_DIR <- .get_script_dir()

# Lightweight absolute-path check that handles both Unix (/foo) and Windows
# (C:\foo) — avoids dragging in R.utils for one predicate.
.is_abs <- function(p) grepl("^([A-Za-z]:)?[/\\\\]", p)
.resolve <- function(p) {
  if (length(p) == 0 || !nzchar(p)) return(p)
  if (.is_abs(p)) p else normalizePath(file.path(SCRIPT_DIR, p), mustWork = FALSE)
}

INPUT_DIR   <- .resolve(INPUT_DIR)
OUTPUT_DIR  <- .resolve(OUTPUT_DIR)
HELPERS_DIR <- .resolve(HELPERS_DIR)


# ── USE ACTIVE CONDA ENV'S R LIBRARY ─────────────────────────────────────────
# Install/load packages from the active conda env (CONDA_PREFIX) rather than a
# system or user library, so chris_conda owns its own package set.
.lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library")
dir.create(.lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(.lib, .libPaths()))


# ── PACKAGE CHECK + AUTO-INSTALL ─────────────────────────────────────────────
# Any required CRAN package not already installed is installed from CRAN
# before the script proceeds. `parallel` is base R and ships with the R
# install, so it's never an install target.
required_pkgs <- c("tidyverse", "data.table", "mgcv", "parallel", "emmeans",
                   "jsonlite", "gratia", "MASS")
cran_mirror   <- "https://cloud.r-project.org"

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                       logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  # Base R packages can't be installed via install.packages() — exclude them.
  installable <- setdiff(missing_pkgs, c("parallel"))
  if (length(installable) > 0) {
    cat("Installing missing R packages from CRAN: ",
        paste(installable, collapse = ", "), "\n", sep = "")
    install.packages(installable, repos = cran_mirror)
    # Re-check after the install pass.
    missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                           logical(1), quietly = TRUE)]
  }
  if (length(missing_pkgs) > 0) {
    stop("Still missing after install attempt: ",
         paste(missing_pkgs, collapse = ", "),
         "\nInstall manually with install.packages() and re-run.")
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(mgcv)
  library(parallel)
  library(emmeans)
  library(jsonlite)
})

options(width = 320)
options(warn = 1)

`%||%` <- function(a, b) if (is.null(a)) b else a


# ── SOURCE CORRECTION HELPERS ────────────────────────────────────────────────
.corrections_path   <- file.path(HELPERS_DIR, "corrections.R")
.tree_metadata_path <- file.path(HELPERS_DIR, "tree_metadata.R")
if (!file.exists(.corrections_path) || !file.exists(.tree_metadata_path)) {
  stop("Could not find correction helper scripts at:\n  ", HELPERS_DIR,
       "\nSet HELPERS_DIR (top of this script) to point at App/R/scripts/.")
}
source(.corrections_path)
source(.tree_metadata_path)


# ── READ DESIGN INPUTS ───────────────────────────────────────────────────────
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Input dir:   ", INPUT_DIR,  "\n")
cat("Output dir:  ", OUTPUT_DIR, "\n\n")

design_path <- file.path(INPUT_DIR, "design.csv")
if (!file.exists(design_path)) stop("Missing design.csv at: ", design_path)
design_df <- read.csv(design_path, stringsAsFactors = FALSE)
if (!all(c("variable", "reference") %in% colnames(design_df))) {
  stop("design.csv must have columns: variable, reference")
}
var_names_display <- as.character(design_df$variable)
ref_values        <- as.character(design_df$reference)
n_vars            <- length(var_names_display)
# R-safe names (e.g. spaces -> dots) for use as data-frame column names.
var_names         <- make.names(var_names_display, unique = TRUE)
ref_map           <- setNames(ref_values, var_names)
# Canonical "+"-joined reference condition (matches the app's convention).
ref_condition     <- paste(ref_values, collapse = "+")

cat("Design variables: ", paste(var_names_display, collapse = ", "), "\n")
cat("Reference levels: ", paste(ref_values,        collapse = ", "), "\n")
cat("Ref condition:    ", ref_condition, "\n\n")

time_frame_path <- file.path(INPUT_DIR, "time_frame.csv")
if (!file.exists(time_frame_path)) stop("Missing time_frame.csv at: ", time_frame_path)
time_frame <- read.csv(time_frame_path, stringsAsFactors = FALSE)
if (!all(c("Group", "Phase", "End") %in% colnames(time_frame))) {
  stop("time_frame.csv must have columns: Group, Phase, End  ",
       "(Start is optional and ignored; phase assignment is End-only ",
       "to match the app's pd.merge_asof(direction='forward') semantics).")
}
if ("Start" %in% colnames(time_frame)) {
  cat("Note: 'Start' column present in time_frame.csv but ignored. ",
      "Phase assignment uses End-only semantics matching the app.\n", sep = "")
}
time_frame$Group <- as.integer(time_frame$Group)
time_frame$End   <- as.integer(time_frame$End)
cat("Time frame:\n")
print(time_frame, row.names = FALSE)
cat("\n")


# ── DISCOVER PLATES ──────────────────────────────────────────────────────────
# Recursively scan INPUT_DIR for any directory containing a plate_layout.csv.
# Each such directory is one plate (Option A: replicates are treated as
# independent plates, each with its own `plate` factor level in the BAM
# model). The plate's name is its path relative to INPUT_DIR with separators
# normalized to underscores, so a hierarchy like
#   Plate 1/Replicate 1/plate_layout.csv
#   Plate 1/Replicate 2/plate_layout.csv
#   Plate 2/Replicate 1/plate_layout.csv
# yields three distinct plate factor levels. Per-replicate layouts are always
# read independently — no inheritance from parent dirs, no caching.
all_dirs   <- list.dirs(INPUT_DIR, recursive = TRUE, full.names = TRUE)
plate_dirs <- all_dirs[file.exists(file.path(all_dirs, "plate_layout.csv"))]
# Don't scan our own output dir if it happens to live under INPUT_DIR.
plate_dirs <- plate_dirs[normalizePath(plate_dirs, mustWork = FALSE) !=
                         normalizePath(OUTPUT_DIR, mustWork = FALSE)]
if (length(plate_dirs) == 0) {
  stop("No plate_layout.csv found anywhere under:\n  ", INPUT_DIR)
}

input_norm  <- normalizePath(INPUT_DIR, mustWork = FALSE)
plate_names <- vapply(plate_dirs, function(d) {
  d_norm <- normalizePath(d, mustWork = FALSE)
  rel    <- substring(d_norm, nchar(input_norm) + 2L)
  if (nchar(rel) == 0L) rel <- basename(d_norm)
  gsub("[/\\\\]", "_", rel)
}, character(1), USE.NAMES = FALSE)

cat("Plates:\n")
for (i in seq_along(plate_dirs)) {
  cat("  ", plate_names[i], "\n", sep = "")
}
cat("\n")


# ── PARSE A PER-PLATE LAYOUT ─────────────────────────────────────────────────
# Layout shape: first column = row letter (A, B, C, ...), remaining columns =
# well-column numbers (1, 2, 3, ...). Each interior cell is a "+"-joined
# Condition string (or "empty" / "" for unused wells).
parse_plate_layout <- function(layout_path, expected_n_vars) {
  raw <- read.csv(layout_path, stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = c("", "NA"))
  if (ncol(raw) < 2) {
    stop("plate_layout.csv at ", layout_path,
         " must have a row-letter column followed by one column per well-column.")
  }
  row_letters <- as.character(raw[[1]])
  grid        <- as.matrix(raw[, -1, drop = FALSE])
  col_numbers <- colnames(raw)[-1]

  # Validate "+"-split arity for every non-empty cell.
  flat <- as.character(grid)
  is_empty <- is.na(flat) | flat == "" | tolower(flat) == "empty"
  non_empty <- flat[!is_empty]
  if (length(non_empty) > 0) {
    splits <- strsplit(non_empty, "+", fixed = TRUE)
    bad <- which(lengths(splits) != expected_n_vars)
    if (length(bad) > 0) {
      stop(sprintf(
        "%s: %d cell(s) do not split into %d '+'-separated parts.\n  Examples: %s",
        layout_path, length(bad), expected_n_vars,
        paste(unique(head(non_empty[bad], 5)), collapse = ", ")))
    }
  }

  # Build well -> Condition lookup.
  well_to_cond <- character(0)
  for (r in seq_along(row_letters)) {
    for (c in seq_along(col_numbers)) {
      well <- paste0(row_letters[r], col_numbers[c])
      v <- grid[r, c]
      if (is.na(v) || v == "" || tolower(v) == "empty") next
      well_to_cond[[well]] <- as.character(v)
    }
  }
  list(
    well_to_cond = well_to_cond,
    row_letters  = row_letters,
    col_numbers  = col_numbers,
    rows         = length(row_letters),
    cols         = length(col_numbers)
  )
}


# ── LOAD ONE PLATE'S RAW CSV(S) AND BIN TO 1-SEC ────────────────────────────
load_plate <- function(plate_dir, plate_name, plate_idx) {
  cat("  [", plate_idx, "] ", plate_name, " ... ", sep = "")

  layout <- parse_plate_layout(file.path(plate_dir, "plate_layout.csv"), n_vars)

  # All .csv / .xls / .xlsx data files in THIS plate dir, minus plate_layout
  # itself. Non-recursive — otherwise a parent dir (e.g. "Plate 1/") that
  # happens to have its own plate_layout.csv would suck in data files from
  # every replicate underneath. Zebrabox emits tab-separated text with a .xls
  # extension, so fread's auto-detection handles both.
  data_files <- list.files(plate_dir, pattern = "\\.(csv|xls|xlsx)$",
                           full.names = TRUE, recursive = FALSE,
                           ignore.case = TRUE)
  data_files <- data_files[!grepl("^plate_layout\\.",
                                   basename(data_files), ignore.case = TRUE)]
  if (length(data_files) == 0) {
    cat("no data files (.csv / .xls / .xlsx) found, skipping.\n")
    return(NULL)
  }

  dfs <- vector("list", length(data_files))
  for (i in seq_along(data_files)) {
    # Zebrabox .xls files are tab-separated and inconsistently shaped:
    # headers may be short (6 cols) or padded out to data2..data16 (22 cols),
    # and data rows can have a trailing tab or vary by 1 column. `fill = TRUE`
    # lets fread tolerate the ragged rows; explicit sep = "\t" avoids it
    # mis-detecting on the short-header files. data1 is read as character
    # because marker rows (SOUND / BACK_LIGHT / etc.) put strings there.
    tmp <- fread(data_files[i], na.strings = c("NA", ""),
                 sep = "\t", fill = TRUE, showProgress = FALSE,
                 colClasses = list(character = "data1"))
    keep_cols <- intersect(c("time", "location", "data1"), names(tmp))
    if (!all(c("time", "location", "data1") %in% keep_cols)) {
      warning(data_files[i], " missing required column(s) — skipped.")
      next
    }
    tmp <- tmp[, ..keep_cols]
    # Mirror Python (data_loader.py L661): drop time<=0 rows before binning,
    # so the time_sec=0 bin contains the same samples in both pipelines.
    tmp <- tmp[time > 0]
    # Mirror Python (data_loader.py L681): coerce data1 to integer; non-numeric
    # marker strings (SOUND / BACK_LIGHT / TOP_LIGHT / "Not identify") become
    # NA via suppressWarnings — equivalent to pd.to_numeric(errors='coerce').
    # DO NOT drop NA rows: bins containing only markers must still emit a row
    # with pixel_diff = 0, which sum(na.rm = TRUE) produces — matching pandas'
    # default sum(skipna = True).
    tmp[, data1 := suppressWarnings(as.integer(data1))]
    dfs[[i]] <- tmp
  }
  dfs <- dfs[!sapply(dfs, is.null)]
  if (length(dfs) == 0) { cat("no usable rows, skipping.\n"); return(NULL) }
  raw <- rbindlist(dfs)

  # Bin to one row per (location, integer second).
  raw[, time_sec := as.integer(time / 1e6)]
  binned <- raw[, .(pixel_diff = sum(data1, na.rm = TRUE)),
                by = .(location, time_sec)]

  # Reconstruct well coordinate ("A1", "B7", ...) from the integer location
  # index using the plate's row/column dimensions from its layout.
  binned[, loc_id := suppressWarnings(as.integer(str_extract(location, "\\d+")))]
  binned <- binned[!is.na(loc_id)]
  cols_p <- layout$cols
  rows_p <- layout$rows
  binned[, loc_coord := paste0(
    layout$row_letters[((loc_id - 1) %/% cols_p) + 1],
    layout$col_numbers[((loc_id - 1) %% cols_p) + 1]
  )]
  # Drop locations beyond the layout's well count (shouldn't happen with a
  # well-formed plate but defensive).
  binned <- binned[loc_id <= rows_p * cols_p]

  # Attach Condition from the per-plate layout.
  binned[, Condition := unname(layout$well_to_cond[loc_coord])]
  binned <- binned[!is.na(Condition) & Condition != ""]

  # Match the app exactly:
  #   plate     = integer plate index (1, 2, ...) — bam_widget.py L242 reads
  #               the integer `plate` column the data loader assigns at
  #               data_loader.py L670 (plate_combined['plate'] = plate_idx).
  #   animal_id = paste(plate, location, sep = "_"), where `location` is the
  #               raw Zebrabox "loc_<n>" string (NOT the reconstructed A1/B7
  #               loc_coord). Mirrors bam_widget.py L242:
  #                 animal_id = df['plate'].astype(str) + "_" + df['location'].astype(str)
  #               and the app R's TweedieAR1 BAM.R L271 reconstruction
  #                 animal_id = paste(plate, location, sep = "_").
  # plate_name (the path-relative directory name) is no longer used for any
  # model-facing column — it remains the per-plate console label only.
  binned[, plate     := plate_idx]
  binned[, animal_id := paste(plate_idx, location, sep = "_")]

  cat(nrow(binned), " rows.\n", sep = "")
  binned[]
}

cat("── Aggregating plate data ──\n")
plate_tabs <- Map(load_plate, plate_dirs, plate_names, seq_along(plate_dirs))
plate_tabs <- plate_tabs[!sapply(plate_tabs, is.null)]
if (length(plate_tabs) == 0) stop("All plates yielded zero usable rows.")
full_df <- rbindlist(plate_tabs)
cat("Aggregated total: ", nrow(full_df), " rows across ",
    length(plate_tabs), " plate(s).\n", sep = "")
cat("Data time_sec range: [", min(full_df$time_sec), ", ",
    max(full_df$time_sec), "] seconds\n\n", sep = "")


# ── JOIN PHASE / GROUP FROM time_frame ──────────────────────────────────────
# Mirrors the app's pd.merge_asof(direction='forward') on End
# (data_loader.py L700-706): for each row at time_sec t, assign the phase
# whose End is the smallest value >= t. Equivalently, phase i covers
# (End_{i-1}, End_i], with phase 1 covering [0, End_1]; Start is not used.
# Rows with t > max(End) get NA and are dropped.
time_frame <- time_frame[order(time_frame$End), ]
phase_idx  <- findInterval(full_df$time_sec - 1L, time_frame$End) + 1L
phase_idx[phase_idx > nrow(time_frame)] <- NA_integer_
full_df[, `:=`(
  Phase = time_frame$Phase[phase_idx],
  Group = time_frame$Group[phase_idx]
)]
n_before <- nrow(full_df)
full_df  <- full_df[!is.na(Phase)]
cat("Phase join: kept ", nrow(full_df), " / ", n_before,
    " rows.\n\n", sep = "")

if (nrow(full_df) == 0) {
  stop("No data rows fell within any time_frame.csv window.\n",
       "  Data time_sec range:    [",
       min(rbindlist(plate_tabs)$time_sec), ", ",
       max(rbindlist(plate_tabs)$time_sec), "]\n",
       "  time_frame.csv covers:  [0, ", max(time_frame$End), "]\n",
       "Edit time_frame.csv so its End values cover this experiment's ",
       "phase boundaries (in seconds).")
}


# ── SPLIT Condition INTO PER-VARIABLE COLUMNS ───────────────────────────────
# Forward direction (variable values -> Condition string) is the user-facing
# convention used in plate_layout.csv. We reverse it here by splitting each
# Condition on "+" in the same order declared in design.csv.
cond_parts <- strsplit(full_df$Condition, "+", fixed = TRUE)
for (vi in seq_along(var_names)) {
  full_df[[var_names[vi]]] <- vapply(cond_parts, `[`, character(1), vi)
}

# Persist the aggregated long-format frame for inspection / re-use.
fwrite(full_df, file.path(OUTPUT_DIR, "aggregated_data.csv"))
cat("Saved: aggregated_data.csv  (", nrow(full_df), " rows)\n\n", sep = "")


# =============================================================================
# BAM MODEL
# =============================================================================

# Coerce factors, relevel against the reference declared in design.csv,
# compute time_in_group (per-Group elapsed time), and mark AR(1) start events
# at animal boundaries AND phase-group boundaries within an animal.
# Mutate full_df in place (it is already a data.table) and alias as gam_df.
# Avoids the deep copy that as_tibble() + a dplyr chain would create on top
# of the existing data.table — peak memory during construction was the prior
# copy of the whole frame.
setDT(full_df)
full_df[, `:=`(
  plate           = as.factor(plate),
  animal_id       = as.factor(animal_id),
  Condition_Combo = as.factor(Condition),
  Group_factor    = as.factor(Group)
)]
full_df[, time_in_group := as.integer(time_sec - min(time_sec)), by = Group]
setorder(full_df, animal_id, Group, time_in_group)
full_df[, start_event := (seq_len(.N) == 1L) |
                        (Group != shift(Group, fill = -1L)),
        by = animal_id]
gam_df <- full_df

for (v in var_names) {
  gam_df[[v]] <- as.factor(gam_df[[v]])
  if (!ref_map[[v]] %in% levels(gam_df[[v]])) {
    stop("Reference level '", ref_map[[v]], "' for variable '", v,
         "' not found in the data. Levels present: ",
         paste(levels(gam_df[[v]]), collapse = ", "))
  }
  gam_df[[v]] <- relevel(gam_df[[v]], ref = ref_map[[v]])
}
if (!ref_condition %in% levels(gam_df$Condition_Combo)) {
  stop("Reference condition '", ref_condition,
       "' not found among Conditions in the data. Levels present: ",
       paste(levels(gam_df$Condition_Combo), collapse = ", "))
}
gam_df$Condition_Combo <- relevel(gam_df$Condition_Combo, ref = ref_condition)

n_groups <- length(unique(gam_df$Group))
n_plates <- length(unique(gam_df$plate))
# When there is only one phase-group in time_frame.csv, Group_factor becomes
# a 1-level factor — including it in the parametric formula then triggers
# "contrasts not defined for 0 degrees of freedom" inside bam(). Drop it
# from the formula and from emmeans specs in that case (this mirrors the
# Arch A path in TweedieAR1 BAM.R).
use_group_factor <- n_groups >= 2L
single_group_id  <- if (!use_group_factor) {
  as.integer(unique(gam_df$Group))[1]
} else NA_integer_
cat("Phase groups: ", n_groups,
    if (!use_group_factor) "  (Group_factor will be dropped from the formula)" else "",
    "\n", sep = "")
cat("Plates:       ", n_plates, "\n")
cat("Animals:      ", length(unique(gam_df$animal_id)), "\n")
cat("Conditions:   ", paste(levels(gam_df$Condition_Combo), collapse = ", "), "\n\n")


# ── BUILD THE FORMULA ───────────────────────────────────────────────────────
# Mirrors the unified-model formula from TweedieAR1 BAM.R (Arch C, sz random
# basis - the production-winning configuration per ABLATION_METHODOLOGY.md).
#
# Parametric:  Var1 * Var2 * ... [* Group_factor]   (full factorial)
# Smooths:
#   s(time_in_group, k = k_start)
#     - shared within-phase trajectory baseline
#   s(time_in_group, by = ..., k = k_start)
#     - per-(condition [× phase-group]) deviation from baseline
#   s(time_in_group, plate, bs = 'sz')      (omitted when only one plate)
#   s(time_in_group, animal_id, bs = 'sz')  - per-animal trajectory
k_start <- 30

parametric_str <- if (use_group_factor) {
  paste(c(var_names, "Group_factor"), collapse = " * ")
} else {
  paste(var_names, collapse = " * ")
}
by_term <- if (use_group_factor) {
  "interaction(Condition_Combo, Group_factor)"
} else {
  "Condition_Combo"
}
plate_smooth  <- if (n_plates >= 2) " + s(time_in_group, plate, bs = 'sz')" else ""
animal_smooth <- " + s(time_in_group, animal_id, bs = 'sz')"

formula_str <- paste0(
  "pixel_diff ~ ", parametric_str,
  " + s(time_in_group, k = ", k_start, ")",
  " + s(time_in_group, by = ", by_term, ", k = ", k_start, ")",
  plate_smooth, animal_smooth
)
cat("BAM formula:\n  ", formula_str, "\n\n", sep = "")
bam_formula <- as.formula(formula_str)

usable_cores <- max(1L, detectCores() - 1L)


# ── MEMORY ESTIMATE ─────────────────────────────────────────────────────────
# Back-of-envelope estimate of bam()'s design matrix size, so we can flag
# OOM-risk runs BEFORE the long pilot fit instead of after. Estimates are
# approximate — mgcv applies sum-to-zero constraints and other reductions
# internally, and the discrete method's actual peak depends on the unique-
# row compression. As a rule of thumb the true peak is 3-5x the X-matrix
# size; we report 4x.
.n_param <- prod(vapply(var_names,
                        function(v) nlevels(gam_df[[v]]), integer(1)))
if (use_group_factor) .n_param <- .n_param * nlevels(gam_df$Group_factor)

.n_baseline   <- k_start
.n_by_smooth  <- k_start * nlevels(gam_df$Condition_Combo)
if (use_group_factor) .n_by_smooth <- .n_by_smooth * nlevels(gam_df$Group_factor)
.n_plate_sm   <- if (n_plates >= 2L) k_start * n_plates else 0L
.n_animal_sm  <- k_start * length(unique(gam_df$animal_id))

.p_total <- .n_param + .n_baseline + .n_by_smooth + .n_plate_sm + .n_animal_sm
.n_rows  <- nrow(gam_df)

.X_mb         <- as.numeric(.n_rows) * as.numeric(.p_total) * 8 / 1e6
.peak_est_mb  <- 4 * .X_mb

cat("── BAM memory estimate ──\n")
cat("  Rows (n):              ", format(.n_rows,  big.mark = ","), "\n", sep = "")
cat("  Columns (p, approx):   ", format(.p_total, big.mark = ","), "\n", sep = "")
cat("    parametric:          ", format(.n_param,      big.mark = ","), "\n", sep = "")
cat("    baseline smooth:     ", format(.n_baseline,   big.mark = ","), "\n", sep = "")
cat("    by-condition smooth: ", format(.n_by_smooth,  big.mark = ","), "\n", sep = "")
cat("    plate smooth:        ", format(.n_plate_sm,   big.mark = ","), "\n", sep = "")
cat("    animal smooth:       ", format(.n_animal_sm,  big.mark = ","), "\n", sep = "")
cat("  Design matrix X:       ~", format(round(.X_mb), big.mark = ","), " MB\n", sep = "")
cat("  Estimated bam() peak:  ~", format(round(.peak_est_mb), big.mark = ","),
    " MB (~4x X)\n", sep = "")
if (.peak_est_mb > 2000) {
  warning("BAM peak memory estimate is ~", round(.peak_est_mb / 1000, 1),
          " GB. If Pass 1 OOMs, lower k_start (currently ", k_start,
          ") or simplify the formula. Conditions = ",
          nlevels(gam_df$Condition_Combo),
          " — every extra Condition_Combo level multiplies the by-smooth's ",
          k_start, " basis columns.")
}
cat("\n")


# ── TRIM gam_df + FREE INTERMEDIATES ────────────────────────────────────────
# Drop columns the BAM formula does not reference. bam() copies its `data`
# argument internally; trimming reduces that copy. aggregated_data.csv already
# persists the untrimmed frame (line above) so Condition / loc_coord / Phase /
# Group / location / loc_id remain available on disk for inspection.
bam_cols <- c("pixel_diff", "time_in_group",
              "plate", "animal_id", "Condition_Combo",
              "start_event", var_names)
if (use_group_factor) bam_cols <- c(bam_cols, "Group_factor")
gam_df <- gam_df[, bam_cols, with = FALSE]

# Free the per-plate list and the per-row condition split — neither is needed
# again. full_df is the same object as gam_df after the in-place build above,
# so removing its name does not free memory but does prevent accidental reuse.
rm(plate_tabs, cond_parts, full_df); gc(verbose = FALSE)


# ── PEAK-RSS HELPER ─────────────────────────────────────────────────────────
# Reads the OS-level peak working set for this R process. Captures everything
# resident — R heap, mgcv C buffers, BLAS scratch — so it reports the number
# that decides OOM. Windows path uses PowerShell's PeakWorkingSet64; Linux
# reads VmHWM from /proc/self/status; macOS has no kernel-tracked peak, so we
# fall back to the current RSS (which under-reports the true high-water mark).
peak_rss_mb <- function() {
  if (.Platform$OS.type == "windows") {
    out <- tryCatch(
      system2("powershell",
              c("-NoProfile", "-Command",
                sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid())),
              stdout = TRUE, stderr = FALSE),
      error = function(e) NA_character_)
    if (length(out) == 0 || is.na(out[1])) return(NA_real_)
    as.numeric(out[1]) / 1e6
  } else if (file.exists("/proc/self/status")) {
    hwm <- grep("^VmHWM:",
                readLines("/proc/self/status"), value = TRUE)
    if (length(hwm) == 0) return(NA_real_)
    as.numeric(sub(".*?(\\d+).*", "\\1", hwm)) / 1024
  } else {
    out <- tryCatch(
      system2("ps", c("-o", "rss=", "-p", Sys.getpid()), stdout = TRUE),
      error = function(e) NA_character_)
    if (length(out) == 0 || is.na(out[1])) return(NA_real_)
    as.numeric(trimws(out[1])) / 1024
  }
}
.fmt_mb <- function(x) if (is.na(x)) "  n/a" else format(round(x), big.mark = ",")
.report_rss <- function(label, before, after) {
  delta <- after - before
  cat(sprintf("  %-22s before: %s MB   after: %s MB   delta: %s MB\n",
              label, .fmt_mb(before), .fmt_mb(after),
              if (is.na(delta)) "n/a"
              else sprintf("%s%s", if (delta >= 0) "+" else "", .fmt_mb(delta))))
}


# ── PASS 1: PILOT FIT (NO AR) → ESTIMATE RHO ────────────────────────────────
cat("Pass 1: pilot fit to estimate AR(1) rho ...\n")
fallback_rho <- 0.20
.rss_p1_before <- peak_rss_mb()
optimal_rho <- tryCatch({
  pilot <- bam(
    formula  = bam_formula, data = gam_df,
    family   = tw(),
    select   = TRUE, method = "fREML",
    discrete = TRUE, nthreads = usable_cores,
    chunk.size = 10000
  )
  rho_est <- acf(resid(pilot), plot = FALSE)$acf[2]
  rho_est <- max(min(rho_est, 0.95), -0.95)
  rm(pilot); gc(verbose = FALSE)
  rho_est
}, error = function(e) {
  warning("Pass 1 failed (", conditionMessage(e), ") — using fallback rho = ",
          fallback_rho)
  fallback_rho
})
.rss_p1_after <- peak_rss_mb()
.report_rss("Pass 1 peak RSS", .rss_p1_before, .rss_p1_after)
cat(sprintf("  Estimated rho = %.4f\n\n", optimal_rho))


# ── PASS 2: REFIT WITH AR(1) ────────────────────────────────────────────────
cat("Pass 2: refit with AR(1) correction ...\n")
.rss_p2_before <- peak_rss_mb()
final_model <- bam(
  formula  = bam_formula, data = gam_df,
  family   = tw(),
  rho      = optimal_rho, AR.start = gam_df$start_event,
  select   = TRUE, method = "fREML",
  discrete = TRUE, nthreads = usable_cores,
  chunk.size = 10000
)
.rss_p2_after <- peak_rss_mb()
.report_rss("Pass 2 peak RSS", .rss_p2_before, .rss_p2_after)
cat("  Done.\n\n")


# =============================================================================
# CONTRASTS
# =============================================================================
# emmeans is always run with `adjust = "none"` here — multiple-testing
# correction is applied post-hoc via apply_correction() so the same raw
# p-values can be re-corrected under a different strategy without refitting.

# Reference-grid size can grow with (vars × Group) — bump rg.limit generously.
factor_sizes <- vapply(c(var_names, "Group_factor"),
                       function(v) nlevels(gam_df[[v]]), integer(1))
emm_options(rg.limit = max(50000L, as.integer(prod(factor_sizes) * 2L)))

all_contrasts_list <- list()


# ── Per-variable effects ────────────────────────────────────────────────────
# For each variable V, contrast each non-reference level vs the reference,
# conditioning on (all other variables [× Group_factor]). The Test_Family
# tag is "<V>_Effect".
for (vi in seq_along(var_names)) {
  focal      <- var_names[vi]
  other_vars <- var_names[-vi]
  ref_val    <- ref_map[[focal]]

  cond_vars <- if (use_group_factor) c(other_vars, "Group_factor") else other_vars
  spec_str  <- if (length(cond_vars) > 0) {
    paste0("~ ", focal, " | ", paste(cond_vars, collapse = " + "))
  } else {
    paste0("~ ", focal)
  }
  emm_base <- emmeans(final_model, specs = as.formula(spec_str))
  emm_c    <- contrast(emm_base, method = "trt.vs.ctrl",
                       ref = ref_val, adjust = "none")

  df_c <- as.data.frame(emm_c) %>%
    mutate(
      Test_Family      = paste0(focal, "_Effect"),
      Passed_Contrasts = as.character(contrast)
    )
  if (use_group_factor) {
    df_c$Group <- as.integer(as.character(df_c$Group_factor))
    df_c <- df_c %>% select(-Group_factor)
  } else {
    df_c$Group <- single_group_id
  }

  if (length(other_vars) > 0) {
    split_cols <- as.data.frame(df_c[, other_vars, drop = FALSE])
    df_c$Split_By <- apply(split_cols, 1, function(r) paste(r, collapse = " + "))
  } else {
    # Match the app's sentinel — tree_metadata.R compares Split_By literally
    # against a tree spec's split_by value, so "All" must match "All", not NA.
    df_c$Split_By <- "All"
  }

  df_c <- df_c %>% select(-any_of(other_vars))
  all_contrasts_list[[paste0("var_", focal)]] <- df_c
}


# ── Full_Interaction: all-pairs vs experimental-vs-control ──────────────────
# IMPORTANT: build the spec from the DESIGN variables (Var1 * Var2 * ...),
# NOT from Condition_Combo. Condition_Combo only appears inside a smooth
# term (the by-interaction smooth), so emmeans strips it from the reference
# grid and errors with "No variable named Condition_Combo in the reference
# grid". Building from the design variables yields the same grid cells; we
# reconstruct the "+"-joined condition labels from the grid afterwards.
inter_spec <- if (use_group_factor) {
  as.formula(paste0("~ ", paste(var_names, collapse = " * "), " | Group_factor"))
} else {
  as.formula(paste0("~ ", paste(var_names, collapse = " * ")))
}
emm_inter_base <- emmeans(final_model, specs = inter_spec)

# Reconstruct condition strings from the (within-group) grid. Rows here match
# the ordering that emmeans's contrast weight vectors expect.
full_inter_grid <- emm_inter_base@grid
within_grid <- if (use_group_factor) {
  first_group <- levels(gam_df$Group_factor)[1]
  full_inter_grid[as.character(full_inter_grid$Group_factor) == first_group,
                  var_names, drop = FALSE]
} else {
  full_inter_grid[, var_names, drop = FALSE]
}
cond_strs <- apply(within_grid, 1,
                   function(r) paste(as.character(r), collapse = "+"))

custom_contrasts <- list()
if (ALL_PAIRS) {
  cat("Full_Interaction: all pairwise contrasts (ALL_PAIRS = TRUE).\n")
  pairs <- combn(seq_along(cond_strs), 2)
  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1, k]; j <- pairs[2, k]
    v <- numeric(length(cond_strs)); v[i] <- 1; v[j] <- -1
    custom_contrasts[[paste(cond_strs[i], "-", cond_strs[j])]] <- v
  }
} else {
  cat("Full_Interaction: experimental-vs-control against '", ref_condition,
      "' (ALL_PAIRS = FALSE).\n", sep = "")
  ref_idx <- which(cond_strs == ref_condition)
  if (length(ref_idx) != 1L) {
    stop("Reference condition '", ref_condition,
         "' not found in the emmeans grid. Grid cells: ",
         paste(cond_strs, collapse = ", "))
  }
  for (i in seq_along(cond_strs)) {
    if (i == ref_idx) next
    v <- numeric(length(cond_strs))
    v[i] <- 1; v[ref_idx] <- -1
    custom_contrasts[[paste(cond_strs[i], "-", cond_strs[ref_idx])]] <- v
  }
}

emm_inter <- if (use_group_factor) {
  contrast(emm_inter_base, method = custom_contrasts,
           by = "Group_factor", adjust = "none")
} else {
  contrast(emm_inter_base, method = custom_contrasts, adjust = "none")
}

df_inter <- as.data.frame(emm_inter) %>%
  mutate(
    Test_Family      = "Full_Interaction",
    Passed_Contrasts = as.character(contrast),
    # Match the app's sentinel — tree_metadata.R compares Split_By literally
    # against a tree spec's split_by value, so "None" must match "None", not NA.
    Split_By         = "None"
  )
if (use_group_factor) {
  df_inter$Group <- as.integer(as.character(df_inter$Group_factor))
  df_inter <- df_inter %>% select(-Group_factor)
} else {
  df_inter$Group <- single_group_id
}
all_contrasts_list[["interaction"]] <- df_inter


# =============================================================================
# CORRECTION + OUTPUT
# =============================================================================

master_results <- bind_rows(all_contrasts_list) %>%
  mutate(
    raw_pvalue = p.value,
    # Tweedie's log link means estimates are on the log scale — divide by log(2)
    # to express them as log2 ratios (matches the BAM script's convention).
    estimate   = estimate / log(2),
    SE         = SE       / log(2)
  ) %>%
  select(Test_Family, Group, Split_By, Passed_Contrasts,
         estimate, SE, raw_pvalue) %>%
  arrange(Test_Family, Group, Passed_Contrasts)

cat("\nApplying correction:\n")
str(CORRECTION_SPEC)

# Tree strategies need tree_id / tree_level / tree_parent_id BEFORE correction.
if (identical(CORRECTION_SPEC$strategy, "tree")) {
  tree_spec <- CORRECTION_SPEC$tree
  if (is.null(tree_spec) || length(tree_spec) == 0) {
    cat("Tree strategy with no spec — using default rescue tree.\n")
    tree_spec <- default_rescue_tree_spec(var_names, ref_map)
  }
  master_results <- assign_tree_metadata(master_results, tree_spec, var_names)
}

final_master <- apply_correction(master_results, CORRECTION_SPEC) %>%
  select(Test_Family, Group, Split_By, Passed_Contrasts,
         estimate, SE, raw_pvalue, adjusted_pvalue, correction_method,
         tidyselect::any_of(c("tree_id", "tree_level",
                              "tree_parent_id", "tree_status")))

cat("\n── Master results ──\n")
print(as.data.frame(final_master), row.names = FALSE)

write_csv(final_master, file.path(OUTPUT_DIR, "master_results.csv"))
cat("\nSaved: master_results.csv\n")


# ── SUMMARY TABLE ───────────────────────────────────────────────────────────
summary_table <- final_master %>%
  group_by(Test_Family, Group) %>%
  summarise(
    n_total     = n(),
    n_sig       = sum(adjusted_pvalue < CORRECTION_SPEC$threshold, na.rm = TRUE),
    mean_effect = mean(abs(estimate), na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  pivot_wider(
    names_from  = Group,
    values_from = c(n_sig, mean_effect),
    names_glue  = "G{Group}_{.value}"
  )

write_csv(summary_table, file.path(OUTPUT_DIR, "summary_statistics_by_group.csv"))
cat("Saved: summary_statistics_by_group.csv\n")


# =============================================================================
# POSTERIOR EQUIVALENCE TESTING (Bayesian)
# =============================================================================
# For each experimental-vs-reference Condition_Combo pair × phase group, draws
# coefficient samples from the BAM's large-sample Gaussian posterior, builds
# the population-level Δη(t) trajectory (excluding the random-effect 'sz'
# smooths), and computes M = max_t |Δη(t)| / log(2). Reports posterior
# median, 95% CI, and Pr(M < delta). Mirrors the app's TweedieAR1 BAM.R PE
# block; differences from the app:
#   * pairs are always experimental-vs-reference (no contrast_selection.json)
#   * random_basis is hardcoded to 'sz' (the standalone formula's choice)

.posterior_equivalence_pair <- function(model, full_df, lhs, rhs, var_names,
                                        group_id, delta_log2, n_draws,
                                        n_time = 200, seed = 42) {
  cond_levels <- levels(full_df$Condition_Combo)
  if (!(lhs %in% cond_levels) || !(rhs %in% cond_levels)) return(NULL)

  group_str <- as.character(group_id)
  if (use_group_factor) {
    group_levels <- levels(full_df$Group_factor)
    if (!(group_str %in% group_levels)) return(NULL)
    group_rows <- full_df[as.character(full_df$Group_factor) == group_str, ]
  } else {
    group_rows <- full_df
  }
  if (nrow(group_rows) == 0) return(NULL)

  parse_cond <- function(cond_str) {
    parts <- strsplit(cond_str, "+", fixed = TRUE)[[1]]
    if (length(parts) != length(var_names)) {
      stop("Cannot parse '", cond_str, "': expected ", length(var_names),
           " parts, got ", length(parts))
    }
    setNames(as.list(parts), var_names)
  }
  lhs_vars <- parse_cond(lhs)
  rhs_vars <- parse_cond(rhs)

  time_grid <- seq(0, max(group_rows$time_in_group), length.out = n_time)

  build_nd <- function(cond_str, vars_list) {
    nd <- data.frame(
      time_in_group = time_grid
    )
    for (v in var_names) {
      nd[[v]] <- factor(vars_list[[v]], levels = levels(full_df[[v]]))
    }
    nd$Condition_Combo <- factor(cond_str, levels = cond_levels)
    if (use_group_factor) {
      nd$Group_factor <- factor(group_str, levels = group_levels)
    }
    # animal_id / plate must be present for predict() but are excluded from
    # the linear predictor below — anchor to any real level in the model.
    nd$animal_id <- factor(full_df$animal_id[1],
                           levels = levels(full_df$animal_id))
    nd$plate <- factor(full_df$plate[1],
                       levels = levels(full_df$plate))
    nd
  }
  nd_lhs <- build_nd(lhs, lhs_vars)
  nd_rhs <- build_nd(rhs, rhs_vars)

  # Standalone hardcodes bs='sz' for animal (and plate when n_plates >= 2),
  # so exclude exactly those terms to get population-level draws.
  excl <- "s(time_in_group,animal_id)"
  if (n_plates >= 2L) excl <- c(excl, "s(time_in_group,plate)")

  X_lhs  <- predict(model, newdata = nd_lhs, type = "lpmatrix", exclude = excl)
  X_rhs  <- predict(model, newdata = nd_rhs, type = "lpmatrix", exclude = excl)
  X_diff <- X_lhs - X_rhs

  # Bayesian (large-sample Gaussian) posterior of β.
  set.seed(seed)
  Cv <- vcov(model, unconditional = TRUE)
  beta_samples <- MASS::mvrnorm(n_draws, mu = coef(model), Sigma = Cv)

  diff_samples <- X_diff %*% t(beta_samples)
  M_log2 <- apply(abs(diff_samples), 2, max) / log(2)

  list(
    M_log2   = M_log2,
    median   = median(M_log2),
    lower95  = unname(quantile(M_log2, 0.025)),
    upper95  = unname(quantile(M_log2, 0.975)),
    pr_equiv = mean(M_log2 < delta_log2)
  )
}

if (isTRUE(POSTERIOR_EQUIVALENCE$enabled)) {
  delta_log2 <- as.numeric(POSTERIOR_EQUIVALENCE$delta   %||% 1.0)
  n_draws    <- as.integer(POSTERIOR_EQUIVALENCE$n_draws %||% 10000)

  # Option b: experimental-vs-reference pairs only.
  non_ref       <- setdiff(levels(gam_df$Condition_Combo), ref_condition)
  pe_kept_pairs <- lapply(non_ref, function(x) list(x, ref_condition))

  pe_group_ids <- if (use_group_factor) {
    sort(as.integer(levels(gam_df$Group_factor)))
  } else {
    single_group_id
  }

  cat("\n========================================\n")
  cat("Posterior Equivalence Testing (Bayesian)\n")
  cat(sprintf("  delta            = %.3f (log_2 fold change)\n", delta_log2))
  cat(sprintf("  posterior draws  = %d\n", n_draws))
  cat(sprintf("  pairs            = %d (experimental-vs-reference)\n",
              length(pe_kept_pairs)))
  cat(sprintf("  groups           = %d\n", length(pe_group_ids)))
  cat("========================================\n")

  .rss_pe_before <- peak_rss_mb()
  pe_summary_rows <- list()
  pe_draws_named  <- list()

  for (g in pe_group_ids) {
    g_str <- as.character(g)
    cat(sprintf("Group %s: posterior equivalence ...\n", g_str))
    for (i in seq_along(pe_kept_pairs)) {
      pair <- pe_kept_pairs[[i]]
      lhs  <- pair[[1]]; rhs <- pair[[2]]
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
  .rss_pe_after <- peak_rss_mb()
  .report_rss("PE peak RSS", .rss_pe_before, .rss_pe_after)

  if (length(pe_summary_rows) > 0) {
    pe_summary_df <- bind_rows(pe_summary_rows) %>%
      arrange(Group, pair) %>%
      as_tibble()
    cat("\n── Posterior Equivalence Summary ──\n")
    print(as.data.frame(pe_summary_df), row.names = FALSE)
    write_csv(pe_summary_df,
              file.path(OUTPUT_DIR, "posterior_equivalence_summary.csv"))
    cat("\nSaved: posterior_equivalence_summary.csv\n")

    draws_df <- bind_rows(lapply(seq_along(pe_draws_named), function(j) {
      key <- names(pe_draws_named)[j]
      m   <- regmatches(key, regexec("^g(.+?)__(.+?)__(.+)$", key))[[1]]
      data.frame(
        Group  = as.integer(m[2]),
        pair   = paste(m[3], "vs", m[4]),
        M_log2 = pe_draws_named[[j]],
        stringsAsFactors = FALSE
      )
    }))
    write_csv(draws_df,
              file.path(OUTPUT_DIR, "posterior_equivalence_draws.csv"))
    cat("Saved: posterior_equivalence_draws.csv  (", nrow(draws_df),
        " draws total)\n", sep = "")
  } else {
    cat("No posterior-equivalence results were produced (all pairs skipped).\n")
  }
}

cat("\nDone.\n")
