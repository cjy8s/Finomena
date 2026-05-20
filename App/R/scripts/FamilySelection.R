# ── FamilySelection.R ──────────────────────────────────────────────────────────
# Optional distribution family selection test for Finomena BAM analysis.
# Compares Tweedie, Gamma(log), and Negative Binomial using the full production
# BAM formula (by-condition smooth, plate factor smooth, per-animal factor
# smooth) fit per phase group. AR1 rho is fixed at 0.25 to halve compute (vs the
# production two-pass rho optimisation); this is the only deviation from the
# downstream model. Family ranking is robust to small rho mis-specification.
#
# Outputs (to output_dir):
#   family_selection_results.csv   — per-family aggregated AIC/BIC/dispersion table
#   family_selection_by_group.csv  — per-group × per-family breakdown
#   family_selection_winner.csv    — winning family name and shift value
#
# Called from the app as:
#   Rscript "FamilySelection.R" <input_csv> <output_dir> <var_names> <ref_values> <ref_condition>
# ──────────────────────────────────────────────────────────────────────────────

# ── Package check ──────────────────────────────────────────────────────────────
cat("R library path:", paste(.libPaths(), collapse = "\n               "), "\n")
required_pkgs <- c("tidyverse", "data.table", "mgcv", "parallel")
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

# ── Argument parsing ───────────────────────────────────────────────────────────
args           <- commandArgs(trailingOnly = TRUE)
input_csv      <- args[1]
output_dir     <- if (length(args) >= 2 && !is.na(args[2])) args[2] else dirname(input_csv)
var_names_raw  <- if (length(args) >= 3 && !is.na(args[3]) && nchar(args[3]) > 0) args[3] else "Genotype,Drug"
ref_values_raw <- if (length(args) >= 4 && !is.na(args[4]) && nchar(args[4]) > 0) args[4] else "WT,DMSO"
ref_condition  <- if (length(args) >= 5 && !is.na(args[5]) && nchar(args[5]) > 0) args[5] else "WT+DMSO"

if (is.na(input_csv) || !file.exists(input_csv)) {
  stop("ERROR: input_csv argument missing or file not found: ", input_csv)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Input CSV:  ", input_csv,  "\n")
cat("Output dir: ", output_dir, "\n")

# ── Variable parsing ────────────────────────────────────────────────────────────
var_names_raw_vec <- str_split(var_names_raw, ",")[[1]]
ref_values        <- str_split(ref_values_raw, ",")[[1]]
n_vars            <- length(var_names_raw_vec)
var_names         <- make.names(var_names_raw_vec, unique = TRUE)
ref_map           <- setNames(ref_values[1:n_vars], var_names)

cat("Variables:  ", paste(var_names, collapse = ", "), "\n")
cat("References: ", paste(ref_map, collapse = ", "), "\n\n")

# ── Data loading ───────────────────────────────────────────────────────────────
full_df <- read_csv(input_csv, show_col_types = FALSE) %>%
  rename(pxl_diff = pixel_diff) %>%
  mutate(
    Condition = as.character(Condition),
    Phase     = as.character(Phase),
    Group     = as.integer(Group),
    time_sec  = as.integer(time_sec),
    pxl_diff  = as.numeric(pxl_diff),
    plate     = as.character(plate),
    location  = as.character(location)
  ) %>%
  filter(!is.na(Condition), Condition != "")

# ── Data prep ──────────────────────────────────────────────────────────────────
# plate is converted to a factor so it can serve as a level for the
# s(time_in_group, plate, bs='fs') smooth, matching production BAM.
gam_df <- full_df %>%
  mutate(
    animal_id       = as.factor(paste(plate, location, sep = "_")),
    plate           = as.factor(plate),
    Condition_Combo = as.factor(Condition)
  ) %>%
  arrange(animal_id, time_sec) %>%
  group_by(animal_id) %>%
  mutate(
    time_in_group = time_sec - min(time_sec),
    start_event   = row_number() == 1
  ) %>%
  ungroup()

if (!all(var_names_raw_vec %in% colnames(gam_df))) {
  gam_df <- gam_df %>%
    separate(Condition, into = var_names, sep = "\\+", remove = FALSE)
} else {
  for (i in seq_along(var_names)) {
    if (var_names_raw_vec[i] != var_names[i] && var_names_raw_vec[i] %in% colnames(gam_df)) {
      gam_df <- gam_df %>% rename(!!var_names[i] := !!var_names_raw_vec[i])
    }
  }
}

for (v in var_names) {
  gam_df[[v]] <- relevel(as.factor(gam_df[[v]]), ref = ref_map[[v]])
}
gam_df$Condition_Combo <- relevel(as.factor(gam_df$Condition_Combo), ref = ref_condition)

# ── Observed zero proportion ────────────────────────────────────────────────────
zero_prop_obs <- mean(gam_df$pxl_diff == 0, na.rm = TRUE)
cat(sprintf("Observed zero proportion: %.2f%%\n", zero_prop_obs * 100))

# Shift so Gamma can run: shift minimum value up to 1 if any zeros/negatives
min_val   <- min(gam_df$pxl_diff, na.rm = TRUE)
shift_val <- if (min_val <= 0) (1 - min_val) else 0
cat(sprintf("Shift for Gamma family:   %.4f\n\n", shift_val))

# ── Full production-style formula ──────────────────────────────────────────────
# Mirrors the BAM script so the family comparison is on the model that will
# actually be used downstream. Rho is fixed (not optimised) for speed.
k_start        <- 30
fixed_rho      <- 0.25
parametric_str <- paste(var_names, collapse = " * ")
test_formula   <- as.formula(paste0(
  "pxl_diff_test ~ ", parametric_str,
  " + s(time_in_group, k = k_start)",
  " + s(time_in_group, by = Condition_Combo, k = k_start)",
  " + s(time_in_group, plate, bs = 'fs', m = 1)",
  " + s(time_in_group, animal_id, bs = 'fs', m = 1)"
))
cat("Test formula:", deparse(test_formula, width.cutoff = 200), "\n")
cat(sprintf("Fixed AR1 rho: %.2f (production optimises this per group)\n\n", fixed_rho))

usable_cores <- max(1, detectCores() - 1)

# ── Candidate families ──────────────────────────────────────────────────────────
candidates <- list(
  Tweedie     = list(family = tw(),                shift = 0),
  Gamma       = list(family = Gamma(link = "log"), shift = shift_val),
  NegBinomial = list(family = nb(),                shift = 0)
)

# ── Per-group × per-family fits ────────────────────────────────────────────────
# For each phase group, we set reference levels and AR1 start markers exactly
# the way the production BAM does, then fit each candidate family with the full
# formula. Metrics are collected per (group, family) and then aggregated.
unique_groups   <- sort(unique(gam_df$Group))
per_group_rows  <- list()   # detailed per-group metrics
fit_failures    <- list()   # family -> error message if any group failed

for (g in unique_groups) {
  cat("========================================\n")
  cat("Group:", g, "\n")

  group_data <- gam_df %>%
    filter(Group == g) %>%
    mutate(time_in_group = time_sec - min(time_sec)) %>%
    arrange(animal_id, time_in_group) %>%
    group_by(animal_id) %>%
    mutate(start_event = row_number() == 1) %>%
    ungroup()

  for (v in var_names) {
    group_data[[v]] <- relevel(as.factor(group_data[[v]]), ref = ref_map[[v]])
  }
  group_data$Condition_Combo <- relevel(as.factor(group_data$Condition_Combo), ref = ref_condition)
  group_data$plate           <- as.factor(group_data$plate)
  group_data$animal_id       <- as.factor(group_data$animal_id)

  group_zero_obs <- mean(group_data$pxl_diff == 0, na.rm = TRUE)
  group_n        <- nrow(group_data)

  for (fname in names(candidates)) {
    cat(sprintf("  Family: %-12s ", fname))
    cand <- candidates[[fname]]
    group_data$pxl_diff_test <- group_data$pxl_diff + cand$shift

    tryCatch({
      m <- bam(
        test_formula,
        data     = group_data,
        family   = cand$family,
        rho      = fixed_rho,
        AR.start = group_data$start_event,
        select   = TRUE,
        method   = "fREML",
        discrete = TRUE,
        nthreads = usable_cores
      )

      aic_val  <- AIC(m)
      bic_val  <- BIC(m)
      p_resid  <- residuals(m, type = "pearson")
      disp     <- sum(p_resid^2, na.rm = TRUE) / m$df.residual
      dev_expl <- summary(m)$dev.expl * 100

      pred_zero <- switch(fname,
        Tweedie = {
          p_hat  <- m$family$getTheta(TRUE)
          phi    <- m$sig2
          mu_hat <- fitted(m)
          mean(exp(-mu_hat^(2 - p_hat) / (phi * (2 - p_hat))), na.rm = TRUE)
        },
        NegBinomial = {
          theta  <- m$family$getTheta(TRUE)
          mu_hat <- fitted(m)
          mean((theta / (theta + mu_hat))^theta, na.rm = TRUE)
        },
        Gamma = NA_real_
      )

      cat(sprintf("AIC=%.1f  BIC=%.1f  disp=%.3f  dev=%.1f%%\n",
                  aic_val, bic_val, disp, dev_expl))

      per_group_rows[[length(per_group_rows) + 1]] <- data.frame(
        Group         = g,
        Family        = fname,
        n_obs         = group_n,
        AIC           = aic_val,
        BIC           = bic_val,
        Dispersion    = disp,
        Dev_Explained = dev_expl,
        Pred_Zero     = pred_zero,
        Obs_Zero      = group_zero_obs,
        Shift_Applied = cand$shift,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      msg <- conditionMessage(e)
      cat("FAILED:", msg, "\n")
      fit_failures[[fname]] <<- msg
      per_group_rows[[length(per_group_rows) + 1]] <<- data.frame(
        Group         = g,
        Family        = fname,
        n_obs         = group_n,
        AIC           = NA_real_,
        BIC           = NA_real_,
        Dispersion    = NA_real_,
        Dev_Explained = NA_real_,
        Pred_Zero     = NA_real_,
        Obs_Zero      = group_zero_obs,
        Shift_Applied = cand$shift,
        stringsAsFactors = FALSE
      )
    })
  }
}

per_group_df <- bind_rows(per_group_rows)
write_csv(per_group_df, file.path(output_dir, "family_selection_by_group.csv"))

# ── Aggregate per-family across groups ─────────────────────────────────────────
# AIC and BIC are summed because each group is an independent fit and the
# information criteria are additive in log-likelihood + penalty. Dispersion and
# Dev_Explained are weighted by n_obs for an interpretable single value. Zero
# match is reported as the pooled observed vs n_obs-weighted predicted rate.
results <- per_group_df %>%
  group_by(Family) %>%
  summarise(
    any_failed     = any(is.na(AIC)),
    AIC_sum        = sum(AIC, na.rm = TRUE),
    BIC_sum        = sum(BIC, na.rm = TRUE),
    disp_w         = sum(Dispersion * n_obs, na.rm = TRUE) / sum(n_obs[!is.na(Dispersion)]),
    dev_w          = sum(Dev_Explained * n_obs, na.rm = TRUE) / sum(n_obs[!is.na(Dev_Explained)]),
    pred_zero_w    = if (all(is.na(Pred_Zero))) NA_real_
                     else sum(Pred_Zero * n_obs, na.rm = TRUE) / sum(n_obs[!is.na(Pred_Zero)]),
    Shift_Applied  = first(Shift_Applied),
    .groups        = "drop"
  ) %>%
  mutate(
    Zero_Match = case_when(
      Family == "Gamma" ~ sprintf("N/A (shift of %.4f applied; Gamma requires >0)", Shift_Applied),
      is.na(pred_zero_w) ~ "FAILED — see by-group CSV",
      TRUE ~ sprintf("Obs: %.1f%% | Pred: %.1f%%", zero_prop_obs * 100, pred_zero_w * 100)
    ),
    AIC = ifelse(any_failed, Inf, round(AIC_sum, 1)),
    BIC = ifelse(any_failed, Inf, round(BIC_sum, 1)),
    Dispersion    = round(disp_w, 3),
    Dev_Explained = round(dev_w, 1)
  ) %>%
  select(Family, AIC, BIC, Dispersion, Dev_Explained, Zero_Match, Shift_Applied)

results_df <- results

# ── Select winner ────────────────────────────────────────────────────────────────
valid <- results_df %>% filter(is.finite(AIC))

if (nrow(valid) == 0) {
  stop("All family tests failed. Check the log above for details.")
}

winner_row   <- valid %>% slice_min(AIC, n = 1)
winner_name  <- winner_row$Family
winner_shift <- winner_row$Shift_Applied

cat("\n========================================\n")
cat("FAMILY SELECTION RESULTS:\n")
print(results_df, row.names = FALSE)
cat(sprintf("\nRECOMMENDED FAMILY: %s  (lowest AIC = %.1f)\n", winner_name, winner_row$AIC))
cat(sprintf("SHIFT VALUE:        %.4f\n", winner_shift))

# ── Save outputs ─────────────────────────────────────────────────────────────────
write_csv(results_df, file.path(output_dir, "family_selection_results.csv"))
write_csv(
  data.frame(
    family         = winner_name,
    shift          = winner_shift,
    zero_prop_obs  = round(zero_prop_obs * 100, 2)
  ),
  file.path(output_dir, "family_selection_winner.csv")
)

cat("\nFamily selection complete. Results saved to:", output_dir, "\n")
