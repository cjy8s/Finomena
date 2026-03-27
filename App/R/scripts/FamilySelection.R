# ── FamilySelection.R ──────────────────────────────────────────────────────────
# Optional distribution family selection test for Finomena GAMM analysis.
# Compares Tweedie, Gamma(log), and Negative Binomial using a simplified model
# with a fixed rho = 0.25 (assumed AR1 correlation) for speed.
#
# Outputs (to output_dir):
#   family_selection_results.csv  — AIC/BIC/dispersion/deviance table
#   family_selection_winner.csv   — winning family name and shift value
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
    "Open the app, go to the GAMM Analysis tab, and click 'Install R Packages'."
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
gam_df <- full_df %>%
  mutate(
    animal_id       = as.factor(paste(plate, location, sep = "_")),
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

# ── Simplified test formula ─────────────────────────────────────────────────────
# Omits by-condition and factor smooths to keep fitting fast; sufficient for
# comparing distributional assumptions.
parametric_str <- paste(var_names, collapse = " * ")
test_formula   <- as.formula(paste0(
  "pxl_diff_test ~ ", parametric_str, " + s(time_in_group, k = 10)"
))
cat("Test formula:", deparse(test_formula), "\n\n")

usable_cores <- max(1, detectCores() - 1)

# ── Candidate families ──────────────────────────────────────────────────────────
candidates <- list(
  Tweedie     = list(family = tw(),                shift = 0),
  Gamma       = list(family = Gamma(link = "log"), shift = shift_val),
  NegBinomial = list(family = nb(),                shift = 0)
)

results <- list()

for (fname in names(candidates)) {
  cat("========================================\n")
  cat("Testing family:", fname, "\n")

  cand                  <- candidates[[fname]]
  gam_df$pxl_diff_test  <- gam_df$pxl_diff + cand$shift

  tryCatch({
    m <- bam(
      test_formula,
      data     = gam_df,
      family   = cand$family,
      rho      = 0.25,
      AR.start = gam_df$start_event,
      method   = "fREML",
      discrete = TRUE,
      nthreads = usable_cores
    )

    aic_val  <- AIC(m)
    bic_val  <- BIC(m)

    # Dispersion: Pearson chi-sq / residual df — well-fitted model gives ~1.0
    p_resid <- residuals(m, type = "pearson")
    disp    <- sum(p_resid^2, na.rm = TRUE) / m$df.residual

    # Deviance explained (%)
    dev_expl <- summary(m)$dev.expl * 100

    # Zero-proportion matching: compare observed vs model-predicted
    zero_match_str <- switch(fname,
      Tweedie = {
        # Compound Poisson-Gamma: P(Y=0) = exp(-mu^(2-p) / (phi*(2-p)))
        p_hat  <- m$family$getTheta(TRUE)   # Tweedie power parameter (1 < p < 2)
        phi    <- m$sig2                    # scale/dispersion parameter
        mu_hat <- fitted(m)
        pred_zero <- mean(exp(-mu_hat^(2 - p_hat) / (phi * (2 - p_hat))), na.rm = TRUE)
        sprintf("Obs: %.1f%% | Pred: %.1f%%", zero_prop_obs * 100, pred_zero * 100)
      },
      NegBinomial = {
        # NB2: P(Y=0) = (theta / (theta + mu))^theta
        theta  <- m$family$getTheta(TRUE)   # size/overdispersion parameter
        mu_hat <- fitted(m)
        pred_zero <- mean((theta / (theta + mu_hat))^theta, na.rm = TRUE)
        sprintf("Obs: %.1f%% | Pred: %.1f%%", zero_prop_obs * 100, pred_zero * 100)
      },
      Gamma = sprintf("N/A (shift of %.4f applied; Gamma requires >0)", cand$shift)
    )

    cat(sprintf(
      "  AIC:          %.1f\n  BIC:          %.1f\n  Dispersion:   %.3f\n  Dev. Expl.:   %.1f%%\n  Zero match:   %s\n",
      aic_val, bic_val, disp, dev_expl, zero_match_str
    ))

    results[[fname]] <- data.frame(
      Family        = fname,
      AIC           = round(aic_val, 1),
      BIC           = round(bic_val, 1),
      Dispersion    = round(disp, 3),
      Dev_Explained = round(dev_expl, 1),
      Zero_Match    = zero_match_str,
      Shift_Applied = cand$shift,
      stringsAsFactors = FALSE
    )

  }, error = function(e) {
    msg <- conditionMessage(e)
    cat("  FAILED:", msg, "\n")
    results[[fname]] <<- data.frame(
      Family        = fname,
      AIC           = Inf,
      BIC           = Inf,
      Dispersion    = NA_real_,
      Dev_Explained = NA_real_,
      Zero_Match    = paste("FAILED:", msg),
      Shift_Applied = cand$shift,
      stringsAsFactors = FALSE
    )
  })
}

# ── Select winner ────────────────────────────────────────────────────────────────
results_df <- bind_rows(results)
valid      <- results_df %>% filter(is.finite(AIC))

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
