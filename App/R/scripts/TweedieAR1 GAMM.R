# ── PACKAGE CHECK ─────────────────────────────────────────────────────────────
# Packages are managed by the Python app (installed into App/R/library/).
# R_LIBS is set by the app before launching this script, so R finds them there.
cat("R library path:", paste(.libPaths(), collapse = "\n               "), "\n")
required_pkgs <- c("tidyverse", "data.table", "mgcv", "parallel", "emmeans")
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
library(emmeans)

# ── ARGUMENT PARSING ──────────────────────────────────────────────────────────
# Called from the app as:
#  Rscript "TweedieAR1 GAMM.R" <input_csv> <output_dir> <var_names> <ref_values> <ref_condition> <global_corr> <contrast_adj> <roles> <family_name> <shift_val>
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
global_correction <- if (length(args) >= 6 && !is.na(args[6]) && nchar(args[6]) > 0) args[6] else "BH"
contrast_adjust   <- if (length(args) >= 7 && !is.na(args[7]) && nchar(args[7]) > 0) args[7] else "dunnett"
roles_raw         <- if (length(args) >= 8 && !is.na(args[8]) && nchar(args[8]) > 0) args[8] else ""
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

cat("Input CSV:        ", input_csv, "\n")
cat("Output dir:       ", output_dir, "\n")
cat("Variables (orig): ", paste(var_names_display, collapse = ", "), "\n")
cat("Variables (R):    ", paste(var_names, collapse = ", "), "\n")
cat("Reference values: ", paste(ref_map, collapse = ", "), "\n")
cat("Ref Condition:    ", ref_condition, "\n")
cat("Global correction:", global_correction, "\n")
cat("Contrast adjust:  ", contrast_adjust, "\n")
cat("Roles:            ", roles_raw, "\n")
cat("Family:           ", family_name, "\n")
cat("Shift applied:    ", shift_val, "\n")

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
gam_df <- full_df %>%
  mutate(
    animal_id = as.factor(paste(plate, location, sep = "_")),
    # Combined condition factor for interaction smooths
    Condition_Combo = as.factor(Condition)
  )

# Convert each variable to factor
for (v in var_names) {
  gam_df[[v]] <- as.factor(gam_df[[v]])
}

# Apply shift to pxl_diff if requested (required when using Gamma family)
if (shift_val != 0) {
  gam_df <- gam_df %>% mutate(pxl_diff = pxl_diff + shift_val)
  cat(sprintf("\nShift of %.4f applied to pxl_diff (family: %s).\n", shift_val, family_name))
}

# ---------------------------------------------------------
# 2. BUILD DYNAMIC GAMM FORMULA
# ---------------------------------------------------------
# Parametric: all main effects + all interactions  (Var1 * Var2 * ... * VarN)
# Smooth: s(time_in_group) + s(time_in_group, by=Condition_Combo) + s(time_in_group, animal_id, bs="fs")
parametric_str <- paste(var_names, collapse = " * ")
formula_str <- paste0(
  "pxl_diff ~ ", parametric_str,
  " + s(time_in_group, k = k_start)",
  " + s(time_in_group, by = Condition_Combo, k = k_start)",
  " + s(time_in_group, animal_id, bs = 'fs', m = 1)"
)
cat("\nGAMM formula:\n  ", formula_str, "\n\n")

models_by_group <- list()
unique_groups <- sort(unique(gam_df$Group))

# Define k value before edf calculation
k_start = 30

# Parse the formula
gamm_formula <- as.formula(formula_str)

for (g in unique_groups) {
  cat("\n========================================\n")
  cat("Processing Group:", g, "...\n")

  # A. Subset data and setup AR1 boundaries
  group_data <- gam_df %>%
    filter(Group == g) %>%
    mutate(time_in_group = time_sec - min(time_sec)) %>%
    arrange(animal_id, time_in_group) %>%
    group_by(animal_id) %>%
    mutate(start_event = row_number() == 1) %>%
    ungroup()

  # Ensure variables are factors and set reference levels
  for (v in var_names) {
    group_data[[v]] <- as.factor(group_data[[v]])
    group_data[[v]] <- relevel(group_data[[v]], ref = ref_map[[v]])
  }

  # Combined condition factor with reference level
  group_data$Condition_Combo <- as.factor(group_data$Condition_Combo)
  group_data$Condition_Combo <- relevel(group_data$Condition_Combo, ref = ref_condition)

  # ---------------------------------------------------------
  # B. PASS 1: FIND THE OPTIMAL RHO
  # ---------------------------------------------------------
  cat("  -> Fitting initial model to estimate autocorrelation (rho)...\n")

  model_no_ar <- bam(
    formula = gamm_formula,
    data = group_data,
    family = chosen_family,
    select = TRUE,
    method = "fREML",
    discrete = TRUE,
    nthreads = usable_cores
  )

  # Extract residuals and calculate the lag-1 correlation securely within each animal
  group_data$res <- resid(model_no_ar)

  optimal_rho <- group_data %>%
    group_by(animal_id) %>%
    summarise(
      animal_rho = cor(res, lag(res), use = "pairwise.complete.obs")
    ) %>%
    pull(animal_rho) %>%
    mean(na.rm = TRUE)

  cat("  -> Calculated Optimal Rho:", round(optimal_rho, 3), "\n")

  # ---------------------------------------------------------
  # C. PASS 2: FIT THE FINAL MODEL
  # ---------------------------------------------------------
  cat("  -> Fitting final model with AR1 correction...\n")

  final_model <- bam(
    formula = gamm_formula,
    data = group_data,
    family = chosen_family,
    rho = optimal_rho,
    AR.start = group_data$start_event,
    select = TRUE,
    method = "fREML",
    discrete = TRUE,
    nthreads = usable_cores
  )

  # Store the final model in the list
  models_by_group[[as.character(g)]] <- final_model

  cat("Successfully fitted Group", g, "!\n")
}

# ---------------------------------------------------------
# 3. GLOBAL MULTIPLE COMPARISONS: ACROSS ALL PHASE GROUPS
# ---------------------------------------------------------
# Dynamically generate per-variable contrasts and the full interaction contrast

all_contrasts_list <- list()

for (g in as.character(sort(as.integer(names(models_by_group))))) {

  target_model <- models_by_group[[g]]

  # ── Per-variable effects (each variable vs reference, split by other variables) ──
  for (vi in seq_along(var_names)) {
    focal_var   <- var_names[vi]
    other_vars  <- var_names[-vi]
    ref_val     <- ref_map[[focal_var]]

    # Build emmeans spec: ~ FocalVar | OtherVar1 + OtherVar2 + ...
    if (length(other_vars) > 0) {
      spec_str <- paste0("~ ", focal_var, " | ", paste(other_vars, collapse = " + "))
    } else {
      spec_str <- paste0("~ ", focal_var)
    }
    spec_formula <- as.formula(spec_str)

    emm_base <- emmeans(target_model, specs = spec_formula)
    emm_contrasts <- contrast(emm_base, method = "trt.vs.ctrl", ref = ref_val, adjust = contrast_adjust)

    df_contrasts <- as.data.frame(emm_contrasts) %>%
      mutate(
        Group = g,
        Test_Family = paste0(focal_var, "_Effect"),
        Tested_Level = as.character(contrast)
      )

    # Rename the "other vars" columns to a generic Split_By for stacking
    if (length(other_vars) > 0) {
      # Create a combined Split_By column from other variable columns
      split_cols <- intersect(other_vars, colnames(df_contrasts))
      if (length(split_cols) > 0) {
        df_contrasts$Split_By <- apply(df_contrasts[, split_cols, drop = FALSE], 1,
                                        function(x) paste(x, collapse = " + "))
        df_contrasts <- df_contrasts %>% select(-all_of(split_cols))
      } else {
        df_contrasts$Split_By <- "All"
      }
    } else {
      df_contrasts$Split_By <- "All"
    }

    all_contrasts_list[[paste0(g, "_", focal_var)]] <- df_contrasts
  }

  # ── Full interaction: all pairwise condition comparisons ──
  # Build the full interaction spec: ~ Var1 * Var2 * ... * VarN
  inter_spec_str <- paste0("~ ", paste(var_names, collapse = " * "))
  inter_spec <- as.formula(inter_spec_str)

  emm_inter_base <- emmeans(target_model, specs = inter_spec)
  emm_inter_contrasts <- contrast(emm_inter_base, method = "pairwise")

  df_inter <- as.data.frame(emm_inter_contrasts) %>%
    mutate(
      Group = g,
      Test_Family = "Full_Interaction",
      Tested_Level = as.character(contrast),
      Split_By = "None"
    )

  all_contrasts_list[[paste0(g, "_interaction")]] <- df_inter
}

# ---------------------------------------------------------
# COMBINE AND APPLY GLOBAL PENALTY
# ---------------------------------------------------------
master_results_df <- bind_rows(all_contrasts_list)

final_master_table <- master_results_df %>%
  group_by(Test_Family) %>%
  mutate(
    Global_FDR_pvalue = p.adjust(p.value, method = global_correction),
    Group = as.integer(Group)   # convert from character so sorting is numerical
  ) %>%
  ungroup() %>%
  select(Test_Family, Group, Split_By, Tested_Level, estimate, SE, p.value, Global_FDR_pvalue) %>%
  arrange(Test_Family, Group)

cat("\n--- CLEAN MASTER TABLE ---\n")
print(final_master_table, n=nrow(final_master_table))

# =============================================================================
# Visualizations
# =============================================================================

# Labels that reflect the chosen global correction method
correction_label     <- if (global_correction == "BH") "FDR-corrected p-value" else "FWER-corrected p-value"
correction_sig_label <- if (global_correction == "BH") "FDR < 0.05" else "FWER < 0.05"

# ---------------------------------------------------------
# FOREST PLOTS: One per variable effect + one for interaction
# ---------------------------------------------------------

# Helper: annotate condition names within a comparison label (e.g. "WT DMSO - KO DMSO")
# Condition names in emmeans output use spaces between variable levels (not "+")
# so we need to match against the "+" form from our role_map
annotate_comparison <- function(label) {
  if (length(role_map) == 0) return(label)
  # emmeans pairwise labels separate conditions with " - "
  parts <- str_split(label, " - ", n = 2)[[1]]
  annotated <- sapply(parts, function(p) {
    p <- trimws(p)
    # Try matching: emmeans uses space-separated, role_map uses "+" separated
    cond_plus <- gsub(" ", "+", p)
    role <- role_map[cond_plus]
    if (!is.na(role)) paste0(p, " [", role, "]") else p
  }, USE.NAMES = FALSE)
  paste(annotated, collapse = " - ")
}

forest_plot_data <- final_master_table %>%
  mutate(
    CI_lower = estimate - 1.96 * SE,
    CI_upper = estimate + 1.96 * SE,
    sig_star = case_when(
      Global_FDR_pvalue < 0.001 ~ "***",
      Global_FDR_pvalue < 0.01 ~ "**",
      Global_FDR_pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Tested_Level_Annotated = sapply(Tested_Level, annotate_comparison),
    Tested_Level_Clean = str_replace_all(Tested_Level_Annotated, " - ", "\nvs\n"),
    Group_Label = paste("Group", Group),
    sig_category = case_when(
      Global_FDR_pvalue < 0.05 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

# Per-variable forest plots
for (vi in seq_along(var_names)) {
  focal_var    <- var_names[vi]
  display_var  <- var_display_map[[focal_var]]
  ref_val      <- ref_map[[focal_var]]
  family_name  <- paste0(focal_var, "_Effect")

  plot_data <- forest_plot_data %>% filter(Test_Family == family_name)
  if (nrow(plot_data) == 0) next

  p <- plot_data %>%
    ggplot(aes(x = estimate, y = interaction(Group, Split_By), color = sig_category)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper), width = 0.3, linewidth = 0.8, orientation = "y") +
    geom_point(size = 3) +
    geom_text(aes(label = sig_star, x = CI_upper), hjust = -0.5, size = 4, show.legend = FALSE) +
    scale_color_manual(
      values = c("Significant" = "#D55E00", "Not significant" = "gray60"),
      name = correction_label
    ) +
    facet_grid(Split_By ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = paste0(display_var, " Effect (vs. ", ref_val, ")"),
      x = "Effect Size (log scale pixel difference)",
      y = NULL,
      caption = "Significance: * p < 0.05   ** p < 0.01   *** p < 0.001"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.y = element_text(angle = 0, face = "bold"),
      axis.text.y = element_text(size = 10),
      plot.caption = element_text(hjust = 1, face = "italic", size = 10)
    )

  fname <- paste0("forest_plot_", tolower(focal_var), "_effect.png")
  ggsave(fname, p, width = 10, height = 8, dpi = 300)
  cat("Saved:", fname, "\n")
}

# Interaction forest plot
inter_plot_data <- forest_plot_data %>% filter(Test_Family == "Full_Interaction")
if (nrow(inter_plot_data) > 0) {
  p3_interaction <- inter_plot_data %>%
    ggplot(aes(x = estimate, y = factor(Group), color = sig_category)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper), width = 0.3, linewidth = 0.8, orientation = "y") +
    geom_point(size = 3) +
    geom_text(aes(label = sig_star, x = CI_upper), hjust = -0.5, size = 3.5, show.legend = FALSE) +
    scale_color_manual(
      values = c("Significant" = "#D55E00", "Not significant" = "gray60"),
      name = correction_label
    ) +
    facet_wrap(~Tested_Level_Annotated, ncol = 3, scales = "free_x") +
    labs(
      title = paste0("Full Interaction (", paste(var_names_display, collapse = " x "), ")"),
      x = "Effect Size (log scale pixel difference)",
      y = "Group",
      caption = "Significance: * p < 0.05   ** p < 0.01   *** p < 0.001"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 9, face = "bold"),
      plot.caption = element_text(hjust = 1, face = "italic", size = 9)
    )

  ggsave("forest_plot_interactions.png", p3_interaction, width = 14, height = 6, dpi = 300)
  cat("Saved: forest_plot_interactions.png\n")
}

# ---------------------------------------------------------
# HEATMAP OF SIGNIFICANCE ACROSS GROUPS
# ---------------------------------------------------------

heatmap_data_diverging <- final_master_table %>%
  mutate(
    Tested_Level_Annotated = sapply(Tested_Level, annotate_comparison),
    signed_sig = -log10(Global_FDR_pvalue) * sign(estimate),
    signed_sig_capped = pmax(pmin(signed_sig, 5), -5),
    row_label = case_when(
      Split_By != "None" ~ paste(Split_By, ":", str_trunc(Tested_Level_Annotated, width = 50)),
      TRUE ~ str_trunc(Tested_Level_Annotated, width = 55)
    )
  )

p_heatmap_diverging_v1 <- heatmap_data_diverging %>%
  ggplot(aes(x = factor(Group, levels = sort(unique(Group))), y = row_label, fill = signed_sig_capped)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(abs(signed_sig) > 1.3,  # p < 0.05
                               sprintf("%.2f", estimate),
                               "")),
            size = 2.3, color = "gray20", fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC",      # Blue = negative effect
    mid = "gray95",       # Gray = not significant
    high = "#B2182B",     # Red = positive effect
    midpoint = 0,
    name = expression("Signed -log"[10]*"(p)"),
    breaks = c(-5, -2, 0, 2, 5),
    labels = c("Decrease\np<0.00001", "p<0.01", "ns", "p<0.01", "Increase\np<0.00001"),
    limits = c(-5, 5)
  ) +
  facet_grid(Test_Family ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Effect Direction and Significance Across Groups",
    subtitle = "Blue = decreased activity | Red = increased activity | Color intensity = significance",
    x = "Group (Experimental Phase)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, face = "bold", size = 10),
    legend.position = "right",
    panel.grid = element_blank(),
    legend.key.height = unit(1.5, "cm")
  )

ggsave("heatmap_diverging_v1.png", p_heatmap_diverging_v1, width = 12, height = 14, dpi = 300)

# ---------------------------------------------------------
# RESCUE LINE PLOTS (only for 2-variable designs with rescue roles)
# ---------------------------------------------------------
# This plot is most meaningful for the classic Genotype+Drug rescue paradigm.
# It is generated only when exactly 2 variables are present and the reference
# condition can form meaningful rescue comparisons.

if (n_vars == 2) {
  cat("Generating rescue assessment plot (2-variable design)...\n")

  rescue_context_data <- final_master_table %>%
    filter(Test_Family == "Full_Interaction") %>%
    mutate(
      Tested_Level_Annotated = sapply(Tested_Level, annotate_comparison),
      CI_lower = estimate - 1.96 * SE,
      CI_upper = estimate + 1.96 * SE,
      Group = as.numeric(Group),
      is_sig = Global_FDR_pvalue < 0.05
    )

  if (nrow(rescue_context_data) > 0) {
    p_rescue_context_faceted <- rescue_context_data %>%
      ggplot(aes(x = Group, y = estimate)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, group = Tested_Level),
                  alpha = 0.2, fill = "gray70") +
      geom_line(aes(group = Tested_Level), linewidth = 1.2, color = "gray30") +
      geom_point(aes(color = is_sig), size = 3.5) +
      scale_color_manual(
        values = c("TRUE" = "#D55E00", "FALSE" = "gray60"),
        labels = c("ns", correction_sig_label),
        name = NULL
      ) +
      scale_x_continuous(breaks = sort(unique(rescue_context_data$Group))) +
      facet_wrap(~Tested_Level_Annotated, ncol = 2, scales = "free_y") +
      labs(
        title = "Effect Size Rescue Assessment",
        subtitle = "Each panel shows a different comparison. Shaded areas = 95% CI",
        x = "Group (Experimental Phase)",
        y = "Effect Size (+/- 95% CI)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        panel.spacing = unit(1, "lines")
      )

    ggsave("rescue_assessment_context_faceted.png", p_rescue_context_faceted, width = 12, height = 8, dpi = 300)
    cat("Saved: rescue_assessment_context_faceted.png\n")
  }
} else {
  cat("Skipping rescue assessment plot (designed for 2-variable experiments).\n")
}

# ---------------------------------------------------------
# BONUS: CREATE A SUMMARY TABLE FOR QUICK REFERENCE
# ---------------------------------------------------------

summary_table <- final_master_table %>%
  group_by(Test_Family, Group) %>%
  summarise(
    n_total = n(),
    n_sig = sum(Global_FDR_pvalue < 0.05),
    mean_effect = mean(abs(estimate)),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(n_sig, mean_effect),
    names_glue = "G{Group}_{.value}"
  )

write_csv(summary_table, "summary_statistics_by_group.csv")
cat("\nAnalysis complete. All outputs saved to:", output_dir, "\n")
