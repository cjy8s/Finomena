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
#  Rscript "TweedieAR1 GAMM.R" --args <input_csv> <output_dir>
args <- commandArgs(trailingOnly = TRUE)
input_csv  <- args[1]
output_dir <- if (length(args) >= 2 && !is.na(args[2])) args[2] else dirname(input_csv)

if (is.na(input_csv) || !file.exists(input_csv)) {
  stop("ERROR: input_csv argument missing or file not found: ", input_csv)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# ── ANALYSIS PARAMETERS (from Python app args[3]–args[7]) ─────────────────────
ref_genotype      <- if (length(args) >= 3 && !is.na(args[3]) && nchar(args[3]) > 0) args[3] else "WT"
ref_drug          <- if (length(args) >= 4 && !is.na(args[4]) && nchar(args[4]) > 0) args[4] else "DMSO"
ref_condition     <- if (length(args) >= 5 && !is.na(args[5]) && nchar(args[5]) > 0) args[5] else "WT+DMSO"
global_correction <- if (length(args) >= 6 && !is.na(args[6]) && nchar(args[6]) > 0) args[6] else "BH"
contrast_adjust   <- if (length(args) >= 7 && !is.na(args[7]) && nchar(args[7]) > 0) args[7] else "dunnett"

cat("Input CSV:        ", input_csv, "\n")
cat("Output dir:       ", output_dir, "\n")
cat("Ref Genotype:     ", ref_genotype, "\n")
cat("Ref Drug:         ", ref_drug, "\n")
cat("Ref Condition:    ", ref_condition, "\n")
cat("Global correction:", global_correction, "\n")
cat("Contrast adjust:  ", contrast_adjust, "\n")

# Find how many cores the computer has, and leave 1 free so the computer doesn't freeze
usable_cores <- max(1, detectCores() - 1)

# ── DATA LOADING (replaces all preprocessing and file reading) ─────────────────
# The Python app exports a CSV with columns:
#   time_sec, location, loc_coord, pixel_diff, Condition, Phase, Group, animal_id, plate
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
  filter(!is.na(Condition), Condition != "") %>%
  separate(Condition, into = c("Genotype", "Drug"), sep = "\\+", remove = FALSE)

# ---------------------------------------------------------
# 1. FINAL DATA PREP FOR MGCV
# ---------------------------------------------------------
gam_df <- full_df %>%
  # A. Convert categorical variables to R Factors (required for mgcv to understand them)
  mutate(
    # Create a truly unique animal ID (Plate 1 Location A1 is different from Plate 2 Location A1)
    animal_id = as.factor(paste(plate, location, sep = "_")),
    Genotype = as.factor(Genotype),
    Drug = as.factor(Drug),
    # We create a combined factor for the interaction smooths
    Geno_Drug = as.factor(Condition) 
  )

# ---------------------------------------------------------
# 2. SETUP THE MODELING LOOP
# ---------------------------------------------------------
# Create an empty list to store the fitted models so you can analyze them later
models_by_group <- list()

# Extract the unique group numbers (e.g., 1 through 8)
unique_groups <- sort(unique(gam_df$Group))

# Define k value before edf calculation
k_start = 30

# Define the GAMM formula
gamm_formula <- pxl_diff ~ 
  Genotype + Drug + Genotype:Drug +
  s(time_in_group, k = k_start) +
  s(time_in_group, by = Geno_Drug, k = k_start) +
  s(time_in_group, animal_id, bs = "fs", m = 1)

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
  
  # Ensure the variables are treated as factors
  group_data$Genotype <- as.factor(group_data$Genotype)
  group_data$Drug <- as.factor(group_data$Drug)
  
  group_data$Genotype <- relevel(group_data$Genotype, ref = ref_genotype)
  group_data$Drug <- relevel(group_data$Drug, ref = ref_drug)

  # Combined Genotype_Drug factor,
  group_data$Geno_Drug <- as.factor(group_data$Geno_Drug)
  group_data$Geno_Drug <- relevel(group_data$Geno_Drug, ref = ref_condition)
  
  # ---------------------------------------------------------
  # B. PASS 1: FIND THE OPTIMAL RHO
  # ---------------------------------------------------------
  cat("  -> Fitting initial model to estimate autocorrelation (rho)...\n")
  
  model_no_ar <- bam(
    formula = gamm_formula,
    data = group_data,
    family = tw(),
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
    family = tw(),
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
drug_contrasts_list <- list()
geno_contrasts_list <- list()
interaction_contrasts_list <- list()

for (g in as.character(sort(as.integer(names(models_by_group))))) {
  
  target_model <- models_by_group[[g]]
  
  # ---------------------------------------------------------
  # QUESTION 1: The Drug Reference (Drug vs. DMSO)
  # ---------------------------------------------------------
  # Step 1: Get the base means
  emm_drug_base <- emmeans(target_model, specs = ~ Drug | Genotype)
  # Step 2: Explicitly isolate ONLY the contrasts
  emm_drug_contrasts <- contrast(emm_drug_base, method = "trt.vs.ctrl", ref = ref_drug, adjust = contrast_adjust)
  
  df_drug <- as.data.frame(emm_drug_contrasts) %>%
    mutate(
      Group = g, 
      Test_Family = "Drug_Effect",
      Tested_Level = as.character(contrast) # Copy the contrast string cleanly
    ) %>%
    rename(Split_By = Genotype) # Rename the grouping column to stack perfectly
  
  drug_contrasts_list[[g]] <- df_drug
  
  # ---------------------------------------------------------
  # QUESTION 2: The Genotype Reference (KO vs. WT)
  # ---------------------------------------------------------
  emm_geno_base <- emmeans(target_model, specs = ~ Genotype | Drug)
  emm_geno_contrasts <- contrast(emm_geno_base, method = "trt.vs.ctrl", ref = ref_genotype, adjust = contrast_adjust)
  
  df_geno <- as.data.frame(emm_geno_contrasts) %>%
    mutate(
      Group = g, 
      Test_Family = "Genotype_Effect",
      Tested_Level = as.character(contrast)
    ) %>%
    rename(Split_By = Drug)
  
  geno_contrasts_list[[g]] <- df_geno
  
  # ---------------------------------------------------------
  # QUESTION 3: The Grand Interaction Reference (All Pairwise)
  # ---------------------------------------------------------
  emm_inter_base <- emmeans(target_model, specs = ~ Genotype * Drug)
  emm_inter_contrasts <- contrast(emm_inter_base, method = "pairwise")
  
  df_inter <- as.data.frame(emm_inter_contrasts) %>%
    # I REMOVED the filter() completely. You will now get all 6 combinations.
    mutate(
      Group = g, 
      Test_Family = "Geno_Drug_Interaction",
      Tested_Level = as.character(contrast),
      Split_By = "None"
    )
  
  interaction_contrasts_list[[g]] <- df_inter
}

# ---------------------------------------------------------
# COMBINE AND APPLY GLOBAL PENALTY
# ---------------------------------------------------------
master_results_df <- bind_rows(
  bind_rows(drug_contrasts_list),
  bind_rows(geno_contrasts_list),
  bind_rows(interaction_contrasts_list)
)

final_master_table <- master_results_df %>%
  group_by(Test_Family) %>%
  mutate(
    Global_FDR_pvalue = p.adjust(p.value, method = global_correction),
    Group = as.integer(Group)   # convert from character so sorting is numerical
  ) %>%
  ungroup() %>%
  # We can drop the redundant 'contrast' column here since Tested_Level has the same info
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
# VISUALIZATION 1: FOREST PLOT (Effect Sizes with Confidence Intervals)
# ---------------------------------------------------------
# This is THE most important visualization for your paper/presentation

forest_plot_data <- final_master_table %>%
  mutate(
    # Calculate confidence intervals
    CI_lower = estimate - 1.96 * SE,
    CI_upper = estimate + 1.96 * SE,
    # Significance markers
    sig_star = case_when(
      Global_FDR_pvalue < 0.001 ~ "***",
      Global_FDR_pvalue < 0.01 ~ "**",
      Global_FDR_pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # Create cleaner labels
    Tested_Level_Clean = str_replace_all(Tested_Level, " - ", "\nvs\n"),
    Group_Label = paste("Group", Group),
    # More informative significance category
    sig_category = case_when(
      Global_FDR_pvalue < 0.05 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

# Separate plot for each Test_Family
p1_drug <- forest_plot_data %>%
  filter(Test_Family == "Drug_Effect") %>%
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
    title = "Drug Effect: Compound - DMSO",
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

p2_geno <- forest_plot_data %>%
  filter(Test_Family == "Genotype_Effect") %>%
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
    title = "Genotype Effect: KO - WT",
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
    plot.caption = element_text(hjust = 1, face = "italic", size = 10)
  )

# For interaction effects, we'll make a more compact version
p3_interaction <- forest_plot_data %>%
  filter(Test_Family == "Geno_Drug_Interaction") %>%
  # Only show the key comparisons for clarity
  filter(str_detect(Tested_Level, "WT DMSO|KO DMSO|WT compound - KO compound")) %>%
  ggplot(aes(x = estimate, y = factor(Group), color = sig_category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper), width = 0.3, linewidth = 0.8, orientation = "y") +
  geom_point(size = 3) +
  geom_text(aes(label = sig_star, x = CI_upper), hjust = -0.5, size = 3.5, show.legend = FALSE) +
  scale_color_manual(
    values = c("Significant" = "#D55E00", "Not significant" = "gray60"),
    name = correction_label
  ) +
  facet_wrap(~Tested_Level, ncol = 3, scales = "free_x") +
  labs(
    title = "Key Genotype × Drug Interactions",
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

# Save individual plots
ggsave("forest_plot_drug_effect.pdf", p1_drug, width = 10, height = 8)
ggsave("forest_plot_genotype_effect.pdf", p2_geno, width = 10, height = 8)
ggsave("forest_plot_interactions.pdf", p3_interaction, width = 14, height = 6)

# ---------------------------------------------------------
# HEATMAP OF SIGNIFICANCE ACROSS GROUPS
# ---------------------------------------------------------
# Shows pattern of significant effects across your 8 experimental phases

heatmap_data_diverging <- final_master_table %>%
  mutate(
    # Signed significance: negative log10(p) with sign of effect
    signed_sig = -log10(Global_FDR_pvalue) * sign(estimate),
    # Cap for visualization
    signed_sig_capped = pmax(pmin(signed_sig, 5), -5),
    # Row labels
    row_label = case_when(
      Test_Family == "Drug_Effect" ~ paste(Split_By, ":", str_remove(Tested_Level, "compound - ")),
      Test_Family == "Genotype_Effect" ~ paste(Split_By, ":", str_remove(Tested_Level, "KO - ")),
      TRUE ~ str_trunc(Tested_Level, width = 40)
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
    name = expression("Signed -log"[10]*"(p)"),  # Proper subscript using expression()
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

ggsave("heatmap_diverging_v1.pdf", p_heatmap_diverging_v1, width = 12, height = 14)

# ---------------------------------------------------------
# RESCUE LINE PLOTS
# ---------------------------------------------------------

rescue_context_data <- final_master_table %>%
  filter(
    Test_Family == "Geno_Drug_Interaction",
    str_detect(Tested_Level, "WT DMSO") | str_detect(Tested_Level, "KO compound")
  ) %>%
  mutate(
    CI_lower = estimate - 1.96 * SE,
    CI_upper = estimate + 1.96 * SE,
    Group = as.numeric(Group),
    is_sig = Global_FDR_pvalue < 0.05,
    # Create cleaner labels
    Comparison = case_when(
      str_detect(Tested_Level, "WT DMSO - KO compound") ~ "WT DMSO - KO compound\n(Rescue target)",
      str_detect(Tested_Level, "WT DMSO - KO DMSO") ~ "WT DMSO - KO DMSO\n(Disease phenotype)",
      str_detect(Tested_Level, "WT DMSO - WT compound") ~ "WT DMSO - WT compound\n(Drug effect in WT)",
      str_detect(Tested_Level, "KO DMSO - KO compound") ~ "KO DMSO - KO compound\n(Drug effect in KO)",
      TRUE ~ Tested_Level
    )
  )
p_rescue_context_faceted <- rescue_context_data %>%
  ggplot(aes(x = Group, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Ribbon for CI - keep it gray/neutral since line color shows significance
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, group = Comparison), 
              alpha = 0.2, fill = "gray70") +
  # Line WITHOUT color mapping - ensures continuity
  geom_line(aes(group = Comparison), linewidth = 1.2, color = "gray30") +
  # Points WITH color mapping - shows which groups are significant
  geom_point(aes(color = is_sig), size = 3.5) +
  scale_color_manual(
    values = c("TRUE" = "#D55E00", "FALSE" = "gray60"),
    labels = c("ns", correction_sig_label),
    name = NULL
  ) +
  scale_x_continuous(breaks = sort(unique(rescue_context_data$Group))) +
  facet_wrap(~Comparison, ncol = 2, scales = "free_y") +
  labs(
    title = "Effect Size Rescue Assessment",
    subtitle = "Each panel shows a different comparison. Shaded areas = 95% CI",
    x = "Group (Experimental Phase)",
    y = "Effect Size (± 95% CI)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(1, "lines")
  )

ggsave("rescue_assessment_context_faceted.pdf", p_rescue_context_faceted, width = 12, height = 8)

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
