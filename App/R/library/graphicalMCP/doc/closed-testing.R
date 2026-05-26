## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  fig.align = "center",
  collapse = TRUE,
  comment = "#>",
  cache.lazy = FALSE
)

## ----setup, message = FALSE, warning = FALSE----------------------------------
library(gt)
library(gMCP)
library(lrstat)
library(graphicalMCP)

## ----create-graph, fig.dim=c(5, 3.4)------------------------------------------
hypotheses <- c(0.5, 0.5, 0, 0, 0, 0)

epsilon <- 1e-5
transitions <- rbind(
  c(0, 0.5, 0.25, 0, 0.25, 0),
  c(0.5, 0, 0, 0.25, 0, 0.25),
  c(0, 0, 0, 0, 1, 0),
  c(epsilon, 0, 0, 0, 0, 1 - epsilon),
  c(0, epsilon, 1 - epsilon, 0, 0, 0),
  c(0, 0, 0, 1, 0, 0)
)

hyp_names <- c("H1", "H2", "H3", "H4", "H5", "H6")
g <- graph_create(hypotheses, transitions, hyp_names)

plot_layout <- rbind(
  c(0.15, 0.5),
  c(0.65, 0.5),
  c(0, 0),
  c(0.5, 0),
  c(0.3, 0),
  c(0.8, 0)
)

plot(
  g,
  layout = plot_layout,
  eps = epsilon,
  edge_curves = c(pairs = 0.8),
  vertex.size = 35
)

## ----graph-test-bonferroni----------------------------------------------------
p_values <- c(0.015, 0.013, 0.01, 0.007, 0.1, 0.0124)
test_results <- graph_test_closure(
  g,
  p = p_values,
  alpha = 0.025,
  verbose = TRUE,
  test_values = TRUE
)

test_results$outputs$adjusted_p

test_results$outputs$rejected
# Same testing results as
# 'graph_test_shortcut(g, p = p_values, alpha = 0.025)$outputs$rejected'

## ----graph-test-verbose-------------------------------------------------------
test_results_verbose <-
  graph_test_closure(
    g,
    p = p_values,
    alpha = 0.025,
    verbose = TRUE
  )

head(test_results_verbose$details$results)

## ----graph-test-test_values---------------------------------------------------
test_results_test_values <-
  graph_test_closure(
    g,
    test_types = "b",
    p = p_values,
    alpha = 0.025,
    test_values = TRUE
  )

head(test_results_test_values$test_values$results)

## ----graph-test-parametric----------------------------------------------------
corr_12 <- matrix(0.5, nrow = 2, ncol = 2)
diag(corr_12) <- 1

test_results_parametric <-
  graph_test_closure(
    g,
    p = p_values,
    alpha = 0.025,
    test_groups = list(c(1, 2), 3:6),
    test_types = c("parametric", "bonferroni"),
    test_corr = list(corr_12, NA),
    test_values = TRUE
  )

test_results_parametric$outputs$adjusted_p
test_results_parametric$outputs$rejected

head(test_results_parametric$test_values$results)

## ----graph-test-parametric-simes----------------------------------------------
test_results_parametric_simes <-
  graph_test_closure(
    g,
    p = p_values, alpha = 0.025,
    test_groups = list(c(1, 2), c(3, 5), c(4, 6)),
    test_types = c("parametric", "simes", "simes"),
    test_corr = list(corr_12, NA, NA),
    test_values = TRUE
  )

test_results_parametric_simes$outputs$adjusted_p
test_results_parametric_simes$outputs$rejected

head(test_results_parametric_simes$test_values$results)

## ----power-calculate-primary--------------------------------------------------
alpha <- 0.025
prop <- c(0.3, 0.181, 0.181)
sample_size <- rep(200, 3)

unpooled_variance <-
  prop[-1] * (1 - prop[-1]) / sample_size[-1] +
  prop[1] * (1 - prop[1]) / sample_size[1]

noncentrality_parameter_primary <-
  -(prop[-1] - prop[1]) / sqrt(unpooled_variance)

marginal_power_primary <- pnorm(
  qnorm(alpha, lower.tail = FALSE),
  mean = noncentrality_parameter_primary,
  sd = 1,
  lower.tail = FALSE
)

names(marginal_power_primary) <- c("H1", "H2")
marginal_power_primary

## ----power-calculate-secondary------------------------------------------------
mean_change_se1 <- c(5, 7.5, 8.25)
sd <- rep(10, 3)
variance <- sd[-1]^2 / sample_size[-1] + sd[1]^2 / sample_size[1]

noncentrality_parameter_se1 <-
  (mean_change_se1[-1] - mean_change_se1[1]) /
    sqrt(variance)

marginal_power_se1 <- pnorm(
  qnorm(alpha, lower.tail = FALSE),
  mean = noncentrality_parameter_se1,
  sd = 1,
  lower.tail = FALSE
)

names(marginal_power_se1) <- c("H3", "H4")
marginal_power_se1

mean_change_se2 <- c(6, 8, 9)

noncentrality_parameter_se2 <-
  (mean_change_se2[-1] - mean_change_se2[1]) /
    sqrt(variance)

marginal_power_se2 <- pnorm(
  qnorm(alpha, lower.tail = FALSE),
  mean = noncentrality_parameter_se2,
  sd = 1,
  lower.tail = FALSE
)

names(marginal_power_se2) <- c("H5", "H6")
marginal_power_se2

## ----power-calculate-correlation----------------------------------------------
corr <- matrix(0, nrow = 6, ncol = 6)

corr[1, 2] <-
  corr[3, 4] <-
  corr[5, 6] <-
  sqrt(
    sample_size[2] / sum(sample_size[1:2]) *
      sample_size[3] / sum(sample_size[c(1, 3)])
  )

rho <- 0.5

corr[1, 3] <-
  corr[1, 5] <-
  corr[2, 4] <-
  corr[2, 6] <-
  corr[3, 5] <-
  corr[4, 6] <-
  rho

corr[1, 4] <- corr[1, 6] <- corr[2, 3] <- corr[2, 5] <- corr[1, 2] * rho

corr[3, 6] <- corr[1, 3] * corr[1, 6]
corr[4, 5] <- corr[1, 4] * corr[1, 6]

corr <- corr + t(corr)
diag(corr) <- 1
colnames(corr) <- hyp_names
rownames(corr) <- hyp_names
corr

## ----power-calculate-success--------------------------------------------------
success_fns <- list(
  # Probability to reject H1
  H1 = function(x) x[1],

  # Expected number of rejections
  `Expected no. of rejections` =
    function(x) x[1] + x[2] + x[3] + x[4] + x[5] + x[6],

  # Probability to reject at least one hypothesis
  `AtLeast1` = function(x) x[1] | x[2] | x[3] | x[4] | x[5] | x[6],

  # Probability to reject all hypotheses
  `All` = function(x) x[1] & x[2] & x[3] & x[4] & x[5] & x[6],

  # Probability to reject both H1 and H2
  `H1andH2` = function(x) x[1] & x[2],

  # Probability to reject all hypotheses for the low dose or the high dose
  `(H1andH3andH5)or(H2andH4andH6)` <-
    function(x) (x[1] & x[3] & x[5]) | (x[2] & x[4] & x[6])
)

## ----power-calculate-calculate------------------------------------------------
set.seed(1234)
power_bonferroni <- graph_calculate_power(
  g,
  alpha = 0.025,
  sim_corr = corr,
  sim_n = 1e5,
  power_marginal = c(
    marginal_power_primary,
    marginal_power_se1,
    marginal_power_se2
  ),
  sim_success = success_fns
)

round(power_bonferroni$power$power_local, 3)

set.seed(1234)
power_parametric <- graph_calculate_power(
  g,
  alpha = 0.025,
  sim_corr = corr,
  sim_n = 1e5,
  power_marginal = c(
    marginal_power_primary,
    marginal_power_se1,
    marginal_power_se2
  ),
  test_groups = list(c(1, 2), 3:6),
  test_types = c("parametric", "bonferroni"),
  test_corr = list(corr_12, NA),
  sim_success = success_fns
)

round(power_parametric$power$power_local, 3)

set.seed(1234)
power_parametric_simes <- graph_calculate_power(
  g,
  alpha = 0.025,
  sim_corr = corr,
  sim_n = 1e5,
  power_marginal = c(
    marginal_power_primary,
    marginal_power_se1,
    marginal_power_se2
  ),
  test_groups = list(c(1, 2), c(3, 5), c(4, 6)),
  test_types = c("parametric", "simes", "simes"),
  test_corr = list(corr_12, NA, NA),
  sim_success = success_fns
)

round(power_parametric_simes$power$power_local, 3)

## ----power-calculate-calculate_verbose----------------------------------------
set.seed(1234)
power_verbose_output_parametric_simes <- graph_calculate_power(
  g,
  alpha = 0.025,
  sim_corr = corr,
  sim_n = 1e5,
  power_marginal = c(
    marginal_power_primary,
    marginal_power_se1,
    marginal_power_se2
  ),
  test_groups = list(c(1, 2), c(3, 5), c(4, 6)),
  test_types = c("parametric", "simes", "simes"),
  test_corr = list(corr_12, NA, NA),
  sim_success = success_fns,
  verbose = TRUE
)

head(power_verbose_output_parametric_simes$details$p_sim, 10)

print(power_verbose_output_parametric_simes,
  indent = 4,
  precision = 6
)

