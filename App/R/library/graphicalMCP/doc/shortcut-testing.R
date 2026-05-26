## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  fig.align = "center",
  collapse = TRUE,
  comment = "#>",
  cache.lazy = FALSE
)

## ----setup, message = FALSE, warning = FALSE----------------------------------
library(gt)
library(graphicalMCP)

## ----create-graph, fig.dim=c(3, 3)--------------------------------------------
hypotheses <- c(0.5, 0.5, 0, 0)
transitions <- rbind(
  c(0, 0.5, 0.5, 0),
  c(0.5, 0, 0, 0.5),
  c(0, 1, 0, 0),
  c(1, 0, 0, 0)
)
hyp_names <- c("H1", "H2", "H3", "H4")

g <- graph_create(hypotheses, transitions, hyp_names)

plot(g, vertex.size = 60)

## ----graph-test---------------------------------------------------------------
p_values <- c(0.018, 0.01, 0.105, 0.006)
test_results <- graph_test_shortcut(g, p = p_values, alpha = 0.025)

test_results$outputs$adjusted_p # Adjusted p-values
test_results$outputs$rejected # Rejections

## ----graph-test-verbose-------------------------------------------------------
test_results$outputs$graph # Final graph after H1, H2 and H4 rejected (as NA's)
test_results$outputs$graph$hypotheses # Hypothesis weights of the final graph
test_results$outputs$graph$transitions # Transition weights of the final graph

test_results_verbose <- graph_test_shortcut(
  g,
  p = p_values,
  alpha = 0.025,
  verbose = TRUE
)

# Intermediate graph after H1 and H2 rejected
test_results_verbose$details$results[[3]]

## ----graph-test-order---------------------------------------------------------
# Obtain all valid orders of rejections
orders <- graph_rejection_orderings(test_results)$valid_orderings
orders

# Intermediate graphs following the order of H2 and H4
graph_update(g, delete = orders[[2]])$intermediate_graphs[[3]]

## ----graph-test-test_values---------------------------------------------------
test_results_test_values <- graph_test_shortcut(
  g,
  p = p_values,
  alpha = 0.025,
  test_values = TRUE
)

test_results_test_values$test_values$results

## ----power-calculate-primary--------------------------------------------------
alpha <- 0.025
prop <- c(0.3, 0.181, 0.181)
sample_size <- rep(200, 3)

unpooled_variance <-
  prop[-1] * (1 - prop[-1]) / sample_size[-1] +
  prop[1] * (1 - prop[1]) / sample_size[1]

noncentrality_parameter_primary <-
  -(prop[-1] - prop[1]) / sqrt(unpooled_variance)

power_marginal_primary <- pnorm(
  qnorm(alpha, lower.tail = FALSE),
  mean = noncentrality_parameter_primary,
  sd = 1,
  lower.tail = FALSE
)

names(power_marginal_primary) <- c("H1", "H2")
power_marginal_primary

## ----power-calculate-secondary------------------------------------------------
mean_change <- c(5, 7.5, 8.25)
sd <- rep(10, 3)
variance <- sd[-1]^2 / sample_size[-1] + sd[1]^2 / sample_size[1]

noncentrality_parameter_secondary <-
  (mean_change[-1] - mean_change[1]) / sqrt(variance)

power_marginal_secondary <- pnorm(
  qnorm(alpha, lower.tail = FALSE),
  mean = noncentrality_parameter_secondary,
  sd = 1,
  lower.tail = FALSE
)

names(power_marginal_secondary) <- c("H3", "H4")
power_marginal_secondary

## ----power-calculate-correlation----------------------------------------------
corr <- matrix(0, nrow = 4, ncol = 4)

corr[1, 2] <-
  corr[3, 4] <-
  sqrt(
    sample_size[2] / sum(sample_size[1:2]) *
      sample_size[3] / sum(sample_size[c(1, 3)])
  )

rho <- 0.5
corr[1, 3] <- corr[2, 4] <- rho
corr[1, 4] <- corr[2, 3] <- corr[1, 2] * rho
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
  `Expected no. of rejections` = function(x) x[1] + x[2] + x[3] + x[4],

  # Probability to reject at least one hypothesis
  `AtLeast1` = function(x) x[1] | x[2] | x[3] | x[4],

  # Probability to reject all hypotheses
  `All` = function(x) x[1] & x[2] & x[3] & x[4],

  # Probability to reject both H1 and H2
  `H1andH2` = function(x) x[1] & x[2],

  # Probability to reject both hypotheses for the low dose or the high dose
  `(H1andH3)or(H2andH4)` = function(x) (x[1] & x[3]) | (x[2] & x[4])
)

## ----power-calculate-calculate------------------------------------------------
set.seed(1234)
power_output <- graph_calculate_power(
  g,
  alpha = 0.025,
  sim_corr = corr,
  sim_n = 1e5,
  power_marginal = c(power_marginal_primary, power_marginal_secondary),
  sim_success = success_fns
)

power_output$power

## ----power-calculate-calculate_verbose----------------------------------------
set.seed(1234)
power_verbose_output <- graph_calculate_power(
  g,
  alpha = 0.025,
  sim_corr = corr,
  sim_n = 1e5,
  power_marginal = c(power_marginal_primary, power_marginal_secondary),
  sim_success = success_fns,
  verbose = TRUE
)

head(power_verbose_output$details$p_sim, 10)

print(power_verbose_output, indent = 4, precision = 6)

## ----power-gMCP, eval=FALSE, echo=FALSE---------------------------------------
#  ## Compared with gMCP
#  
#  library(gMCP)
#  
#  g_gMCP <- as_graphMCP(g)
#  
#  set.seed(1234)
#  result <- calcPower(
#    graph = g_gMCP, mean = c(
#      noncentrality_parameter_primary,
#      noncentrality_parameter_secondary
#    ), f = success_fns, corr.sim = corr, alpha = 0.025,
#    n.sim = 1e+05
#  )
#  
#  # Local power
#  output$power$power_local
#  result$LocalPower
#  # User-defined power
#  as.numeric(output$power$power_success)
#  c(result$func1, result$func2, result$func3, result$func4, result$func5, result$func6)

