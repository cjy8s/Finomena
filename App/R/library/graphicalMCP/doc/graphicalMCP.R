## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)

## ----create-graph, fig.dim=c(3, 3)--------------------------------------------
library(graphicalMCP)
# A graph of two primary hypotheses (H1 and H2) and two secondary hypotheses (H3
# and H4)
hypotheses <- c(0.5, 0.5, 0, 0)
transitions <- rbind(
  c(0, 0, 1, 0),
  c(0, 0, 0, 1),
  c(0, 1, 0, 0),
  c(1, 0, 0, 0)
)
hyp_names <- c("H1", "H2", "H3", "H4")
example_graph <- graph_create(hypotheses, transitions, hyp_names)
example_graph

## ----plot-graph, fig.dim=c(3, 3)----------------------------------------------
plot(example_graph, vertex.size = 60)

## ----update-graph-------------------------------------------------------------
updated_example <- graph_update(
  example_graph,
  delete = c(TRUE, TRUE, FALSE, TRUE)
)

updated_example

## ----plot-updated-graph, fig.dim=c(3, 3)--------------------------------------
plot(updated_example, vertex.size = 60)

## ----test-graph-shortcut------------------------------------------------------
test_results <- graph_test_shortcut(
  example_graph,
  p = c(0.01, 0.005, 0.03, 0.01),
  alpha = 0.025
)
test_results$outputs$rejected

## ----test-graph---------------------------------------------------------------
test_results_closed <- graph_test_closure(
  example_graph,
  p = c(0.01, 0.005, 0.03, 0.01),
  alpha = 0.025,
  test_types = "bonferroni",
  test_groups = list(1:4)
)
test_results_closed$outputs$rejected

## ----power--------------------------------------------------------------------
set.seed(1234)
power_results <- graph_calculate_power(
  example_graph,
  sim_n = 1e6,
  power_marginal = c(0.9, 0.9, 0.8, 0.8)
)
power_results$power$power_local

