## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  fig.align = "center",
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(3, 3)
)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(graphicalMCP)

## ----bonferroni---------------------------------------------------------------
set.seed(1234)
alpha <- 0.025
m <- 3
bonferroni_graph <- bonferroni(m)
# transitions <- matrix(0, m, m)
# bonferroni_graph <- graph_create(rep(1 / m, m), transitions)

plot(
  bonferroni_graph,
  layout = igraph::layout_in_circle(
    as_igraph(bonferroni_graph),
    order = c(2, 1, 3)
  ),
  vertex.size = 70
)

p_values <- runif(m, 0, alpha)

test_results <-
  graph_test_shortcut(
    bonferroni_graph,
    p = p_values,
    alpha = alpha
  )

test_results$outputs$rejected

## ----weighted-bonferroni------------------------------------------------------
set.seed(1234)
alpha <- 0.025
hypotheses <- c(0.5, 0.3, 0.2)
weighted_bonferroni_graph <- bonferroni_weighted(hypotheses)
# m <- length(hypotheses)
# transitions <- matrix(0, m, m)
# weighted_bonferroni_graph <- graph_create(hypotheses, transitions)

plot(
  weighted_bonferroni_graph,
  layout = igraph::layout_in_circle(
    as_igraph(bonferroni_graph),
    order = c(2, 1, 3)
  ),
  vertex.size = 70
)

p_values <- runif(m, 0, alpha)

test_results <-
  graph_test_shortcut(
    weighted_bonferroni_graph,
    p = p_values,
    alpha = alpha
  )

test_results$outputs$rejected

## ----bonferroni-holm----------------------------------------------------------
set.seed(1234)
alpha <- 0.025
m <- 3
holm_graph <- bonferroni_holm(m)
# transitions <- matrix(1 / (m - 1), m, m)
# diag(transitions) <- 0
# holm_graph <- graph_create(rep(1 / m, m), transitions)

plot(holm_graph,
  layout = igraph::layout_in_circle(
    as_igraph(holm_graph),
    order = c(2, 1, 3)
  ),
  vertex.size = 70
)

p_values <- runif(m, 0, alpha)

test_results <- graph_test_shortcut(holm_graph, p = p_values, alpha = alpha)

test_results$outputs$rejected

## ----weighted-bonferroni-holm-------------------------------------------------
set.seed(1234)
alpha <- 0.025
hypotheses <- c(0.5, 0.3, 0.2)
weighted_holm_graph <- bonferroni_holm_weighted(hypotheses)
# m <- length(hypotheses)
# transitions <- matrix(1 / (m - 1), m, m)
# diag(transitions) <- 0
# weighted_holm_graph <- graph_create(hypotheses, transitions)

plot(weighted_holm_graph,
  layout = igraph::layout_in_circle(
    as_igraph(holm_graph),
    order = c(2, 1, 3)
  ),
  vertex.size = 70
)

p_values <- runif(m, 0, alpha)

test_results <- graph_test_shortcut(weighted_holm_graph, p = p_values, alpha = alpha)

test_results$outputs$rejected

## ----fixed-sequence-----------------------------------------------------------
set.seed(1234)
alpha <- 0.025
m <- 3
fixed_sequence_graph <- fixed_sequence(m)
# transitions <- rbind(
#   c(0, 1, 0),
#   c(0, 0, 1),
#   c(0, 0, 0)
# )
# fixed_sequence_graph <- graph_create(c(1, 0, 0), transitions)

plot(fixed_sequence_graph, nrow = 1, asp = 0.5, vertex.size = 50)

p_values <- runif(m, 0, alpha)

test_results <-
  graph_test_shortcut(
    fixed_sequence_graph,
    p = p_values,
    alpha = alpha
  )

test_results$outputs$rejected

## ----fallback-----------------------------------------------------------------
set.seed(1234)
alpha <- 0.025
hypotheses <- c(0.5, 0.3, 0.2)
m <- length(hypotheses)
fallback_graph <- fallback(hypotheses)
# transitions <- rbind(
#   c(0, 1, 0),
#   c(0, 0, 1),
#   c(0, 0, 0)
# )
# fallback_graph <- graph_create(hypotheses, transitions)

plot(fallback_graph, nrow = 1, asp = 0.5, vertex.size = 50)

p_values <- runif(m, 0, alpha)

test_results <-
  graph_test_shortcut(
    fallback_graph,
    p = p_values,
    alpha = alpha
  )

test_results$outputs$rejected

## ----fallback-improved--------------------------------------------------------
set.seed(1234)
alpha <- 0.025
hypotheses <- c(0.5, 0.3, 0.2)
m <- length(hypotheses)
fallback_improved_1_graph <- fallback_improved_1(hypotheses)
# transitions <- rbind(
#   c(0, 1, 0),
#   c(0, 0, 1),
#   c(hypotheses[seq_len(m - 1)] / sum(hypotheses[seq_len(m - 1)]), 0)
# )
# fallback_improved_1_graph <- graph_create(hypotheses, transitions)
plot_layout <- rbind(
  c(0, 0.8),
  c(0.8, 0.8),
  c(0.4, 0)
)

plot(
  fallback_improved_1_graph,
  layout = plot_layout,
  asp = 0.5,
  edge_curves = c(pairs = 0.5),
  vertex.size = 70
)

epsilon <- 0.0001
fallback_improved_2_graph <- fallback_improved_2(hypotheses, epsilon)
# transitions <- rbind(
#   c(0, 1, 0),
#   c(1 - epsilon, 0, epsilon),
#   c(1, 0, 0)
# )
# fallback_improved_2_graph <- graph_create(hypotheses, transitions)
plot_layout <- rbind(
  c(0, 0.8),
  c(0.8, 0.8),
  c(0.4, 0)
)

plot(
  fallback_improved_2_graph,
  layout = plot_layout,
  eps = 0.0001,
  asp = 0.5,
  edge_curves = c(pairs = 0.5),
  vertex.size = 70
)

p_values <- runif(m, 0, alpha)

test_results <-
  graph_test_shortcut(
    fallback_improved_2_graph,
    p = p_values,
    alpha = alpha
  )

test_results$outputs$rejected

## ----serial-gatekeeping-------------------------------------------------------
set.seed(1234)
alpha <- 0.025
m <- 3
epsilon <- 0.0001

transitions <- rbind(
  c(0, 1, 0),
  c(1 - epsilon, 0, epsilon),
  c(0, 0, 0)
)

serial_gatekeeping_graph <- graph_create(c(0.5, 0.5, 0), transitions)

plot_layout <- rbind(
  c(0, 0.8),
  c(0.8, 0.8),
  c(0.4, 0)
)

plot(
  serial_gatekeeping_graph,
  layout = plot_layout,
  eps = 0.0001,
  asp = 0.5,
  edge_curves = c(pairs = 0.5),
  vertex.size = 70
)

p_values <- runif(m, 0, alpha)

test_results <-
  graph_test_shortcut(
    serial_gatekeeping_graph,
    p = p_values,
    alpha = alpha
  )

test_results$outputs$rejected

## ----parallel-gatekeeping-----------------------------------------------------
set.seed(1234)
alpha <- 0.025
m <- 4

transitions <- rbind(
  c(0, 0, 0.5, 0.5),
  c(0, 0, 0.5, 0.5),
  c(0, 0, 0, 1),
  c(0, 0, 1, 0)
)

parallel_gatekeeping_graph <- graph_create(c(0.5, 0.5, 0, 0), transitions)

plot(parallel_gatekeeping_graph, vertex.size = 70)

p_values <- runif(m, 0, alpha)

test_results <-
  graph_test_shortcut(
    parallel_gatekeeping_graph,
    p = p_values,
    alpha = alpha
  )

test_results$outputs$rejected

## ----parallel-gatekeeping-improved--------------------------------------------
set.seed(1234)
alpha <- 0.025
m <- 4
epsilon <- 0.0001

transitions <- rbind(
  c(0, 0, 0.5, 0.5),
  c(0, 0, 0.5, 0.5),
  c(epsilon, 0, 0, 1 - epsilon),
  c(0, epsilon, 1 - epsilon, 0)
)

parallel_gatekeeping_improved_graph <-
  graph_create(c(0.5, 0.5, 0, 0), transitions)

plot(parallel_gatekeeping_improved_graph, eps = 0.0001, vertex.size = 70)

p_values <- runif(m, 0, alpha)

test_results <-
  graph_test_shortcut(
    parallel_gatekeeping_improved_graph,
    p = p_values,
    alpha = alpha
  )

test_results$outputs$rejected

## ----simple-successive-1------------------------------------------------------
set.seed(1234)
alpha <- 0.025
m <- 4
simple_successive_graph <- simple_successive_1()
# transitions <- rbind(
#   c(0, 0, 1, 0),
#   c(0, 0, 0, 1),
#   c(0, 1, 0, 0),
#   c(1, 0, 0, 0)
# )
# simple_successive_graph <- graph_create(c(0.5, 0.5, 0, 0), transitions)

plot(simple_successive_graph, layout = "grid", nrow = 2, vertex.size = 70)

p_values <- runif(m, 0, alpha)

test_results <-
  graph_test_shortcut(
    simple_successive_graph,
    p = p_values,
    alpha = alpha
  )

test_results$outputs$rejected

## ----simple-successive-var----------------------------------------------------
set.seed(1234)
alpha <- 0.025
m <- 4

successive_var <- simple_successive_var <- function(gamma) {
  graph_create(
    c(0.5, 0.5, 0, 0),
    rbind(
      c(0, gamma, 1 - gamma, 0),
      c(gamma, 0, 0, 1 - gamma),
      c(0, 1, 0, 0),
      c(1, 0, 0, 0)
    )
  )
}

successive_var_graph <- successive_var(0.5)
plot(successive_var_graph, layout = "grid", nrow = 2, vertex.size = 70)
p_values <- runif(m, 0, alpha)

test_results <-
  graph_test_shortcut(
    successive_var_graph,
    p = p_values,
    alpha = alpha
  )

test_results$outputs$rejected

## ----hochberg-----------------------------------------------------------------
set.seed(1234)
alpha <- 0.025
m <- 3
hochberg_graph <- hochberg(m)

plot(
  hochberg_graph,
  layout = igraph::layout_in_circle(
    as_igraph(hochberg_graph),
    order = c(2, 1, 3)
  ),
  vertex.size = 70
)

p_values <- runif(m, 0, alpha)

test_results <- graph_test_closure(
  hochberg_graph,
  p = p_values,
  alpha = alpha,
  test_types = "hochberg"
)

test_results$outputs$rejected

## ----hommel-------------------------------------------------------------------
set.seed(1234)
alpha <- 0.025
m <- 3
hommel_graph <- hommel(m)

plot(
  hommel_graph,
  layout = igraph::layout_in_circle(
    as_igraph(hommel_graph),
    order = c(2, 1, 3)
  ),
  vertex.size = 70
)

p_values <- runif(m, 0, alpha)

test_results <- graph_test_closure(
  hommel_graph,
  p = p_values,
  alpha = alpha,
  test_types = "simes"
)

test_results$outputs$rejected

## ----sidak--------------------------------------------------------------------
set.seed(1234)
alpha <- 0.025
m <- 3
sidak_graph <- sidak(m)

plot(
  sidak_graph,
  layout = igraph::layout_in_circle(
    as_igraph(bonferroni_graph),
    order = c(2, 1, 3)
  ),
  vertex.size = 70
)

p_values <- runif(m, 0, alpha)
corr <- diag(m)

test_results <- graph_test_closure(
  sidak_graph,
  p = p_values, alpha = alpha,
  test_types = "parametric",
  test_corr = list(corr)
)

test_results$outputs$rejected

## ----dunnett-test-------------------------------------------------------------
set.seed(1234)
alpha <- 0.025
m <- 3
dunnett_graph <- dunnett_single_step(m)

plot(
  dunnett_graph,
  layout = igraph::layout_in_circle(
    as_igraph(dunnett_graph),
    order = c(2, 1, 3)
  ),
  vertex.size = 70
)

p_values <- runif(m, 0, alpha)
corr <- matrix(0.5, m, m)
diag(corr) <- 1

test_results <- graph_test_closure(
  dunnett_graph,
  p = p_values, alpha = alpha,
  test_types = "parametric",
  test_corr = list(corr)
)

test_results$outputs$rejected

## ----weighted-dunnett-test----------------------------------------------------
set.seed(1234)
alpha <- 0.025
hypotheses <- c(0.5, 0.3, 0.2)
dunnett_graph <- dunnett_single_step_weighted(hypotheses)

plot(
  dunnett_graph,
  layout = igraph::layout_in_circle(
    as_igraph(dunnett_graph),
    order = c(2, 1, 3)
  ),
  vertex.size = 70
)

p_values <- runif(m, 0, alpha)
corr <- matrix(0.5, m, m)
diag(corr) <- 1

test_results <- graph_test_closure(
  dunnett_graph,
  p = p_values, alpha = alpha,
  test_types = "parametric",
  test_corr = list(corr)
)

test_results$outputs$rejected

## ----dunnett-procedure--------------------------------------------------------
set.seed(1234)
alpha <- 0.025
hypotheses <- c(0.5, 0.3, 0.2)
dunnett_graph <- dunnett_closure_weighted(hypotheses)

plot(
  dunnett_graph,
  layout = igraph::layout_in_circle(
    as_igraph(dunnett_graph),
    order = c(2, 1, 3)
  ),
  vertex.size = 70
)

p_values <- runif(m, 0, alpha)
corr <- matrix(0.5, m, m)
diag(corr) <- 1

test_results <- graph_test_closure(
  dunnett_graph,
  p = p_values, alpha = alpha,
  test_types = "parametric",
  test_corr = list(corr)
)

test_results$outputs$rejected

