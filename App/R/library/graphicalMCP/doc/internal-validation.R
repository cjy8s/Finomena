## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  fig.align = "center",
  collapse = TRUE,
  comment = "#>",
  cache.lazy = FALSE
)

## ----setup, message = FALSE, warning = FALSE----------------------------------
library(graphicalMCP)

## ----bonferroni-results-------------------------------------------------------
out <- read.csv(here::here("vignettes/internal-validation_bonferroni.csv"))
# Matching power using the adjusted p-value approach
all.equal(out$adjusted_p, rep(TRUE, nrow(out)))
# Matching power using the adjusted significance level approach
all.equal(out$adjusted_significance_level, rep(TRUE, nrow(out)))

## ----hochberg-results---------------------------------------------------------
out <- read.csv(here::here("vignettes/internal-validation_hochberg.csv"))
# Matching power using the adjusted p-value approach
all.equal(out$adjusted_p, rep(TRUE, nrow(out)))
# Matching power using the adjusted significance level approach
all.equal(out$adjusted_significance_level, rep(TRUE, nrow(out)))

## ----simes-results------------------------------------------------------------
out <- read.csv(here::here("vignettes/internal-validation_simes.csv"))
# Matching power using the adjusted p-value approach
all.equal(out$adjusted_p, rep(TRUE, nrow(out)))
# Matching power using the adjusted significance level approach
all.equal(out$adjusted_significance_level, rep(TRUE, nrow(out)))

## ----parametric-results-------------------------------------------------------
out <- read.csv(here::here("vignettes/internal-validation_parametric.csv"))
# Matching power using the adjusted p-value approach
all.equal(out$adjusted_p, rep(TRUE, nrow(out)))
# Matching power using the adjusted significance level approach
all.equal(out$adjusted_significance_level, rep(TRUE, nrow(out)))

## ----mixed-results------------------------------------------------------------
out <- read.csv(here::here("vignettes/internal-validation_mixed.csv"))
# Matching power using the adjusted p-value approach
all.equal(out$adjusted_p, rep(TRUE, nrow(out)))
# Matching power using the adjusted significance level approach
all.equal(out$adjusted_significance_level, rep(TRUE, nrow(out)))

## ----parametric-mixed-results-------------------------------------------------
out <- read.csv(here::here("vignettes/internal-validation_parametric-mixed.csv"))
# Matching power using the adjusted p-value approach
all.equal(out$adjusted_p, rep(TRUE, nrow(out)))
# Matching power using the adjusted significance level approach
all.equal(out$adjusted_significance_level, rep(TRUE, nrow(out)))

## ----bonferroni, include = FALSE, eval = FALSE--------------------------------
#  n_graph <- 1000
#  m <- 6
#  alpha <- 0.025
#  n_sim <- 1e2
#  sim_corr <- matrix(0.5, m, m)
#  diag(sim_corr) <- 1
#  out_adjusted_p <- rep(NA, n_graph)
#  out_adjusted_significance_level <- rep(NA, n_graph)
#  for (i in 1:n_graph) {
#    set.seed(1234 + i - 1)
#    graph <- random_graph(m)
#    marginal_power <- runif(m, 0.5, 0.9)
#    results_calculate_power <- graph_calculate_power(
#      graph,
#      alpha = alpha,
#      power_marginal = marginal_power,
#      sim_corr = sim_corr,
#      sim_n = n_sim,
#      verbose = TRUE
#    )
#    p_sim <- results_calculate_power$details$p_sim
#    adjusted_p_results <- p_sim
#    adjusted_significance_level_results <- p_sim
#    for (j in 1:n_sim) {
#      temp <- graph_test_shortcut(
#        graph,
#        p_sim[j, ],
#        alpha = alpha,
#        test_values = TRUE
#      )
#      adjusted_p_results[j, ] <- temp$outputs$rejected
#      temp_adjusted_significance_level <- temp$test_values$results$Inequality_holds
#      names(temp_adjusted_significance_level) <- temp$test_values$results$Hypothesis
#      adjusted_significance_level_results[j, ] <-
#        temp_adjusted_significance_level[names(temp$outputs$rejected)]
#    }
#    results_adjusted_p <- colMeans(adjusted_p_results)
#    names(results_adjusted_p) <- names(temp$outputs$rejected)
#    # Matching power using the adjusted p-value approach
#    out_adjusted_p[i] <- all.equal(
#      results_calculate_power$power$power_local,
#      results_adjusted_p,
#      tolerance = 1 / n_sim
#    )
#  
#    results_significance_level <- colMeans(adjusted_significance_level_results)
#    names(results_significance_level) <- names(temp$outputs$rejected)
#    # Matching power using the adjusted significance level approach
#    out_adjusted_significance_level[i] <- all.equal(
#      results_calculate_power$power$power_local,
#      results_significance_level,
#      tolerance = 1 / n_sim
#    )
#  }
#  out <- cbind(out_adjusted_p, out_adjusted_significance_level)
#  colnames(out) <- c("adjusted_p", "adjusted_significance_level")
#  write.csv(
#    out,
#    here::here("vignettes/internal-validation_bonferroni.csv"),
#    row.names = FALSE
#  )

## ----hochberg, include = FALSE, eval = FALSE----------------------------------
#  n_graph <- 1000
#  m <- 4
#  alpha <- 0.025
#  n_sim <- 1e2
#  sim_corr <- matrix(0.5, m, m)
#  diag(sim_corr) <- 1
#  out_adjusted_p <- rep(NA, n_graph)
#  out_adjusted_significance_level <- rep(NA, n_graph)
#  for (i in 1:n_graph) {
#    set.seed(1234 + i - 1)
#    groups <- sample(1:m)
#    test_groups <- list(groups[1:(m / 2)], groups[(m / 2 + 1):m])
#    graph <- random_graph(m)
#    graph$hypotheses[groups[1:(m / 2)]] <- sum(graph$hypotheses[groups[1:(m / 2)]]) / 3
#    graph$hypotheses[groups[(m / 2 + 1):m]] <- sum(graph$hypotheses[groups[(m / 2 + 1):m]]) / 3
#    graph$transitions <- matrix(1 / (m - 1), nrow = m, ncol = m)
#    diag(graph$transitions) <- 0
#    marginal_power <- runif(m, 0.5, 0.9)
#    results_calculate_power <- graph_calculate_power(
#      graph,
#      alpha = alpha,
#      power_marginal = marginal_power,
#      sim_corr = sim_corr,
#      sim_n = n_sim,
#      test_types = c("hochberg", "hochberg"),
#      test_groups = test_groups,
#      verbose = TRUE
#    )
#    p_sim <- results_calculate_power$details$p_sim
#    adjusted_p_results <- p_sim
#    adjusted_significance_level_results <- p_sim
#    for (j in 1:n_sim) {
#      temp <- graph_test_closure(
#        graph,
#        p_sim[j, ],
#        alpha = alpha,
#        test_types = c("hochberg", "hochberg"),
#        test_groups = test_groups,
#        test_values = TRUE
#      )
#      adjusted_p_results[j, ] <- temp$outputs$rejected
#      intersection <- unique(temp$test_values$results$Intersection)
#      intersection_rej <- intersection
#      for (k in 1:length(intersection)) {
#        intersection_rej[k] <- any(subset(
#          temp$test_values$results,
#          Intersection == intersection[k]
#        )$Inequality_holds)
#      }
#      for (k in 1:m) {
#        id <- substr(intersection, start = k, stop = k) == "1"
#        adjusted_significance_level_results[j, k] <- suppressWarnings(all(intersection_rej[id]))
#      }
#    }
#    results_adjusted_p <- colMeans(adjusted_p_results)
#    names(results_adjusted_p) <- names(temp$outputs$rejected)
#    # Matching power using the adjusted p-value approach
#    out_adjusted_p[i] <- all.equal(
#      results_calculate_power$power$power_local,
#      results_adjusted_p,
#      tolerance = 1 / n_sim
#    )
#  
#    results_significance_level <- colMeans(adjusted_significance_level_results)
#    names(results_significance_level) <- names(temp$outputs$rejected)
#    # Matching power using the adjusted significance level approach
#    out_adjusted_significance_level[i] <- all.equal(
#      results_calculate_power$power$power_local,
#      results_significance_level,
#      tolerance = 1 / n_sim
#    )
#  }
#  out <- cbind(out_adjusted_p, out_adjusted_significance_level)
#  colnames(out) <- c("adjusted_p", "adjusted_significance_level")
#  write.csv(
#    out,
#    here::here("vignettes/internal-validation_hochberg.csv"),
#    row.names = FALSE
#  )

## ----simes, include = FALSE, eval = FALSE-------------------------------------
#  n_graph <- 1000
#  m <- 4
#  alpha <- 0.025
#  n_sim <- 1e2
#  sim_corr <- matrix(0.5, m, m)
#  diag(sim_corr) <- 1
#  out_adjusted_p <- rep(NA, n_graph)
#  out_adjusted_significance_level <- rep(NA, n_graph)
#  for (i in 1:n_graph) {
#    set.seed(1234 + i - 1)
#    groups <- sample(1:m)
#    test_groups <- list(groups[1:(m / 2)], groups[(m / 2 + 1):m])
#    graph <- random_graph(m)
#    marginal_power <- runif(m, 0.5, 0.9)
#    results_calculate_power <- graph_calculate_power(
#      graph,
#      alpha = alpha,
#      power_marginal = marginal_power,
#      sim_corr = sim_corr,
#      sim_n = n_sim,
#      test_types = c("simes", "simes"),
#      test_groups = test_groups,
#      verbose = TRUE
#    )
#    p_sim <- results_calculate_power$details$p_sim
#    adjusted_p_results <- p_sim
#    adjusted_significance_level_results <- p_sim
#    for (j in 1:n_sim) {
#      temp <- graph_test_closure(
#        graph,
#        p_sim[j, ],
#        alpha = alpha,
#        test_types = c("simes", "simes"),
#        test_groups = test_groups,
#        test_values = TRUE
#      )
#      adjusted_p_results[j, ] <- temp$outputs$rejected
#      intersection <- unique(temp$test_values$results$Intersection)
#      intersection_rej <- intersection
#      for (k in 1:length(intersection)) {
#        intersection_rej[k] <- any(subset(
#          temp$test_values$results,
#          Intersection == intersection[k]
#        )$Inequality_holds)
#      }
#      for (k in 1:m) {
#        id <- substr(intersection, start = k, stop = k) == "1"
#        adjusted_significance_level_results[j, k] <- suppressWarnings(all(intersection_rej[id]))
#      }
#    }
#    results_adjusted_p <- colMeans(adjusted_p_results)
#    names(results_adjusted_p) <- names(temp$outputs$rejected)
#    # Matching power using the adjusted p-value approach
#    out_adjusted_p[i] <- all.equal(
#      results_calculate_power$power$power_local,
#      results_adjusted_p,
#      tolerance = 1 / n_sim
#    )
#  
#    results_significance_level <- colMeans(adjusted_significance_level_results)
#    names(results_significance_level) <- names(temp$outputs$rejected)
#    # Matching power using the adjusted significance level approach
#    out_adjusted_significance_level[i] <- all.equal(
#      results_calculate_power$power$power_local,
#      results_significance_level,
#      tolerance = 1 / n_sim
#    )
#  }
#  out <- cbind(out_adjusted_p, out_adjusted_significance_level)
#  colnames(out) <- c("adjusted_p", "adjusted_significance_level")
#  write.csv(
#    out,
#    here::here("vignettes/internal-validation_simes.csv"),
#    row.names = FALSE
#  )

## ----parametric, include = FALSE, eval = FALSE--------------------------------
#  n_graph <- 1000
#  m <- 4
#  alpha <- 0.025
#  n_sim <- 1e2
#  sim_corr <- matrix(0.5, m, m)
#  diag(sim_corr) <- 1
#  out_adjusted_p <- rep(NA, n_graph)
#  out_adjusted_significance_level <- rep(NA, n_graph)
#  for (i in 1:n_graph) {
#    set.seed(1234 + i - 1)
#    groups <- sample(1:m)
#    test_groups <- list(groups[1:(m / 2)], groups[(m / 2 + 1):m])
#    test_corr <- list(
#      sim_corr[groups[1:(m / 2)], groups[1:(m / 2)]],
#      sim_corr[groups[(m / 2 + 1):m], groups[(m / 2 + 1):m]]
#    )
#    graph <- random_graph(m)
#    marginal_power <- runif(m, 0.5, 0.9)
#    results_calculate_power <- graph_calculate_power(
#      graph,
#      alpha = alpha,
#      power_marginal = marginal_power,
#      sim_corr = sim_corr,
#      sim_n = n_sim,
#      test_types = c("parametric", "parametric"),
#      test_groups = test_groups,
#      test_corr = test_corr,
#      verbose = TRUE
#    )
#    p_sim <- results_calculate_power$details$p_sim
#    adjusted_p_results <- p_sim
#    adjusted_significance_level_results <- p_sim
#    for (j in 1:n_sim) {
#      temp <- graph_test_closure(
#        graph,
#        p_sim[j, ],
#        alpha = alpha,
#        test_types = c("parametric", "parametric"),
#        test_groups = test_groups,
#        test_corr = test_corr,
#        test_values = TRUE
#      )
#      adjusted_p_results[j, ] <- temp$outputs$rejected
#      intersection <- unique(temp$test_values$results$Intersection)
#      intersection_rej <- intersection
#      for (k in 1:length(intersection)) {
#        intersection_rej[k] <- any(subset(
#          temp$test_values$results,
#          Intersection == intersection[k]
#        )$Inequality_holds)
#      }
#      for (k in 1:m) {
#        id <- substr(intersection, start = k, stop = k) == "1"
#        adjusted_significance_level_results[j, k] <- suppressWarnings(all(intersection_rej[id]))
#      }
#    }
#    results_adjusted_p <- colMeans(adjusted_p_results)
#    names(results_adjusted_p) <- names(temp$outputs$rejected)
#    # Matching power using the adjusted p-value approach
#    out_adjusted_p[i] <- all.equal(
#      results_calculate_power$power$power_local,
#      results_adjusted_p,
#      tolerance = 1 / n_sim
#    )
#  
#    results_significance_level <- colMeans(adjusted_significance_level_results)
#    names(results_significance_level) <- names(temp$outputs$rejected)
#    # Matching power using the adjusted significance level approach
#    out_adjusted_significance_level[i] <- all.equal(
#      results_calculate_power$power$power_local,
#      results_significance_level,
#      tolerance = 1 / n_sim
#    )
#  }
#  out <- cbind(out_adjusted_p, out_adjusted_significance_level)
#  colnames(out) <- c("adjusted_p", "adjusted_significance_level")
#  write.csv(
#    out,
#    here::here("vignettes/internal-validation_parametric.csv"),
#    row.names = FALSE
#  )

## ----mixed, include = FALSE, eval = FALSE-------------------------------------
#  n_graph <- 1000
#  m <- 4
#  alpha <- 0.025
#  n_sim <- 1e2
#  sim_corr <- matrix(0.5, m, m)
#  diag(sim_corr) <- 1
#  out_adjusted_p <- rep(NA, n_graph)
#  out_adjusted_significance_level <- rep(NA, n_graph)
#  for (i in 1:n_graph) {
#    set.seed(1234 + i - 1)
#    groups <- sample(1:m)
#    test_groups <- list(groups[1:(m / 2)], groups[(m / 2 + 1):m])
#    test_corr <- list(
#      sim_corr[groups[1:(m / 2)], groups[1:(m / 2)]],
#      sim_corr[groups[(m / 2 + 1):m], groups[(m / 2 + 1):m]]
#    )
#    test_types <- sample(c("bonferroni", "hochberg", "simes"), 2)
#    graph <- random_graph(m)
#    graph$hypotheses[groups[1:(m / 2)]] <- sum(graph$hypotheses[groups[1:(m / 2)]]) / 3
#    graph$hypotheses[groups[(m / 2 + 1):m]] <- sum(graph$hypotheses[groups[(m / 2 + 1):m]]) / 3
#    graph$transitions <- matrix(1 / (m - 1), nrow = m, ncol = m)
#    diag(graph$transitions) <- 0
#    marginal_power <- runif(m, 0.5, 0.9)
#    results_calculate_power <- graph_calculate_power(
#      graph,
#      alpha = alpha,
#      power_marginal = marginal_power,
#      sim_corr = sim_corr,
#      sim_n = n_sim,
#      test_types = test_types,
#      test_groups = test_groups,
#      verbose = TRUE
#    )
#    p_sim <- results_calculate_power$details$p_sim
#    adjusted_p_results <- p_sim
#    adjusted_significance_level_results <- p_sim
#    for (j in 1:n_sim) {
#      temp <- graph_test_closure(
#        graph,
#        p_sim[j, ],
#        alpha = alpha,
#        test_types = test_types,
#        test_groups = test_groups,
#        test_values = TRUE
#      )
#      adjusted_p_results[j, ] <- temp$outputs$rejected
#      intersection <- unique(temp$test_values$results$Intersection)
#      intersection_rej <- intersection
#      for (k in 1:length(intersection)) {
#        intersection_rej[k] <- any(subset(
#          temp$test_values$results,
#          Intersection == intersection[k]
#        )$Inequality_holds)
#      }
#      for (k in 1:m) {
#        id <- substr(intersection, start = k, stop = k) == "1"
#        adjusted_significance_level_results[j, k] <- suppressWarnings(all(intersection_rej[id]))
#      }
#    }
#    results_adjusted_p <- colMeans(adjusted_p_results)
#    names(results_adjusted_p) <- names(temp$outputs$rejected)
#    # Matching power using the adjusted p-value approach
#    out_adjusted_p[i] <- all.equal(
#      results_calculate_power$power$power_local,
#      results_adjusted_p,
#      tolerance = 1 / n_sim
#    )
#  
#    results_significance_level <- colMeans(adjusted_significance_level_results)
#    names(results_significance_level) <- names(temp$outputs$rejected)
#    # Matching power using the adjusted significance level approach
#    out_adjusted_significance_level[i] <- all.equal(
#      results_calculate_power$power$power_local,
#      results_significance_level,
#      tolerance = 1 / n_sim
#    )
#  }
#  out <- cbind(out_adjusted_p, out_adjusted_significance_level)
#  colnames(out) <- c("adjusted_p", "adjusted_significance_level")
#  write.csv(
#    out,
#    here::here("vignettes/internal-validation_mixed.csv"),
#    row.names = FALSE
#  )

## ----parametric-mixed, include = FALSE, eval = FALSE--------------------------
#  n_graph <- 1000
#  m <- 4
#  alpha <- 0.025
#  n_sim <- 1e2
#  sim_corr <- matrix(0.5, m, m)
#  diag(sim_corr) <- 1
#  out_adjusted_p <- rep(NA, n_graph)
#  out_adjusted_significance_level <- rep(NA, n_graph)
#  for (i in 1:n_graph) {
#    set.seed(1234 + i - 1)
#    groups <- sample(1:m)
#    test_groups <- list(groups[1:(m / 2)], groups[(m / 2 + 1):m])
#    test_corr <- list(sim_corr[groups[1:(m / 2)], groups[1:(m / 2)]], NA)
#    test_types <- c(
#      "parametric",
#      sample(c("bonferroni", "hochberg", "simes"), 1)
#    )
#    graph <- random_graph(m)
#    graph$hypotheses[groups[1:(m / 2)]] <- sum(graph$hypotheses[groups[1:(m / 2)]]) / 3
#    graph$hypotheses[groups[(m / 2 + 1):m]] <- sum(graph$hypotheses[groups[(m / 2 + 1):m]]) / 3
#    graph$transitions <- matrix(1 / (m - 1), nrow = m, ncol = m)
#    diag(graph$transitions) <- 0
#    marginal_power <- runif(m, 0.5, 0.9)
#    results_calculate_power <- graph_calculate_power(
#      graph,
#      alpha = alpha,
#      power_marginal = marginal_power,
#      sim_corr = sim_corr,
#      sim_n = n_sim,
#      test_types = test_types,
#      test_groups = test_groups,
#      test_corr = test_corr,
#      verbose = TRUE
#    )
#    p_sim <- results_calculate_power$details$p_sim
#    adjusted_p_results <- p_sim
#    adjusted_significance_level_results <- p_sim
#    for (j in 1:n_sim) {
#      temp <- graph_test_closure(
#        graph,
#        p_sim[j, ],
#        alpha = alpha,
#        test_types = test_types,
#        test_groups = test_groups,
#        test_corr = test_corr,
#        test_values = TRUE
#      )
#      adjusted_p_results[j, ] <- temp$outputs$rejected
#      intersection <- unique(temp$test_values$results$Intersection)
#      intersection_rej <- intersection
#      for (k in 1:length(intersection)) {
#        intersection_rej[k] <- any(subset(
#          temp$test_values$results,
#          Intersection == intersection[k]
#        )$Inequality_holds)
#      }
#      for (k in 1:m) {
#        id <- substr(intersection, start = k, stop = k) == "1"
#        adjusted_significance_level_results[j, k] <- suppressWarnings(all(intersection_rej[id]))
#      }
#    }
#    results_adjusted_p <- colMeans(adjusted_p_results)
#    names(results_adjusted_p) <- names(temp$outputs$rejected)
#    # Matching power using the adjusted p-value approach
#    out_adjusted_p[i] <- all.equal(
#      results_calculate_power$power$power_local,
#      results_adjusted_p,
#      tolerance = 1 / n_sim
#    )
#  
#    results_significance_level <- colMeans(adjusted_significance_level_results)
#    names(results_significance_level) <- names(temp$outputs$rejected)
#    # Matching power using the adjusted significance level approach
#    out_adjusted_significance_level[i] <- all.equal(
#      results_calculate_power$power$power_local,
#      results_significance_level,
#      tolerance = 1 / n_sim
#    )
#  }
#  out <- cbind(out_adjusted_p, out_adjusted_significance_level)
#  colnames(out) <- c("adjusted_p", "adjusted_significance_level")
#  write.csv(
#    out,
#    here::here("vignettes/internal-validation_parametric-mixed.csv"),
#    row.names = FALSE
#  )

