# Adaptive Endpoint Selection Analysis - Journal Submission Version
# Minimal, clean implementation focusing on 4 key scenarios
# Simulates adaptive endpoint selection vs bonferroni fixed designs in rare disease trials

# Setup ----

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Parameters ----

ALPHA <- 0.025
RHO <- 0.5
CONTROL_MEANS <- c(0, 0, 0)
ENDPOINT_SD <- c(7, 4, 6)
N_ITERATIONS <- 10000
INFO_FRACTIONS <- c(0.25, 0.3, 0.5, 0.6, 0.75, 0.8)

# Effect size scenarios (standardized effects)
SCENARIOS <- list(
  "FVC only" = c(3.5, 0, 0),
  "FVC + 6MWT" = c(3.5, 1, 0),
  "All endpoints" = c(3.5, 2, 3),
  "Null" = c(0, 0, 0)
)

SCENARIO_LABELS <- c(
  "FVC=0.5  6MWT=0  MIP=0",
  "FVC=0.5  6MWT=0.25  MIP=0",
  "FVC=0.5  6MWT=0.5  MIP=0.5",
  "FVC=0  6MWT=0  MIP=0"
)

# Core Functions ----

create_covariance_matrix <- function(sd_vector, correlation) {
  n <- length(sd_vector)
  sigma <- outer(sd_vector, sd_vector) * correlation
  diag(sigma) <- sd_vector^2
  return(sigma)
}

calculate_z_statistic <- function(treatment, control) {
  n1 <- length(treatment)
  n2 <- length(control)
  pooled_var <- (var(treatment) * (n1-1) + var(control) * (n2-1)) / (n1+n2-2)
  pooled_se <- sqrt(pooled_var * (1/n1 + 1/n2))
  return((mean(treatment) - mean(control)) / pooled_se)
}

select_optimal_endpoint <- function(z_statistics, info_fraction, alpha = ALPHA) {
  adjusted_z <- (qnorm(1-alpha) - z_statistics * sqrt(info_fraction)) / sqrt(1-info_fraction)
  return(which.min(adjusted_z))
}

# Main Simulation Function ----

run_adaptive_simulation <- function(treatment_means, info_fraction, n_iter = N_ITERATIONS) {

  # Setup
  sigma <- create_covariance_matrix(ENDPOINT_SD, RHO)
  n_trt <- 120
  n_ctrl <- 60

  # Results storage
  results <- list(
    selection_counts = rep(0, 3),
    bias_sum = rep(0, 3),
    label_bias_sum = rep(0, 3),
    label_count = rep(0, 3),  # Count successful studies per endpoint
    power_count = 0,
    bonferroni_power = rep(0, 3),
    bonferroni_bias = rep(0, 3),
    bonferroni_label_bias = rep(0, 3)
  )

  set.seed(227)

  for(i in 1:n_iter) {
    # Generate data
    trt_data <- rmvnorm(n_trt, treatment_means, sigma)
    ctrl_data <- rmvnorm(n_ctrl, CONTROL_MEANS, sigma)

    # Fixed design analysis (Bonferroni)
    z_fixed <- sapply(1:3, function(j) calculate_z_statistic(trt_data[,j], ctrl_data[,j]))
    p_fixed <- 1 - pnorm(z_fixed)
    effect_sizes <- z_fixed * sqrt(1/n_trt + 1/n_ctrl)

    # Bonferroni correction
    bonf_significant <- p_fixed < (ALPHA / 3)
    results$bonferroni_power <- results$bonferroni_power + bonf_significant

    # Calculate bias for Bonferroni (same way as adaptive method)
    true_effects <- treatment_means / ENDPOINT_SD
    bonf_bias <- effect_sizes - true_effects

    # All endpoints get bias calculated if tested
    results$bonferroni_bias <- results$bonferroni_bias + bonf_bias

    # Label bias only for significant endpoints
    results$bonferroni_label_bias <- results$bonferroni_label_bias + bonf_bias * bonf_significant

    # Adaptive design - random cohort (focus on this for main analysis)
    stage1_ids_trt <- sample(n_trt, floor(n_trt * info_fraction))
    stage1_ids_ctrl <- sample(n_ctrl, floor(n_ctrl * info_fraction))

    # Stage 1 data (information cohort)
    trt_s1 <- trt_data[stage1_ids_trt, ]
    ctrl_s1 <- ctrl_data[stage1_ids_ctrl, ]

    # Stage 2 data (confirmatory cohort)
    trt_s2 <- trt_data[-stage1_ids_trt, ]
    ctrl_s2 <- ctrl_data[-stage1_ids_ctrl, ]

    # Calculate stage 1 test statistics
    z_stage1 <- sapply(1:3, function(j) calculate_z_statistic(trt_s1[,j], ctrl_s1[,j]))

    # Select primary endpoint
    selected_ep <- select_optimal_endpoint(z_stage1, info_fraction)
    results$selection_counts[selected_ep] <- results$selection_counts[selected_ep] + 1

    # Calculate final test statistic for selected endpoint
    z_stage2 <- calculate_z_statistic(trt_s2[,selected_ep], ctrl_s2[,selected_ep])

    # Partition hypothesis testing (full implementation matching original)
    p_stage1 <- 1 - pnorm(z_stage1)

    # H0{1,2,3} - global null
    min_p_global <- min(min(p_stage1) * 3,
                       min(p_stage1[p_stage1 != min(p_stage1)]) * 2,
                       max(p_stage1),
                       1 - 1e-5)
    h0_123 <- sqrt(info_fraction) * qnorm(1 - min_p_global) +
              sqrt(1 - info_fraction) * z_stage2

    # Conditional hypothesis tests based on selected endpoint
    if (selected_ep %in% c(2, 3)) {
      min_p_23 <- min(min(p_stage1[c(2, 3)]) * 2, max(p_stage1[c(2, 3)]), 1 - 1e-5)
      h0_23 <- sqrt(info_fraction) * qnorm(1 - min_p_23) +
               sqrt(1 - info_fraction) * z_stage2
    } else {
      h0_23 <- 10  # Large value ensures rejection
    }

    if (selected_ep %in% c(1, 3)) {
      min_p_13 <- min(min(p_stage1[c(1, 3)]) * 2, max(p_stage1[c(1, 3)]), 1 - 1e-5)
      h0_13 <- sqrt(info_fraction) * qnorm(1 - min_p_13) +
               sqrt(1 - info_fraction) * z_stage2
    } else {
      h0_13 <- 10
    }

    if (selected_ep %in% c(1, 2)) {
      min_p_12 <- min(min(p_stage1[c(1, 2)]) * 2, max(p_stage1[c(1, 2)]), 1 - 1e-5)
      h0_12 <- sqrt(info_fraction) * qnorm(1 - min_p_12) +
               sqrt(1 - info_fraction) * z_stage2
    } else {
      h0_12 <- 10
    }

    # All partition hypotheses must be rejected for success
    critical_value <- qnorm(1 - ALPHA)
    partition_tests <- c(h0_123, h0_23, h0_13, h0_12)
    hypothesis_rejected <- all(partition_tests > critical_value)

    # Calculate effect and bias for selected endpoint
    selected_effect <- z_fixed[selected_ep] * sqrt(1/n_trt + 1/n_ctrl)
    true_effect <- treatment_means[selected_ep] / ENDPOINT_SD[selected_ep]
    endpoint_bias <- selected_effect - true_effect

    # Always accumulate bias for selected endpoint
    results$bias_sum[selected_ep] <- results$bias_sum[selected_ep] + endpoint_bias

    if(hypothesis_rejected) {
      results$power_count <- results$power_count + 1
      # Label bias only for successful studies
      results$label_bias_sum[selected_ep] <- results$label_bias_sum[selected_ep] + endpoint_bias
      results$label_count[selected_ep] <- results$label_count[selected_ep] + 1
    }
  }

  # Calculate final metrics
  selection_probs <- results$selection_counts / n_iter
  mean_bias <- results$bias_sum / pmax(results$selection_counts, 1)
  label_bias <- results$label_bias_sum / pmax(results$label_count, 1)
  power <- results$power_count / n_iter
  bonf_power <- results$bonferroni_power / n_iter
  bonf_bias <- results$bonferroni_bias / n_iter
  bonf_label_bias <- results$bonferroni_label_bias / pmax(results$bonferroni_power, 1)

  return(list(
    selection_probs = selection_probs,
    mean_bias = mean_bias,
    label_bias = label_bias,
    power = power,
    bonferroni_power = bonf_power,
    bonferroni_bias = bonf_bias,
    bonferroni_label_bias = bonf_label_bias
  ))
}

# Run Simulations ----

run_all_simulations <- function() {
  results <- expand.grid(
    scenario = names(SCENARIOS),
    info_fraction = INFO_FRACTIONS,
    stringsAsFactors = FALSE
  )

  results$adaptive_bias_fvc <- NA
  results$adaptive_bias_6mwt <- NA
  results$adaptive_bias_mip <- NA
  results$adaptive_label_bias_fvc <- NA
  results$adaptive_label_bias_6mwt <- NA
  results$adaptive_label_bias_mip <- NA
  results$bonf_bias_fvc <- NA
  results$bonf_bias_6mwt <- NA
  results$bonf_bias_mip <- NA
  results$bonf_label_bias_fvc <- NA
  results$bonf_label_bias_6mwt <- NA
  results$bonf_label_bias_mip <- NA
  results$adaptive_power <- NA
  results$bonf_power_any <- NA

  cat("Running simulations...\n")
  for(i in 1:nrow(results)) {
    cat("Scenario:", results$scenario[i], "IF:", results$info_fraction[i], "\n")

    sim_result <- run_adaptive_simulation(
      treatment_means = SCENARIOS[[results$scenario[i]]],
      info_fraction = results$info_fraction[i]
    )

    results$adaptive_bias_fvc[i] <- sim_result$mean_bias[1]
    results$adaptive_bias_6mwt[i] <- sim_result$mean_bias[2]
    results$adaptive_bias_mip[i] <- sim_result$mean_bias[3]
    results$adaptive_label_bias_fvc[i] <- sim_result$label_bias[1]
    results$adaptive_label_bias_6mwt[i] <- sim_result$label_bias[2]
    results$adaptive_label_bias_mip[i] <- sim_result$label_bias[3]
    results$bonf_bias_fvc[i] <- sim_result$bonferroni_bias[1]
    results$bonf_bias_6mwt[i] <- sim_result$bonferroni_bias[2]
    results$bonf_bias_mip[i] <- sim_result$bonferroni_bias[3]
    results$bonf_label_bias_fvc[i] <- sim_result$bonferroni_label_bias[1]
    results$bonf_label_bias_6mwt[i] <- sim_result$bonferroni_label_bias[2]
    results$bonf_label_bias_mip[i] <- sim_result$bonferroni_label_bias[3]
    results$adaptive_power[i] <- sim_result$power
    results$bonf_power_any[i] <- max(sim_result$bonferroni_power)
  }

  return(results)
}

# Plotting Functions ----

create_publication_plots <- function(sim_results) {

  # Prepare data for bias plots
  bias_data <- sim_results |>
    select(scenario, info_fraction,
           adaptive_bias_fvc, adaptive_bias_6mwt, adaptive_bias_mip,
           bonf_bias_fvc, bonf_bias_6mwt, bonf_bias_mip) |>
    pivot_longer(cols = -c(scenario, info_fraction),
                 names_to = "method_endpoint", values_to = "bias") |>
    mutate(
      method = ifelse(grepl("^adaptive", method_endpoint), "Adaptive", "Bonferroni"),
      endpoint = case_when(
        grepl("fvc$", method_endpoint) ~ "FVC",
        grepl("6mwt$", method_endpoint) ~ "6MWT",
        grepl("mip$", method_endpoint) ~ "MIP"
      ),
      endpoint = factor(endpoint, levels = c("FVC", "6MWT", "MIP")),
      method = factor(method, levels = c("Adaptive", "Bonferroni"))
    )

  # Label bias data
  label_bias_data <- sim_results |>
    select(scenario, info_fraction,
           adaptive_label_bias_fvc, adaptive_label_bias_6mwt, adaptive_label_bias_mip,
           bonf_label_bias_fvc, bonf_label_bias_6mwt, bonf_label_bias_mip) |>
    pivot_longer(cols = -c(scenario, info_fraction),
                 names_to = "method_endpoint", values_to = "label_bias") |>
    mutate(
      method = ifelse(grepl("^adaptive", method_endpoint), "Adaptive", "Bonferroni"),
      endpoint = case_when(
        grepl("fvc$", method_endpoint) ~ "FVC",
        grepl("6mwt$", method_endpoint) ~ "6MWT",
        grepl("mip$", method_endpoint) ~ "MIP"
      ),
      endpoint = factor(endpoint, levels = c("FVC", "6MWT", "MIP")),
      method = factor(method, levels = c("Adaptive", "Bonferroni"))
    )

  # Create plots
  bias_plot <- ggplot(bias_data, aes(x = factor(info_fraction), y = bias,
                                    color = endpoint, linetype = method)) +
    geom_line(aes(group = interaction(method, endpoint)), size = 1) +
    facet_wrap(~ scenario, nrow = 1,
               labeller = labeller(scenario = setNames(SCENARIO_LABELS, names(SCENARIOS)))) +
    theme_bw() +
    labs(x = "Information Fraction", y = "Mean Bias",
         color = "Endpoint", linetype = "Method") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(-0.01, 0.5))

  label_bias_plot <- ggplot(label_bias_data, aes(x = factor(info_fraction), y = label_bias,
                                                 color = endpoint, linetype = method)) +
    geom_line(aes(group = interaction(method, endpoint)), size = 1) +
    facet_wrap(~ scenario, nrow = 1,
               labeller = labeller(scenario = setNames(SCENARIO_LABELS, names(SCENARIOS)))) +
    theme_bw() +
    labs(x = "Information Fraction", y = "Mean Label Bias",
         color = "Endpoint", linetype = "Method") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(0, 0.5))

  # Combine plots
  combined_figure <- ggarrange(bias_plot, label_bias_plot,
                              common.legend = TRUE, legend = "bottom",
                              labels = c("A", "B"), ncol = 1)

  return(list(
    combined = combined_figure,
    bias = bias_plot,
    label_bias = label_bias_plot
  ))
}

# Main Execution ----

main_analysis <- function() {
  cat("Starting adaptive endpoint selection analysis...\n")

  # Run simulations
  simulation_results <- run_all_simulations()

  # Create plots
  plots <- create_publication_plots(simulation_results)

  # Save results
  save(simulation_results, file = "adaptive_endpoint_results.RData")

  # Save main figure
  ggsave("adaptive_endpoint_main_figure.pdf", plots$combined,
         width = 12, height = 8, dpi = 300)

  cat("Analysis complete! Results saved.\n")
  return(list(results = simulation_results, plots = plots))
}

# Uncomment to run:
final_results <- main_analysis()
