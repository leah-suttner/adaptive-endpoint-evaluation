# Random Cohort Variability Analysis - Journal Submission Version
# Evaluates endpoint selection variability due to information cohort sampling
# Companion analysis to adaptive endpoint selection study

# Setup ----

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr)

# Parameters ----

ALPHA <- 0.025
RHO <- 0.5
CONTROL_MEANS <- c(0, 0, 0)
ENDPOINT_SD <- c(7, 4, 6)
N_DATASETS <- 1000
N_RESAMPLES <- 1000
INFO_FRACTIONS <- c(0.25, 0.3, 0.5, 0.6, 0.75, 0.8)

# Effect size scenarios (reduced to 4 key scenarios)
SCENARIOS <- rbind(
  c(3.5, 0, 0),    # FVC only
  c(3.5, 1, 0),    # FVC + 6MWT
  c(3.5, 2, 3),    # All endpoints
  c(0, 0, 0)       # Null
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

# Random Cohort Simulation Function ----

simulate_random_cohort_variability <- function(treatment_means, info_fraction, n_resamples = N_RESAMPLES) {

  # Setup
  sigma <- create_covariance_matrix(ENDPOINT_SD, RHO)
  n_trt <- 120
  n_ctrl <- 60

  # Generate single dataset (this stays fixed across resamples)
  trt_data <- rmvnorm(n_trt, treatment_means, sigma)
  ctrl_data <- rmvnorm(n_ctrl, CONTROL_MEANS, sigma)

  # Calculate fixed design test statistics (for reference)
  z_fixed <- sapply(1:3, function(j) calculate_z_statistic(trt_data[,j], ctrl_data[,j]))

  # Storage for resample results
  results <- matrix(0, nrow = n_resamples, ncol = 7)
  # Columns: [H0_123, H0_23, H0_13, H0_12, z_selected, selected_endpoint, effect_size]

  for(i in 1:n_resamples) {
    # Random information cohort sampling
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

    # Calculate stage 2 test statistics
    z_stage2 <- sapply(1:3, function(j) calculate_z_statistic(trt_s2[,j], ctrl_s2[,j]))

    # Select primary endpoint based on stage 1 data
    selected_ep <- select_optimal_endpoint(z_stage1, info_fraction)

    # Calculate p-values for stage 1
    p_stage1 <- 1 - pnorm(z_stage1)

    # Partition hypothesis test statistics
    # H0{1,2,3} - global null
    min_p_adjusted <- min(min(p_stage1) * 3,
                         min(p_stage1[p_stage1 != min(p_stage1)]) * 2,
                         max(p_stage1),
                         1 - 1e-5)
    h0_123 <- sqrt(info_fraction) * qnorm(1 - min_p_adjusted) +
              sqrt(1 - info_fraction) * z_stage2[selected_ep]

    # Conditional hypothesis tests based on selected endpoint
    if (selected_ep %in% c(2, 3)) {
      # H0{2,3}
      min_p_23 <- min(min(p_stage1[c(2, 3)]) * 2, max(p_stage1[c(2, 3)]), 1 - 1e-5)
      h0_23 <- sqrt(info_fraction) * qnorm(1 - min_p_23) +
               sqrt(1 - info_fraction) * z_stage2[selected_ep]
    } else {
      h0_23 <- 10  # Large value ensures rejection
    }

    if (selected_ep %in% c(1, 3)) {
      # H0{1,3}
      min_p_13 <- min(min(p_stage1[c(1, 3)]) * 2, max(p_stage1[c(1, 3)]), 1 - 1e-5)
      h0_13 <- sqrt(info_fraction) * qnorm(1 - min_p_13) +
               sqrt(1 - info_fraction) * z_stage2[selected_ep]
    } else {
      h0_13 <- 10
    }

    if (selected_ep %in% c(1, 2)) {
      # H0{1,2}
      min_p_12 <- min(min(p_stage1[c(1, 2)]) * 2, max(p_stage1[c(1, 2)]), 1 - 1e-5)
      h0_12 <- sqrt(info_fraction) * qnorm(1 - min_p_12) +
               sqrt(1 - info_fraction) * z_stage2[selected_ep]
    } else {
      h0_12 <- 10
    }

    # Store results
    effect_size <- z_fixed[selected_ep] * sqrt(1/n_trt + 1/n_ctrl)
    results[i, ] <- c(h0_123, h0_23, h0_13, h0_12, z_fixed[selected_ep], selected_ep, effect_size)
  }

  # Calculate endpoint selection probabilities
  selection_counts <- table(factor(results[, 6], levels = 1:3))
  selection_probs <- as.numeric(selection_counts / n_resamples)

  # Calculate power (all partition hypotheses rejected)
  critical_value <- qnorm(1 - ALPHA)
  power_by_endpoint <- sapply(1:3, function(ep) {
    ep_results <- results[results[, 6] == ep, 1:4, drop = FALSE]
    if (nrow(ep_results) > 0) {
      sum(rowSums(ep_results > critical_value) == 4) / nrow(ep_results)
    } else {
      0
    }
  })

  overall_power <- sum(sapply(1:3, function(ep) {
    ep_results <- results[results[, 6] == ep, 1:4, drop = FALSE]
    if (nrow(ep_results) > 0) {
      sum(rowSums(ep_results > critical_value) == 4)
    } else {
      0
    }
  })) / n_resamples

  return(list(
    selection_probs = selection_probs,
    power_by_endpoint = power_by_endpoint,
    overall_power = overall_power,
    all_results = results
  ))
}

# Run Variability Analysis ----

run_variability_analysis <- function() {

  n_scenarios <- nrow(SCENARIOS)
  n_fractions <- length(INFO_FRACTIONS)
  total_combinations <- n_scenarios * n_fractions

  # Storage arrays for results across datasets
  prob_select_array <- array(0, dim = c(N_DATASETS, total_combinations, 3))
  power_array <- array(0, dim = c(N_DATASETS, total_combinations, 1))

  set.seed(1018)

  cat("Running random cohort variability analysis...\n")
  cat("Total datasets:", N_DATASETS, "x", N_RESAMPLES, "resamples each\n")
  cat("Total scenarios:", n_scenarios, "x", n_fractions, "information fractions =", total_combinations, "combinations per dataset\n")

  start_time <- Sys.time()

  for (dataset_idx in 1:N_DATASETS) {
    if (dataset_idx %% 50 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      remaining_datasets <- N_DATASETS - dataset_idx
      est_time_remaining <- (elapsed / dataset_idx) * remaining_datasets
      cat("Dataset:", dataset_idx, "/", N_DATASETS,
          "| Elapsed:", round(elapsed, 1), "min",
          "| Est. remaining:", round(est_time_remaining, 1), "min\n")
    }

    combination_idx <- 0

    for (frac_idx in seq_along(INFO_FRACTIONS)) {
      for (scen_idx in seq_len(n_scenarios)) {
        combination_idx <- combination_idx + 1

        # Progress for each scenario (like main simulation)
        if (dataset_idx <= 5 || dataset_idx %% 100 == 0) {
          cat("Dataset", dataset_idx, "- Scenario:", names(SCENARIOS)[scen_idx],
              "IF:", INFO_FRACTIONS[frac_idx], "\n")
        }

        # Run simulation for this scenario/fraction combination
        sim_result <- simulate_random_cohort_variability(
          treatment_means = SCENARIOS[scen_idx, ],
          info_fraction = INFO_FRACTIONS[frac_idx]
        )

        # Store results
        prob_select_array[dataset_idx, combination_idx, ] <- sim_result$selection_probs
        power_array[dataset_idx, combination_idx, 1] <- sim_result$overall_power
      }
    }
  }

  total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat("Variability analysis complete! Total time:", round(total_elapsed, 1), "minutes\n")

  return(list(
    prob_select = prob_select_array,
    power = power_array,
    scenarios = SCENARIOS[1:4, ],  # Only first 4 scenarios
    info_fractions = INFO_FRACTIONS,
    scenario_labels = SCENARIO_LABELS
  ))
}

# Plotting Functions ----

create_variability_plots <- function(variability_results) {

  # Prepare selection probability data
  prob_select_long <- melt(variability_results$prob_select)
  names(prob_select_long) <- c("Dataset", "Setting", "Endpoint", "SelectionProb")

  # Add descriptive labels
  prob_select_long <- prob_select_long |>
    mutate(
      endpoint_name = case_when(
        Endpoint == 1 ~ "FVC",
        Endpoint == 2 ~ "6MWT",
        Endpoint == 3 ~ "MIP"
      ),
      endpoint_name = factor(endpoint_name, levels = c("FVC", "6MWT", "MIP")),
      info_fraction = case_when(
        Setting %in% 1:4 ~ 0.25,
        Setting %in% 5:8 ~ 0.30,
        Setting %in% 9:12 ~ 0.50,
        Setting %in% 13:16 ~ 0.60,
        Setting %in% 17:20 ~ 0.75,
        Setting %in% 21:24 ~ 0.80
      ),
      scenario_num = ((Setting - 1) %% 4) + 1,
      # Create compact labels with less whitespace
      scenario_label = case_when(
        scenario_num == 1 ~ "FVC=0.5 6MWT=0 MIP=0",
        scenario_num == 2 ~ "FVC=0.5 6MWT=0.25 MIP=0",
        scenario_num == 3 ~ "FVC=0.5 6MWT=0.5 MIP=0.5",
        scenario_num == 4 ~ "FVC=0 6MWT=0 MIP=0"
      ),
      scenario_label = factor(scenario_label, levels = c(
        "FVC=0.5 6MWT=0 MIP=0",
        "FVC=0.5 6MWT=0.25 MIP=0",
        "FVC=0.5 6MWT=0.5 MIP=0.5",
        "FVC=0 6MWT=0 MIP=0"
      ))
    )

  # Selection variability boxplot
  selection_plot <- ggplot(prob_select_long,
                          aes(x = factor(info_fraction), y = SelectionProb, color = endpoint_name)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ scenario_label, nrow = 1) +
    labs(
      x = "Information Cohort Fraction (ICF) (%)",
      y = "Endpoint Selection",
      color = "Endpoint"
    ) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent) +
    scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%")) +
    theme(
      text = element_text(size = 8),
      strip.text = element_text(size = 7),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      legend.margin = margin(t = -5, unit = "pt"),
      axis.title.x = element_text(margin = margin(t = 5, unit = "pt"))
    )

  # Power variability data
  power_long <- melt(variability_results$power[,,1])
  names(power_long) <- c("Dataset", "Setting", "Power")

  power_long <- power_long |>
    mutate(
      info_fraction = case_when(
        Setting %in% 1:4 ~ 0.25,
        Setting %in% 5:8 ~ 0.30,
        Setting %in% 9:12 ~ 0.50,
        Setting %in% 13:16 ~ 0.60,
        Setting %in% 17:20 ~ 0.75,
        Setting %in% 21:24 ~ 0.80
      ),
      scenario_num = ((Setting - 1) %% 4) + 1,
      # Create compact labels with less whitespace
      scenario_label = case_when(
        scenario_num == 1 ~ "FVC=0.5 6MWT=0 MIP=0",
        scenario_num == 2 ~ "FVC=0.5 6MWT=0.25 MIP=0",
        scenario_num == 3 ~ "FVC=0.5 6MWT=0.5 MIP=0.5",
        scenario_num == 4 ~ "FVC=0 6MWT=0 MIP=0"
      ),
      scenario_label = factor(scenario_label, levels = c(
        "FVC=0.5 6MWT=0 MIP=0",
        "FVC=0.5 6MWT=0.25 MIP=0",
        "FVC=0.5 6MWT=0.5 MIP=0.5",
        "FVC=0 6MWT=0 MIP=0"
      )),
      # Calculate variance for power (binomial: p*(1-p))
      power_variance = Power * (1 - Power)
    )

  # Power variability boxplot
  power_plot <- ggplot(power_long,
                      aes(x = factor(info_fraction), y = power_variance)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ scenario_label, nrow = 1) +
    theme_bw() +
    labs(
      x = "Information Cohort Fraction (ICF) (%)",
      y = "Success Variability"
    ) +
    scale_x_discrete(labels = function(x) paste0(as.numeric(x) * 100, "%")) +
    theme(
      text = element_text(size = 8),
      strip.text = element_text(size = 7),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )

  return(list(
    selection_plot = selection_plot,
    power_plot = power_plot
  ))
}

# Main Execution ----

main_variability_analysis <- function() {
  cat("Starting random cohort variability analysis...\n")

  # Run analysis
  variability_results <- run_variability_analysis()

  # Create plots using improved plotting functions
  plots <- create_variability_plots(variability_results)

  # Save results
  save(variability_results, file = "random_cohort_variability_results.RData")

  # Save figures with improved formatting (standalone with facet labels)
  ggsave("endpoint_selection_variability.pdf", plots$selection_plot,
         width = 180, height = 70, units = "mm", dpi = 300)

  ggsave("study_success_variability.pdf", plots$power_plot,
         width = 180, height = 70, units = "mm", dpi = 300)

  # Display plots for preview
  print("Selection Variability Plot:")
  print(plots$selection_plot)

  print("Power Variability Plot:")
  print(plots$power_plot)

  cat("\nVariability analysis complete! Results and figures saved.\n")
  cat("Figures saved:\n")
  cat("- endpoint_selection_variability.pdf\n")
  cat("- study_success_variability.pdf\n")

  # Summary Statistics
  cat("\nGenerating summary statistics...\n")

  # Prepare data for summary statistics
  prob_select_long <- melt(variability_results$prob_select)
  names(prob_select_long) <- c("Dataset", "Setting", "Endpoint", "SelectionProb")

  prob_select_long <- prob_select_long |>
    mutate(
      endpoint_name = case_when(
        Endpoint == 1 ~ "FVC",
        Endpoint == 2 ~ "6MWT",
        Endpoint == 3 ~ "MIP"
      ),
      endpoint_name = factor(endpoint_name, levels = c("FVC", "6MWT", "MIP")),
      scenario_num = ((Setting - 1) %% 4) + 1,
      scenario_label = case_when(
        scenario_num == 1 ~ "FVC=0.5 6MWT=0 MIP=0",
        scenario_num == 2 ~ "FVC=0.5 6MWT=0.25 MIP=0",
        scenario_num == 3 ~ "FVC=0.5 6MWT=0.5 MIP=0.5",
        scenario_num == 4 ~ "FVC=0 6MWT=0 MIP=0"
      )
    )

  cat("\nSelection variability by scenario:\n")
  selection_summary <- prob_select_long |>
    group_by(scenario_label, endpoint_name) |>
    summarise(
      mean_selection = mean(SelectionProb),
      sd_selection = sd(SelectionProb),
      .groups = "drop"
    )
  print(selection_summary)

  # Power summary
  power_long <- melt(variability_results$power[,,1])
  names(power_long) <- c("Dataset", "Setting", "Power")

  power_long <- power_long |>
    mutate(
      info_fraction = case_when(
        Setting %in% 1:4 ~ 0.25,
        Setting %in% 5:8 ~ 0.30,
        Setting %in% 9:12 ~ 0.50,
        Setting %in% 13:16 ~ 0.60,
        Setting %in% 17:20 ~ 0.75,
        Setting %in% 21:24 ~ 0.80
      ),
      scenario_num = ((Setting - 1) %% 4) + 1,
      scenario_label = case_when(
        scenario_num == 1 ~ "FVC=0.5 6MWT=0 MIP=0",
        scenario_num == 2 ~ "FVC=0.5 6MWT=0.25 MIP=0",
        scenario_num == 3 ~ "FVC=0.5 6MWT=0.5 MIP=0.5",
        scenario_num == 4 ~ "FVC=0 6MWT=0 MIP=0"
      ),
      power_variance = Power * (1 - Power)
    )

  cat("\nPower variability by scenario:\n")
  power_summary <- power_long |>
    group_by(scenario_label, info_fraction) |>
    summarise(
      mean_power = mean(Power),
      sd_power = sd(Power),
      mean_variance = mean(power_variance),
      .groups = "drop"
    )
  print(power_summary)

  return(list(
    results = variability_results,
    plots = plots,
    summaries = list(selection = selection_summary, power = power_summary)
  ))
}

# Uncomment to run analysis
variability_analysis <- main_variability_analysis()