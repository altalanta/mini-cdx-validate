suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(boot)
})

source("R/utils.R")

load_analytical_data <- function(path = "data/synthetic/analytical_runs.csv") {
  readr::read_csv(path, show_col_types = FALSE)
}

load_truth_data <- function(path = "data/synthetic/truth_reference.csv") {
  readr::read_csv(path, show_col_types = FALSE)
}

# Compute LOD values using probit regression with bootstrap CIs
generate_lod_summary <- function(data, targets = c(0.5, 0.95), n_boot = 300) {
  data %>%
    filter(nominal_vaf > 0) %>%
    mutate(log_nominal = log10(nominal_vaf)) %>%
    group_by(variant_id, gene) %>%
    group_modify(~ compute_variant_lod(.x, targets = targets, n_boot = n_boot)) %>%
    ungroup()
}

compute_variant_lod <- function(df, targets, n_boot) {
  model <- glm(detected ~ log_nominal, family = binomial(link = "probit"), data = df)
  coefs <- coef(model)
  slope <- coefs[["log_nominal"]]
  intercept <- coefs[["(Intercept)"]]

  estimate_lod <- function(int, slp, target_prob) {
    10 ^ ((qnorm(target_prob) - int) / slp)
  }

  point_estimates <- tibble(
    target = targets,
    estimate = estimate_lod(intercept, slope, targets)
  )

  boot_fun <- function(data, indices) {
    sampled <- data[indices, , drop = FALSE]
    boot_model <- try(glm(detected ~ log_nominal, family = binomial(link = "probit"), data = sampled), silent = TRUE)
    if (inherits(boot_model, "try-error")) {
      return(rep(NA_real_, length(targets)))
    }
    boot_coefs <- coef(boot_model)
    sapply(targets, function(prob) estimate_lod(boot_coefs[["(Intercept)"]], boot_coefs[["log_nominal"]], prob))
  }

  boot_res <- boot::boot(df, statistic = boot_fun, R = n_boot)
  boot_df <- as_tibble(boot_res$t, .name_repair = "minimal")
  names(boot_df) <- paste0("lod", targets * 100)

  ci_bounds <- map2_dfr(names(boot_df), point_estimates$target, function(col_nm, tgt) {
    vals <- boot_df[[col_nm]]
    tibble(
      target = tgt,
      lower = quantile(vals, probs = 0.025, na.rm = TRUE),
      upper = quantile(vals, probs = 0.975, na.rm = TRUE)
    )
  })

  summary <- point_estimates %>%
    left_join(ci_bounds, by = "target") %>%
    mutate(metric = paste0("LOD", target * 100)) %>%
    select(metric, estimate, lower, upper)

  pivot_wider(summary, names_from = metric, values_from = c(estimate, lower, upper))
}

# Precision summaries and pseudo variance components
summarise_precision <- function(data) {
  precision <- data %>%
    filter(nominal_vaf > 0) %>%
    group_by(variant_id, gene, nominal_vaf) %>%
    summarise(
      mean_measured = mean(measured_vaf),
      sd_measured = sd(measured_vaf),
      cv_pct = 100 * sd_measured / mean_measured,
      n = n(),
      .groups = "drop"
    )

  variance_components <- data %>%
    filter(nominal_vaf > 0) %>%
    group_by(variant_id, gene) %>%
    group_modify(~ {
      mod <- aov(measured_vaf ~ site + run_id + operator, data = .x)
      tidy(mod) %>%
        transmute(
          component = term,
          mean_sq = meansq,
          pct_contrib = mean_sq / sum(mean_sq, na.rm = TRUE)
        )
    }) %>%
    ungroup()

  list(precision = precision, variance_components = variance_components)
}

# Accuracy against truth reference at a target VAF threshold
summarise_accuracy <- function(data, truth, threshold = 5) {
  joined <- data %>%
    left_join(truth, by = c("sample_id", "variant_id", "gene"))

  positives <- joined %>% filter(truth_vaf >= threshold)
  negatives <- joined %>% filter(truth_vaf < threshold)

  ppa <- if (nrow(positives) > 0) sum(positives$detected == 1) else 0
  npa <- if (nrow(negatives) > 0) sum(negatives$detected == 0) else 0

  ppa_ci <- if (nrow(positives) > 0) wilson_ci(ppa, nrow(positives)) else tibble(estimate = NA_real_, lower = NA_real_, upper = NA_real_)
  npa_ci <- if (nrow(negatives) > 0) wilson_ci(npa, nrow(negatives)) else tibble(estimate = NA_real_, lower = NA_real_, upper = NA_real_)

  tibble(
    metric = c("PPA", "NPA"),
    estimate = c(ppa_ci$estimate[1], npa_ci$estimate[1]),
    lower = c(ppa_ci$lower[1], npa_ci$lower[1]),
    upper = c(ppa_ci$upper[1], npa_ci$upper[1]),
    denominator = c(nrow(positives), nrow(negatives))
  )
}

# Contamination detection at unexpected <0.3% VAF hits
summarise_contamination <- function(data, truth) {
  joined <- data %>%
    left_join(truth, by = c("sample_id", "variant_id", "gene"))

  contamination <- joined %>%
    filter(truth_present == 0, measured_vaf < 0.3)

  per_run <- contamination %>%
    group_by(run_id) %>%
    summarise(
      contamination_events = sum(detected == 1),
      total_blanks = n(),
      rate = ifelse(total_blanks > 0, contamination_events / total_blanks, NA_real_),
      .groups = "drop"
    )

  list(events = contamination, summary = per_run)
}

run_analytical_validation <- function(output_dir = "results/analytical", threshold = 5, n_boot = 300) {
  ensure_dir(file.path(output_dir, "dummy"))
  data <- load_analytical_data()
  truth <- load_truth_data()

  lod <- generate_lod_summary(data, targets = c(0.5, 0.95), n_boot = n_boot)
  precision <- summarise_precision(data)
  accuracy <- summarise_accuracy(data, truth, threshold = threshold)
  contamination <- summarise_contamination(data, truth)

  readr::write_csv(lod, file.path(output_dir, "lod_summary.csv"))
  readr::write_csv(precision$precision, file.path(output_dir, "precision_summary.csv"))
  readr::write_csv(precision$variance_components, file.path(output_dir, "variance_components.csv"))
  readr::write_csv(accuracy, file.path(output_dir, "accuracy_summary.csv"))
  readr::write_csv(contamination$summary, file.path(output_dir, "contamination_by_run.csv"))
  readr::write_csv(contamination$events, file.path(output_dir, "contamination_events.csv"))

  list(
    lod = lod,
    precision = precision,
    accuracy = accuracy,
    contamination = contamination
  )
}

if (sys.nframe() == 0) {
  set_analysis_seed(42)
  invisible(run_analytical_validation())
}
