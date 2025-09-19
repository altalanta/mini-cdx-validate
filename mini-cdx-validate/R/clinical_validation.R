suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(pROC)
  library(boot)
  library(survival)
})

source("R/utils.R")

load_clinical_data <- function(path = "data/synthetic/clinical_cohort.csv") {
  readr::read_csv(path, show_col_types = FALSE)
}

calc_classification_metrics <- function(data) {
  tp <- sum(data$biomarker_call == 1 & data$response == 1)
  fp <- sum(data$biomarker_call == 1 & data$response == 0)
  tn <- sum(data$biomarker_call == 0 & data$response == 0)
  fn <- sum(data$biomarker_call == 0 & data$response == 1)

  total_pos <- tp + fn
  total_neg <- tn + fp
  total_test_pos <- tp + fp
  total_test_neg <- tn + fn

  sens_ci <- if (total_pos > 0) wilson_ci(tp, total_pos) else tibble(estimate = NA_real_, lower = NA_real_, upper = NA_real_)
  spec_ci <- if (total_neg > 0) wilson_ci(tn, total_neg) else tibble(estimate = NA_real_, lower = NA_real_, upper = NA_real_)
  ppv_ci <- if (total_test_pos > 0) wilson_ci(tp, total_test_pos) else tibble(estimate = NA_real_, lower = NA_real_, upper = NA_real_)
  npv_ci <- if (total_test_neg > 0) wilson_ci(tn, total_test_neg) else tibble(estimate = NA_real_, lower = NA_real_, upper = NA_real_)

  tibble(
    metric = c("Sensitivity", "Specificity", "PPV", "NPV"),
    estimate = c(sens_ci$estimate[1], spec_ci$estimate[1], ppv_ci$estimate[1], npv_ci$estimate[1]),
    lower = c(sens_ci$lower[1], spec_ci$lower[1], ppv_ci$lower[1], npv_ci$lower[1]),
    upper = c(sens_ci$upper[1], spec_ci$upper[1], ppv_ci$upper[1], npv_ci$upper[1]),
    numerator = c(tp, tn, tp, tn),
    denominator = c(total_pos, total_neg, total_test_pos, total_test_neg)
  )
}

fit_response_model <- function(data) {
  glm(response ~ biomarker_call + tmb, data = data, family = binomial())
}

compute_roc_pr <- function(data, model) {
  data <- data %>% mutate(score = predict(model, newdata = data, type = "link"))
  roc_obj <- pROC::roc(response ~ score, data = data, quiet = TRUE, direction = "<")
  roc_df <- tibble(
    threshold = roc_obj$thresholds,
    sensitivity = roc_obj$sensitivities,
    specificity = roc_obj$specificities
  )

  positives <- sum(data$response == 1)
  negatives <- sum(data$response == 0)
  ordered <- data %>% arrange(desc(score))
  tp <- cumsum(ordered$response == 1)
  fp <- cumsum(ordered$response == 0)
  precision <- tp / pmax(tp + fp, 1)
  recall <- tp / pmax(positives, 1)
  pr_df <- tibble(
    threshold_rank = seq_len(nrow(ordered)),
    precision = precision,
    recall = recall
  ) %>%
    bind_rows(tibble(
      threshold_rank = 0,
      precision = sum(data$response == 1) / nrow(data),
      recall = 0
    ), .) %>%
    arrange(threshold_rank)

  list(data = data, roc_df = roc_df, pr_df = pr_df, roc_obj = roc_obj)
}

bootstrap_auc_ppv <- function(data, model, n_boot = 500) {
  boot_fun <- function(dat, indices) {
    sampled <- dat[indices, ]
    fit <- glm(response ~ biomarker_call + tmb, data = sampled, family = binomial())
    score <- predict(fit, newdata = sampled, type = "link")
    roc_val <- try(pROC::roc(sampled$response, score, quiet = TRUE, direction = "<"), silent = TRUE)
    auc_val <- if (inherits(roc_val, "try-error")) NA_real_ else as.numeric(roc_val$auc)

    tp <- sum(sampled$biomarker_call == 1 & sampled$response == 1)
    fp <- sum(sampled$biomarker_call == 1 & sampled$response == 0)
    ppv_val <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_

    c(auc_val, ppv_val)
  }

  boot::boot(data, statistic = boot_fun, R = n_boot)
}

subgroup_forest <- function(data) {
  tmb_cut <- median(data$tmb)
  data <- data %>% mutate(tmb_group = if_else(tmb >= tmb_cut, "High TMB", "Low TMB"))

  tumor_res <- data %>%
    group_by(tumor_type) %>%
    group_modify(function(.x, .y) {
      subgroup_or(.x, paste("Tumor:", .y$tumor_type))
    })

  tmb_res <- data %>%
    group_by(tmb_group) %>%
    group_modify(function(.x, .y) {
      subgroup_or(.x, .y$tmb_group)
    })

  bind_rows(tumor_res, tmb_res)
}

subgroup_or <- function(df, label) {
  if (length(unique(df$biomarker_call)) < 2 || length(unique(df$response)) < 2) {
    return(tibble(subgroup = label, n = nrow(df), or = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  fit <- glm(response ~ biomarker_call, data = df, family = binomial())
  tidy_fit <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  biomarker_row <- tidy_fit %>% filter(term == "biomarker_call")
  tibble(
    subgroup = label,
    n = nrow(df),
    or = biomarker_row$estimate,
    lower = biomarker_row$conf.low,
    upper = biomarker_row$conf.high
  )
}

km_curves <- function(data) {
  surv_obj <- survival::Surv(time = data$pfs_time, event = data$pfs_event)
  fit <- survival::survfit(surv_obj ~ biomarker_call, data = data)
  tidy_fit <- broom::tidy(fit) %>%
    mutate(biomarker_call = if_else(str_detect(strata, "biomarker_call=1"), 1L, 0L))
  tidy_fit
}

run_clinical_validation <- function(output_dir = "results/clinical", n_boot = 500) {
  ensure_dir(file.path(output_dir, "dummy"))
  data <- load_clinical_data()
  metrics <- calc_classification_metrics(data)
  model <- fit_response_model(data)
  roc_pr <- compute_roc_pr(data, model)
  boot_res <- bootstrap_auc_ppv(roc_pr$data, model, n_boot = n_boot)
  boot_df <- as_tibble(boot_res$t, .name_repair = "minimal")
  names(boot_df) <- c("auc", "ppv")

  auc_ci <- quantile(boot_df$auc, probs = c(0.025, 0.975), na.rm = TRUE)
  ppv_ci <- quantile(boot_df$ppv, probs = c(0.025, 0.975), na.rm = TRUE)

  bootstrap_summary <- tibble(
    metric = c("AUC", "PPV"),
    estimate = c(as.numeric(roc_pr$roc_obj$auc), metrics$estimate[metrics$metric == "PPV"]),
    lower = c(auc_ci[[1]], ppv_ci[[1]]),
    upper = c(auc_ci[[2]], ppv_ci[[2]])
  )

  forest <- subgroup_forest(roc_pr$data)
  km <- km_curves(data)

  readr::write_csv(metrics, file.path(output_dir, "classification_metrics.csv"))
  readr::write_csv(roc_pr$roc_df, file.path(output_dir, "roc_curve.csv"))
  readr::write_csv(roc_pr$pr_df, file.path(output_dir, "pr_curve.csv"))
  readr::write_csv(forest, file.path(output_dir, "subgroup_forest.csv"))
  readr::write_csv(km, file.path(output_dir, "km_curve.csv"))
  readr::write_csv(bootstrap_summary, file.path(output_dir, "bootstrap_summary.csv"))

  list(
    metrics = metrics,
    roc = roc_pr,
    bootstrap = bootstrap_summary,
    forest = forest,
    km = km
  )
}

if (sys.nframe() == 0) {
  set_analysis_seed(42)
  invisible(run_clinical_validation())
}
