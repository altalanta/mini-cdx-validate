suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

source("R/utils.R")

plot_lod_curves <- function(data, output = "reports/figures/lod_curves.png") {
  ensure_dir(output)
  curve_df <- data %>%
    filter(nominal_vaf > 0) %>%
    group_by(variant_id, gene) %>%
    group_modify(~ {
      model <- glm(detected ~ log10(nominal_vaf), data = .x, family = binomial(link = "probit"))
      grid <- tibble(nominal_vaf = seq(min(.x$nominal_vaf), max(.x$nominal_vaf), length.out = 100))
      pred <- predict(model, newdata = mutate(grid, log_nominal = log10(nominal_vaf)), type = "link")
      tibble(nominal_vaf = grid$nominal_vaf, detection_prob = pnorm(pred))
    })

  p <- ggplot() +
    geom_point(data = data, aes(x = nominal_vaf, y = detected, color = gene), alpha = 0.3, position = position_jitter(height = 0.03)) +
    geom_line(data = curve_df, aes(x = nominal_vaf, y = detection_prob, color = gene), linewidth = 1) +
    scale_x_log10() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Nominal VAF (%)", y = "Detection probability", color = "Gene", title = "Probit LOD curves") +
    theme_minimal(base_size = 12)

  ggsave(output, p, width = 7, height = 5, dpi = 300)
  output
}

plot_precision_cv <- function(precision_summary, output = "reports/figures/precision_cv.png") {
  ensure_dir(output)
  p <- precision_summary %>%
    ggplot(aes(x = nominal_vaf, y = cv_pct, color = gene)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    labs(x = "Nominal VAF (%)", y = "%CV", title = "Precision across dilution series", color = "Gene") +
    theme_minimal(base_size = 12)
  ggsave(output, p, width = 7, height = 5, dpi = 300)
  output
}

plot_roc_pr <- function(roc_df, pr_df, output_prefix = "reports/figures/clinical") {
  ensure_dir(paste0(output_prefix, "_tmp.png"))
  roc_plot <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(color = "#1b9e77", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = "1 - Specificity", y = "Sensitivity", title = "ROC Curve") +
    theme_minimal(base_size = 12)

  pr_plot <- ggplot(pr_df, aes(x = recall, y = precision)) +
    geom_line(color = "#d95f02", linewidth = 1) +
    labs(x = "Recall", y = "Precision", title = "Precision-Recall Curve") +
    theme_minimal(base_size = 12)

  roc_path <- paste0(output_prefix, "_roc.png")
  pr_path <- paste0(output_prefix, "_pr.png")
  ggsave(roc_path, roc_plot, width = 6, height = 5, dpi = 300)
  ggsave(pr_path, pr_plot, width = 6, height = 5, dpi = 300)
  c(roc = roc_path, pr = pr_path)
}

plot_subgroup_forest <- function(forest_df, output = "reports/figures/subgroup_forest.png") {
  ensure_dir(output)
  df <- forest_df %>% mutate(subgroup = forcats::fct_rev(forcats::as_factor(subgroup)))
  p <- ggplot(df, aes(x = subgroup, y = or, ymin = lower, ymax = upper)) +
    geom_pointrange(color = "#7570b3") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    coord_flip() +
    labs(x = NULL, y = "Odds ratio (biomarker positive vs negative)", title = "ORR subgroup odds ratios") +
    theme_minimal(base_size = 12)
  ggsave(output, p, width = 7, height = 5, dpi = 300)
  output
}

plot_km_curve <- function(km_df, output = "reports/figures/km_curve.png") {
  ensure_dir(output)
  p <- ggplot(km_df, aes(x = time, y = estimate, color = factor(biomarker_call))) +
    geom_step() +
    labs(x = "Months", y = "PFS probability", color = "Biomarker", title = "Kaplan-Meier PFS by biomarker call") +
    theme_minimal(base_size = 12)
  ggsave(output, p, width = 7, height = 5, dpi = 300)
  output
}
