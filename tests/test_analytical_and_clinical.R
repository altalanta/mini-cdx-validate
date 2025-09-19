pkg_root <- normalizePath(file.path(".."))
old_wd <- setwd(pkg_root)
on.exit(setwd(old_wd), add = TRUE)

set_analysis_seed(42)

if (!file.exists("data/synthetic/analytical_runs.csv") ||
    !file.exists("data/synthetic/clinical_cohort.csv")) {
  source("R/generate_synthetic_data.R")
}

analytical_results <- run_analytical_validation(n_boot = 120)
clinical_results <- run_clinical_validation(n_boot = 200)

lod_summary <- analytical_results$lod
precision_df <- analytical_results$precision$precision
accuracy_df <- analytical_results$accuracy
bootstrap_df <- clinical_results$bootstrap
roc_auc <- as.numeric(clinical_results$roc$roc_obj$auc)

sens_est <- clinical_results$metrics$estimate[clinical_results$metrics$metric == "Sensitivity"]
spec_est <- clinical_results$metrics$estimate[clinical_results$metrics$metric == "Specificity"]

prevalence_grid <- seq(0.1, 0.5, by = 0.1)
ppv_curve <- (sens_est * prevalence_grid) /
  (sens_est * prevalence_grid + (1 - spec_est) * (1 - prevalence_grid))

if (length(sens_est) == 0 || length(spec_est) == 0) {
  stop("Clinical metrics failed to compute sensitivity/specificity estimates")
}

test_that("LOD95 exceeds LOD50 for each variant", {
  expect_true(all(lod_summary$`estimate_LOD95` > lod_summary$`estimate_LOD50`))
})

test_that("Precision improves at higher VAF", {
  low_cv <- stats::median(precision_df$cv_pct[precision_df$nominal_vaf <= 0.5], na.rm = TRUE)
  high_cv <- stats::median(precision_df$cv_pct[precision_df$nominal_vaf >= 5], na.rm = TRUE)
  expect_lt(high_cv, low_cv)
})

test_that("Wilson confidence intervals remain within [0,1]", {
  expect_true(all(accuracy_df$lower >= 0 & accuracy_df$lower <= 1))
  expect_true(all(accuracy_df$upper >= 0 & accuracy_df$upper <= 1))
})

test_that("Bootstrap confidence intervals are numeric", {
  expect_true(all(is.finite(bootstrap_df$lower)))
  expect_true(all(is.finite(bootstrap_df$upper)))
})

test_that("ROC AUC lies between 0.5 and 1", {
  expect_gt(roc_auc, 0.5)
  expect_lt(roc_auc, 1)
})

test_that("PPV increases monotonically with prevalence", {
  expect_true(all(diff(ppv_curve) >= -1e-8))
})
