suppressPackageStartupMessages({
  library(tidyverse)
})

source("R/utils.R")

required_positives <- function(true_sensitivity = 0.80, target_lower = 0.70, conf_level = 0.95, max_n = 500) {
  for (n in seq_len(max_n)) {
    successes <- round(true_sensitivity * n)
    ci <- wilson_ci(successes, n, conf_level)
    if (ci$lower >= target_lower) {
      return(n)
    }
  }
  stop("Increase max_n; requirements not met within search range")
}

sample_size_table <- function(true_sensitivity = 0.80, target_lower = 0.70, conf_level = 0.95, span = 5) {
  n_star <- required_positives(true_sensitivity, target_lower, conf_level)
  grid <- seq(pmax(5, n_star - span), n_star + span)
  tibble(n_positives = grid) %>%
    mutate(
      expected_hits = round(true_sensitivity * n_positives),
      ci = map(n_positives, ~ wilson_ci(round(true_sensitivity * .x), .x, conf_level)),
      lower = map_dbl(ci, ~ .x$lower),
      upper = map_dbl(ci, ~ .x$upper),
      meets_target = lower >= target_lower
    ) %>%
    select(-ci)
}

run_sample_size_calc <- function(output_dir = "results/sample_size") {
  ensure_dir(file.path(output_dir, "dummy"))
  table <- sample_size_table()
  readr::write_csv(table, file.path(output_dir, "sensitivity_sample_size.csv"))
  table
}

if (sys.nframe() == 0) {
  run_sample_size_calc()
}
