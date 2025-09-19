# Utility functions shared across analytical and clinical scripts

#' Ensure deterministic behavior across scripts
#' @param seed integer seed value
set_analysis_seed <- function(seed = 42) {
  set.seed(seed)
  invisible(seed)
}

#' Logistic transform helpers
logit <- function(p) {
  stopifnot(is.numeric(p))
  log(p / (1 - p))
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

#' Wilson score confidence interval for a binomial proportion
#' @param success number of successes
#' @param total number of trials
#' @param conf_level confidence level (default 0.95)
#' @return tibble with estimate, lower, upper
wilson_ci <- function(success, total, conf_level = 0.95) {
  if (any(total <= 0)) {
    stop("Total trials must be > 0")
  }
  z <- qnorm(1 - (1 - conf_level) / 2)
  p_hat <- success / total
  denom <- 1 + z^2 / total
  center <- (p_hat + z^2 / (2 * total)) / denom
  half_width <- z * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * total)) / total) / denom
  tibble::tibble(
    estimate = p_hat,
    lower = pmax(0, center - half_width),
    upper = pmin(1, center + half_width)
  )
}

#' Convenience wrapper to format a confidence interval
format_ci <- function(lower, upper, digits = 3) {
  paste0("[", formatC(lower, format = "f", digits = digits), ", ",
         formatC(upper, format = "f", digits = digits), "]")
}

#' Simple helper to ensure directories exist before writing
ensure_dir <- function(path) {
  dir <- dirname(path)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  path
}
