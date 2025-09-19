#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
})

source("R/utils.R")
set_analysis_seed(42)

# Parameters for the dilution series
variants <- tribble(
  ~variant_id, ~gene,
  "VAR1", "EGFR",
  "VAR2", "KRAS",
  "VAR3", "BRAF"
)
nominal_levels <- c(0, 0.1, 0.25, 0.5, 1, 2, 5, 10)
sites <- c("SiteA", "SiteB", "SiteC")
run_numbers <- 1:3
replicates <- 1:3
operators <- c("Operator1", "Operator2")

analytical_runs <- expand_grid(
  variants,
  nominal_vaf = nominal_levels,
  site = sites,
  run_number = run_numbers,
  replicate = replicates
) %>%
  mutate(
    operator = if_else(replicate %% 2 == 1, operators[1], operators[2]),
    run_id = sprintf("%s_Run%02d", site, run_number),
    depth = pmax(200, round(rnorm(n(), mean = 800, sd = 150))),
    qc_prob = plogis(-2 + 0.7 * (nominal_vaf < 0.5) + 1.2 * (nominal_vaf == 0)),
    qc_flag = rbinom(n(), 1, qc_prob)
  )

# Log-normal noise around the nominal value keeps CV stable on the log scale
noise <- rnorm(nrow(analytical_runs), mean = 0, sd = 0.2)
measured_vaf <- if_else(
  analytical_runs$nominal_vaf == 0,
  exp(rnorm(nrow(analytical_runs), log(0.02 + runif(nrow(analytical_runs), 0, 0.02)), 0.25)),
  analytical_runs$nominal_vaf * exp(noise)
)

# Inject a light contamination signal among blank wells
is_blank <- analytical_runs$nominal_vaf == 0
contam_index <- sample(which(is_blank), size = ceiling(0.07 * sum(is_blank)))
measured_vaf[contam_index] <- runif(length(contam_index), min = 0.03, max = 0.25)

# Detection probability depends on VAF, depth, and QC status
log_vaf <- if_else(analytical_runs$nominal_vaf > 0, log10(analytical_runs$nominal_vaf), -3)
logit_p <- -2.2 + 1.6 * log_vaf + 0.001 * (analytical_runs$depth - 600) - 1.1 * analytical_runs$qc_flag
prob_detect <- inv_logit(logit_p)

analytical_runs <- analytical_runs %>%
  mutate(
    measured_vaf = pmax(measured_vaf, 0.01),
    detected = rbinom(n(), 1, pmin(pmax(prob_detect, 0.001), 0.999)),
    sample_id = sprintf("S%05d", row_number()),
    qc_flag = as.integer(qc_flag)
  ) %>%
  select(sample_id, variant_id, gene, nominal_vaf, measured_vaf, depth,
         site, run_id, operator, replicate, detected, qc_flag)

# Truth reference links sample and variant to expected VAF and panel
truth_reference <- analytical_runs %>%
  distinct(sample_id, variant_id, gene, nominal_vaf) %>%
  mutate(
    truth_present = as.integer(nominal_vaf > 0),
    truth_vaf = nominal_vaf,
    reference_panel = if_else(variant_id == "VAR1", "GIAB_mock", "consensus")
  ) %>%
  select(sample_id, variant_id, gene, truth_present, truth_vaf, reference_panel)

# Clinical cohort of ~350 subjects with 20% biomarker prevalence
n_subjects <- 350
subject_df <- tibble(
  subject_id = sprintf("SUBJ%03d", seq_len(n_subjects)),
  tumor_type = sample(c("NSCLC", "CRC", "Other"), n_subjects, replace = TRUE, prob = c(0.45, 0.35, 0.20)),
  age = pmin(89, pmax(25, round(rnorm(n_subjects, mean = 63, sd = 9)))),
  sex = sample(c("F", "M"), n_subjects, replace = TRUE)
)

# TMB distribution skewed right; tie biomarker prevalence to higher TMB
subject_df <- subject_df %>%
  mutate(
    tmb = round(rgamma(n_subjects, shape = 6, scale = 2), 1),
    tmb_centered = tmb - median(tmb),
    biomarker_prob = plogis(-1.4 + 0.06 * tmb_centered + 0.3 * (tumor_type == "NSCLC")),
    biomarker_call = rbinom(n_subjects, 1, biomarker_prob)
  )

# Clinical response model assumes higher odds when biomarker positive
response_lp <- -1.1 + 1.0 * subject_df$biomarker_call + 0.02 * subject_df$tmb_centered + 0.3 * (subject_df$tumor_type == "NSCLC")
response <- rbinom(n_subjects, 1, inv_logit(response_lp))

# Progression-free survival with modest benefit for biomarker positives
base_hazard <- 0.12
hazard <- base_hazard * exp(-0.35 * subject_df$biomarker_call)
pfs_time <- rexp(n_subjects, rate = hazard)

# Administrative censoring at ~18 months
censor_time <- runif(n_subjects, min = 12, max = 20)
observed_time <- pmin(pfs_time, censor_time)
pfs_event <- as.integer(pfs_time <= censor_time)

clinical_cohort <- subject_df %>%
  mutate(
    response = response,
    pfs_time = round(observed_time, 2),
    pfs_event = pfs_event
  ) %>%
  select(subject_id, tumor_type, age, sex, tmb, biomarker_call, response, pfs_time, pfs_event)

# Write outputs
analytical_path <- ensure_dir("data/synthetic/analytical_runs.csv")
truth_path <- ensure_dir("data/synthetic/truth_reference.csv")
clinical_path <- ensure_dir("data/synthetic/clinical_cohort.csv")

readr::write_csv(analytical_runs, analytical_path)
readr::write_csv(truth_reference, truth_path)
readr::write_csv(clinical_cohort, clinical_path)

message("Synthetic data written to data/synthetic/")
