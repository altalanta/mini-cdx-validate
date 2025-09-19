# mini-cdx-validate

Compact, reproducible showcase of analytical and clinical validation workflows for a synthetic next-generation sequencing (NGS) companion diagnostic. All data are deterministically simulated with `set.seed(42)` to stay FDA-aligned while remaining reviewable in a small footprint.

## Quick Start

1. Install core R packages:
   ```r
   install.packages(c("tidyverse", "broom", "pROC", "boot", "rmarkdown", "testthat", "survival", "markdown"))
   ```
2. Generate data, run analyses, and render reports:
   ```sh
   make all
   ```

Artifacts appear in `results/` and `reports/`. Rendered HTML summaries: `reports/CDx_Analytical_Validation.html`, `reports/CDx_Clinical_Validation.html`.

## Repository Map

| Component | Description |
|-----------|-------------|
| `R/generate_synthetic_data.R` | Creates analytical and clinical synthetic datasets (dilution series, truth table, cohort). |
| `R/analytical_validation.R` | LOD (probit), precision, accuracy, contamination summaries with bootstrap CIs. |
| `R/clinical_validation.R` | Clinical metrics, ROC/PR, subgroup forest plot, KM curves, bootstrap inference. |
| `R/sample_size_calc.R` | Wilson-based binomial power exploration for sensitivity targets. |
| `R/figures.R` | ggplot helpers for LOD, precision, ROC/PR, forest, and KM visuals. |
| `sas/clinical_metrics.sas` | PROC FREQ/LOGISTIC reproduction of key clinical metrics. |
| `reports/*.Rmd` | Executive-style analytical and clinical validation reports (HTML outputs in `reports/`). |
| `sap/*.md` | Protocol, SAP, and data dictionary excerpts. |
| `regulatory/*.md` | Mock Q-Sub package and PMA statistical summary. |
| `tests/` | `testthat` regression tests spanning analytical and clinical checks. |
| `Makefile` | Orchestrates data generation, analyses, reporting, and cleanup. |
| `.github/workflows/ci.yml` | CI pipeline running R CMD check, tests, and report renders. |

## Limitations

- All results rely on synthetic data; no patient data are included.
- Statistical choices (e.g., Wilson intervals, bootstrap percents) are illustrative; adjust to match sponsor SOPs before real submissions.
- SAS script provides corroboration of key clinical metrics but is not wired into Make targets.

## Useful Follow-Ups

- Review the rendered HTML reports and compare against the regulatory summaries.
- Inspect `results/` CSV outputs for traceability of reported figures.
- Run `make clean` to remove generated outputs if you need a fresh build.
