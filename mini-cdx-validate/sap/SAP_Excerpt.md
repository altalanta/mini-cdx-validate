# Statistical Analysis Plan Excerpt

## Analysis Populations

- **Analytical evaluation set:** All synthetic runs meeting QC acceptance criteria, including blank controls for contamination assessment.
- **Clinical full analysis set (FAS):** All subjects with valid biomarker, ORR, and PFS measurements.
- **Clinical biomarker-positive/negative subsets:** Defined by assay rules using the 5% VAF decision threshold and integrated TMB adjustment.

## Handling of Missing Data

- Analytical runs with missing detection calls are imputed as failures for LOD and accuracy summaries.
- Clinical endpoints with missing ORR or PFS data are excluded from respective analyses; sensitivity analyses are not required given deterministic data generation.
- Missing covariates (none expected in synthetic data) would be imputed using median (numeric) or mode (categorical) values.

## Decision Thresholds and Estimands

- Analytical decision threshold: 5% VAF for defining truth-positive events.
- Clinical biomarker-positive call: Logistic predictor combining variant status and TMB with deterministic coefficients.
- Primary clinical estimand: Odds ratio for ORR comparing biomarker-positive vs biomarker-negative subjects in the FAS.

## Confidence Interval Methods

- Binomial metrics (sensitivity, specificity, PPA, NPA, PPV, NPV): Wilson score intervals.
- LOD estimates: Bootstrap percentile intervals (B = 300).
- Odds ratios: Wald intervals based on logistic regression.
- PFS: Greenwood variance for Kaplan–Meier curves.

## Multiplicity

- No formal multiplicity adjustment is planned; subgroup estimates are descriptive and interpreted in the context of the primary estimand.

## Resampling Strategy

- Bootstrap with 500 replicates for ROC AUC and PPV confidence intervals.
- Stratified resampling not required because biomarker prevalence is fixed deterministically at ~20%.
