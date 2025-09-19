# Study Protocol Excerpt

## Objective

Demonstrate analytical and clinical performance of the synthetic NGS companion diagnostic (CDx) that detects actionable variants in EGFR, KRAS, and BRAF to guide targeted therapy selection.

## Analytical Endpoints

- Limit of detection (LOD50/LOD95) for each representative variant across dilution series (0.1%–10% VAF).
- Intra- and inter-run precision measured by %CV at prespecified nominal VAF levels.
- Positive/negative percent agreement (PPA/NPA) versus a mock Genome in a Bottle (GIAB) truth set at the 5% VAF decision threshold.
- Contamination monitoring rate defined by unexpected calls at VAF <0.3% in blank controls.

## Clinical Endpoints

- Sensitivity, specificity, PPV, and NPV of biomarker calls relative to overall response rate (ORR).
- Area under the ROC curve (AUC) for predicting ORR using biomarker status and tumour mutational burden (TMB).
- PFS comparison by biomarker status using Kaplan–Meier methodology.
- Subgroup ORR odds ratios for tumour type and TMB strata.

## Design Overview

- Analytical: Three sites, three runs per site, two operators, and three replicates per nominal VAF level for each variant, leveraging a standardised synthetic reference panel.
- Clinical: Observational dataset of ~350 oncology subjects with NSCLC, CRC, or other tumours, reflecting ~20% biomarker prevalence.
- All simulations executed with `set.seed(42)` to ensure deterministic reproducibility.

## Acceptance Criteria

- LOD95 ≤1% VAF for all monitored variants.
- %CV ≤15% at ≥5% VAF and ≤30% below 1% VAF.
- PPA and NPA ≥90% with 95% lower confidence bounds ≥85%.
- Contamination rate <5% per run in blank controls.
- Clinical sensitivity ≥75%, specificity ≥70%, and PPV ≥60% with supporting ROC AUC ≥0.70.
