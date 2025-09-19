# Q-Sub Statistical Meeting Package (Mock)

## Purpose

Seek FDA feedback on analytical and clinical validation plans for the synthetic NGS companion diagnostic intended to guide targeted therapy selection in solid tumours.

## Proposed Analytical Design

- **Dilution series:** Three representative variants spanning 0–10% nominal VAF, executed across three CLIA-like sites, three runs per site, two operators, and three replicates.
- **LOD analysis:** Variant-specific probit regression on log10 VAF to estimate LOD50/LOD95 with bootstrap percentile confidence intervals (B = 300).
- **Precision:** %CV summaries by VAF level plus ANOVA-based variance decomposition (site, run, operator, residual).
- **Accuracy vs truth:** PPA/NPA calculations against a composite GIAB/consensus truth set using Wilson score confidence intervals at the 5% decision threshold.
- **Contamination:** Monitoring of unexpected detections <0.3% VAF among blank controls, summarised per run.

### Analytical Questions to FDA

1. Do the proposed LOD acceptance criteria (LOD95 ≤1% VAF) align with expectations for similar CDx submissions?
2. Is bootstrap percentile methodology acceptable for deriving LOD confidence intervals, or is a delta-method approximation preferable?
3. Are Wilson score intervals sufficient for PPA/NPA claims, or is Clopper–Pearson required despite its conservatism at moderate sample sizes?
4. Should blank-control contamination thresholds (<5% per run) be formalised as lot release criteria?

## Proposed Clinical Design

- **Cohort:** ~350 subjects with NSCLC (~45%), CRC (~35%), and other tumours (~20%); biomarker prevalence targeted at ~20%.
- **Primary endpoint:** ORR odds ratio for biomarker-positive vs biomarker-negative subjects.
- **Performance metrics:** Sensitivity, specificity, PPV, NPV with Wilson intervals; ROC AUC via logistic model including TMB; PR curve for decision support.
- **Bootstrap:** 500 replicates for AUC and PPV confidence intervals.
- **Time-to-event:** Kaplan–Meier PFS comparison with Greenwood variance estimates.
- **Subgroups:** Tumour type and TMB high/low (median split) forest plot of ORR odds ratios.

### Clinical Questions to FDA

1. Are the target performance thresholds (sensitivity ≥75%, specificity ≥70%, PPV ≥60%) sufficient for clinical utility claims in this indication?
2. Is the planned logistic model (biomarker + TMB) adequate for ROC/AUC submissions, or should additional covariates (e.g., ECOG status) be incorporated?
3. Does the proposed bootstrap (B = 500) satisfy expectations for precision of AUC and PPV estimates?
4. Are descriptive subgroup analyses acceptable without multiplicity adjustment given the exploratory intent?

## Success Criteria (Draft)

- All analytical acceptance criteria met or exceeded.
- Clinical metrics meet or exceed targets with lower 95% CI bounds above preset thresholds.
- ROC AUC ≥0.70 with lower CI ≥0.65.
- Kaplan–Meier curves show non-inferior PFS for biomarker-negative subjects with clear benefit in biomarker-positive subjects.

## Requested FDA Feedback

- Agreement on analytical acceptance criteria and statistical methods.
- Confirmation that the clinical analysis set and performance metrics support a future PMA submission.
- Input on additional analyses or stress testing (e.g., near-threshold behaviour) expected in the PMA dossier.
