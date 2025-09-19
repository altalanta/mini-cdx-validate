# PMA Statistical Summary (Mock)

## Executive Summary

The synthetic companion diagnostic met all pre-specified analytical and clinical validation criteria. Analytical performance demonstrated robust sensitivity (LOD95 ≤0.6% VAF) and precision, while clinical analyses showed consistent benefit for biomarker-positive subjects across tumour types with favourable predictive metrics.

## Analytical Validation Highlights

- **LOD:** All three sentinel variants achieved LOD95 estimates between 0.42% and 0.58% VAF, comfortably below the 1% target (see `reports/CDx_Analytical_Validation.html`, Table “LOD50 and LOD95 estimates”).
- **Precision:** %CV decreased from ~18% at 0.1% VAF to <8% at ≥5% VAF, with site/run/operator components each contributing <20% of total variance.
- **Accuracy:** PPA and NPA both exceeded 92%, with Wilson lower bounds above 0.88 (refer to analytical report accuracy table).
- **Contamination:** Blank control hit rates remained <3% per run; contamination events are tabulated in the analytical report appendix.

## Clinical Validation Highlights

- **Classification metrics:** Sensitivity 82% (95% CI: 75–88), specificity 74% (95% CI: 68–80), PPV 64% (95% CI: 58–71), NPV 88% (95% CI: 83–92). Detailed tables reside in `reports/CDx_Clinical_Validation.html`.
- **ROC / PR:** Logistic predictor (biomarker + TMB) yielded an AUC of 0.79 (bootstrap 95% CI: 0.74–0.83) and a favourable precision-recall profile at 20% prevalence.
- **Subgroups:** Forest plot shows ORR odds ratios >1 across NSCLC, CRC, Other, and both TMB strata, without any strata crossing unity.
- **PFS:** Kaplan–Meier curves confirmed longer PFS in biomarker-positive subjects (median 14.6 vs 11.2 months) with maintained safety in biomarker-negative population.

## Risk–Benefit Considerations

- Analytical evidence supports reliable detection down to 0.5% VAF with tight precision, minimising the risk of false negatives near clinical decision boundaries.
- Clinical metrics demonstrate meaningful enrichment of responders among biomarker-positive patients while maintaining acceptable specificity to avoid overtreatment.
- Survival trends corroborate ORR findings and suggest durable benefit for biomarker-positive individuals.
- Residual risk is primarily associated with near-threshold samples; ongoing monitoring will include control charting of blank-run contamination to detect drift.

## Planned Post-Approval Commitments

- Annual proficiency testing across all sites using the synthetic reference panel.
- Real-world evidence collection to monitor biomarker prevalence and confirm PPV in broader populations.
- Periodic review of near-threshold calls to ensure decision rules remain calibrated.
