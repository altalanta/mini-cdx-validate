# Data Dictionary

## analytical_runs.csv

| Column         | Type    | Description |
|----------------|---------|-------------|
| sample_id      | character | Unique synthetic sample identifier per run replicate. |
| variant_id     | character | Representative variant code (VAR1–VAR3). |
| gene           | character | Gene symbol (EGFR, KRAS, BRAF). |
| nominal_vaf    | numeric | Nominal VAF (%) used in the dilution series (0–10). |
| measured_vaf   | numeric | Observed VAF (%) with log-normal noise and contamination injections. |
| depth          | integer | Sequencing depth (reads) for the assay replicate. |
| site           | character | Testing site (SiteA–SiteC). |
| run_id         | character | Unique run identifier combining site and run number. |
| operator       | character | Operator performing the run (Operator1/Operator2). |
| replicate      | integer | Replicate number within the run (1–3). |
| detected       | integer | Detection outcome (1 = variant detected, 0 = not detected). |
| qc_flag        | integer | QC status (0 = pass, 1 = flagged). |

## truth_reference.csv

| Column          | Type    | Description |
|-----------------|---------|-------------|
| sample_id       | character | Links to analytical_runs sample_id. |
| variant_id      | character | Variant code as in analytical_runs. |
| gene            | character | Matching gene symbol. |
| truth_present   | integer | 1 if variant present in truth set (>0 VAF), 0 for blanks. |
| truth_vaf       | numeric | Truth VAF (%) used to define PPA/NPA thresholding. |
| reference_panel | character | Source of mock reference (GIAB_mock or consensus). |

## clinical_cohort.csv

| Column         | Type    | Description |
|----------------|---------|-------------|
| subject_id     | character | Unique subject identifier. |
| tumor_type     | character | Tumour type (NSCLC, CRC, Other). |
| age            | integer | Age in years (25–89). |
| sex            | character | Biological sex (F/M). |
| tmb            | numeric | Tumour mutational burden (mut/Mb). |
| biomarker_call | integer | 1 = biomarker positive, 0 = negative per CDx algorithm. |
| response       | integer | ORR status (1 = responder, 0 = non-responder). |
| pfs_time       | numeric | PFS time in months (censored at ~20 months). |
| pfs_event      | integer | Event indicator (1 = progression/death, 0 = censored). |
