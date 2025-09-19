/* Clinical metrics for synthetic cohort */
%let root = ..;

proc import datafile="&root./data/synthetic/clinical_cohort.csv"
    out=work.clinical_cohort
    dbms=csv
    replace;
    guessingrows=max;
run;

proc freq data=work.clinical_cohort;
    tables biomarker_call*response / nocol norow nopercent out=work.classification_counts;
run;

proc logistic data=work.clinical_cohort;
    class tumor_type (ref='Other') / param=ref;
    model response(event='1') = biomarker_call tmb tumor_type;
    ods output OddsRatios=work.or_table;
run;

proc print data=work.classification_counts label;
    title "Sensitivity and Specificity Components";
run;

proc print data=work.or_table label;
    title "Logistic Regression Odds Ratios";
run;
