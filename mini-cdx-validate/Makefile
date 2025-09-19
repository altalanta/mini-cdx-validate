.PHONY: all data analytical clinical clean

all: data analytical clinical

data:
	Rscript R/generate_synthetic_data.R

analytical: data
	Rscript R/analytical_validation.R
	Rscript R/render_report.R reports/CDx_Analytical_Validation.Rmd reports/CDx_Analytical_Validation.html

clinical: data
	Rscript R/clinical_validation.R
	Rscript R/render_report.R reports/CDx_Clinical_Validation.Rmd reports/CDx_Clinical_Validation.html

clean:
	rm -rf results/*
	rm -f reports/*.html
	rm -f reports/figures/*.png
