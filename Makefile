R := Rscript --vanilla
SRC_DIR := src
CFG_DIR := cfg
DAT_DIR := dat
TMP_DIR := tmp
OUT_DIR := out

RAW_INFANT_DIR := $(DAT_DIR)/10-nchs-us_cohort_linked_infant_deaths_births
RAW_FETUS_DIR := $(DAT_DIR)/10-nchs-us_fetal_deaths
COD_TABLE := $(DAT_DIR)/10-cod-list/cod.csv

DOWNLOAD_STAMP := $(TMP_DIR)/.downloaded-nchs-data.stamp

TMP_OUTPUTS := \
	$(TMP_DIR)/20-fetoinfant.qs \
	$(TMP_DIR)/21-fetoinfant.qs \
	$(TMP_DIR)/50-competing_risks_model_fits.qs \
	$(TMP_DIR)/51-lifetables_and_parametric_fits_by_cod.qs \
	$(TMP_DIR)/52-competing_risk_model_parameter_tables.qs \
	$(TMP_DIR)/total_vs_cod_hazard_alignment_check.svg \
	$(DOWNLOAD_STAMP)

OUT_OUTPUTS := \
	$(OUT_DIR)/30-fetoinfant_lifetables.qs \
	$(OUT_DIR)/31-descriptives.csv \
	$(OUT_DIR)/40-lifetable_decomposition.qs \
	$(OUT_DIR)/40-gestational_age_pattern.svg \
	$(OUT_DIR)/40-perinatal_popdynamics.svg \
	$(OUT_DIR)/52-partab_strata.csv \
	$(OUT_DIR)/53-competing_risks_statistics.qs \
	$(OUT_DIR)/53-bhump_by_cod.csv \
	$(OUT_DIR)/54-parametric_decompositions.qs \
	$(OUT_DIR)/55-backwards_extrapolation.svg \
	$(OUT_DIR)/56-hzrd-origineducation.svg \
	$(OUT_DIR)/56-compression.svg \
	$(OUT_DIR)/60-birthhump_cod_joint.qs \
	$(OUT_DIR)/60-birthhump_cod_joint.svg \
	$(OUT_DIR)/60-birthhump_cod_separate.qs \
	$(OUT_DIR)/60-birthhump_cod_separate.svg \
	$(OUT_DIR)/61-overall_hazard_and_fit.qs \
	$(OUT_DIR)/61-overall_hazard_and_fit.svg \
	$(OUT_DIR)/62-hazards_by_social_strata.qs \
	$(OUT_DIR)/62-hazards_by_social_strata.svg \
	$(OUT_DIR)/63-hazards_by_cod.qs \
	$(OUT_DIR)/63-hazards_by_cod.svg

.PHONY: all download clean help

all: $(OUT_OUTPUTS)

help:
	@printf '%s\n' \
		'Targets:' \
		'  make / make all   Build all analysis outputs from prepared raw data.' \
		'  make download     Run the optional raw-data download script.' \
		'  make clean        Remove generated tmp/ and out/ artifacts from this pipeline.' \
		'  make <artifact>   Build a specific target, e.g. out/63-hazards_by_cod.svg.'

download: $(DOWNLOAD_STAMP)

clean:
	rm -f $(TMP_OUTPUTS) $(OUT_OUTPUTS)

$(TMP_DIR) $(OUT_DIR):
	mkdir -p $@

$(DOWNLOAD_STAMP): \
	$(SRC_DIR)/10-download_nchs_data_on_us_births_fetal_and_infant_deaths.R \
	$(CFG_DIR)/codebook-us_cohort_infant_births_deaths_minimal.yaml \
	$(CFG_DIR)/codebook-us_fetal_deaths_minimal.yaml | $(TMP_DIR)
	$(R) $<
	touch $@

$(TMP_DIR)/20-fetoinfant.qs: \
	$(SRC_DIR)/20-prepare_us_fetoinfants.R \
	$(SRC_DIR)/00-codebook.R \
	$(CFG_DIR)/codebook-us_cohort_infant_births_deaths_minimal.yaml \
	$(CFG_DIR)/codebook-us_fetal_deaths_minimal.yaml \
	$(RAW_INFANT_DIR) \
	$(RAW_FETUS_DIR) | $(TMP_DIR)
	$(R) $<

$(TMP_DIR)/21-fetoinfant.qs: \
	$(SRC_DIR)/21-derive_research_microdata.R \
	$(TMP_DIR)/20-fetoinfant.qs \
	$(COD_TABLE) | $(TMP_DIR)
	$(R) $<

$(OUT_DIR)/30-fetoinfant_lifetables.qs: \
	$(SRC_DIR)/30-fetoinfant_life_table_aggregation.R \
	$(TMP_DIR)/21-fetoinfant.qs \
	$(SRC_DIR)/00-fnct-feto_infant_lt.R \
	$(CFG_DIR)/config.yaml | $(OUT_DIR)
	$(R) $<
	
$(OUT_DIR)/31-descriptives.csv: \
	$(SRC_DIR)/31-descriptive_tables.R \
	$(OUT_DIR)/30-fetoinfant_lifetables.qs \
	$(SRC_DIR)/00-fnct-feto_infant_lt.R \
	$(CFG_DIR)/config.yaml | $(OUT_DIR)
	$(R) $<

$(OUT_DIR)/40-lifetable_decomposition.qs \
$(OUT_DIR)/40-gestational_age_pattern.svg \
$(OUT_DIR)/40-perinatal_popdynamics.svg &: \
	$(SRC_DIR)/40-life_table_decomposition_analysis.R \
	$(OUT_DIR)/30-fetoinfant_lifetables.qs \
	$(SRC_DIR)/00-figure_specifications.R \
	$(SRC_DIR)/00-fnct-feto_infant_lt.R | $(OUT_DIR)
	$(R) $<

$(TMP_DIR)/50-competing_risks_model_fits.qs: \
	$(SRC_DIR)/50-fit_competing_risks_model.R \
	$(OUT_DIR)/30-fetoinfant_lifetables.qs \
	$(SRC_DIR)/00-figure_specifications.R \
	$(SRC_DIR)/00-fnct-feto_infant_lt.R \
	$(SRC_DIR)/00-fnct-parametric_survival_model.R \
	$(CFG_DIR)/config.yaml | $(TMP_DIR)
	$(R) $<

$(TMP_DIR)/51-lifetables_and_parametric_fits_by_cod.qs \
$(TMP_DIR)/total_vs_cod_hazard_alignment_check.svg &: \
	$(SRC_DIR)/51-combine_lifetables_and_parametric_fits_by_cod.R \
	$(TMP_DIR)/50-competing_risks_model_fits.qs \
	$(SRC_DIR)/00-figure_specifications.R \
	$(CFG_DIR)/config.yaml | $(TMP_DIR)
	$(R) $<

$(TMP_DIR)/52-competing_risk_model_parameter_tables.qs \
$(OUT_DIR)/52-partab_strata.csv &: \
	$(SRC_DIR)/52-parameter_tables.R \
	$(TMP_DIR)/50-competing_risks_model_fits.qs \
	$(CFG_DIR)/config.yaml | $(TMP_DIR) $(OUT_DIR)
	$(R) $<

$(OUT_DIR)/53-competing_risks_statistics.qs \
$(OUT_DIR)/53-bhump_by_cod.csv &: \
	$(SRC_DIR)/53-competing_risks_inference.R \
	$(TMP_DIR)/50-competing_risks_model_fits.qs \
	$(SRC_DIR)/00-fnct-feto_infant_lt.R \
	$(SRC_DIR)/00-fnct-parametric_survival_model.R \
	$(CFG_DIR)/config.yaml | $(OUT_DIR)
	$(R) $<

$(OUT_DIR)/54-parametric_decompositions.qs: \
	$(SRC_DIR)/54-parametric_decomposition.R \
	$(TMP_DIR)/50-competing_risks_model_fits.qs \
	$(SRC_DIR)/00-fnct-parametric_survival_model.R \
	$(CFG_DIR)/config.yaml | $(OUT_DIR)
	$(R) $<

$(OUT_DIR)/55-backwards_extrapolation.svg: \
	$(SRC_DIR)/55-backwards_extrapolation.R \
	$(TMP_DIR)/50-competing_risks_model_fits.qs \
	$(SRC_DIR)/00-fnct-parametric_survival_model.R \
	$(SRC_DIR)/00-figure_specifications.R \
	$(CFG_DIR)/config.yaml | $(OUT_DIR)
	$(R) $<

$(OUT_DIR)/56-hzrd-origineducation.svg \
$(OUT_DIR)/56-compression.svg &: \
	$(SRC_DIR)/56-mortality_compression.R \
	$(OUT_DIR)/30-fetoinfant_lifetables.qs \
	$(TMP_DIR)/50-competing_risks_model_fits.qs \
	$(SRC_DIR)/00-fnct-parametric_survival_model.R \
	$(SRC_DIR)/00-fnct-feto_infant_lt.R \
	$(SRC_DIR)/00-figure_specifications.R \
	$(CFG_DIR)/config.yaml | $(OUT_DIR)
	$(R) $<

$(OUT_DIR)/60-birthhump_cod_joint.qs \
$(OUT_DIR)/60-birthhump_cod_joint.svg \
$(OUT_DIR)/60-birthhump_cod_separate.qs \
$(OUT_DIR)/60-birthhump_cod_separate.svg &: \
	$(SRC_DIR)/60-plot_birthhump_by_cod.R \
	$(TMP_DIR)/51-lifetables_and_parametric_fits_by_cod.qs \
	$(SRC_DIR)/00-figure_specifications.R \
	$(CFG_DIR)/config.yaml | $(OUT_DIR)
	$(R) $<

$(OUT_DIR)/61-overall_hazard_and_fit.qs \
$(OUT_DIR)/61-overall_hazard_and_fit.svg &: \
	$(SRC_DIR)/61-plot_total_hazard.R \
	$(OUT_DIR)/30-fetoinfant_lifetables.qs \
	$(TMP_DIR)/50-competing_risks_model_fits.qs \
	$(SRC_DIR)/00-figure_specifications.R \
	$(CFG_DIR)/config.yaml | $(OUT_DIR)
	$(R) $<

$(OUT_DIR)/62-hazards_by_social_strata.qs \
$(OUT_DIR)/62-hazards_by_social_strata.svg &: \
	$(SRC_DIR)/62-plot_hazards_by_stratum.R \
	$(TMP_DIR)/50-competing_risks_model_fits.qs \
	$(SRC_DIR)/00-figure_specifications.R \
	$(SRC_DIR)/00-fnct-parametric_survival_model.R \
	$(CFG_DIR)/config.yaml | $(OUT_DIR)
	$(R) $<

$(OUT_DIR)/63-hazards_by_cod.qs \
$(OUT_DIR)/63-hazards_by_cod.svg &: \
	$(SRC_DIR)/63-plot_hazards_by_cod.R \
	$(TMP_DIR)/51-lifetables_and_parametric_fits_by_cod.qs \
	$(OUT_DIR)/53-competing_risks_statistics.qs \
	$(SRC_DIR)/00-figure_specifications.R \
	$(CFG_DIR)/config.yaml | $(OUT_DIR)
	$(R) $<
