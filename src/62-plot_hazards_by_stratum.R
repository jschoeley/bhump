# Plot hazards over population strata

# Init ------------------------------------------------------------

library(tidyverse)

paths <- list()
paths$input <- list(
  competing_risk_model_fits = 'tmp/50-competing_risks_model_fits.rds',
  figure_specs = 'src/00-figure_specifications.R',
  parametric_functions = 'src/00-fnct-parametric_survival_model.R',
  config = 'cfg/config.yaml'
)
paths$output <- list(
  figures = 'out'
)

# figure specs
source(paths$input$figure_specs)
source(paths$input$parametric_functions)

config <- yaml::read_yaml(paths$input$config)

# constants
cnst <-
  list(
    gestage_brk = seq(24, 77, by = 4),
    lifetable_breaks = 24:77,
    left_truncation_gestage = 24,
    right_censoring_gestage = 77
  )

# Load data -------------------------------------------------------

fit <- readRDS(paths$input$competing_risk_model_fits)

# Plot hazards by social strata -----------------------------------

hazards_by_social_strata <- list()

hazards_by_social_strata$data <- list(
  sex = fit$sex |> select(stratum, pred_summary, lifetable),
  cohort = fit$cohort |> select(stratum, pred_summary, lifetable),
  origin = fit$origin |> select(stratum, pred_summary, lifetable),
  education = fit$education |> select(stratum, pred_summary, lifetable)
)

# by sex
hzrd_sex <- PlotHazards(
  hazards_by_social_strata$data$sex,
  ylim_hx = c(0.6, 80), ylim_Fx = c(0, 1500), ar = 0.7,
  notitle = TRUE
)

# by cohort
hzrd_cohort <- PlotHazards(
  hazards_by_social_strata$data$cohort,
  ylim_hx = c(0.6, 80), ylim_Fx = c(0, 1500), ar = 0.7,
  notitle = TRUE
)

# by origin
hzrd_origin <- PlotHazards(
  hazards_by_social_strata$data$origin,
  ylim_hx = c(0.6, 80), ylim_Fx = c(0, 1500), ar = 0.7,
  notitle = TRUE
)

# by education
hzrd_education <- PlotHazards(
  hazards_by_social_strata$data$education
  |> filter(stratum != 'Unknown'),
  ylim_hx = c(0.6, 80), ylim_Fx = c(0, 1500), notitle = TRUE, ar = 0.7,
)

hazards_by_social_strata$plot <- cowplot::plot_grid(
  hzrd_sex, hzrd_origin, hzrd_cohort, hzrd_education,
  nrow = 4, align = 'hv', labels = 'AUTO', axis = 'l'
)
hazards_by_social_strata$plot

# Export ----------------------------------------------------------

saveRDS(hazards_by_social_strata$data, 'out/62-hazards_by_social_strata.rds')
fig_spec$ExportPDF(
  hazards_by_social_strata$plot,
  filename = '62-hazards_by_social_strata',
  path = 'out',
  width = fig_spec$width,
  height = fig_spec$width*1.4
)
