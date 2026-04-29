# Plot hazards over population strata

# Init --------------------------------------------------------------------

here::i_am('src/62-plot_hazards_by_stratum.R'); setwd(here::here())

library(qs2)
library(tidyverse)

paths <- list()
paths$input <- list(
  competing_risk_model_fits.qs = 'tmp/50-competing_risks_model_fits.qs',
  figure_specs.R = 'src/00-figure_specifications.R',
  parametric_functions.R = 'src/00-fnct-parametric_survival_model.R',
  config.yaml = 'cfg/config.yaml'
)
paths$output <- list(
  hazards_by_social_strata.qs = 'out/62-hazards_by_social_strata.qs',
  hazards_by_social_strata.svg = 'out/62-hazards_by_social_strata.svg'
)

# figure specs
source(paths$input$figure_specs.R)
source(paths$input$parametric_functions.R)

config <- yaml::read_yaml(paths$input$config.yaml)

# constants
cnst <-
  list(
    gestage_brk = seq(24, 77, by = 4),
    lifetable_breaks = 24:77,
    left_truncation_gestage = 24,
    right_censoring_gestage = 77
  )

# Input -------------------------------------------------------------------

fit <- qs_read(paths$input$competing_risk_model_fits.qs)

# Plot hazards by social strata -------------------------------------------

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

# Export ------------------------------------------------------------------

qs_save(hazards_by_social_strata$data, paths$output$hazards_by_social_strata.qs)
fig_spec$ExportSVG(
  hazards_by_social_strata$plot,
  paths$output$hazards_by_social_strata.svg,
  width = fig_spec$width,
  height = fig_spec$width*1.4
)
