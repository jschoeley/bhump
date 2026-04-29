# Analyze mortality compression

# Init --------------------------------------------------------------------

here::i_am('src/56-mortality_compression.R'); setwd(here::here())

library(qs2)
library(tidyverse)

paths <- list()
paths$input <- list(
  parametric_functions.R = 'src/00-fnct-parametric_survival_model.R',
  lifetable_functions.R = 'src/00-fnct-feto_infant_lt.R',
  fetoinfant_lifetables.qs = 'out/30-fetoinfant_lifetables.qs',
  competing_risk_model_fits.qs = 'tmp/50-competing_risks_model_fits.qs',
  figure_specs.R = 'src/00-figure_specifications.R',
  config.yaml = 'cfg/config.yaml'
)
paths$output <- list(
  #overall_hazard_and_fit.qs = 'out/61-overall_hazard_and_fit.qs',
  #overall_hazard_and_fit.svg = 'out/61-overall_hazard_and_fit.svg'
  hzrd_origineducation.svg = 'out/56-hzrd-origineducation.svg',
  compression.svg = 'out/56-compression.svg'
)

# figure specs
source(paths$input$figure_specs.R)
source(paths$input$parametric_functions.R)
source(paths$input$lifetable_functions.R)

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

filt <- qs_read(paths$input$fetoinfant_lifetables.qs)

# Mortality compression towards birth -----------------------------

fit <- list()
fig <- list()

fit$origineducation <-
  filt$origineducation |>
  filter(!grepl('Unknown', stratum)) |>
  FitFetoinfantSurvival(
    control = ControlFitFetoinfantSurvival(simulate = FALSE)
  )
fig$hzrd_origineducation <- PlotHazards(
  fit$origineducation, colors = c(
    rep('black',7),
    'red', # white academic
    rep('black',1),
    'blue', # black primary
    rep('black',6)
  ), legend = FALSE)
fig$hzrd_origineducation


fig$compression <-
  tibble(
    F76 = ProbFetoInfantDeath(fit$origineducation, x = 77) |>
      filter(!grepl('Unknown', stratum)) |> pull(avg_total_Fx),
    rho = BirthHumpDeaths(fit$origineducation, x = 77) |>
      filter(!grepl('Unknown', stratum)) |> pull(avg_p_birth)
  ) |>
  ggplot(aes(x = F76*100,y = rho*100)) +
  geom_smooth(method = 'lm', se = FALSE, color = 'grey40') +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  labs(y = '% feto-infant deaths contributed by birth-hump',
       x = '% probability of feto infant death') +
  fig_spec$MyGGplotTheme()
fig$compression

# Export ------------------------------------------------------------------

fig_spec$ExportSVG(
  fig$hzrd_origineducation,
  paths$output$hzrd_origineducation.svg,
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

fig_spec$ExportSVG(
  fig$compression,
  paths$output$compression.svg,
  width = fig_spec$width*0.8,
  height = fig_spec$width*0.6
)
