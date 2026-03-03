# A latent-competing risks model of the
# feto-infant mortality trajectory over
# age of gestation

# Init ------------------------------------------------------------

set.seed(1987)

library(tidyverse)

paths <- list()
paths$input <- list(
  figure_specs = 'src/00-figure_specifications.R',
  fetoinfant_lifetables = 'out/30-fetoinfant_lifetables.RData',
  lifetable_functions = 'src/00-fnct-feto_infant_lt.R',
  parametric_functions = 'src/00-fnct-parametric_survival_model.R',
  config = 'cfg/config.yaml'
)
paths$output <- list(
  competing_risk_model_fits = 'tmp/50-competing_risks_model_fits.rds'
)

# figure specs
source(paths$input$figure_specs)
# fetoinfant lifetable functions
source(paths$input$lifetable_functions)
# fetoinfant parametric functions
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

# model fits
fit <- list()

# Data ------------------------------------------------------------

# a list of feto-infant lifetable as FILT objects
load(paths$input$fetoinfant_lifetables)

# Fit feto-infant survival by social strata -----------------------

# total
fit$total14 <-
  filt$total14 |> FitFetoinfantSurvival(
    control = ControlFitFetoinfantSurvival(method = 'deoptim')
  )
PlotHazards(fit$total14, components = FALSE, colors = 'black')

# by sex
fit$sex <-
  filt$sex14 |>
  filter(stratum != 'Unknown') |>
  FitFetoinfantSurvival(
    control = ControlFitFetoinfantSurvival(method = 'deoptim')
  )
PlotHazards(fit$sex)

# by cohort
fit$cohort <-
  filt$cohort |>
  FitFetoinfantSurvival(
    control = ControlFitFetoinfantSurvival(method = 'deoptim'))
PlotHazards(fit$cohort)

# by origin
fit$origin <-
  filt$origin14 |>
  filter(
    stratum %in% c('Hispanic', 'Non-Hispanic White', 'Non-Hispanic Black')
  ) |> 
  FitFetoinfantSurvival(
    control = ControlFitFetoinfantSurvival(method = 'deoptim')
  )
PlotHazards(fit$origin)

fit$education <-
  filt$education14 |>
  filter(stratum != 'Unknown') |>
  FitFetoinfantSurvival(
    control = ControlFitFetoinfantSurvival(method = 'deoptim')
  )
PlotHazards(fit$education)

# Fit feto-infant survival by cause of death ----------------------

# by cause of death
fit$pcml_complications <-
  FitFetoinfantSurvival(
    filt$pcml_complications,
    control = ControlFitFetoinfantSurvival(
      method = 'deoptim',
      DEoptim_control = DEoptim.control(strategy = 3)
    )
  )
PlotHazards(fit$pcml_complications)

fit$congenital_malformations <-
  FitFetoinfantSurvival(
    filt$congenital_malformations,
    control = ControlFitFetoinfantSurvival(method = 'deoptim')
  )
PlotHazards(fit$congenital_malformations)

fit$maternal_complications <-
  FitFetoinfantSurvival(
    filt$maternal_complications,
    control = ControlFitFetoinfantSurvival(method = 'deoptim')
  )
PlotHazards(fit$maternal_complications)

fit$infections_parasites_toxins <-
  FitFetoinfantSurvival(
    filt$infections_parasites_toxins,
    control = ControlFitFetoinfantSurvival(
      method = 'deoptim', hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$infections_parasites_toxins)

# we fit a restricted model without birth hump
fit$prematurity <-
  FitFetoinfantSurvival(
    filt$prematurity,
    control = ControlFitFetoinfantSurvival(
      method = 'optim', simulate = TRUE, hessian_inverse = 'choleskypivot',
      model = 'flexible1',  lambda1 = 1e6, lambda2 = 10,
      # exclude birth hump parameters as they don't contribute to fit
      exclude_from_hessian_inverse = c(5, 7, 8)
    )
  )
PlotHazards(fit$prematurity)

# we fit a restricted model without birth hump
fit$accidents_and_violence <-
  FitFetoinfantSurvival(
    filt$accidents_and_violence,
    control = ControlFitFetoinfantSurvival(
      method = 'optim', simulate = TRUE,
      model = 'flexible2',
      hessian_inverse = 'choleskypivot',
      # exclude birth hump parameters as they don't contribute to fit
      exclude_from_hessian_inverse = c(5, 7, 8)
    )
  )
PlotHazards(fit$accidents_and_violence)

# we fit a restricted model without birth hump
fit$sids <-
  FitFetoinfantSurvival(
    filt$sids,
    control = ControlFitFetoinfantSurvival(
      simulate = TRUE,
      model = 'flexible2',
      hessian_inverse = 'choleskypivot',
      # exclude birth hump parameters as they don't contribute to fit
      exclude_from_hessian_inverse = c(5, 7, 8)
    )
  )
PlotHazards(fit$sids)

fit$unspecific_stillbirth <-
  FitFetoinfantSurvival(
    filt$unspecific_stillbirth,
    control = ControlFitFetoinfantSurvival(
      method = 'deoptim', model = 'basic'
    )
  )
PlotHazards(fit$unspecific_stillbirth)

fit$other <-
  FitFetoinfantSurvival(
    filt$other,
    control = ControlFitFetoinfantSurvival(
      method = 'deoptim', model = 'basic'
    )
  )
PlotHazards(fit$other)

# Export ----------------------------------------------------------

saveRDS(fit, paths$output$competing_risk_model_fits)
