# A latent-competing risks model of the
# feto-infant mortality trajectory over
# age of gestation

# Init --------------------------------------------------------------------

here::i_am('src/50-fit_competing_risks_model.R'); setwd(here::here())

set.seed(1987)

library(qs2)
library(tidyverse)

paths <- list()
paths$input <- list(
  figure_specs.R = 'src/00-figure_specifications.R',
  fetoinfant_lifetables.qs = 'out/30-fetoinfant_lifetables.qs',
  lifetable_functions.R = 'src/00-fnct-feto_infant_lt.R',
  parametric_functions.R = 'src/00-fnct-parametric_survival_model.R',
  config.yaml = 'cfg/config.yaml'
)
paths$output <- list(
  competing_risk_model_fits.qs = 'tmp/50-competing_risks_model_fits.qs'
)

# figure specs
source(paths$input$figure_specs.R)
# fetoinfant lifetable functions
source(paths$input$lifetable_functions.R)
# fetoinfant parametric functions
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

# model fits
fit <- list()

# Data --------------------------------------------------------------------

# a list of feto-infant lifetable as FILT objects
filt <- qs_read(paths$input$fetoinfant_lifetables.qs)

# Fit feto-infant survival by social strata -------------------------------

# total
fit$total14 <-
  filt$total14 |> FitFetoinfantSurvival(
    control = ControlFitFetoinfantSurvival(hessian_inverse = 'choleskypivot')
  )
PlotHazards(fit$total14, components = FALSE, colors = 'black')

# by sex
fit$sex <-
  filt$sex14 |>
  filter(stratum != 'Unknown') |>
  FitFetoinfantSurvival(
    control = ControlFitFetoinfantSurvival(hessian_inverse = 'choleskypivot')
  )
PlotHazards(fit$sex)

# by cohort
fit$cohort <-
  filt$cohort |>
  FitFetoinfantSurvival(
    control = ControlFitFetoinfantSurvival(hessian_inverse = 'choleskypivot'))
PlotHazards(fit$cohort)

# by origin
fit$origin <-
  filt$origin14 |>
  filter(
    stratum %in% c('Hispanic', 'Non-Hispanic White', 'Non-Hispanic Black')
  ) |> 
  FitFetoinfantSurvival(
    control = ControlFitFetoinfantSurvival(hessian_inverse = 'choleskypivot')
  )
PlotHazards(fit$origin)

fit$education <-
  filt$education14 |>
  filter(stratum != 'Unknown') |>
  FitFetoinfantSurvival(
    control = ControlFitFetoinfantSurvival(
      hessian_inverse = 'choleskypivot',
      zeta_range = c(36, 40)-cnst$left_truncation_gestage
    )
  )
PlotHazards(fit$education)

# Fit feto-infant survival by cause of death ------------------------------

# by cause of death
fit$pcm <-
  FitFetoinfantSurvival(
    filt$pcm,
    control = ControlFitFetoinfantSurvival(hessian_inverse = 'choleskypivot')
  )
PlotHazards(fit$pcm)

fit$labor <-
  FitFetoinfantSurvival(
    filt$labor,
    control = ControlFitFetoinfantSurvival(hessian_inverse = 'choleskypivot')
  )
PlotHazards(fit$labor)

fit$congenital <-
  FitFetoinfantSurvival(
    filt$congenital,
    control = ControlFitFetoinfantSurvival(hessian_inverse = 'choleskypivot')
  )
PlotHazards(fit$congenital)

fit$maternal <-
  FitFetoinfantSurvival(
    filt$maternal,
    control = ControlFitFetoinfantSurvival(hessian_inverse = 'choleskypivot')
  )
PlotHazards(fit$maternal)

fit$convulsions <-
  FitFetoinfantSurvival(
    filt$convulsions,
    control = ControlFitFetoinfantSurvival(
      zeta_range = c(39, 43)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$convulsions)

fit$sepsis <-
  FitFetoinfantSurvival(
    filt$sepsis,
    control = ControlFitFetoinfantSurvival(
      model = 'flexible1',
      zeta_range = c(41, 44)-cnst$left_truncation_gestage,lambda2 = 1e1,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$sepsis)

fit$hypoxia <-
  FitFetoinfantSurvival(
    filt$hypoxia,
    control = ControlFitFetoinfantSurvival(
      zeta_range = c(38, 42)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$hypoxia)

fit$respiratory <-
  FitFetoinfantSurvival(
    filt$respiratory,
    control = ControlFitFetoinfantSurvival(
      model = 'flexible1',
      zeta_range = c(37, 41)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$respiratory)

fit$prematurity <-
  FitFetoinfantSurvival(
    filt$prematurity,
    control = ControlFitFetoinfantSurvival(
      model = 'flexible1', lambda2 = 1,
      zeta_range = c(38, 40)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$prematurity)

# we fit a restricted model without birth hump
fit$sids <-
  FitFetoinfantSurvival(
    filt$sids,
    control = ControlFitFetoinfantSurvival(
      model = 'flexible2', lambda1 = 1e3, lambda2 = 1e1,
      zeta_range = c(40, 45)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot',
    )
  )
PlotHazards(fit$sids)

fit$unspecific <-
  FitFetoinfantSurvival(
    filt$unspecific,
    control = ControlFitFetoinfantSurvival(
      zeta_range = c(37, 41)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$unspecific)

fit$otherspecific <-
  FitFetoinfantSurvival(
    filt$otherspecific,
    control = ControlFitFetoinfantSurvival(
      model = 'flexible1', lambda2 = 10,
      zeta_range = c(39, 41)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$otherspecific)

fit$unknown <-
  FitFetoinfantSurvival(
    filt$unknown,
    control = ControlFitFetoinfantSurvival(
      zeta_range = c(37, 41)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$unknown)

# Export ------------------------------------------------------------------

qs_save(fit, paths$output$competing_risk_model_fits.qs)
