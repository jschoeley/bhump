# A latent-competing risks model of the
# feto-infant mortality trajectory over
# age of gestation

# Init ------------------------------------------------------------

set.seed(1987)

library(qs2)
library(tidyverse)

paths <- list()
paths$input <- list(
  figure_specs = 'src/00-figure_specifications.R',
  fetoinfant_lifetables = 'out/30-fetoinfant_lifetables.qs',
  lifetable_functions = 'src/00-fnct-feto_infant_lt.R',
  parametric_functions = 'src/00-fnct-parametric_survival_model.R',
  config = 'cfg/config.yaml'
)
paths$output <- list(
  competing_risk_model_fits = 'tmp/50-competing_risks_model_fits.qs'
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
filt <- qs_read(paths$input$fetoinfant_lifetables)

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
fit$PlacentaCordMembrane <-
  FitFetoinfantSurvival(
    filt$PlacentaCordMembrane,
    control = ControlFitFetoinfantSurvival(
      method = 'optim', hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$PlacentaCordMembrane)

fit$LabourBirth <-
  FitFetoinfantSurvival(
    filt$LabourBirth,
    control = ControlFitFetoinfantSurvival(
      method = 'optim', hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$LabourBirth)

fit$Malformations <-
  FitFetoinfantSurvival(
    filt$Malformations,
    control = ControlFitFetoinfantSurvival(
      method = 'optim', hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$Malformations)

fit$Maternal <-
  FitFetoinfantSurvival(
    filt$Maternal,
    control = ControlFitFetoinfantSurvival(
      method = 'optim', hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$Maternal)

fit$CerebralConvulsions <-
  FitFetoinfantSurvival(
    filt$CerebralConvulsions,
    control = ControlFitFetoinfantSurvival(
      method = 'optim', hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$CerebralConvulsions)

fit$Infections <-
  FitFetoinfantSurvival(
    filt$Infections,
    control = ControlFitFetoinfantSurvival(
      method = 'optim', hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$Infections)

fit$HypoxiaAsphyxie <-
  FitFetoinfantSurvival(
    filt$HypoxiaAsphyxie,
    control = ControlFitFetoinfantSurvival(
      method = 'deoptim', hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$HypoxiaAsphyxie)

fit$RespiratoryCardio <-
  FitFetoinfantSurvival(
    filt$RespiratoryCardio,
    control = ControlFitFetoinfantSurvival(
      model = 'flexible1',
      method = 'optim', hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$RespiratoryCardio)

fit$FetalGrowthPremature <-
  FitFetoinfantSurvival(
    filt$FetalGrowthPremature,
    control = ControlFitFetoinfantSurvival(
      model = 'flexible1',
      method = 'optim', hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$FetalGrowthPremature)

# we fit a restricted model without birth hump
fit$SID <-
  FitFetoinfantSurvival(
    filt$SID,
    control = ControlFitFetoinfantSurvival(
      model = 'flexible2',
      method = 'optim', hessian_inverse = 'choleskypivot',
      lambda1 = 1e3,
      # exclude birth hump parameters as they don't contribute to fit
      exclude_from_hessian_inverse = c(5, 7, 8)
    )
  )
PlotHazards(fit$SID)

fit$UnspecificStillbirth <-
  FitFetoinfantSurvival(
    filt$UnspecificStillbirth,
    control = ControlFitFetoinfantSurvival(
      method = 'optim', hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$UnspecificStillbirth)

fit$Other <-
  FitFetoinfantSurvival(
    filt$Other,
    control = ControlFitFetoinfantSurvival(
      method = 'optim', hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$Other)

fit$Unknown <-
  FitFetoinfantSurvival(
    filt$Unknown,
    control = ControlFitFetoinfantSurvival(
      method = 'optim', hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$Other)

# Export ----------------------------------------------------------

qs_save(fit, paths$output$competing_risk_model_fits)
