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

# Fit feto-infant survival by cause of death ----------------------

# by cause of death
fit$PlacentaCordMembrane <-
  FitFetoinfantSurvival(
    filt$PlacentaCordMembrane,
    control = ControlFitFetoinfantSurvival(hessian_inverse = 'choleskypivot')
  )
PlotHazards(fit$PlacentaCordMembrane)

fit$LabourBirth <-
  FitFetoinfantSurvival(
    filt$LabourBirth,
    control = ControlFitFetoinfantSurvival(hessian_inverse = 'choleskypivot')
  )
PlotHazards(fit$LabourBirth)

fit$Malformations <-
  FitFetoinfantSurvival(
    filt$Malformations,
    control = ControlFitFetoinfantSurvival(hessian_inverse = 'choleskypivot')
  )
PlotHazards(fit$Malformations)

fit$Maternal <-
  FitFetoinfantSurvival(
    filt$Maternal,
    control = ControlFitFetoinfantSurvival(hessian_inverse = 'choleskypivot')
  )
PlotHazards(fit$Maternal)

fit$CerebralConvulsions <-
  FitFetoinfantSurvival(
    filt$CerebralConvulsions,
    control = ControlFitFetoinfantSurvival(
      zeta_range = c(39, 43)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$CerebralConvulsions)

fit$Infections <-
  FitFetoinfantSurvival(
    filt$Infections,
    control = ControlFitFetoinfantSurvival(
      zeta_range = c(38, 42)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$Infections)

fit$HypoxiaAsphyxie <-
  FitFetoinfantSurvival(
    filt$HypoxiaAsphyxie,
    control = ControlFitFetoinfantSurvival(
      zeta_range = c(38, 42)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$HypoxiaAsphyxie)

fit$RespiratoryCardio <-
  FitFetoinfantSurvival(
    filt$RespiratoryCardio,
    control = ControlFitFetoinfantSurvival(
      model = 'flexible1',
      zeta_range = c(37, 41)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$RespiratoryCardio)

fit$FetalGrowthPremature <-
  FitFetoinfantSurvival(
    filt$FetalGrowthPremature,
    control = ControlFitFetoinfantSurvival(
      model = 'flexible1',
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$FetalGrowthPremature)

# we fit a restricted model without birth hump
fit$SID <-
  FitFetoinfantSurvival(
    filt$SID,
    control = ControlFitFetoinfantSurvival(
      model = 'flexible2',
      zeta_range = c(40, 45)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot',
    )
  )
PlotHazards(fit$SID)

fit$UnspecificStillbirth <-
  FitFetoinfantSurvival(
    filt$UnspecificStillbirth,
    control = ControlFitFetoinfantSurvival(
      zeta_range = c(37, 41)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$UnspecificStillbirth)

fit$Other <-
  FitFetoinfantSurvival(
    filt$Other,
    control = ControlFitFetoinfantSurvival(
      zeta_range = c(37, 41)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$Other)

fit$Unknown <-
  FitFetoinfantSurvival(
    filt$Unknown,
    control = ControlFitFetoinfantSurvival(
      zeta_range = c(37, 41)-cnst$left_truncation_gestage,
      hessian_inverse = 'choleskypivot'
    )
  )
PlotHazards(fit$Unknown)


# Export ----------------------------------------------------------

qs_save(fit, paths$output$competing_risk_model_fits)
