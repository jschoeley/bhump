# Inference from fitted competing risks model

# Init ------------------------------------------------------------

library(qs2)
library(tidyverse)

paths <- list()
paths$input <- list(
  competing_risk_model_fits = 'tmp/50-competing_risks_model_fits.qs',
  lifetable_functions = 'src/00-fnct-feto_infant_lt.R',
  parametric_functions = 'src/00-fnct-parametric_survival_model.R',
  config = 'cfg/config.yaml'
)
paths$output <- list(
  competing_risks_statistics = 'out/53-competing_risks_statistics.qs'
)

# fetoinfant lifetable functions
source(paths$input$lifetable_functions)
# fetoinfant parametric functions
source(paths$input$parametric_functions)

config <- yaml::read_yaml(paths$input$config)

# constants
cnst <-
  list(
    # censoring age, end of analysis time
    right_censoring_gestage = 77
  )

# tables
tab <- list()

# Load data -------------------------------------------------------

fit <- qs_read(paths$input$competing_risk_model_fits)

# Share of deaths due to birth hump -----------------------------

# the share of feto-infant deaths in the year following fetal viability
# which are contributed by the birth hump component

tab$rho <- list(
  # by sex
  sex = BirthHumpDeaths(fit$sex, x = cnst$right_censoring_gestage),
  # by cohort
  cohort = BirthHumpDeaths(fit$cohort, x = cnst$right_censoring_gestage),
  # by origin
  origin = BirthHumpDeaths(fit$origin, x = cnst$right_censoring_gestage),
  # by education
  education = BirthHumpDeaths(fit$education, x = cnst$right_censoring_gestage)
)

# by cause of death
tab$rho_by_cod <- list(
  BirthHumpDeaths(fit$PlacentaCordMembrane, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$LabourBirth, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$Malformations, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$Maternal, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$CerebralConvulsions, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$Infections, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$HypoxiaAsphyxie, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$RespiratoryCardio, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$FetalGrowthPremature, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$SID, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$UnspecificStillbirth, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$Other, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$Unknown, x = cnst$right_censoring_gestage)
)
names(tab$rho_by_cod) <- config$cod_lookup$key
tab$rho_by_cod <- bind_rows(tab$rho_by_cod, .id = 'cod')

# Probability of Fetoinfant death ---------------------------------

# the probability of a feto-infant death in the year following fetal
# viability

tab$Fx_by_stratum <- list(
  # by sex
  sex = ProbFetoInfantDeath(fit$sex),
  # by cohort
  cohort = ProbFetoInfantDeath(fit$cohort),
  # by origin
  origin = ProbFetoInfantDeath(fit$origin),
  # by education
  education = ProbFetoInfantDeath(fit$education)
)

# by cause of death
tab$Fx_by_cod <- list(
  ProbFetoInfantDeath(fit$PlacentaCordMembrane, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$LabourBirth, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$Malformations, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$Maternal, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$CerebralConvulsions, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$Infections, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$HypoxiaAsphyxie, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$RespiratoryCardio, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$FetalGrowthPremature, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$SID, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$UnspecificStillbirth, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$Other, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$Unknown, x = cnst$right_censoring_gestage)
)
names(tab$Fx_by_cod) <- config$cod_lookup$key
tab$Fx_by_cod <- bind_rows(tab$Fx_by_cod, .id = 'cod')

# COD composition of birth hump -----------------------------------

# the probability of dying of a specific cause if the death is
# due to the transitional risk

# cause-of-death composition of birth-hump
tab$bhump_by_cod <-
  tab$Fx_by_cod %>%
  mutate(
    p_birth = avg_birth_Fx / sum(avg_birth_Fx)*100,
    p_ontogen = avg_ontogen_Fx / sum(avg_ontogen_Fx)*100,
    p_total = avg_total_Fx / sum(avg_total_Fx)*100
  ) %>%
  arrange(p_birth) %>%
  select(cod, p_birth, p_ontogen, p_total)

# Export ----------------------------------------------------------

qs_save(tab, paths$output$competing_risks_statistics)
