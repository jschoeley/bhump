# Inference from fitted competing risks model

# Init --------------------------------------------------------------------

library(qs2)
library(tidyverse)

paths <- list()
paths$input <- list(
  competing_risk_model_fits.qs = 'tmp/50-competing_risks_model_fits.qs',
  lifetable_functions.R = 'src/00-fnct-feto_infant_lt.R',
  parametric_functions.R = 'src/00-fnct-parametric_survival_model.R',
  config.yaml = 'cfg/config.yaml'
)
paths$output <- list(
  competing_risks_statistics.qs = 'out/53-competing_risks_statistics.qs',
  bhump_by_cod.csv = 'out/53-bhump_by_cod.csv'
)

# fetoinfant lifetable functions
source(paths$input$lifetable_functions.R)
# fetoinfant parametric functions
source(paths$input$parametric_functions.R)

config <- yaml::read_yaml(paths$input$config.yaml)

# constants
cnst <-
  list(
    # censoring age, end of analysis time
    right_censoring_gestage = 77
  )

# tables
tab <- list()

# Load data ---------------------------------------------------------------

fit <- qs_read(paths$input$competing_risk_model_fits.qs)

# Share of deaths due to birth hump ---------------------------------------

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
  BirthHumpDeaths(fit$pcm, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$labor, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$congenital, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$maternal, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$convulsions, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$sepsis, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$hypoxia, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$respiratory, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$prematurity, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$sids, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$unspecific, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$otherspecific, x = cnst$right_censoring_gestage),
  BirthHumpDeaths(fit$unknown, x = cnst$right_censoring_gestage)
)
names(tab$rho_by_cod) <- config$cod_lookup$key
tab$rho_by_cod <- bind_rows(tab$rho_by_cod, .id = 'cod')

# Probability of Fetoinfant death -----------------------------------------

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
  ProbFetoInfantDeath(fit$pcm, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$labor, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$congenital, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$maternal, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$convulsions, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$sepsis, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$hypoxia, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$respiratory, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$prematurity, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$sids, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$unspecific, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$otherspecific, x = cnst$right_censoring_gestage),
  ProbFetoInfantDeath(fit$unknown, x = cnst$right_censoring_gestage)
)
names(tab$Fx_by_cod) <- config$cod_lookup$key
tab$Fx_by_cod <- bind_rows(tab$Fx_by_cod, .id = 'cod')

# COD composition of birth hump -------------------------------------------

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
  mutate(cod = factor(
    cod,
    levels = config$cod_lookup$key,
    labels = config$cod_lookup$shortlabel
  )) |>
  #arrange(-p_birth) %>%
  select(cod, p_birth, p_ontogen, p_total)

# Export ------------------------------------------------------------------

qs_save(tab, paths$output$competing_risks_statistics.qs)
write_csv(
  mutate(tab$bhump_by_cod, across(where(is.numeric), ~ round(.x, 2))),
  paths$output$bhump_by_cod.csv
)
