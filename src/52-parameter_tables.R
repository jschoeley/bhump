# Export parameter tables of competing risks model

# Init ------------------------------------------------------------

library(qs2)
library(tidyverse)

paths <- list()
paths$input <- list(
  config = 'cfg/config.yaml',
  competing_risk_model_fits = 'tmp/50-competing_risks_model_fits.qs'
)
paths$output <- list(
  competing_risk_model_parameter_tables = 'tmp/52-competing_risks_model_parameter_tables.qs'
)

config <- yaml::read_yaml(paths$input$config)

# parameter tables
partab <- list()

# Load data -------------------------------------------------------

fit <- qs_read(paths$input$competing_risk_model_fits)

# Parameter tables --------------------------------------------------------

PrintParameterTable <- function (filt_fit) {
  
  tab_of_pars <-
    filt_fit %>%
    unnest_legacy(par_rescaled_summary) %>%
    transmute(
      name,
      stratum,
      avg = avg,
      ci025 = ci025,
      ci975 = ci975
    ) %>%
    mutate_at(
      c('avg', 'ci025', 'ci975'),
      ~ formatC(., format = 'e', digits = 1)
    )
  
  return(tab_of_pars)
  
}

# by sex
partab$sex <- PrintParameterTable(fit$sex)
# by cohort
partab$cohort <- PrintParameterTable(fit$cohort)
# by age
partab$origin <- PrintParameterTable(fit$origin)
# by education
partab$education <- PrintParameterTable(fit$education)

# Range of rate of ontogenescence ---------------------------------

do.call(rbind, partab) |>
  filter(name == 'beta1') |>
  pull(avg) |>
  range()

# Range of level of feto-infant mortality -------------------------

do.call(rbind, partab) |>
  filter(name == 'alpha1') |>
  pull(avg) |>
  range()

# Export ----------------------------------------------------------

qs_save(partab, paths$output$competing_risk_model_parameter_tables)
