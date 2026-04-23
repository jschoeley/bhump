# Aggregate microdata on fetal and infant death to feto-infant lifetable

# Init --------------------------------------------------------------------

library(qs2)
library(tidyverse)
library(yaml)

paths <- list()
paths$input <- list(
  fetoinfant.qs = 'tmp/21-fetoinfant.qs',
  lifetable_functions.R = 'src/00-fnct-feto_infant_lt.R',
  config.yaml = 'cfg/config.yaml'
)
paths$output <- list(
  fetoinfant_lt.qs = 'out/30-fetoinfant_lifetables.qs'
)

# multistate aggregation of events and exposure times
source(paths$input$lifetable_functions.R)

cnst <- list()
# aggregate over each single week LMP from 24 to 77
cnst$lifetable_breaks <- 24:77

config <- read_yaml(paths$input$config.yaml)

# Load data ---------------------------------------------------------------

# individual level event history data on feto-infant survival
fetoinfant_event_histories <- qs_read(paths$input$fetoinfant.qs)

# Recode categorical strata -----------------------------------------------

fetoinfant_event_histories <-
  fetoinfant_event_histories |>
  mutate(
    maternal_origin = case_when(
      race_and_hispanic_orig_of_mother %in% c(
        'Mexican', 'Puerto Rican', 'Cuban',
        'Central or South American',
        'Other and unknown Hispanic', 'Hispanic'
      ) ~ 'Hispanic',
      race_and_hispanic_orig_of_mother %in% c(
        'Non-Hispanic other races',
        'Non-Hispanic AIAN (only)',
        'Non-Hispanic Asian (only)',
        'Non-Hispanic NHOPI (only)',
        'Non-Hispanic more than one race'
      ) ~ 'Other',
      race_and_hispanic_orig_of_mother %in% c(
        'Non-Hispanic Black (only)',
        'Non-Hispanic Black'
      ) ~ 'Non-Hispanic Black',
      race_and_hispanic_orig_of_mother %in% c(
        'Non-Hispanic White (only)',
        'Non-Hispanic White'
      ) ~ 'Non-Hispanic White',
      race_and_hispanic_orig_of_mother ==
        'Origin unknown or not stated' ~ 'Unknown',
      is.na(race_and_hispanic_orig_of_mother) ~ 'Unknown'
    ),
    maternal_education = case_when(
      education_of_mother %in% c(1:2) ~ 'Primary',
      education_of_mother %in% c(3:4) ~ 'High school',
      education_of_mother %in% c(5:6) ~ 'Associate',
      education_of_mother %in% c(6:8) ~ 'Bachelor, Master, Doctorate',
      education_of_mother %in% c(9, ' ') ~ 'Unknown',
      is.na(education_of_mother) ~ 'Unknown'
    )
  )

# Aggregate observations into multistate lt -------------------------------

filt <- list()

# total (cohort 2014)
filt$total14 <-
  FILT(
    fetoinfant_event_histories %>%
      filter(date_of_conception_y == 2014),
    breaks = cnst$lifetable_breaks
  )
filt$total14

# by sex (cohort 2014)
filt$sex14 <-
  FILT(
    fetoinfant_event_histories %>%
      filter(date_of_conception_y == 2014) %>%
      mutate(sex = forcats::fct_na_value_to_level(sex, level = 'Unknown')),
    breaks = cnst$lifetable_breaks,
    stratum = sex
  )
filt$sex14

# by cohort
filt$cohort <-
  FILT(
    fetoinfant_event_histories,
    breaks = cnst$lifetable_breaks,
    stratum = date_of_conception_y
  )
filt$cohort

# by origin of mother (cohort 2014)
filt$origin14 <-
  fetoinfant_event_histories %>%
  filter(date_of_conception_y == 2014) %>%
  FILT(
    df = .,
    breaks = cnst$lifetable_breaks,
    stratum = maternal_origin
  )
filt$origin14

# by maternal education
filt$education14 <-
  fetoinfant_event_histories %>%
  filter(date_of_conception_y == 2014) %>%
  FILT(
    df = .,
    breaks = cnst$lifetable_breaks,
    stratum = maternal_education
  )
filt$education14

# by maternal origin x education
filt$origineducation14 <-
  fetoinfant_event_histories %>%
  filter(date_of_conception_y == 2014) %>%
  mutate(
    origineduc = interaction(maternal_origin,
                             maternal_education)
  ) %>%
  FILT(
    df = .,
    breaks = cnst$lifetable_breaks,
    stratum = origineduc
  )
filt$origineducation14

# Cause of death ----------------------------------------------------------

# cohort 2014
map(config$cod_lookup$key, ~{
  cat(.x, '\n')
  filt[[.x]] <<-
    fetoinfant_event_histories %>%
    filter(date_of_conception_y == 2014) %>%
    FILT(
      df = .,
      d_out = 'destination2', death_state_name = .x,
      breaks = cnst$lifetable_breaks
    )  
})

# Export ------------------------------------------------------------------

qs_save(filt, file = paths$output$fetoinfant_lt.qs)
