# Export descriptive tables

# Init --------------------------------------------------------------------

here::i_am('src/31-descriptive_tables.R'); setwd(here::here())

library(qs2)
library(tidyverse)

paths <- list()
paths$input <- list(
  fnct_feto_infant_lt.R = 'src/00-fnct-feto_infant_lt.R',
  fetoinfant_lifetables.qs = 'out/30-fetoinfant_lifetables.qs',
  config.yaml = 'cfg/config.yaml'
)
paths$output <- list(
  descriptives.csv = 'out/31-descriptives.csv'
)

config <- yaml::read_yaml(paths$input$config.yaml)
source(paths$input$fnct_feto_infant_lt.R)

# constants
cnst <- list()

# Load data ---------------------------------------------------------------

filt <- qs_read(paths$input$fetoinfant_lifetables.qs)

# Create descriptive table ------------------------------------------------

descriptives <- rbind(
  bind_cols(
    var = " ",
    FILTCohortSize(filt$total) |> 
      left_join(FILTExposures(filt$total)) |>
      left_join(FILTTransitionsStratifiedTotal(filt$total)) |>
      select(stratum, N, E, D_F, D_I, B, C)    
  ),
  bind_cols(
    var = "Sex",
    FILTCohortSize(filt$sex14) |> 
      left_join(FILTExposures(filt$sex14)) |>
      left_join(FILTTransitionsStratifiedTotal(filt$sex14)) |>
      select(stratum, N, E, D_F, D_I, B, C)
  ),
  bind_cols(
    var = "Cohort",
    FILTCohortSize(filt$cohort) |> 
      left_join(FILTExposures(filt$cohort)) |>
      left_join(FILTTransitionsStratifiedTotal(filt$cohort)) |>
      select(stratum, N, E, D_F, D_I, B, C)
  ),
  bind_cols(
    var = "Origin",
    FILTCohortSize(filt$origin14) |> 
      left_join(FILTExposures(filt$origin14)) |>
      left_join(FILTTransitionsStratifiedTotal(filt$origin14)) |>
      select(stratum, N, E, D_F, D_I, B, C)
  ),
  bind_cols(
    var = "Education",
    FILTCohortSize(filt$education14) |> 
      left_join(FILTExposures(filt$education14)) |>
      left_join(FILTTransitionsStratifiedTotal(filt$education14)) |>
      select(stratum, N, E, D_F, D_I, B, C)
  )
)

# Export ------------------------------------------------------------------

write_csv(
  mutate(descriptives, across(where(is.numeric), ~ round(.x, 2))),
  paths$output$descriptives.csv
)
