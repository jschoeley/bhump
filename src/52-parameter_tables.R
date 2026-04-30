# Export parameter tables of competing risks model

# Init --------------------------------------------------------------------

here::i_am('src/52-parameter_tables.R'); setwd(here::here())

library(qs2)
library(tidyverse)
library(gt)

paths <- list()
paths$input <- list(
  config.yaml = 'cfg/config.yaml',
  competing_risk_model_fits.qs = 'tmp/50-competing_risks_model_fits.qs'
)
paths$output <- list(
  stratapara.qs = 'tmp/52-stratapara.qs',
  stratapara.csv = 'out/52-stratapara.csv',
  stratapara.tex = 'out/52-stratapara.tex',
  codpara.qs = 'tmp/52-codpara.qs',
  codpara.csv = 'out/52-codpara.csv',
  codpara.tex = 'out/52-codpara.tex'
)

config <- yaml::read_yaml(paths$input$config.yaml)

# Load data ---------------------------------------------------------------

fit <- qs_read(paths$input$competing_risk_model_fits.qs)

# Parameter tables --------------------------------------------------------

PrintParameterTable <- function (filt_fit) {
  
  tab_of_pars <-
    filt_fit %>%
    unnest_legacy(par_rescaled_summary) %>%
    transmute(
      name,
      stratum = as.character(stratum),
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

# Parameter tables across population strata -------------------------------

stratapara <- list()

# total
stratapara$total <- PrintParameterTable(fit$total14)
# by sex
stratapara$sex <- PrintParameterTable(fit$sex)
# by cohort
stratapara$cohort <- PrintParameterTable(fit$cohort)
# by age
stratapara$origin <- PrintParameterTable(fit$origin)
# by education
stratapara$education <- PrintParameterTable(fit$education)

# turn into wide format table
stratapara.csv <-
  bind_rows(stratapara, .id = 'var') |>
  mutate(value = paste0(avg, ' (', ci025, ', ', ci975, ')')) |>
  select(var, name, stratum, value) |>
  pivot_wider(names_from = name, id_cols = c(var, stratum), values_from = value) |>
  select(-alpha2, -beta2)

# latex format
stratapara.tex <-
  stratapara.csv |>
  group_by(var) |>
  # variable lable only in first row
  mutate(var = c(var[1], rep('',n()-1))) |>
  ungroup() |>
  gt() |>
  tab_header(
    title = "Model parameters for competing-risks survival fit to feto-infant mortality across population strata."
  )

# Parameter tables by cod -------------------------------------------------

codpara <- list(
  pcm = PrintParameterTable(fit$pcm),
  labor = PrintParameterTable(fit$labor),
  congenital = PrintParameterTable(fit$congenital),
  maternal = PrintParameterTable(fit$maternal),
  convulsions = PrintParameterTable(fit$convulsions),
  sepsis = PrintParameterTable(fit$sepsis),
  hypoxia = PrintParameterTable(fit$hypoxia),
  respiratory = PrintParameterTable(fit$respiratory),
  prematurity = PrintParameterTable(fit$prematurity),
  sids = PrintParameterTable(fit$sids),
  otherspecific = PrintParameterTable(fit$otherspecific),
  unspecific = PrintParameterTable(fit$unspecific),
  unknown = PrintParameterTable(fit$unknown)
)

codpara.csv <-
  codpara |>
  bind_rows(.id = 'var') |>
  mutate(value = paste0(avg, ' (', ci025, ', ', ci975, ')')) |>
  select(var, name, value) |>
  pivot_wider(names_from = name, id_cols = c(var), values_from = value) |>
  mutate(
    alpha2 = ifelse(beta1==beta2, '', alpha2),
    beta2 = ifelse(beta1==beta2, '', beta2)
  )

# latex format
codpara.tex <-
  codpara.csv |>
  gt() |>
  tab_header(
    title = "Model parameters for competing-risks survival fit to cause specific feto-infant life tables."
  )

# Export ------------------------------------------------------------------

qs_save(stratapara, paths$output$stratapara.qs)
write_csv(stratapara.csv, paths$output$stratapara.csv)
gtsave(stratapara.tex, paths$output$stratapara.tex)

qs_save(codpara, paths$output$codpara.qs)
write_csv(codpara.csv, paths$output$codpara.csv)
gtsave(codpara.tex, paths$output$codpara.tex)
