# Decompose differences in feto-infant mortality between groups
# into level, ontogenescence, and transition components

# The decomposition is performed via the horiuchi method. After
# fitting the competing risks hazard model we evaluate the predicted
# difference in overall feto-infant probability of death from
# week 24 through 76 for two groups. We then attribute this difference
# to differences in the fitted model parameters, distinguishing between
# level, ontogenescent (slope), and birth-hump contributions.

# Init --------------------------------------------------------------------

here::i_am('src/54-parametric_decomposition.R'); setwd(here::here())

library(tidyverse)
library(qs2)
library(readr)
library(gt)

paths <- list()
paths$input <- list(
  competing_risk_model_fits.qs = 'tmp/50-competing_risks_model_fits.qs',
  parametric_decomposition.R = 'src/00-fnct-parametric_decomposition.R',
  parametric_survival_model.R = 'src/00-fnct-parametric_survival_model.R',
  config.yaml = 'cfg/config.yaml'
)
paths$output <- list(
  parametric_decompositions.qs =
    'out/54-parametric_decompositions.qs',
  parametric_decompositions.csv =
    'out/54-parametric_decompositions.csv',
  parametric_decompositions.tex =
    'out/54-parametric_decompositions.tex'
)

config <- yaml::read_yaml(paths$input$config.yaml)

# fetoinfant parametric functions
source(paths$input$parametric_decomposition.R)
source(paths$input$parametric_survival_model.R)

# constants
cnst <- list()

# Load data ---------------------------------------------------------------

# a list of feto-infant lifetable as FILT objects
fit <- qs_read(paths$input$competing_risk_model_fits.qs)

# Parametric decomposition of mortality differences -----------------------

decomp <- list(
  sex = list(),
  cohort = list(),
  origin = list(),
  education = list()
)

# by sex
decomp$sex$female_male <- DecomposeFetoInfantDeaths(fit$sex, Female, Male)
decomp$sex$tabl <-
  bind_cols(
    `Name` = decomp$sex$female_male$tabl[,1][[1]],
    `Male-Female` = decomp$sex$female_male$tabl[,2][[1]]
  )

# by cohort
decomp$cohort$`89vs99` <-
  DecomposeFetoInfantDeaths(fit$cohort, `1989`, `1999`)
decomp$cohort$`99vs09` <-
  DecomposeFetoInfantDeaths(fit$cohort, `1999`, `2009`)
decomp$cohort$`09vs14` <-
  DecomposeFetoInfantDeaths(fit$cohort, `2009`, `2014`)
decomp$cohort$`89vs14` <-
  DecomposeFetoInfantDeaths(fit$cohort, `1989`, `2014`)
decomp$cohort$tabl <-
  bind_cols(
    `Name` = decomp$cohort$`89vs99`$tabl[,1][[1]],
    `1989-1999` = decomp$cohort$`89vs99`$tabl[,-1][[1]],
    `1999-2009` = decomp$cohort$`99vs09`$tabl[,-1][[1]],
    `2009-2014` = decomp$cohort$`09vs14`$tabl[,-1][[1]],
    `1989-2014` = decomp$cohort$`89vs14`$tabl[,-1][[1]]
  )

# by origin
decomp$origin$white_black <-
  DecomposeFetoInfantDeaths(fit$origin, `Non-Hispanic White`, `Non-Hispanic Black`)
decomp$origin$white_hispanic <-
  DecomposeFetoInfantDeaths(fit$origin, `Non-Hispanic White`, `Hispanic`)
decomp$origin$black_hispanic <-
  DecomposeFetoInfantDeaths(fit$origin, `Non-Hispanic Black`, `Hispanic`)
decomp$origin$tabl <-
  bind_cols(
    `Name` = decomp$origin$white_black$tabl[,1][[1]],
    `White-Black` = decomp$origin$white_black$tabl[,-1][[1]],
    `White-Hispanic` = decomp$origin$white_hispanic$tabl[,-1][[1]],
    `Black-Hispanic` = decomp$origin$black_hispanic$tabl[,-1][[1]]
  )

# by education
decomp$education$primary_academic <-
  DecomposeFetoInfantDeaths(fit$education, `Primary`, `Bachelor, Master, Doctorate`)
decomp$education$highschool_academic <-
  DecomposeFetoInfantDeaths(fit$education, `High school`, `Bachelor, Master, Doctorate`)
decomp$education$associate_academic <-
  DecomposeFetoInfantDeaths(fit$education, `Associate`, `Bachelor, Master, Doctorate`)
decomp$education$tabl <-
  bind_cols(
    `Name` = decomp$education$primary_academic$tabl[,1][[1]],
    `Primary-Academic` = decomp$education$primary_academic$tabl[,-1][[1]],
    `Highschool-Academic` = decomp$education$highschool_academic$tabl[,-1][[1]],
    `Associate-Academic` = decomp$education$associate_academic$tabl[,-1][[1]]
  )

# Prepare for export ------------------------------------------------------

parametric_decompositions.csv <- bind_rows(
  sex =
    decomp$sex$tabl[,-1] |> t() |> as.data.frame() |> rownames_to_column(),
  origin =
    decomp$origin$tabl[,-1] |> t() |> as.data.frame() |> rownames_to_column(),
  cohort =
    decomp$cohort$tabl[,-1] |> t() |> as.data.frame() |> rownames_to_column(),
  education =
    decomp$education$tabl[,-1] |> t() |> as.data.frame() |> rownames_to_column(),
  .id = 'Var'
)
colnames(parametric_decompositions.csv) <-
  c('Var', 'Contrast', unlist(decomp$sex$tabl[,1]))


# latex format
parametric_decompositions.tex <-
  parametric_decompositions.csv |>
  group_by(Var) |>
  # variable lable only in first row
  mutate(Var = c(Var[1], rep('',n()-1))) |>
  ungroup() |>
  gt() |>
  tab_header(
    title = "Contributions of model components to differences in one year post-viability feto-infant survival."
  )

# Export ------------------------------------------------------------------

qs_save(decomp, paths$output$parametric_decompositions.qs)
write_csv(parametric_decompositions.csv, paths$output$parametric_decompositions.csv)
gtsave(parametric_decompositions.tex, paths$output$parametric_decompositions.tex)