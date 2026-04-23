# Extrapolate ontogenescent component until conception

# Init --------------------------------------------------------------------

library(tidyverse)
library(qs2)

paths <- list()
paths$input <- list(
  competing_risk_model_fits.qs = 'tmp/50-competing_risks_model_fits.qs',
  parametric_functions.R = 'src/00-fnct-parametric_survival_model.R',
  config.yaml = 'cfg/config.yaml',
  figure_specs.R = 'src/00-figure_specifications.R'
)
paths$output <- list(
  backwards_extrapolation.qs =
    'out/55-backwards_extrapolation.qs',
  backwards_extrapolation.svg =
    'out/55-backwards_extrapolation.svg'
)

config <- yaml::read_yaml(paths$input$config.yaml)

source(paths$input$figure_specs.R)
# fetoinfant parametric functions
source(paths$input$parametric_functions.R)

# constants
cnst <-
  list(
    cod = config$cod_lookup$key
  )

# Input -------------------------------------------------------------------

fit <- qs_read(paths$input$competing_risk_model_fits.qs)

# Backwards extrapolation -------------------------------------------------

p <-
  fit$cohort |>
  filter(stratum == '2014') |>
  unnest(par_summary) |>
  pull(avg)

p2 <- RescaleParameters(p, model = 'basic', split = 14)
p2$zeta <- p2$zeta + 24
p2$alpha1 <- exp(log(p2$alpha1)-24*(-p2$beta1))

odds_of_conception_death <- 1/(1-FetoinfantSurv((76+24), pars = p2, model = 'basic'))

backwards_extrapolation <-
  tibble(
  x = 0:(76+24),
  hx = FetoinfantHzrd(0:(76+24), pars = p2, model = 'basic'),
  Sx = FetoinfantSurv(0:(76+24), pars = p2, model = 'basic')
) |>
  ggplot() +
  aes(x = x, y = Sx) +
  geom_point() +
  annotate(
    'text', x = 75, y = 1, 
    label = paste0('F(100)=1:', round(odds_of_conception_death, 0))
  ) +
  fig_spec$MyGGplotTheme() +
  labs(
    x = 'Weeks of gestation',
    y = 'S(x)'
  )

# Export ------------------------------------------------------------------

fig_spec$ExportSVG(
  backwards_extrapolation,
  paths$output$backwards_extrapolation.svg,
  width = fig_spec$width,
  height = 0.7*fig_spec$width
)
