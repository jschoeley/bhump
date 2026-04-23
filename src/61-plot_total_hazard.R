# Plot total hazard with fitted competing risks model

# Init --------------------------------------------------------------------

library(qs2)
library(tidyverse)

paths <- list()
paths$input <- list(
  fetoinfant_lifetables.qs = 'out/30-fetoinfant_lifetables.qs',
  competing_risk_model_fits.qs = 'tmp/50-competing_risks_model_fits.qs',
  figure_specs.R = 'src/00-figure_specifications.R',
  config.yaml = 'cfg/config.yaml'
)
paths$output <- list(
  overall_hazard_and_fit.qs = 'out/61-overall_hazard_and_fit.qs',
  overall_hazard_and_fit.svg = 'out/61-overall_hazard_and_fit.svg'
)

# figure specs
source(paths$input$figure_specs.R)

config <- yaml::read_yaml(paths$input$config.yaml)

# constants
cnst <- list()

# Input -------------------------------------------------------------------

filt <- qs_read(paths$input$fetoinfant_lifetables.qs)
fit <- qs_read(paths$input$competing_risk_model_fits.qs)

# Plot overall hazard -----------------------------------------------------

overall_hazard_and_fit <- list()
overall_hazard_and_fit$data <-
  list(
    lifetable = filt$total14,
    fit = fit$total14$pred_summary[[1]]
  )

# total
overall_hazard_and_fit$plot <-
  overall_hazard_and_fit$data$lifetable |>
  ggplot(aes(x = x)) +
  geom_col(aes(y = B/3e9), position = position_nudge(x = +0.5)) +
  geom_point(aes(y = D/E), color = 'grey50', position = position_nudge(x = +0.5)) +
  geom_line(
    aes(x = x, y = avg_total_hx),
    data = overall_hazard_and_fit$data$fit
  ) +
  geom_line(
    aes(x = x, y = avg_birth_hx),
    data = overall_hazard_and_fit$data$fit,
    color = 'blue'
  ) +
  geom_line(
    aes(x = x, y = avg_ontogen_hx),
    data = overall_hazard_and_fit$data$fit,
    color = 'red', lty = 2
  ) +
  scale_x_continuous(breaks = seq(24, 76, 2)) +
  scale_y_continuous(labels = ~.x*1e5) +
  coord_cartesian(expand = FALSE, ylim = c(0, 50/1e5)) +
  labs(
    x = 'Week of gestation',
    y = 'Feto-infant deaths per 100,000 weeks at risk'
  ) +
  fig_spec$MyGGplotTheme()
overall_hazard_and_fit$plot

# Export ------------------------------------------------------------------

qs_save(overall_hazard_and_fit$data, paths$output$overall_hazard_and_fit.qs)
fig_spec$ExportSVG(
  overall_hazard_and_fit$plot,
  paths$output$overall_hazard_and_fit.svg,
  width = fig_spec$width,
  height = fig_spec$width*0.7
)
