# Plot the birth hump component by cause of death

# Init --------------------------------------------------------------------

library(qs2)
library(tidyverse)

paths <- list()
paths$input <- list(
  lifetables_and_parametric_fits_by_cod.qs =
    'tmp/51-lifetables_and_parametric_fits_by_cod.qs',
  figure_specs.R = 'src/00-figure_specifications.R',
  config.yaml = 'cfg/config.yaml'
)
paths$output <- list(
  birthhump_cod_separate.qs = 'out/60-birthhump_cod_separate.qs',
  birthhump_cod_separate.svg = 'out/60-birthhump_cod_separate.svg',
  birthhump_cod_joint.qs = 'out/60-birthhump_cod_joint.qs',
  birthhump_cod_joint.svg = 'out/60-birthhump_cod_joint.svg'
)

# figure specs
source(paths$input$figure_specs.R)

config <- yaml::read_yaml(paths$input$config.yaml)

# constants
cnst <-
  list(
    gestage_brk = seq(24, 77, by = 4),
    lifetable_breaks = 24:77,
    left_truncation_gestage = 24,
    right_censoring_gestage = 77
  )

# Input -------------------------------------------------------------------

filt_cod <- qs_read(paths$input$lifetables_and_parametric_fits_by_cod.qs)

# Plot birth hump composition by COD --------------------------------------

birthhump_cod_joint <- list()
birthhump_cod_joint$data <-
  filt_cod |>
  select(cod, x, avg_birth_hx)

birthhump_cod_joint$plot <-
  birthhump_cod_joint$data |>
  ggplot(aes(x = x)) +
  geom_area(aes(y = avg_birth_hx*1e5, fill = cod)) +
  scale_fill_manual(values = rev(config$cod_lookup$color)) +
  scale_x_continuous(breaks = c(24, 40, 76)) +
  scale_y_continuous(expand = c(0,0)) +
  fig_spec$MyGGplotTheme(panel_border = FALSE) +
  labs(y = 'Feto-infant deaths per 100k weeks',
       x = 'Week of gestation',
       fill = 'Cause of death') +
  theme(
    legend.position = c(0.7, 0.5),
    legend.background = element_blank(),
    legend.key.size = unit(0.4, 'cm'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6)
  )
birthhump_cod_joint

# Separately plot birth hump by cod ---------------------------------------

birthhump_cod_separate <- list()
birthhump_cod_separate$data <-
  filt_cod |>
  select(cod, x, avg_birth_hx) |>
  mutate(
    cod = factor(cod, config$cod_lookup$label, config$cod_lookup$label)
  )# |>
  #filter(!cod %in% c('Prematurity', 'Accidents and violence', 'Sudden Infant Death'))

birthhump_cod_separate$plot <-
  birthhump_cod_separate$data |>
  ggplot(aes(x = x)) +
  geom_area(aes(y = avg_birth_hx*1e5, fill = cod)) +
  scale_fill_manual(values = config$cod_lookup$color) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(24, 40, 76)) +
  facet_wrap(~cod, ncol = 3) +
  fig_spec$MyGGplotTheme(panel_border = TRUE) +
  labs(y = 'Feto-infant deaths per 100k weeks',
       x = 'Week of gestation',
       fill = 'Cause of death') +
  guides(fill = 'none')
birthhump_cod_separate$plot

# Export ------------------------------------------------------------------

qs_save(birthhump_cod_joint$data, paths$output$birthhump_cod_joint.qs)
fig_spec$ExportSVG(
  birthhump_cod_joint$plot,
  paths$output$birthhump_cod_joint.svg,
  width = fig_spec$width,
  height = fig_spec$width*0.6
)

qs_save(birthhump_cod_separate$data, paths$output$birthhump_cod_separate.qs)
fig_spec$ExportSVG(
  birthhump_cod_separate$plot,
  paths$output$birthhump_cod_separate.svg,
  width = fig_spec$width,
  height = fig_spec$width*1.4
)
