# Plot hazards over population strata

# Init --------------------------------------------------------------------

here::i_am('src/63-plot_hazards_by_cod.R'); setwd(here::here())

library(qs2)
library(tidyverse)

paths <- list()
paths$input <- list(
  lifetables_and_parametric_fits_by_cod.qs =
    'tmp/51-lifetables_and_parametric_fits_by_cod.qs',
  competing_risks_statistics.qs =
    'out/53-competing_risks_statistics.qs',
  figure_specs.R = 'src/00-figure_specifications.R',
  config.yaml = 'cfg/config.yaml'
)
paths$output <- list(
  hazards_by_cod.svg = 'out/63-hazards_by_cod.svg',
  hazards_by_cod.qs = 'out/63-hazards_by_cod.qs'
)

# figure specs
source(paths$input$figure_specs.R)

config <- yaml::read_yaml(paths$input$config.yaml)

# constants
cnst <- list()

# Load data ---------------------------------------------------------------

filt_cod <- qs_read(paths$input$lifetables_and_parametric_fits_by_cod.qs)
cr <- qs_read(paths$input$competing_risks_statistics.qs)

# Plot hazards by cod -----------------------------------------------------

hazards_by_cod <- list()

hazards_by_cod$data <- list(
  hazards = mutate(
    filt_cod,
    cod = factor(cod, config$cod_lookup$label, config$cod_lookup$shortlabel)
  ),
  rho_by_cod = mutate(
    cr$rho_by_cod,
    cod = factor(cod, config$cod_lookup$key, config$cod_lookup$shortlabel)
  ),
  Fx_by_cod = mutate(
    cr$Fx_by_cod,
    cod = factor(cod, config$cod_lookup$key, config$cod_lookup$shortlabel)
  )
)

hazards_by_cod$plot <-
  hazards_by_cod$data$hazards |>
  ggplot(aes(x = x)) +
  geom_line(aes(y = avg_total_hx*1e5, color = cod), na.rm = FALSE) +
  geom_point(aes(y = m*1e5, color = cod), size = 0.5, alpha = 0.5) +
  geom_text(aes(
    label = paste0('ρ=', formatC(avg_p_birth*100,
                                 digits = 1, format = 'f')),
    color = cod
  ), x = 76, y = 0.7, family = 'sans', parse = FALSE,
  hjust = 1, vjust = 0, size = 2.5,
  data = hazards_by_cod$data$rho_by_cod
  ) +
  geom_text(aes(
    label = paste0('F(77)=1:', formatC(avg_total_iFx,
                                       digits = 0, format = 'f')),
    color = cod
  ), x = 76, y = 0.2, family = 'sans', parse = FALSE,
  hjust = 1, vjust = 0, size = 2.5,
  data = 
    hazards_by_cod$data$Fx_by_cod
  ) +
  scale_color_manual(values = config$cod_lookup$color) +
  scale_y_continuous(expand = c(0,0), trans = 'log10',
                     breaks = c(0.01, 0.1, 1, 10),
                     labels = c('0.01', '0.1', '1', '10')) +
  scale_x_continuous(breaks = c(24, 40, 76)) +
  facet_wrap(~cod, ncol = 3) +
  fig_spec$MyGGplotTheme(panel_border = TRUE) +
  labs(y = 'Feto-infant deaths per 100,000 weeks at risk',
       x = 'Week of gestation',
       fill = 'Cause of death') +
  guides(color = 'none') +
  coord_cartesian(ylim = c(0.001, 15)) +
  geom_vline(
    aes(xintercept = 40),
    lty = 3,
    size = fig_spec$line_size_m
  )

# Export ------------------------------------------------------------------

qs_save(hazards_by_cod$data, paths$output$hazards_by_cod.qs)
fig_spec$ExportSVG(
  hazards_by_cod$plot,
  paths$output$hazards_by_cod.svg,
  width = fig_spec$width,
  height = 1.2*fig_spec$width
)
