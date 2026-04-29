# Combine cause-specific feto-infant life tables with parametric fits

# Init --------------------------------------------------------------------

here::i_am('src/51-combine_lifetables_and_parametric_fits_by_cod.R'); setwd(here::here())

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
  lifetables_and_parametric_fits_by_cod.qs =
    'tmp/51-lifetables_and_parametric_fits_by_cod.qs',
  total_vs_cod_hazard_alignment_check.svg =
    'tmp/total_vs_cod_hazard_alignment_check.svg'
)

# figure specs
source(paths$input$figure_specs.R)

config <- yaml::read_yaml(paths$input$config.yaml)

# constants
cnst <-
  list(
    cod = config$cod_lookup$key
  )

# Load data ---------------------------------------------------------------

fit <- qs_read(paths$input$competing_risk_model_fits.qs)

# Assemble lifetables and fit for cause of death --------------------------

lifetables_and_parametric_fits_by_cod <- map(cnst$cod, ~{
  parametric <-  fit[[.x]] %>%
    unnest_legacy(pred_summary) %>%
    select(stratum, x, starts_with(c('avg', 'q025', 'q975')))
  lifetable <- fit[[.x]] %>%
    unnest_legacy(lifetable) %>%
    select(stratum, x, n, starts_with(c('N', 'E', 'D', 'B', 'C', 'm', 'p', 'S', 'F')))
  full_join(parametric, lifetable, by = c('stratum', 'x')) %>% mutate(cod = .x)
}) %>% bind_rows() %>%
  group_by(cod) %>%
  mutate(
    max = max(avg_birth_hx, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    cod = factor(cod, rev(config$cod_lookup$key),
                 rev(config$cod_lookup$label))
  )

# test total lifetable deaths should be equal to sum of cause-specific deaths
fit$total14 %>% unnest_legacy(lifetable) %>% pull(D) %>% sum()
sum(lifetables_and_parametric_fits_by_cod$D, na.rm = T)

# test if sum of cause specific hazards approximately aligns with
# total life table death rates
total_vs_cod_hazard_alignment_check <-
  lifetables_and_parametric_fits_by_cod %>%
  ggplot(aes(x = x)) +
  geom_area(aes(y = avg_total_hx*1e5, fill = cod)) +
  scale_fill_manual(values = rev(config$cod_lookup$color)) +
  scale_x_continuous(breaks = c(24, 40, 76)) +
  scale_y_continuous(expand = c(0,0)) +
  fig_spec$MyGGplotTheme(panel_border = TRUE) +
  labs(y = 'Feto-infant deaths per 100k weeks',
       x = 'Week of gestation',
       fill = 'Cause of death') +
  theme(legend.position = c(0.7, 0.5), legend.background = element_blank()) +
  geom_step(aes(x = x, y = m*1e5), size = 2,
            data = fit$total14 %>% unnest(lifetable)) +
  geom_line(aes(x = x, y = avg_total_hx*1e5), size = 2,
            data = fit$total14 %>% unnest(pred_summary))

# Export ------------------------------------------------------------------

qs_save(
  lifetables_and_parametric_fits_by_cod,
  paths$output$lifetables_and_parametric_fits_by_cod.qs
)

fig_spec$ExportSVG(
  total_vs_cod_hazard_alignment_check,
  filename = paths$output$total_vs_cod_hazard_alignment_check.svg,
  width = fig_spec$width,
  height = 0.7*fig_spec$width
)
