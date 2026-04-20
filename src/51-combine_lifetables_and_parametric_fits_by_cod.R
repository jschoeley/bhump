# Combine cause-specific feto-infant life tables with parametric fits

# Init ------------------------------------------------------------

library(qs2)
library(tidyverse)

paths <- list()
paths$input <- list(
  fetoinfant_lifetables = 'out/30-fetoinfant_lifetables.qs',
  competing_risk_model_fits = 'tmp/50-competing_risks_model_fits.qs',
  figure_specs = 'src/00-figure_specifications.R',
  config = 'cfg/config.yaml'
)
paths$output <- list(
  lifetables_and_parametric_fits_by_cod =
    'tmp/51-lifetables_and_parametric_fits_by_cod.qs'
)

# figure specs
source(paths$input$figure_specs)

config <- yaml::read_yaml(paths$input$config)

# constants
cnst <-
  list(
    cod = config$cod_lookup$key
  )

# Load data -------------------------------------------------------

fit <- qs_read(paths$input$competing_risk_model_fits)

# Assemble lifetables and fit for cause of death ------------------

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
lifetables_and_parametric_fits_by_cod %>%
  ggplot(aes(x = x)) +
  geom_area(aes(y = avg_total_hx*1e5, fill = cod)) +
  scale_fill_manual(values = rev(config$cod_lookup$color)) +
  scale_x_continuous(breaks = c(24, 40, 76)) +
  scale_y_continuous(expand = c(0,0)) +
  fig_spec$MyGGplotTheme(panel_border = TRUE) +
  labs(y = 'Feto-infant deaths per 100K person-weeks',
       x = 'Week of gestation',
       fill = 'Cause of death') +
  theme(legend.position = c(0.7, 0.5)) +
  geom_step(aes(x = x, y = m*1e5), size = 2,
            data = fit$total14 %>% unnest(lifetable)) +
  geom_line(aes(x = x, y = avg_total_hx*1e5), size = 2,
            data = fit$total14 %>% unnest(pred_summary))

# Export ----------------------------------------------------------

qs_save(
  lifetables_and_parametric_fits_by_cod,
  paths$output$lifetables_and_parametric_fits_by_cod
)
