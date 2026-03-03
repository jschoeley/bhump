
# Mortality compression towards birth -----------------------------

fit$origineducation <-
  filt$origineducation |>
  filter(!grepl('Unknown', stratum)) |>
  FitFetoinfantSurvival(
    control = ControlFitFetoinfantSurvival(simulate = FALSE)
  )
fig$hzrd_origineducation <- PlotHazards(
  fit$origineducation, colors = c(
    rep('black',7),
    'red', # white academic
    rep('black',1),
    'blue', # black primary
    rep('black',6)
  ), legend = FALSE)
fig$hzrd_origineducation

fig_spec$ExportPDF(
  fig$hzrd_origineducation,
  '50-hzrd-origineducation',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

fig$compression <-
  tibble(
    F76 = ProbFetoInfantDeath(fit$origineducation, x = 77) |>
      filter(!grepl('Unknown', stratum)) |> pull(avg_total_Fx),
    rho = BirthHumpDeaths(fit$origineducation, x = 77) |>
      filter(!grepl('Unknown', stratum)) |> pull(avg_p_birth)
  ) |>
  ggplot(aes(x = F76*100,y = rho*100)) +
  geom_smooth(method = 'lm', se = FALSE, color = 'grey40') +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  labs(y = '% feto-infant deaths contributed by birth-hump',
       x = '% probability of feto infant death') +
  fig_spec$MyGGplotTheme()
fig$compression
fig_spec$ExportPDF(
  fig$compression,
  '50-compression',
  'out',
  width = fig_spec$width*0.8,
  height = fig_spec$width*0.6
)
