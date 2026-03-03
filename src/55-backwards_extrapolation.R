# Backwards extrapolation -----------------------------------------

p <-
  fit$cohort |>
  filter(stratum == '2014') |>
  unnest(par_summary) |>
  pull(avg)

p2 <- RescaleParameters(p, model = 'basic', split = 14)
p2$zeta <- p2$zeta + 24
p2$alpha1 <- exp(log(p2$alpha1)-24*(-p2$beta1))

FetoinfantHzrd(0:(76+24), pars = p2, model = 'basic') |>
  plot()
FetoinfantSurv(0:(76+24), pars = p2, model = 'basic') |>
  plot()
