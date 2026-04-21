# Decompose differences in feto-infant mortality between groups
# into level, ontogenescence, and transition components

# The decomposition is performed via the horiuchi method. After
# fitting the competing risks hazard model we evaluate the predicted
# difference in overall feto-infant probability of death from
# week 24 through 76 for two groups. We then attribute this difference
# to differences in the fitted model parameters, distinguishing between
# level, ontogenescent (slope), and birth-hump contributions.

# Init ------------------------------------------------------------

library(tidyverse)
library(qs2)

paths <- list()
paths$input <- list(
  competing_risk_model_fits = 'tmp/50-competing_risks_model_fits.qs',
  parametric_functions = 'src/00-fnct-parametric_survival_model.R',
  config = 'cfg/config.yaml'
)
paths$output <- list(
  parametric_decompositions =
    'out/54-parametric_decompositions.qs'
)

config <- yaml::read_yaml(paths$input$config)

# fetoinfant parametric functions
source(paths$input$parametric_functions)

# constants
cnst <-
  list(
    cod = config$cod_lookup$key
  )

# Load data -------------------------------------------------------

# a list of feto-infant lifetable as FILT objects
fit <- qs_read(paths$input$competing_risk_model_fits)

# Parametric decomposition of mortality differences -------------

DecomposeFetoInfantDeaths <-
  function (filt_fit, pop1, pop2) {
    
    require(DemoDecomp)
    require(rlang)
    
    FetoinfantDeaths <-
      function (pars, control) {
        pars_rescaled <- RescaleParameters(
          pars, control$model, control$split,
          zeta_range = control$zeta_range,
          beta1_range = control$beta1_range,
          beta2_range = control$beta2_range
        )
        (1-FetoinfantSurv(76-24, pars_rescaled))*1e5
      }
    
    control <- filt_fit$control[[1]]
    
    total_difference <-
      filt_fit %>%
      unnest_legacy(par_draws) %>%
      select(draw, stratum, name, value) %>%
      filter(
        stratum %in%
          c(as_name(enquo(pop1)), as_name(enquo(pop2)))
      ) %>%
      spread(stratum, value) %>%
      group_by(draw) %>%
      summarise(
        {{pop1}} := FetoinfantDeaths({{pop1}}, control),
        {{pop2}} := FetoinfantDeaths({{pop2}}, control),
        diff = {{pop2}}-{{pop1}},
        reldiff = diff/{{pop1}}
      ) %>%
      ungroup() %>%
      summarise(
        {{pop1}} := mean({{pop1}}),
        {{pop2}} := mean({{pop2}}),
        diff_avg = mean(diff),
        diff_se = sd(diff),
        diff_q025 = quantile(diff, 0.025),
        diff_q975 = quantile(diff, 0.975),
        reldiff_avg = mean(reldiff),
        reldiff_se = sd(reldiff),
        reldiff_q025 = quantile(reldiff, 0.025),
        reldiff_q975 = quantile(reldiff, 0.975)
      ) %>%
      pivot_longer(everything())
    
    parameter_decomp_draw <-
      filt_fit %>%
      unnest_legacy(par_draws) %>%
      select(draw, stratum, name, value) %>%
      filter(
        stratum %in% c(as_name(enquo(pop1)), as_name(enquo(pop2)))
      ) %>%
      spread(stratum, value) %>%
      group_by(draw) %>%
      mutate(
        contribution =
          horiuchi(
            FetoinfantDeaths,
            pars1 = {{pop1}},
            pars2 = {{pop2}},
            N = 1e2,
            control
          )
      )
    
    parameter_decomp_summary <-
      parameter_decomp_draw %>%
      group_by(name) %>%
      summarise(
        '{{pop1}}_avg' :=
          mean({{pop1}}),
        '{{pop2}}_avg' :=
          mean({{pop2}}),
        '{{pop1}}_se' := sd({{pop1}}),
        '{{pop2}}_se' := sd({{pop2}}),
        '{{pop1}}_ci025' := quantile({{pop1}}, 0.025),
        '{{pop2}}_ci975' := quantile({{pop1}}, 0.975),
        contribution_avg = mean(contribution),
        contribution_se = sd(contribution),
        contribution_ci025 = quantile(contribution, 0.025),
        contribution_ci975 = quantile(contribution, 0.975),
      )
    
    component_decomp <-
      parameter_decomp_draw %>%
      group_by(draw) %>%
      group_modify(~{
        tibble(
          level =
            .x %>%
            filter(grepl('alpha', name)) %>%
            pull('contribution'),
          ontog = .x %>%
            filter(grepl('beta', name)) %>%
            pull('contribution'),
          trans =
            .x %>%
            filter(grepl('gamma|tau|zeta|sigma', name)) %>%
            pull('contribution') %>% sum()
        )  
      }) %>%
      ungroup() %>%
      summarise(
        level_avg = mean(level),
        level_se = sd(level),
        level_ci025 = quantile(level, 0.025),
        level_ci975 = quantile(level, 0.975),
        ontog_avg = mean(ontog),
        ontog_se = sd(ontog),
        ontog_ci025 = quantile(ontog, 0.025),
        ontog_ci975 = quantile(ontog, 0.975),
        trans_avg = mean(trans),
        trans_se = sd(trans),
        trans_ci025 = quantile(trans, 0.025),
        trans_ci975 = quantile(trans, 0.975),
      ) %>%
      pivot_longer(everything())
    
    list(
      diff = total_difference,
      para = parameter_decomp_summary,
      comp = component_decomp
    )
    
  }

decomp <- list(
  sex = list(),
  cohort = list(),
  origin = list(),
  education = list()
)

# by sex
decomp$sex$female_male <- DecomposeFetoInfantDeaths(fit$sex, Female, Male)

# by cohort
decomp$cohort$`89vs99` <-
  DecomposeFetoInfantDeaths(fit$cohort, `1989`, `1999`)
decomp$cohort$`99vs09` <-
  DecomposeFetoInfantDeaths(fit$cohort, `1999`, `2009`)
decomp$cohort$`09vs14` <-
  DecomposeFetoInfantDeaths(fit$cohort, `2009`, `2014`)
decomp$cohort$`89vs14` <-
  DecomposeFetoInfantDeaths(fit$cohort, `1989`, `2014`)

# by origin
decomp$origin$white_black <-
  DecomposeFetoInfantDeaths(fit$origin, `Non-Hispanic White`, `Non-Hispanic Black`)
decomp$origin$white_hispanic <-
  DecomposeFetoInfantDeaths(fit$origin, `Non-Hispanic White`, `Hispanic`)
decomp$origin$black_hispanic <-
  DecomposeFetoInfantDeaths(fit$origin, `Non-Hispanic Black`, `Hispanic`)

# by education
decomp$education$primary_academic <-
  DecomposeFetoInfantDeaths(fit$education, `Primary`, `Bachelor, Master, Doctorate`)

# Export ----------------------------------------------------------

qs_save(decomp, paths$output$parametric_decompositions)
