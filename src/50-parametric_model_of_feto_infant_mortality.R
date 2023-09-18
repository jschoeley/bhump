# A latent-competing risks model of the
# feto-infant mortality trajectory over
# age of gestation

# Init ------------------------------------------------------------

library(tidyverse)
library(maxLik)
library(cowplot)

paths <- list()
paths$input <- list(
  fetoinfant_lifetables = 'out/30-fetoinfant_lifetables.RData',
  figure_specs = 'src/00-figure_specifications.R',
  lifetable_functions = 'src/00-fnct-feto_infant_lt.R',
  parametric_functions = 'src/00-fnct-parametric_survival_model.R',
  config = 'cfg/config.yaml'
)
paths$output <- list(
  fetoinfant = 'tmp/21-fetoinfant.RData'
)

# figure specs
source(paths$input$figure_specs)
# fetoinfant lifetable functions
source(paths$input$lifetable_functions)
# fetoinfant parametric functions
source(paths$input$parametric_functions)

config <- yaml::read_yaml(paths$input$config)

# constants
cnst <-
  list(
    gestage_brk = seq(24, 77, by = 4),
    lifetable_breaks = 24:77,
    left_truncation_gestage = 24,
    right_censoring_gestage = 77
  )

# figures
fig <- list()
# tables
tab <- list()
# model fits
fit <- list()

# Data ------------------------------------------------------------

# a list of feto-infant lifetable as FILT objects
load(paths$input$fetoinfant_lifetables)

# Fit fetoinfant survival ---------------------------------------

# by sex
fit$sex09 <-
  FitFetoinfantSurvival(filt$sex09, simulate = TRUE, model = 'basic',
                        hessian_inverse = 'pseudo')
PlotHazards(fit$sex09)

# by cohort
fit$cohort <-
  FitFetoinfantSurvival(filt$cohort, model = 'basic',
                        hessian_inverse = 'pseudo')
PlotHazards(fit$cohort)

# by origin
fit$origin <-
  FitFetoinfantSurvival(filt$origin09, model = 'basic',
                        hessian_inverse = 'pseudo')
PlotHazards(fit$origin)

# by cause of death
fit$pcml_complications <-
  FitFetoinfantSurvival(filt$pcml_complications, model = 'basic', hessian_inverse = 'cholesky')
PlotHazards(fit$pcml_complications)
fit$untreatable_neoplasms <-
  FitFetoinfantSurvival(filt$untreatable_neoplasms, model = 'basic',
                        hessian_inverse = 'cholesky')
PlotHazards(fit$untreatable_neoplasms)
fit$treatable_neoplasms <-
  FitFetoinfantSurvival(filt$treatable_neoplasms, model = 'basic',
                        hessian_inverse = 'cholesky')
PlotHazards(fit$treatable_neoplasms)
fit$maternal <-
  FitFetoinfantSurvival(filt$maternal, model = 'flexible1',
                        hessian_inverse = 'cholesky', lambda2 = 1)
PlotHazards(fit$maternal)
fit$infections_and_parasites <-
  FitFetoinfantSurvival(
    filt$infections_and_parasites, model = 'basic',
    hessian_inverse = 'cholesky'
  )
PlotHazards(
  fit$infections_and_parasites,
  colors = config$cod_lookup$color[config$cod_lookup$key=='infections_and_parasites']
)
fit$prematurity <-
  FitFetoinfantSurvival(filt$prematurity, model = 'flexible1', lambda2 = 0.84)
PlotHazards(fit$prematurity)
fit$accidents_and_violence <-
  FitFetoinfantSurvival(filt$accidents_and_violence, model = 'flexible2',
                        simulate = TRUE)
PlotHazards(fit$accidents_and_violence)
fit$sids <-
  FitFetoinfantSurvival(filt$sids, simulate = FALSE,
                        model = 'flexible2')
PlotHazards(fit$sids)
fit$unspecific <-
  FitFetoinfantSurvival(filt$unspecific, model = 'basic')
PlotHazards(fit$unspecific)
fit$other <-
  FitFetoinfantSurvival(filt$other, model = 'basic', hessian_inverse = 'cholesky')
PlotHazards(fit$other)


# Plot hazards --------------------------------------------------

# by sex
fig$hzrd_sex09 <- PlotHazards(fit$sex09)
fig_spec$ExportPDF(
  fig$hzrd_sex09,
  '50-hzrd-sex09',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

# by cohort
fig$hzrd_cohort <- PlotHazards(fit$cohort)
fig_spec$ExportPDF(
  fig$hzrd_cohort,
  '50-hzrd-cohort',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

# by origin
fig$hzrd_origin <- PlotHazards(fit$origin)
fig_spec$ExportPDF(
  fig$hzrd_origin,
  '50-hzrd-origin',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

# by cause of death
fig$hzrd_pcml_complications <- PlotHazards(
  fit$pcml_complications, colors =
    config$cod_lookup$color[config$cod_lookup$key == 'pcml_complications']
)
fig_spec$ExportPDF(
  fig$hzrd_pcml_complications,
  '50-hzrd-pcml_complications',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)
fig$hzrd_untreatable_neoplasms <- PlotHazards(
  fit$untreatable_neoplasms, colors =
    config$cod_lookup$color[config$cod_lookup$key == 'untreatable_neoplasms']
)
fig_spec$ExportPDF(
  fig$hzrd_untreatable_neoplasms,
  '50-hzrd-untreatable_neoplasms',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)
fig$hzrd_treatable_neoplasms <- PlotHazards(
  fit$treatable_neoplasms, colors =
    config$cod_lookup$color[config$cod_lookup$key == 'treatable_neoplasms']
)
fig_spec$ExportPDF(
  fig$hzrd_treatable_neoplasms,
  '50-hzrd-treatable_neoplasms',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)
fig$hzrd_maternal <- PlotHazards(
  fit$maternal,
  colors = config$cod_lookup$color[config$cod_lookup$key == 'maternal']
)
fig_spec$ExportPDF(
  fig$hzrd_maternal,
  '50-hzrd-maternal',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)
fig$hzrd_infections_and_parasites <-
  PlotHazards(
    fit$infections_and_parasites,
    colors = config$cod_lookup$color[config$cod_lookup$key == 'infections_and_parasites']
  )
fig_spec$ExportPDF(
  fig$hzrd_infections_and_parasites,
  '50-hzrd-infections_and_parasites',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)
fig$hzrd_prematurity <- PlotHazards(
  fit$prematurity, colors =
    config$cod_lookup$color[config$cod_lookup$key == 'prematurity']
)
fig_spec$ExportPDF(
  fig$hzrd_prematurity,
  '50-hzrd-prematurity',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)
fig$hzrd_accidents_and_violence <- PlotHazards(
  fit$accidents_and_violence, colors =
    config$cod_lookup$color[config$cod_lookup$key == 'accidents_and_violence']
)
fig_spec$ExportPDF(
  fig$hzrd_accidents_and_violence,
  '50-hzrd-accidents_and_violence',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)
fig$hzrd_sids <- PlotHazards(
  fit$sids,colors =
    config$cod_lookup$color[config$cod_lookup$key == 'sids']
)
fig_spec$ExportPDF(
  fig$hzrd_sids,
  '50-hzrd-sids',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)
fig$hzrd_unspecific_stillbirth <- PlotHazards(
  fit$unspecific, colors =
    config$cod_lookup$color[config$cod_lookup$key == 'unspecific']
)
fig_spec$ExportPDF(
  fig$hzrd_unspecific_stillbirth,
  '50-hzrd-unspecific_stillbirth',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)
fig$hzrd_other <- PlotHazards(
  fit$other, colors =
    config$cod_lookup$color[config$cod_lookup$key == 'other']
)
fig_spec$ExportPDF(
  fig$hzrd_other,
  '50-hzrd-other',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

# Share of deaths due to birth hump -----------------------------

# by sex
BirthHumpDeaths(fit$sex09, x = 77)
# by cohort
BirthHumpDeaths(fit$cohort, x = 77)
# by origin
BirthHumpDeaths(fit$origin, x = 77)

# by cause of death
BirthHumpDeaths(fit$pcml_complications, x = 77)
BirthHumpDeaths(fit$untreatable_neoplasms, x = 77)
BirthHumpDeaths(fit$treatable_neoplasms, x = 77)
BirthHumpDeaths(fit$maternal, x = 77)
BirthHumpDeaths(fit$infections_and_parasites, x = 77)
BirthHumpDeaths(fit$prematurity, x = 77)
BirthHumpDeaths(fit$accidents_and_violence, x = 77)
BirthHumpDeaths(fit$sids, x = 77)
BirthHumpDeaths(fit$unspecific, x = 77)
BirthHumpDeaths(fit$other, x = 77)

# Plot transitional component -------------------------------------

cod <- c('maternal', 'treatable_neoplasms','untreatable_neoplasms',
         'accidents_and_violence', 'unspecific',
         'pcml_complications', 'prematurity', 'sids','infections_and_parasites', 'other')

filt_cod <- map(cod, ~{
  fit[[.x]]$pred_summary[[1]] %>%
    mutate(cod = .x)
}) %>% bind_rows() %>%
  group_by(cod) %>%
  mutate(
    max = max(avg_birth_hx, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    cod = factor(cod, rev(config$cod_lookup$key), rev(config$cod_lookup$label))
  )

fig$birthhump_cod_joint <-
  filt_cod %>%
  ggplot(aes(x = x)) +
  geom_area(aes(y = avg_birth_hx*1e5, fill = cod)) +
  scale_fill_manual(values = rev(config$cod_lookup$color)) +
  scale_x_continuous(breaks = c(24, 40, 76)) +
  scale_y_continuous(expand = c(0,0)) +
  fig_spec$MyGGplotTheme(panel_border = TRUE) +
  labs(y = 'Feto-infant deaths per 100K person-weeks',
       x = 'Week of gestation',
       fill = 'Cause of death') +
  theme(legend.position = c(0.7, 0.5))
fig$birthhump_cod_joint

fig_spec$ExportPDF(
  fig$birthhump_cod_joint,
  '50-birthhump_cod_joint',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.6
)

fig$birthhump_cod_separate <-
  filt_cod %>%
  mutate(
    cod = factor(cod, config$cod_lookup$label, config$cod_lookup$label)
  ) %>%
  ggplot(aes(x = x)) +
  geom_area(aes(y = avg_birth_hx*1e5, fill = cod)) +
  scale_fill_manual(values = config$cod_lookup$color) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(24, 40, 76)) +
  facet_wrap(~cod, nrow = 2) +
  fig_spec$MyGGplotTheme(panel_border = TRUE) +
  labs(y = 'Feto-infant deaths per 100K person-weeks',
       x = 'Week of gestation',
       fill = 'Cause of death') +
  guides(fill = 'none')
fig$birthhump_cod_separate

fig_spec$ExportPDF(
  fig$birthhump_cod_separate,
  '50-birthhump_cod_separate',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.6
)

# Probability of Fetoinfant death ---------------------------------

# by sex
ProbFetoInfantDeath(fit$sex09)
# by cohort
ProbFetoInfantDeath(fit$cohort)
# by origin
ProbFetoInfantDeath(fit$origin)

# by cohort
# by cause of death

bind_rows(
  pcml_complications = ProbFetoInfantDeath(fit$pcml_complications),
  untreatable_neoplasms = ProbFetoInfantDeath(fit$untreatable_neoplasms),
  treatable_neoplasms = ProbFetoInfantDeath(fit$treatable_neoplasms),
  maternal = ProbFetoInfantDeath(fit$maternal),
  infections_and_parasites = ProbFetoInfantDeath(fit$infections_and_parasites),
  prematurity = ProbFetoInfantDeath(fit$prematurity),
  accidents_and_violence = ProbFetoInfantDeath(fit$accidents_and_violence),
  sids = ProbFetoInfantDeath(fit$sids),
  unspecific_stillbirth = ProbFetoInfantDeath(fit$unspecific),
  other = ProbFetoInfantDeath(fit$other),
  .id = 'cod'
) %>%
  mutate(p_birth = avg_birth_Fx / sum(avg_birth_Fx)*100) %>%
  arrange(p_birth) %>%
  select(cod, p_birth)

# Parameter tables --------------------------------------------------------

PrintParameterTable <- function (filt_fit) {

    tab_of_pars <-
      filt_fit %>%
      unnest_legacy(par_summary) %>%
      transmute(
        name,
        stratum,
        avg = exp(avg),
        ci025 = exp(ci025),
        ci975 = exp(ci975)
      ) %>%
      mutate_at(
        c('avg', 'ci025', 'ci975'),
        ~ formatC(., format = 'e', digits = 1)
      )
      

    return(tab_of_pars)

}

# by sex
PrintParameterTable(fit$sex09)
# by cohort
PrintParameterTable(fit$cohort)
# by age
PrintParameterTable(fit$origin)

# Parametric decomposition of mortality differences -------------

DecomposeFetoInfantDeaths <-
  function (filt_fit, pop1, pop2) {

    require(DemoDecomp)
    require(rlang)

    FetoinfantDeaths <-
      function (pars) { (1-FetoinfantSurv(76-24, pars))*1e5 }

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
        {{pop1}} := FetoinfantDeaths({{pop1}}),
        {{pop2}} := FetoinfantDeaths({{pop2}}),
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
            N = 1e2
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
          level = .x %>% slice(1) %>% pull('contribution'),
          ontog = .x %>% slice(2) %>% pull('contribution'),
          trans = .x %>% slice(3:5) %>% pull('contribution') %>% sum()
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

# by sex
DecomposeFetoInfantDeaths(fit$sex09, Female, Male)

# by cohort
DecomposeFetoInfantDeaths(fit$cohort, `1989`, `1999`)
DecomposeFetoInfantDeaths(fit$cohort, `1999`, `2009`)

# by origin
DecomposeFetoInfantDeaths(fit$origin, `Non-Hispanic White`, `Non-Hispanic Black`)
DecomposeFetoInfantDeaths(fit$origin, `Non-Hispanic White`, `Hispanic`)
DecomposeFetoInfantDeaths(fit$origin, `Non-Hispanic Black`, `Hispanic`)
