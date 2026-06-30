# Parametric decomposition ------------------------------------------------

ConstructCITabEntry <- function (abc) {
  NumberFormat <- function (x) {
    formatC(x, drop0trailing = FALSE, format = 'f', 1)
  }
  paste0(
    NumberFormat(abc[1]), ' (',
    NumberFormat(abc[2]), ', ',
    NumberFormat(abc[3]), ')'
  )
}

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
            pull('contribution') %>% sum(),
          ontog =
            .x %>%
            filter(grepl('beta', name)) %>%
            pull('contribution') %>% sum(),
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
    
    table_decomp <- data_frame(
      # names
      names = c(
        'Level contribution',
        'Ontogenescent contribution',
        'Birth hump contribution',
        'Delta F(52)',
        'Delta F(52)%'
      ),
      values = c(
        component_decomp |>
          filter(name %in% c('level_avg', 'level_ci025', 'level_ci975')) |>
          pivot_wider() |>
          unlist() |>
          ConstructCITabEntry(),
        component_decomp |>
          filter(name %in% c('ontog_avg', 'ontog_ci025', 'ontog_ci975')) |>
          pivot_wider() |>
          unlist() |>
          ConstructCITabEntry(),
        component_decomp |>
          filter(name %in% c('trans_avg', 'trans_ci025', 'trans_ci975')) |>
          pivot_wider() |>
          unlist() |>
          ConstructCITabEntry(),
        total_difference |>
          filter(name %in% c('diff_avg', 'diff_q025', 'diff_q975')) |>
          pivot_wider() |>
          unlist() |>
          ConstructCITabEntry(),
        total_difference |>
          filter(name %in% c('reldiff_avg', 'reldiff_q025', 'reldiff_q975')) |>
          mutate(value = value*100) |>
          pivot_wider() |>
          unlist() |>
          ConstructCITabEntry()
      )
    )
    
    list(
      diff = total_difference,
      para = parameter_decomp_summary,
      comp = component_decomp,
      tabl = table_decomp
    )
    
  }
