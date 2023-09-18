# Scaling functions -----------------------------------------------

# [a,b] -> [-Inf, Inf]
ScaleLogit <- function(x, a, b) {
  x_norm <- (x - a) / (b - a)
  logit_x <- log(x_norm / (1 - x_norm))
  
  return(logit_x)
}

# [-Inf,Inf] -> [a,b]
ScaleInverseLogit <- function(logit_x, a, b) {
  x_norm <- 1 / (1 + exp(-logit_x))
  # scale x back to [a, b]
  x <- a + x_norm * (b - a)
  
  return(x)
}

RescaleParameters <- function (
    pars, model, split,
    c_range = c(10, 20), b_range = c(-0.7, 0.7), beta_range = c(-0.7, 0.7)
) {
  rescaled_parameters <- switch(
    model,
    basic = function () {
      p <- list()
      p$a1 = exp(pars[1])
      p$b = exp(pars[2])
      p$a2 = exp(pars[3])
      p$c = ScaleInverseLogit(pars[4], c_range[1], c_range[2])
      p$s = exp(pars[5])
      p$alpha = exp(pars[1]-p$c*p$b)
      p$beta = p$b
      return(p)
    },
    flexible1 = function () {
      p <- list()
      p$a1 = exp(pars[1])
      p$b = exp(pars[2])
      p$a2 = exp(pars[5])
      p$c = ScaleInverseLogit(pars[6], c_range[1], c_range[2])
      p$s = exp(pars[7])
      p$alpha = exp(pars[3])
      p$beta = exp(pars[4])
      return(p)
    },
    flexible2 = function () {
      p <- list()
      p$a1 = exp(pars[1])
      p$b = ScaleInverseLogit(pars[2], b_range[1], b_range[2])
      p$a2 = exp(pars[5])
      p$c = ScaleInverseLogit(pars[6], c_range[1], c_range[2])
      p$s = exp(pars[7])
      p$alpha = exp(pars[3])
      p$beta = ScaleInverseLogit(pars[4], beta_range[1], beta_range[2])
      return(p)
    }
  )
  rescaled_parameters()
}

# Parameter Initialization ----------------------------------------

InitializeParameters <- function (m, model = 'basic', split = 14) {
  initial_scaled_parameters <- switch(
    model,
    basic = function () {
      p <- list()
      # log hazard at t=0
      p$a1 = log(ifelse(m[1] == 0, mean(m), m[1]))
      # log relative rate of ontogenescent mortality decline
      #p$b = log(abs(mean(diff(log(m)))))
      p$b = (
        abs(
          log(mean(tail(m, 5), na.rm = TRUE)) -
            log(mean(head(m, 5), na.rm = TRUE))
        )
      ) / (length(m)-2.5 - 2.5)
      p$b = log(p$b)
      if (is.infinite(p$b) | is.na(p$b)) { p$b = log(0.1) }
      # log hazard at split point
      p$alpha = p$a1-exp(p$b)*split
      # log relative rate of decline after split
      p$beta = p$b
      # log hazard contribution of birth hump at peak
      p$a2 = p$a1 - log(2)
      # location of birth component (0 gets rescaled to 39 weeks)
      p$c = 0
      # log spread of birth component
      p$s = log(1)
      return(p)
    },
    flexible1 = function () {
      p <- list()
      # log hazard at t=0
      p$a1 = log(ifelse(m[1] == 0, mean(m), m[1]))
      # log relative rate of ontogenescent mortality decline
      p$b = log(abs(mean(diff(log(m)))))
      if (is.infinite(p$b) | is.na(p$b)) { p$b = log(0.1) }
      # log hazard at split point
      p$alpha = p$a1-exp(p$b)*split
      # log relative rate of decline after split
      p$beta = p$b
      # log hazard contribution of birth hump at peak
      p$a2 = p$a1 - log(2)
      # location of birth component (0 gets rescaled to 39 weeks)
      p$c = 0
      # log spread of birth component
      p$s = log(1)
      return(p)
    },
    flexible2 = function () {
      p <- list()
      # log hazard at t=0
      p$a1 = log(mean(m))
      # ScaledLogit relative rate of ontogenescent mortality decline
      p$b = mean(m)
      # log hazard at split point
      p$alpha = p$a1-p$b*split
      # ScaledLogit relative rate of decline after split
      p$beta = 1e-9
      # log hazard contribution of birth hump at peak
      p$a2 = p$a1 - log(2)
      # location of birth component (0 gets rescaled to 39 weeks)
      p$c = 0
      # log spread of birth component
      p$s = log(1)
      return(p)
    }
  )
  
  init_pars_pre <- initial_scaled_parameters()
  
  init_pars_ordered <- switch(
    model,
    basic = c(
      '1:a1' = init_pars_pre[['a1']], '2:b' = init_pars_pre[['b']],
      '3:a2' = init_pars_pre[['a2']], '4:c' = init_pars_pre[['c']],
      '5:s' = init_pars_pre[['s']]
    ),
    flexible1 = c(
      '1:a1' = init_pars_pre[['a1']], '2:b' = init_pars_pre[['b']],
      '3:alpha' = init_pars_pre[['alpha']], '4:beta' = init_pars_pre[['beta']],
      '5:a2' = init_pars_pre[['a2']], '6:c' = init_pars_pre[['c']], '7:s' = init_pars_pre[['s']]
    ),
    flexible2 = c(
      '1:a1' = init_pars_pre[['a1']], '2:b' = init_pars_pre[['b']],
      '3:alpha' = init_pars_pre[['alpha']], '4:beta' = init_pars_pre[['beta']],
      '5:a2' = init_pars_pre[['a2']], '6:c' = init_pars_pre[['c']], '7:s' = init_pars_pre[['s']]
    )
  )
  
  return(init_pars_ordered)
  
}

# Feto-infant parametric survival -------------------------------

# baseline characteristics
NegGompertzHzrd <-
  function (x, a, b, alpha, beta, split) {
    I1 = ifelse(x < split, 1, 0)
    I2 = ifelse(x >= split, 1, 0)
    a*exp(-b*x)*I1 + alpha*exp(-beta*(x-split))*I2
  }
NegGompertzCumHzrd <-
  function (x, a, b, alpha, beta, split) {
    I1 = ifelse(x < split, 1, 0)
    I2 = ifelse(x >= split, 1, 0)
    Ssplit = ((a - a*exp(-b*split)) / b)
    I1*((a - a*exp(-b*x)) / b) +
      I2*((alpha - alpha*exp(-beta*(x-split))) / beta + Ssplit)
  }
NegGompertzSurv <-
  function (x, a, b, alpha, beta, split) {
    exp(-NegGompertzCumHzrd(x, a, b, alpha, beta, split))
  }

# birth hump characteristics
GaussianHzrd <-
  function (x, a, c, s) {
    a*exp(-(x-c)^2/(2*s^2))
  }
GaussianCumHzrd <-
  function (x, a, c, s) {
    denom <- sqrt(2)*s
    Hx <- a*sqrt(pi/2)*s*(pracma::erf(c/denom) + pracma::erf((x-c)/denom))
    return(Hx)
  }
GaussianSurv <-
  function (x, a, c, s) {
    exp(-GaussianCumHzrd(x, a, c, s))
  }

# feto-infant competing risks survival
FetoinfantSurv <-
  function (x, pars, split, component = 'total', model = 'basic') {
    
    # ontogenescent survival
    ontogen_surv <-
      NegGompertzSurv(
        x = x,
        a = pars[['a1']], b = pars[['b']],
        alpha = pars[['alpha']], beta = pars[['beta']],
        split = pars[['c']]
      )
    # birth survival
    birth_surv <-
      GaussianSurv(
        x = x,
        a = pars[['a2']], c = pars[['c']], s = pars[['s']]
      )
    
    Sx <-
      switch (component,
              'total' = ontogen_surv * birth_surv,
              'ontogen' = ontogen_surv,
              'birth' = birth_surv,
              stop('Component must be one of "total", "ontogen" or "birth"')
      )
    return(Sx)
  }

# feto-infant competing risks hazard
FetoinfantHzrd <-
  function (x, pars, split, component = 'total', model = 'basic') {
    
    # ontogenescent hazard
    ontogen_hzrd <-
      NegGompertzHzrd(
        x = x,
        a = pars[['a1']], b = pars[['b']],
        alpha = pars[['alpha']], beta = pars[['beta']],
        split = pars[['c']]
      )
    
    # birth hazard
    birth_hzrd <-
      GaussianHzrd(
        x = x,
        a = pars[['a2']], c = pars[['c']], s = pars[['s']]
      )
    
    hx <-
      switch (component,
              'total' = ontogen_hzrd + birth_hzrd,
              'ontogen' = ontogen_hzrd,
              'birth' = birth_hzrd,
              stop('Component must be one of "total", "ontogen" or "birth"')
      )
    return(hx)
    
  }

# Objective function --------------------------------------------

# interval censored likelihood
IntervalCensoredLogLike <-
  function (pars, age, width, obsDx, obsCx, SurvFnct,
            lambda1 = 0, lambda2 = 0, lambda3 = 0, split, model = 'basic', ...) {
    
    pars2 <- RescaleParameters(pars, model, split)
    cat(unlist(pars2), '\n')
    
    # predict survival on basis of parameter estimates
    predSurvL <- SurvFnct(x = age, pars = pars2, split = split,
                          component = 'total', model = model)
    predSurvR <- SurvFnct(x = age+width, pars = pars2, split = split,
                          component = 'total', model = model)
    
    loglike <-
      obsDx*log(predSurvL-predSurvR) + obsCx*log(predSurvR)
    
    penalty <- 0
    if (isTRUE(model == 'basic')) {
      # penalize birth hump magnitude
      penalty <- lambda1*exp(pars[3])
    }
    if (isTRUE(model == 'flexible1')) {
      penalty <-
        # penalize birth hump magnitude
        penalty <- lambda1*exp(pars[5]) +
        # penalize discontinuities between the two ontogenescent segments
        lambda2*(pars[3] - (pars[1]-exp(pars[2])*ScaleInverseLogit(pars[6], 10, 20)))^2 +
        # penalize differences in slope between the two ontogenescent segments
        lambda3*(pars[4]-pars[2])^2
    }
    if (isTRUE(model == 'flexible2')) {
      penalty <-
        # penalize birth hump magnitude
        penalty <- lambda1*exp(pars[5]) +
        # penalize discontinuities between the two ontogenescent segments
        lambda2*(pars[3] - (pars[1]-pars[2]*ScaleInverseLogit(pars[6], 10, 20)))^2 +
        # penalize differences in slope between the two ontogenescent segments
        lambda3*(pars[4]-pars[2])^2
    }
    
    ploglike <- loglike - penalty
    
    return(ploglike)
    
  }

# Fit -------------------------------------------------------------

# Fit a parametric model of fetoinfant survival to
# the fetoinfant lifetable assuming left truncation age
# at 24 and right censoring at 77. Derive various predictions.
# Sample parameters from multivariate normal distribution
# derived from hessian in order to calculate CIs.
FitFetoinfantSurvival <-
  function (
    filt, lambda1 = 0, lambda2 = 0, lambda3 = 0, split = 15,
    model = 'basic',
    simulate = TRUE, hessian_inverse = 'cholesky'
  ) {
    
    require(maxLik)
    require(pracma)
    
    stopifnot(any(class(filt) == 'FILT'))
    
    fit_lifetable <-
      filt %>%
      left_join(FILTMortalityRates(.)) %>%
      left_join(FILTSurvival(.)) %>%
      group_by(stratum) %>%
      group_modify(~{
        
        init_pars <- InitializeParameters(.x$m, model, split)
        
        modelfit <-
          maxLik(
            logLik = IntervalCensoredLogLike,
            start  = init_pars,
            method = 'CG',
            # data and arguments to objective function
            age =
              pull(.x, x)-24,
            width =
              pull(.x, n),
            obsDx =
              pull(.x, D),
            obsCx =
              pull(.x, C),
            # hazard function
            SurvFnct =
              FetoinfantSurv,
            lambda1 = lambda1,
            lambda2 = lambda2,
            lambda3 = lambda3,
            split = split,
            model = model,
            # options
            iterlim = 1e4
          )
        
        # 1000 draws from the posterior parameter distribution
        # assuming multivariate normal derived from hessian
        par_draw <-
          expand_grid(
            draw = 1:1000,
            name = names(init_pars)
          )
        if (isTRUE(simulate)) {
          hessian_inverse <- 
            switch(hessian_inverse,
                   cholesky = chol2inv(-modelfit$hessian),
                   pseudo = pinv(-modelfit$hessian))
          par_draw$value <- MASS::mvrnorm(
            n = 1e3,
            mu = modelfit$estimate,
            Sigma = hessian_inverse
          ) %>% t() %>% c()
        } else {
          par_draw$value <- rep(modelfit$estimate, 1e3)
        }
        
        # summarise parameter draws into
        # point estimates and credible intervals
        pars <-
          par_draw %>%
          group_by(name) %>%
          summarise(
            avg = mean(value),
            se = sd(value),
            ci025 = quantile(value, 0.025),
            ci975 = quantile(value, 0.975)
          )
        
        # number of weeks from observation
        # start to end
        omega <-
          cnst$right_censoring_gestage -
          cnst$left_truncation_gestage
        
        # predictions by posterior draw
        pred_draw <-
          par_draw %>%
          group_by(draw) %>%
          group_modify(~{
            pars <- RescaleParameters(.x$value, model, split)
            tibble(
              # weeks since left truncation age
              x =
                seq(0, omega, length.out = 1000),
              # width of age interval
              n = omega/1000,
              # probability of fetoinfant survival until x
              total_lx =
                FetoinfantSurv(
                  x,
                  pars = pars,
                  component = 'total',
                  model = model
                ),
              # probability of fetoinfant death until x
              total_Fx =
                1-total_lx,
              # probability of fetoinfant death until x
              # (one in x won't survive)
              total_iFx =
                1/total_Fx,
              # hazard of fetoinfant death at x
              total_hx =
                FetoinfantHzrd(
                  x,
                  pars = pars,
                  component = 'total',
                  model = model
                ),
              # hazard of birth component at x
              birth_hx =
                FetoinfantHzrd(
                  x,
                  pars = pars,
                  component = 'birth',
                  model = model
                ),
              # hazard of ontogenescent component at x
              ontogen_hx =
                FetoinfantHzrd(
                  x,
                  pars = pars,
                  component = 'ontogen',
                  model = model
                ),
              # cumulative probability of death due
              # to birth component
              birth_Fx =
                cumsum(total_lx*birth_hx*n),
              # cumulative probability of death due
              # to ontogenescent component
              ontogen_Fx =
                cumsum(total_lx*ontogen_hx*n),
              # share of deaths due to birth component
              p_birth =
                birth_Fx/total_Fx
            ) %>%
              mutate(
                # transform back to gestational age
                x = x + cnst$left_truncation_gestage
              )
          }) %>%
          ungroup()
        
        # summarize predictions over draws
        pred <-
          pred_draw %>%
          pivot_longer(
            -c('draw', 'x', 'n')
          ) %>%
          group_by(name, x, n) %>%
          summarise(
            avg = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE),
            q025 = quantile(value, 0.025, na.rm = TRUE),
            q975 = quantile(value, 0.975, na.rm = TRUE)
          ) %>%
          pivot_wider(
            names_from = 'name',
            values_from = c('avg', 'se', 'q025', 'q975')
          ) %>%
          ungroup()
        
        tibble(
          convergence = modelfit$code,
          loglike = maxValue(modelfit),
          par_draws = list(par_draw),
          par_summary = list(pars),
          hessian = list(modelfit$hessian),
          model = list(modelfit),
          pred_draws = list(pred_draw),
          pred_summary = list(pred),
          lifetable = list(.x)
        )
        
      }) %>%
      ungroup()
    
    fit_lifetable
    
  }

# Plotting functions ----------------------------------------------

# This function requires the output of FitFetoinfantSurvival()
# and plots lifetable estimates of fetoinfant survival and mortality
# versus the parametric model estimates.
PlotHazards <-
  function(
    filt_fit,
    ylab,
    xbrk = cnst$gestage_brk,
    xlim =
      c(cnst$left_truncation_gestage,
        cnst$right_censoring_gestage-0.1),
    scaler = 1e5,
    legend = TRUE,
    colors = fig_spec$discrete_colors
  ) {
    
    # aspect ratio
    ar <- 0.85
    
    lifetables <-
      filt_fit %>% unnest_legacy(lifetable)
    predictions <-
      filt_fit %>% unnest_legacy(pred_summary)
    
    plot_hzrd <-
      lifetables %>%
      ggplot() +
      geom_point(
        aes(
          x = x+0.5*n, y = m*scaler, group = stratum,
          color = as.character(stratum)
        ),
        size = fig_spec$point_size_m,
        alpha = 0.2
      ) +
      geom_line(
        aes(
          x = x, y = avg_total_hx*scaler, group = stratum,
          color = as.character(stratum)
        ),
        size = fig_spec$line_size_m,
        data = predictions
      ) +
      geom_ribbon(
        aes(
          x = x,
          ymin = q025_total_hx*scaler,
          ymax = q975_total_hx*scaler,
          group = stratum,
          fill = as.character(stratum)
        ),
        alpha = 0.2,
        data = predictions
      ) +
      geom_vline(
        aes(xintercept = 40),
        lty = 3,
        size = fig_spec$line_size_m
      ) +
      scale_y_continuous(
        paste0('Feto-infant deaths per\n',
               formatC(scaler, format = 'd', big.mark = ','),
               ' person-weeks at risk'),
        breaks = c(seq(2, 10, 2), seq(20, 80, 20)),
        trans = 'log10'
      ) +
      scale_x_continuous(
        'Week of gestation',
        breaks = xbrk,
        limits = xlim
      ) +
      scale_color_manual(
        '',
        values = colors
      ) +
      scale_fill_manual(
        '',
        values = colors
      ) +
      fig_spec$MyGGplotTheme(ar = ar) +
      theme(
        legend.position = ifelse(legend, c(0.8, 0.8), 'none'),
        legend.background = element_blank(),
        legend.box.background = element_blank()
      )
    
    plot_surv <-
      lifetables %>%
      ggplot() +
      geom_vline(
        aes(xintercept = 40),
        lty = 3,
        size = fig_spec$line_size_m
      ) +
      geom_point(
        aes(
          x = x, y = F_empirical*scaler,
          color = as.character(stratum)
        ),
        alpha = 0.2,
        size = fig_spec$point_size_m
      ) +
      geom_ribbon(
        aes(
          x = x,
          ymin = (1-q025_total_lx)*scaler,
          ymax = (1-q975_total_lx)*scaler,
          fill = as.character(stratum),
        ),
        alpha = 0.2,
        data = predictions
      ) +
      geom_line(
        aes(
          x = x, y = (1-avg_total_lx)*scaler,
          color = as.character(stratum),
        ),
        size = fig_spec$line_size_m,
        data = predictions
      ) +
      scale_y_continuous(
        paste0('Cumulative feto-infant deaths\n',
               'out of ', formatC(scaler, format = 'd',
                                  big.mark = ','),
               ' cohort members')
      ) +
      scale_color_manual(values = colors) +
      scale_x_continuous(
        '', breaks = xbrk
      ) +
      scale_fill_manual(
        '',
        values = colors
      ) +
      fig_spec$MyGGplotTheme(ar = ar) +
      theme(
        legend.position = 'none'
      )
    
    cowplot::plot_grid(plot_hzrd, plot_surv, ncol = 2, align = 'h')
    
  }

# Birth hump ------------------------------------------------------

# share of feto-infant deaths until x contributed by birth hump
# component
BirthHumpDeaths <- function (filt_fit, x) {
  
  filt_fit %>%
    unnest_legacy(pred_summary) %>%
    filter(x == {{x}}) %>%
    select(stratum, avg_p_birth,
           se_p_birth,
           q025_p_birth, q975_p_birth) %>%
    ungroup()
}

# probability of feto-infant death
ProbFetoInfantDeath <-
  function (filt_fit, x) {
    
    filt_fit %>%
      unnest_legacy(pred_summary) %>%
      filter(x == 77) %>%
      select(contains(c(
        'stratum',
        'total_Fx', 'total_iFx',
        'birth_Fx', 'birth_iFx'
      )))
    
  }