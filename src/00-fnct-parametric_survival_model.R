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

# Takes a vector of parameters on optimization scale
# – as returned by InitializeParameters – and transforms
# it to a vector of parameters on prediction scale as used in
# the FetoinfantSurv/Hzrd functions
RescaleParameters <- function (
    pars, model, split,
    zeta_range = c(10, 20),
    beta1_range = c(-0.7, 0.7),
    beta2_range = c(-0.7, 0.7)
) {
  rescaled_parameters <- switch(
    model,
    # negative exponential plus birth hump
    basic = function () {
      p <- list()
      p$alpha1 = exp(pars[1])
      p$beta1 = exp(pars[2])
      p$gamma = exp(pars[3])
      p$zeta = ScaleInverseLogit(pars[4], zeta_range[1], zeta_range[2])
      p$sigma = exp(pars[5])
      p$alpha2 = exp(pars[1]-p$zeta*p$beta1)
      p$beta2 = p$beta1
      p$tau = exp(pars[6])
      return(p)
    },
    # segmented two-part negative exponential plus birth hump
    flexible1 = function () {
      p <- list()
      p$alpha1 = exp(pars[1])
      p$beta1 = exp(pars[2])
      p$gamma = exp(pars[5])
      p$zeta = ScaleInverseLogit(pars[6], zeta_range[1], zeta_range[2])
      p$sigma = exp(pars[7])
      p$alpha2 = exp(pars[3])
      p$beta2 = exp(pars[4])
      p$tau = exp(pars[8])
      return(p)
    },
    # segmented two-part exponential (possibly non-monotonous) plus birth hump
    flexible2 = function () {
      p <- list()
      p$alpha1 = exp(pars[1])
      p$beta1 = ScaleInverseLogit(pars[2], beta1_range[1], beta1_range[2])
      p$gamma = exp(pars[5])
      p$zeta = ScaleInverseLogit(pars[6], zeta_range[1], zeta_range[2])
      p$sigma = exp(pars[7])
      p$alpha2 = exp(pars[3])
      p$beta2 = ScaleInverseLogit(pars[4], beta2_range[1], beta2_range[2])
      p$tau = exp(pars[8])
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
      p$alpha1 = log(ifelse(m[1] == 0, mean(m), m[1]))
      # log relative rate of ontogenescent mortality decline
      #p$b = log(abs(mean(diff(log(m)))))
      p$beta1 = (
        abs(
          log(mean(tail(m, 5), na.rm = TRUE)) -
            log(mean(head(m, 5), na.rm = TRUE))
        )
      ) / (length(m)-2.5 - 2.5)
      p$beta1 = log(p$beta1)
      if (is.infinite(p$beta1) | is.na(p$beta1)) { p$beta1 = log(0.1) }
      # log hazard at split point (zeta)
      p$alpha2 = p$alpha1-exp(p$beta1)*split
      # log relative rate of decline after split
      p$beta2 = p$beta1
      # log hazard contribution of birth hump at peak
      p$gamma = p$alpha1 - log(2)
      # location of birth component (0 gets rescaled to 39 weeks)
      p$zeta = log(0.07)
      # log spread of birth component
      p$sigma = log(1.6)
      # birth hump log left-skewness
      p$tau = log(2.7)
      return(p)
    },
    
    flexible1 = function () {
      p <- list()
      # log hazard at t=0
      p$alpha1 = log(ifelse(m[1] == 0, mean(m), m[1]))
      # log relative rate of ontogenescent mortality decline
      p$beta1 = log(abs(mean(diff(log(m)))))
      if (is.infinite(p$beta1) | is.na(p$beta1)) { p$beta1 = log(0.1) }
      # log hazard at split point
      p$alpha2 = p$alpha1-exp(p$beta1)*split
      # log relative rate of decline after split
      p$beta2 = p$beta1
      # log hazard contribution of birth hump at peak
      p$gamma = p$alpha1 - log(2)
      # location of birth component (0 gets rescaled to 39 weeks)
      p$zeta = log(0.07)
      # log spread of birth component
      p$sigma = log(1.6)
      # birth hump log left-skewness
      p$tau = log(2.7)
      return(p)
    },
    
    flexible2 = function () {
      p <- list()
      # log hazard at t=0
      p$alpha1 = log(mean(m))
      # ScaledLogit relative rate of ontogenescent mortality decline
      p$beta1 = mean(m)
      # log hazard at split point
      p$alpha2 = p$alpha1-p$beta1*split
      # ScaledLogit relative rate of decline after split
      p$beta2 = 1e-9
      # log hazard contribution of birth hump at peak
      p$gamma = p$alpha1 - log(2)
      # location of birth component (0 gets rescaled to 39 weeks)
      p$zeta = log(0.07)
      # log spread of birth component
      p$sigma = log(1.6)
      # birth hump log left-skewness
      p$tau = log(2.7)
      return(p)
    }
  
  )
  
  init_pars_pre <- initial_scaled_parameters()
  
  init_pars_ordered <- switch(
    model,
    basic = c(
      '1:alpha1' = init_pars_pre[['alpha1']],
      '2:beta1' = init_pars_pre[['beta1']],
      '3:gamma' = init_pars_pre[['gamma']],
      '4:zeta' = init_pars_pre[['zeta']],
      '5:sigma' = init_pars_pre[['sigma']],
      '6:tau' = init_pars_pre[['tau']]
    ),
    flexible1 = c(
      '1:alpha1' = init_pars_pre[['alpha1']],
      '2:beta1' = init_pars_pre[['beta1']],
      '3:alpha2' = init_pars_pre[['alpha2']],
      '4:beta2' = init_pars_pre[['beta2']],
      '5:gamma' = init_pars_pre[['gamma']],
      '6:zeta' = init_pars_pre[['zeta']],
      '7:sigma' = init_pars_pre[['sigma']],
      '8:tau' = init_pars_pre[['tau']]
    ),
    flexible2 = c(
      '1:alpha1' = init_pars_pre[['alpha1']],
      '2:beta1' = init_pars_pre[['beta1']],
      '3:alpha2' = init_pars_pre[['alpha2']],
      '4:beta2' = init_pars_pre[['beta2']],
      '5:gamma' = init_pars_pre[['gamma']],
      '6:zeta' = init_pars_pre[['zeta']],
      '7:sigma' = init_pars_pre[['sigma']],
      '8:tau' = init_pars_pre[['tau']]
    )
  )
  
  return(init_pars_ordered)
  
}

# Feto-infant parametric survival -------------------------------

# baseline characteristics
NegGompertzHzrd <-
  function (x, alpha1, beta1, alpha2, beta2, zeta) {
    I1 = ifelse(x < zeta, 1, 0)
    I2 = 1-I1
    alpha1*exp(-beta1*x)*I1 + alpha2*exp(-beta2*(x-zeta))*I2
  }
NegGompertzCumHzrd <-
  function (x, alpha1, beta1, alpha2, beta2, zeta) {
    I1 = ifelse(x < zeta, 1, 0)
    I2 = 1-I1
    
    H_zeta = ((alpha1 - alpha1*exp(-beta1*zeta)) / beta1)
    I1*((alpha1 - alpha1*exp(-beta1*x)) / beta1) +
      I2*((alpha2 - alpha2*exp(-beta2*(x-zeta))) / beta2 + H_zeta)
  }
NegGompertzSurv <-
  function (x, alpha1, beta1, alpha2, beta2, zeta) {
    exp(-NegGompertzCumHzrd(x, alpha1, beta1, alpha2, beta2, zeta))
  }

# birth hump characteristics
GaussianHzrd <-
  function(x, gamma, sigma, tau, zeta){
    I1 = ifelse(x < zeta, 1, 0)
    I2 = 1-I1
    gamma * exp( (-(x-zeta)^2) / (sigma + I1*tau) )
  }
GaussianCumHzrd <- #how to change this?
  function (x, gamma, sigma, tau, zeta) {
    I1 = ifelse(x < zeta, 1, 0)
    I2 = 1-I1
    
    A <- gamma*sqrt(pi)
    B <- sqrt(sigma)
    C <- sqrt(sigma+tau)
    AC <- A*C
    
    Erf1 <- pracma::erf(zeta/C)
    Erf2 <- pracma::erf((zeta-x)/B)
    Erf3 <- pracma::erf((zeta-x)/C)
    
    I1*(AC*(Erf1 - Erf3))/2 +  
      I2*(AC*Erf1 - A*B*Erf2)/2
  }
GaussianSurv <-
  function (x, gamma, sigma, tau, zeta) {
    exp(-GaussianCumHzrd(x, gamma, sigma, tau, zeta))
  }

# feto-infant competing risks survival
FetoinfantSurv <-
  function (x, pars, component = 'total', model = 'basic') {
    
    # ontogenescent survival
    ontogen_surv <-
      NegGompertzSurv(
        x = x,
        alpha1 = pars[['alpha1']], beta1 = pars[['beta1']], # what to parse here depends on the beginning of the code
        alpha2 = pars[['alpha2']], beta2 = pars[['beta2']],
        zeta = pars[['zeta']]
      )
    # birth survival
    birth_surv <-
      GaussianSurv(
        x = x,
        gamma = pars[['gamma']], sigma = pars[['sigma']],
        tau = pars[['tau']], zeta = pars[['zeta']]
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
  function (x, pars, component = 'total', model = 'basic') {
    
    # ontogenescent hazard
    ontogen_hzrd <-
      NegGompertzHzrd(
        x = x,
        alpha1 = pars[['alpha1']], beta1 = pars[['beta1']],
        alpha2 = pars[['alpha2']], beta2 = pars[['beta2']],
        zeta = pars[['zeta']]
      )
    
    # birth hazard
    birth_hzrd <-
      GaussianHzrd(
        x = x,
        gamma = pars[['gamma']], sigma = pars[['sigma']],
        tau = pars[['tau']], zeta = pars[['zeta']]
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
            lambda1 = 0, lambda2 = 0, lambda3 = 0, split,
            model = 'basic',
            zeta_range = c(10, 20),
            beta1_range = c(-0.7, 0.7),
            beta2_range = c(-0.7, 0.7),
            llsum = FALSE, llscale = 1, ...) {
    
    # rescale to constrained scale for evaluating survival functions
    pars2 <- RescaleParameters(pars, model, split,
                               zeta_range = zeta_range,
                               beta1_range = beta1_rage,
                               beta2_range = beta2_range)
    cat(unlist(pars2), '\n')
    
    # predict survival on basis of parameter estimates
    predSurvL <- SurvFnct(x = age, pars = pars2,
                          component = 'total', model = model)
    predSurvR <- SurvFnct(x = age+width, pars = pars2,
                          component = 'total', model = model)
    
    loglike <-
      obsDx*log(predSurvL-predSurvR) + obsCx*log(predSurvR)
    
    penalty <- 0
    
    if (isTRUE(model == 'basic')) {
      # penalize birth hump magnitude
      penalty <- lambda1*pars2$gamma
    }
    
    if (isTRUE(model == 'flexible1')) {
      penalty <-
        # penalize birth hump magnitude
        lambda1*pars2$gamma +
        # penalize discontinuities between the two ontogenescent segments on log scale
        # (log-alpha2 - (log-alpha1 - beta1*zeta))^2
        lambda2*(pars[3] - (pars[1]-pars2$beta1*pars2$zeta))^2 +
        # penalize differences in slope between the two ontogenescent segments
        # (log-beta2 - log-beta1)^2
        lambda3*(pars[4]-pars[2])^2
    }
    
    if (isTRUE(model == 'flexible2')) {
      penalty <-
        # penalize birth hump magnitude
        lambda1*pars$gamma +
        # penalize discontinuities between the two ontogenescent segments on log scale
        # (log-alpha2 - (log-alpha1 - beta1*zeta))^2
        lambda2*(pars[3] - (pars[1]-pars2$beta1*pars2$zeta))^2 +
        # penalize differences in slope between the two ontogenescent segments
        lambda3*(pars[4]-pars[2])^2
    }
    
    if (isTRUE(llsum)) { loglike <- sum(loglike) } 
    loglike[is.nan(loglike)] <- Inf
    loglike <- (loglike-penalty)*llscale
    
    return(loglike)
    
  }

# Fit -------------------------------------------------------------

ControlFitFetoinfantSurvival <- function (
    model = 'basic',
    split = 15,
    # fitting options
    lambda1 = 0, lambda2 = 0, lambda3 = 0,
    zeta_range = c(10, 20),
    beta1_range = c(-0.7, 0.7),
    beta2_range = c(-0.7, 0.7),
    method = 'optim',
    DEoptim_control = DEoptim::DEoptim.control(),
    # simulation
    simulate = TRUE, nsim = 1e3, hessian_inverse = 'cholesky',
    exclude_from_hessian_inverse = NULL
) {
  control_pars <- list(
    model = model,
    split = split,
    lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3,
    zeta_range = zeta_range,
    beta1_range = beta1_range,
    beta2_range = beta2_range,
    method = method,
    DEoptim_control = DEoptim_control,
    simulate = simulate, nsim = nsim, hessian_inverse = hessian_inverse,
    exclude_from_hessian_inverse = exclude_from_hessian_inverse
  )
  return(control_pars)
}

# Invert a Square Matrix with options for partial inversion and pseudo-inverses
SquareMatrixInverse <- function (X, inverse = 'cholesky', exclude_from_inverse = NULL) {
  
  # input checks and preparation
  if (!isTRUE(is.matrix(X))) stop('X is not matrix')
  N = NROW(X); M = NCOL(X) 
  if (N!=M) stop('X is not square')
  if (is.null(colnames(X))) {
    colnames(X) <- 1:N
  }
  names_X <- colnames(X)
  
  # exclude entries from matrix if specified
  n_exclude <- length(exclude_from_inverse)
  if (!is.null(exclude_from_inverse)) {
    X_ <- X[-exclude_from_inverse, -exclude_from_inverse]
  } else {
    X_ <- X
  }
  
  # invert (partial) matrix  
  inverse_X <- 
    switch(inverse,
           cholesky = chol2inv(chol(X_)),
           choleskypivot = chol2inv(chol(X_, pivot = TRUE)),
           pseudo = pinv(X_))
  
  # if needed, restore excluded entries
  if (!is.null(exclude_from_inverse)) {
    names_to_exclude <- names_X[exclude_from_inverse]
    N_removed <- length(exclude_from_inverse)
    names_to_keep <- names_X[-exclude_from_inverse]
    # reconstruct full matrix with 0s for previously removed entries
    inverse_X <-
      inverse_X %>%
      cbind(matrix(0, nrow = N-N_removed, ncol = N_removed)) %>%
      rbind(matrix(0, ncol = N, nrow = N_removed))
    colnames(inverse_X) <- c(names_to_keep, names_to_exclude)
    rownames(inverse_X) <- c(names_to_keep, names_to_exclude)
    inverse_X <- inverse_X[names_X,names_X] # restore original order  
  }
  
  return(inverse_X)
}

# Fit a parametric model of fetoinfant survival to
# the fetoinfant lifetable assuming left truncation age
# at 24 and right censoring at 77. Derive various predictions.
# Sample parameters from multivariate normal distribution
# derived from hessian in order to calculate CIs.
FitFetoinfantSurvival <-
  function (
    filt, control = ControlFitFetoinfantSurvival()
  ) {
    
    require(pracma)
    
    stopifnot(any(class(filt) == 'FILT'))
    
    fit_lifetable <-
      filt %>%
      left_join(FILTMortalityRates(.)) %>%
      left_join(FILTSurvival(.)) %>%
      group_by(stratum) %>%
      group_modify(~{
        
        ### initial parameters and box constraints ###
        
        init_pars <-
          InitializeParameters(.x$m, control$model, control$split)
        
        if (identical(control$model, 'basic')) {
          lower = c(-20, -20, -20, -20, -10, -10)
          upper = c( -2,  -1,   1,  10,   3,   5)
        }
        if (identical(control$model, 'flexible1')) {
          lower = c(-20, -20, -20, -20, -20, -10, -10, -10)
          upper = c( -2,  -1,  -2,  -1,   1,  10,   5,   5)
        }
        if (identical(control$model, 'flexible2')) {
          lower = c(-20, -10, -20, -10, -20, -10, -10, -10)
          upper = c( -2,  10,  -2,  10,   1,  10,   5,   5)
        }
        
        ### the fit ###
        
        if (identical(control$method, 'deoptim')) {
          require(DEoptim)
          modelfit <- DEoptim(
            fn = IntervalCensoredLogLike,
            lower = lower,
            upper = upper,
            # data and arguments to objective function
            age =
              pull(.x, x)-24,
            width =
              pull(.x, n),
            obsDx =
              pull(.x, D),
            obsCx =
              pull(.x, C),
            # survival function
            SurvFnct = FetoinfantSurv,
            lambda1 = control$lambda1,
            lambda2 = control$lambda2,
            lambda3 = control$lambda3,
            model = control$model,
            zeta_range = control$zeta_range,
            beta1_range = control$beta1_range,
            beta2_range = control$beta2_range,
            llsum = TRUE,
            llscale = -1,
            control = control$DEoptim_control
          )
          pars_estimate <- modelfit$optim$bestmem
          hessian <- numDeriv::hessian(
            IntervalCensoredLogLike,
            x = pars_estimate,
            age =
              pull(.x, x)-24,
            width =
              pull(.x, n),
            obsDx =
              pull(.x, D),
            obsCx =
              pull(.x, C),
            # survival function
            SurvFnct = FetoinfantSurv,
            lambda1 = control$lambda1,
            lambda2 = control$lambda2,
            lambda3 = control$lambda3,
            model = control$model,
            zeta_range = control$zeta_range,
            llsum = TRUE
          )
        }
        
        if (identical(control$method, 'optim')) {
          require(maxLik)
          modelfit <-
            maxLik(
              logLik = IntervalCensoredLogLike,
              start  = init_pars,
              method = 'BFGS',
              # data and arguments to objective function
              age =
                pull(.x, x)-24,
              width =
                pull(.x, n),
              obsDx =
                pull(.x, D),
              obsCx =
                pull(.x, C),
              # survival function
              SurvFnct = FetoinfantSurv,
              lambda1 = control$lambda1,
              lambda2 = control$lambda2,
              lambda3 = control$lambda3,
              model = control$model,
              zeta_range = control$zeta_range,
              llsum = FALSE,
              llscale = 1,
              # options
              iterlim = 1e4
            )
          pars_estimate <- coef(modelfit)
          hessian <- modelfit$hessian
        }
        
        ### posterior simulation ###
        
        # 1000 draws from the posterior parameter distribution
        # assuming multivariate normal derived from hessian
        par_draw <-
          expand_grid(
            draw = 1:control$nsim,
            name = names(init_pars)
          )
        if (isTRUE(control$simulate)) {
          hessian_inverse <- 
            SquareMatrixInverse(
              -hessian,
              inverse = control$hessian_inverse,
              exclude_from_inverse = control$exclude_from_hessian_inverse
            )
          par_draw$value <- MASS::mvrnorm(
            n = control$nsim,
            mu = pars_estimate,
            Sigma = hessian_inverse
          ) %>% t() %>% c()
        } else {
          par_draw$value <- rep(pars_estimate, control$nsim)
        }
        
        # rescale posterior simulations for use with
        # FetoinfantSurv/Hzrd functions
        par_draw_rescaled <-
          par_draw %>%
          group_by(draw) %>%
          group_modify(~{
            pars <- RescaleParameters(
              .x$value, control$model, control$split,
              zeta_range = control$zeta_range,
              beta1_range = control$beta1_range,
              beta2_range = control$beta2_range
            )
            tibble(name = names(pars), value = unlist(pars))
          })
        
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
        
        pars_rescaled <-
          par_draw_rescaled %>%
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
            pars <- RescaleParameters(
              .x$value, control$model, control$split,
              zeta_range = control$zeta_range,
              beta1_range = control$beta1_range,
              beta2_range = control$beta2_range
            )
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
                  model = control$model
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
                  model = control$model
                ),
              # hazard of birth component at x
              birth_hx =
                FetoinfantHzrd(
                  x,
                  pars = pars,
                  component = 'birth',
                  model = control$model
                ),
              # hazard of ontogenescent component at x
              ontogen_hx =
                FetoinfantHzrd(
                  x,
                  pars = pars,
                  component = 'ontogen',
                  model = control$model
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
        
        ### stratum specific output object ###
        
        tibble(
          # maximum log likelihood
          loglike = list(pars_estimate),
          # posterior parameter simulations on fitting scale
          par_draws = list(par_draw),
          par_summary = list(pars),
          # posterior parameter simulations on evaluation scale
          par_rescaled_draws = list(par_draw_rescaled),
          par_rescaled_summary = list(pars_rescaled),
          # posterior derived statistics simulations
          pred_draws = list(pred_draw),
          pred_summary = list(pred),
          # information on fit
          hessian = list(hessian),
          model = list(modelfit),
          control = list(control),
          par_init = list(init_pars),
          # input data
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
    xbrk = cnst$gestage_brk,
    xlim =
      c(cnst$left_truncation_gestage,
        cnst$right_censoring_gestage-0.1),
    ylim_hx = NULL,
    ylim_Fx = NULL,
    scaler = 1e5,
    legend = TRUE,
    colors = fig_spec$discrete_colors,
    components = FALSE,
    ar = 0.85,
    notitle = FALSE
  ) {
    
    legend_pos <- c(0.8, 0.8)
    if (!isTRUE(legend)) {
      legend_pos <- 'none'
    }
    
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
        breaks = c(1, seq(2, 10, 2), seq(20, 80, 20)),
        trans = 'log10',
        name = NULL
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
        legend.position = legend_pos,
        legend.background = element_blank(),
        legend.box.background = element_blank()
      ) +
      labs(
        subtitle = paste0('Feto-infant deaths per\n',
                          formatC(scaler, format = 'd', big.mark = ','),
                          ' weeks at risk')
      ) +
      coord_cartesian(ylim = ylim_hx)
    
    if (isTRUE(components)) {
      plot_hzrd <-
        plot_hzrd +
        geom_line(
          aes(
            x = x, y = avg_birth_hx*scaler, group = stratum,
            color = as.character(stratum)
          ),
          size = fig_spec$line_size_m,
          data = predictions
        ) +
        geom_line(
          aes(
            x = x, y = avg_ontogen_hx*scaler, group = stratum,
            color = as.character(stratum)
          ),
          size = fig_spec$line_size_m,
          data = predictions
        )
    }
    
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
        name = NULL
      ) +
      scale_color_manual(values = colors) +
      scale_x_continuous(
        '', breaks = xbrk
      ) +
      scale_fill_manual(
        '',
        values = colors
      ) +
      coord_cartesian(ylim = ylim_Fx) +
      fig_spec$MyGGplotTheme(ar = ar) +
      theme(
        legend.position = 'none'
      ) +
      labs(subtitle = paste0('Cumulative feto-infant deaths\n',
                          'out of ', formatC(scaler, format = 'd',
                                             big.mark = ','),
                          ' cohort members')
      )
    
    if (isTRUE(notitle)) {
      plot_surv <- plot_surv + labs(subtitle = NULL)
      plot_hzrd <- plot_hzrd + labs(subtitle = NULL)
    }
    
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
        'birth_Fx', 'birth_iFx',
        'ontogen_Fx', 'ontogen_iFx'
      )))
    
  }
