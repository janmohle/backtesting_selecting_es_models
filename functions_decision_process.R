####################################################################
################### Function calculates VaR and ES #################
####################################################################
VaR_ES_forecast <- function(data_zoo,
                            var_spec,
                            mean_spec,
                            dist_spec,
                            tolerance_lvl,
                            estimate=TRUE,
                            fixed_pars=NULL, # Required if estimate=FALSE
                            seed_opt=NULL,
                            current_loop=NA_integer_,
                            solver_problems_log=FALSE){
  
  # Store last date
  last_date <- tail(index(data_zoo), 1)
  
  if(estimate){
    
    # GARCH Specification
    garchspec <- ugarchspec(variance.model = var_spec,
                            mean.model = mean_spec,
                            distribution.model = dist_spec)
    
    # Fitting GARCH model
    garchfit <- tryCatch(
      {
        ugarchfit(spec = garchspec,
                  data = data_zoo,
                  solver = 'hybrid',
                  solver.control = list(n.restarts = 10,#25                                       # gosolnp: Number of restarts
                                        n.sim = 200,#500                                          # gosolnp: Number of initial parameter simulation -> min objective function is taken for optimization initiation
                                        rseed = if (!is.null(seed_opt)) seed_opt + 111 else NULL,  # gosolnp: Seed for inital parameter simulation
                                        trace = 0))
      }, error = function(e){
        if(solver_problems_log){
          write(paste0('\n1: ', var_spec[['model']], ' ', dist_spec, ' ', paste(mean_spec[['armaOrder']], collapse = ' '),' loop ', current_loop),
                file = 'run_logs/solver_problems.txt', append = TRUE)
        }
        return(NULL)
      }
    )
    
    # Assign NULL to garchfit if solver did not converge
    if(!is.null(garchfit) && garchfit@fit$convergence == 1){
      if(solver_problems_log){
        write(paste0('\n2: ', var_spec[['model']], ' ', dist_spec, ' ', paste(mean_spec[['armaOrder']], collapse = ' '),' loop ', current_loop),
              file = 'run_logs/solver_problems.txt', append = TRUE)
      }
      garchfit <- NULL
    }
    
    # Try nloptr+AUGLAG+PRAXIS if previous fitting with hybrid solver failed
    if(is.null(garchfit)){
      garchfit <- tryCatch(
        {
          ugarchfit(spec = garchspec,
                    data = data_zoo,
                    solver = 'nloptr',
                    solver.control = list(solver = 8,                                                 # AUGLAG+PRAXIS solver
                                          maxeval = 1000,#added
                                          ranseed = if (!is.null(seed_opt)) seed_opt + 111 else NULL,  # Seed for random processes
                                          print_level = 0))
        }, error = function(e){
          if(solver_problems_log){
            write(paste0('\n3: ', var_spec[['model']], ' ', dist_spec, ' ', paste(mean_spec[['armaOrder']], collapse = ' '),' loop ', current_loop),
                  file = 'run_logs/solver_problems.txt', append = TRUE)
          }
          return(NULL)
        }
      )
    }
    
    # Assign NULL to garchfit if solver did not converge
    if(!is.null(garchfit) && garchfit@fit$convergence == 1){
      if(solver_problems_log){
        write(paste0('\n4: ', var_spec[['model']], ' ', dist_spec, ' ', paste(mean_spec[['armaOrder']], collapse = ' '),' loop ', current_loop),
              file = 'run_logs/solver_problems.txt', append = TRUE)
      }
      garchfit <- NULL
    }
    
    # Try gosolnp alone with high number of simulations and tries if previous fitting with hybrid and nloptr solver failed
    if(is.null(garchfit)){
      garchfit <- tryCatch(
        {
          ugarchfit(spec = garchspec,
                    data = data_zoo,
                    solver = 'gosolnp',
                    solver.control = list(n.restarts = 20,#100
                                          n.sim = 500,#2000
                                          rseed = if (!is.null(seed_opt)) seed_opt + 231 else NULL,
                                          trace = 0))
        }, error = function(e){
          if(solver_problems_log){
            write(paste0('\n5: ', var_spec[['model']], ' ', dist_spec, ' ', paste(mean_spec[['armaOrder']], collapse = ' '),' loop ', current_loop),
                  file = 'run_logs/solver_problems.txt', append = TRUE)
          }
          return(NULL)
        }
      )
    }
    
    # Assign NULL to garchfit if solver did not converge
    if(!is.null(garchfit) && garchfit@fit$convergence == 1){
      if(solver_problems_log){
        write(paste0('\n6: ', var_spec[['model']], ' ', dist_spec, ' ', paste(mean_spec[['armaOrder']], collapse = ' '),' loop ', current_loop),
              file = 'run_logs/solver_problems.txt', append = TRUE)
      }
      garchfit <- NULL
    }
    
    # Assign NULL to garchfit if coef contains NAs
    if(!is.null(garchfit) && anyNA(rugarch::coef(garchfit))){
      if(solver_problems_log){
        write(paste0('\n7: ', var_spec[['model']], ' ', dist_spec, ' ', paste(mean_spec[['armaOrder']], collapse = ' '),' loop ', current_loop),
              file = 'run_logs/solver_problems.txt', append = TRUE)
      }
      garchfit <- NULL
    }
    
    # Assign NULL to garchfit if coef contains infinite numbers
    if(!is.null(garchfit) && any(is.infinite(rugarch::coef(garchfit)))){
      if(solver_problems_log){
        write(paste0('\n8: ', var_spec[['model']], ' ', dist_spec, ' ', paste(mean_spec[['armaOrder']], collapse = ' '),' loop ', current_loop),
              file = 'run_logs/solver_problems.txt', append = TRUE)
      }
      garchfit <- NULL
    }
    
  } else {
    
    # GARCH Specification with all parameters fixed
    garchspec <- ugarchspec(variance.model = var_spec,
                            mean.model = mean_spec,
                            distribution.model = dist_spec,
                            fixed.pars = fixed_pars)
    
    # Filter data
    garchfit <- ugarchfilter(spec = garchspec,
                             data = data_zoo)
  }
  
  if(!is.null(garchfit)){
    
    # Extracting coefficients
    coef_fit <- rugarch::coef(garchfit)
    
    # Extracting skewness, shape and lambda parameter (if not there, NA gets assigned)
    skew <- if ('skew' %in% names(coef_fit)) unname(coef_fit['skew']) else NA_real_
    shape <- if ('shape' %in% names(coef_fit)) unname(coef_fit['shape']) else NA_real_
    lambda <- if ('ghlambda' %in% names(coef_fit)) unname(coef_fit['ghlambda']) else NA_real_
    
    # Forecasting
    if(estimate){
      garchforecast <- ugarchforecast(fitORspec = garchfit,
                                      n.ahead = 1)
    } else {
      garchforecast <- ugarchforecast(fitORspec = garchspec,
                                      data = data_zoo,
                                      n.ahead = 1)
    }
    
    # Extracting mu and sigma
    mu <- rugarch::fitted(garchforecast)
    sigma <- rugarch::sigma(garchforecast)
    
    # Calculate VaR
    VaR <- mu + sigma * rugarch::qdist(distribution = dist_spec,
                                       p = tolerance_lvl,
                                       skew = skew,
                                       shape = shape,
                                       lambda = lambda)
    
    # Calculate ES
    ES <- mu + sigma / tolerance_lvl * stats::integrate(f = function(x) rugarch::qdist(distribution = dist_spec,
                                                                                       p = x,
                                                                                       skew = skew,
                                                                                       shape = shape,
                                                                                       lambda = lambda),
                                                        lower = 0,
                                                        upper = tolerance_lvl)[['value']]
  
    # Return results
    results <- tibble(last_date = last_date,
                      mu = as.double(mu),
                      sigma = as.double(sigma),
                      skew = as.double(skew),
                      shape = as.double(shape),
                      lambda = as.double(lambda),
                      dist_spec = dist_spec,
                      VaR = as.double(VaR),
                      ES = as.double(ES))
  } else {
    
    if(solver_problems_log){
      # Record complete fail of estimation
      write(paste0('\n\nComplete fail: ', var_spec[['model']], ' ', dist_spec, ' ', paste(mean_spec[['armaOrder']], collapse = ' '),' loop ', current_loop, '\n'),
            file = 'run_logs/solver_problems.txt', append = TRUE)
    }
      
    results <- tibble(last_date = last_date,
                      mu = NA_real_,
                      sigma = NA_real_,
                      skew = NA_real_,
                      shape = NA_real_,
                      lambda = NA_real_,
                      dist_spec = NA_character_,
                      VaR = NA_real_,
                      ES = NA_real_)
  }
  
  results_zoo <- zoo(x = dplyr::select(results, -'last_date'),
                     order.by = results[['last_date']])
  return(results_zoo)
}

####################################################################
################### Loss function of Fissler et. al. ###############
####################################################################
loss_func <- function(a,
                      r,
                      VaR,
                      ES,
                      g1,
                      g2,
                      G2){
  I <- ifelse(r <= VaR, 1, 0)
  loss <- (I - a) * (g1(VaR) - g1(r)) + (1/a) * g2(ES) * I * (VaR - r) + g2(ES) * (ES - VaR) - G2(ES)
  return(loss)
}


##############################################################################
#############  Automatic lag selection for CC backtest  ######################
##############################################################################
auto_lag <- function(y, max_lag=NULL, na_rm=TRUE, q=2.4) {
  
  if(na_rm){
    y <- y[!is.na(y)]
  } else {
    if(anyNA(y)){
      return(NA)
    }
  }
  
  n <- length(y)
  
  if(is.null(max_lag)){
    max_lag <- floor(sqrt(n))
  }
  
  y <- y - mean(y)
  
  if(all(y == 0)){
    return(1)
  }
  
  auto_cov <- acf(y,lag.max = max_lag, type = 'covariance', plot = FALSE)$acf
  
  rob_var <- sapply(1:max_lag, function(j) {
    mean(y[(1 + j):n]^2 * y[1:(n - j)]^2)
  })
  
  if(any(rob_var == 0)){
    return(1)
  }
  
  auto_cor_rob_2 <- (auto_cov[2:(max_lag + 1)]^2) / rob_var
  Q_p  <- n * cumsum(auto_cor_rob_2)
  
  penalty <- if(sqrt(n) * max(sqrt(abs(auto_cor_rob_2)), na.rm = TRUE) <= sqrt(q * log(n))) log(n) else 2
  
  which.max(Q_p - (1:max_lag) * penalty)
}


####################################################################
################### Backtesting and loss function routine ##########
####################################################################
evaluate_var_es_backtests_and_loss <- function(return_data_tbl, # expects columns: Return and Date
                                               spec_est_lst,
                                               tolerance_lvl,
                                               sign_lvl_backtests,
                                               est_window,
                                               estimate=TRUE,
                                               lags_es_cc=NULL, # if NULL then automatic lag selection
                                               max_lag_es_cc=NULL, # if NULL then max_lag=floor(sqrt(n))
                                               n_boot_de=NULL,
                                               current_seed=NULL,
                                               current_loop=NA_integer_,
                                               run_log=FALSE,
                                               backtests=c('VaR_UC','VaR_CC','UC','CC','A_ESR','S_ESR')){
  
  # Initiate list for loss values
  loss_lst <- list()
  
  # Turn tibble with returns into zoo object
  return_data_zoo <- zoo(x = return_data_tbl[['Return']],
                         order.by = return_data_tbl[['Date']])
  
  # Rolling forecasting over estimation specifications
  for(spec_est_name in names(spec_est_lst)){
    
    # Extract current specification
    spec_est_i <- spec_est_lst[[spec_est_name]]
    
    # Calculate GARCH model
    rolling_VaR_ES <- rollapply(data = return_data_zoo,
                                width = est_window,
                                FUN = function(x) VaR_ES_forecast(data_zoo = x,
                                                                  var_spec = spec_est_i[['var_spec']],
                                                                  mean_spec = spec_est_i[['mean_spec']],
                                                                  dist_spec = spec_est_i[['dist_spec']],
                                                                  tolerance_lvl = tolerance_lvl,
                                                                  estimate = estimate,
                                                                  fixed_pars = spec_est_i[['fixed_pars']],
                                                                  seed_opt = current_seed,
                                                                  current_loop = current_loop,
                                                                  solver_problems_log = run_log),
                                align = 'right',
                                coredata = FALSE)
    
    rolling_VaR_ES_lead <- stats::lag(rolling_VaR_ES,
                                      k = -1)
    rolling_VaR_ES_tbl <- as_tibble(rolling_VaR_ES_lead) %>%
      mutate(across(-c(dist_spec), as.numeric),
             Date = index(rolling_VaR_ES_lead))
    VaR_ES_results_tbl <- inner_join(return_data_tbl, rolling_VaR_ES_tbl, by = 'Date')
    
    # Calculate u
    VaR_ES_results_tbl <- VaR_ES_results_tbl %>%
      mutate(u = rugarch::pdist(distribution = spec_est_i[['dist_spec']],
                                q = Return,
                                mu = mu,
                                sigma = sigma,
                                skew = skew,
                                shape = shape,
                                lambda = lambda))
    
    # Calculate loss (ES must be < 0, otherwise loss = NA)
    VaR_ES_results_tbl <- VaR_ES_results_tbl %>%
      mutate(loss = loss_func(a = tolerance_lvl,
                              r = Return,
                              VaR = VaR,
                              ES = ES,
                              g1 = function(x) 0 * x,
                              g2 = function(x) 1/(-x),
                              G2 = function(x) -log(-x)
                              )
      )
    
    # Remove NAs for backtests
    VaR_ES_results_tbl_noNA <- VaR_ES_results_tbl %>%
      tidyr::drop_na(Return, VaR, ES, u, loss)
    
    if(run_log){
      # Record how often ES >= 0
      es_ge_0_count <- sum(VaR_ES_results_tbl[['ES']] >= 0, na.rm = TRUE)
      
      # Record NAs
      cols_check <- c('Return', 'VaR', 'ES', 'u', 'loss')
      na_counts <- VaR_ES_results_tbl %>%
        dplyr::summarise(dplyr::across(all_of(cols_check), ~sum(is.na(.))))
      
      if(any(na_counts > 0)){
        if(es_ge_0_count > 0){
          write(paste0('\nest_', spec_est_name, ' loop ', current_loop, ': ', paste(names(na_counts), na_counts, collapse = ', '), '\n',
                       'ES >= 0: ', es_ge_0_count, ' times\n'),
                file = 'run_logs/NA_information.txt',
                append = TRUE
          )
        } else {
          write(paste0('\nest_', spec_est_name, ' loop ', current_loop, ': ', paste(names(na_counts), na_counts, collapse = ', '), '\n'),
                file = 'run_logs/NA_information.txt',
                append = TRUE
          )
        }
      }
    }
    
    # Backtests
    pvals <- list()
    for(tst in backtests){
      pvals[[tst]] <- NA_real_
    }
    
    # VaR unconditional coverage test from Kupiec and conditional coverage test from Christofferson
    if(('VaR_UC' %in% backtests) || ('VaR_CC' %in% backtests)){
      
      VaR_uc_cc_res <- tryCatch({rugarch::VaRTest(alpha = tolerance_lvl,
                                                  actual = VaR_ES_results_tbl_noNA[['Return']],
                                                  VaR = VaR_ES_results_tbl_noNA[['VaR']],
                                                  conf.level = (1-sign_lvl_backtests))},
                                error = function(e){return(NULL)})
      
      if(!is.null(VaR_uc_cc_res)){
        if('VaR_UC' %in% backtests){
          pvals[['VaR_UC']] <- VaR_uc_cc_res[['uc.LRp']]
        }
        if('VaR_CC' %in% backtests){
          pvals[['VaR_CC']] <- VaR_uc_cc_res[['cc.LRp']]
        }
      }
    }

    # ES unconditional and conditional backtests
    if(('UC' %in% backtests) || ('CC' %in% backtests)){
      if(!is.null(current_seed)){set.seed(current_seed+1)}
      
      u <- VaR_ES_results_tbl_noNA[['u']]
      CumVio <- (1 / tolerance_lvl) * (tolerance_lvl - u) * ifelse(u <= tolerance_lvl, 1, 0)
      
      ##TEST
      print(if(!is.null(lags_es_cc)) lags_es_cc else auto_lag(y = CumVio,
                                                              max_lag = max_lag_es_cc))
      
      shortfall_de_test_res <- tryCatch({tstests::shortfall_de_test(x = u,
                                                                    alpha = tolerance_lvl,
                                                                    lags = if(!is.null(lags_es_cc)) lags_es_cc else auto_lag(y = CumVio,
                                                                                                                             max_lag = max_lag_es_cc),
                                                                    boot = !is.null(n_boot_de),
                                                                    n_boot = n_boot_de)},
                                        error = function(e){return(NULL)})
      
      if(!is.null(shortfall_de_test_res)){
        if('UC' %in% backtests){
          pvals[['UC']] <- shortfall_de_test_res[['p_value']][1]
        }
        if('CC' %in% backtests){
          pvals[['CC']] <- shortfall_de_test_res[['p_value']][2]
        }
      }
    }

    # Auxiliary ESR backtest
    if('A_ESR' %in% backtests){
      if(!is.null(current_seed)){set.seed(current_seed+2)}
      
      pvals[['A_ESR']] <- tryCatch({esback::esr_backtest(r = VaR_ES_results_tbl_noNA[['Return']],
                                                         q = VaR_ES_results_tbl_noNA[['VaR']],
                                                         e = VaR_ES_results_tbl_noNA[['ES']],
                                                         alpha = tolerance_lvl,
                                                         version = 2)[['pvalue_twosided_asymptotic']]},
                                   error = function(e){return(NA_real_)})
    }
    
    # Strict ESR backtest
    if('S_ESR' %in% backtests){
      if(!is.null(current_seed)){set.seed(current_seed+3)}
      
      pvals[['S_ESR']] <- tryCatch({esback::esr_backtest(r = VaR_ES_results_tbl_noNA[['Return']],
                                                         e = VaR_ES_results_tbl_noNA[['ES']],
                                                         alpha = tolerance_lvl,
                                                         version = 1)[['pvalue_twosided_asymptotic']]},
                                   error = function(e){return(NA_real_)})
    }
    
    pvals_selected <- unlist(pvals)
    
    # Save loss column if backtesting passed
    if(!anyNA(pvals_selected) && all(pvals_selected > sign_lvl_backtests)){
      loss_lst[[spec_est_name]] <- VaR_ES_results_tbl[['loss']] #NAs are also returned!
    }
    
    # Log that a backtest evaluated as NA
    if(run_log && anyNA(pvals_selected)){
      write(paste0('\nBacktest NA in est_', spec_est_name, ' loop ', current_loop, ': ', paste(names(pvals_selected), pvals_selected, collapse = ', '), '\n'),
            file = 'run_logs/NA_information.txt',
            append = TRUE
      )
    }
  }

  # Create tibble of losses
  loss_tbl <- as_tibble(loss_lst)
  
  # Remove NAs
  loss_tbl_noNA <- tidyr::drop_na(loss_tbl)
  
  # Report number of dropped rows due to NAs
  n_rows_dropped <- nrow(loss_tbl) - nrow(loss_tbl_noNA)
  if(run_log && (n_rows_dropped > 0)){
    write(paste0('\nRows removed due to NAs in loss table in loop ', current_loop, ': ', n_rows_dropped, '\n'),
      file = 'run_logs/NA_information.txt',
      append = TRUE
    )
  }
  
  # Return tibble with losses without NAs
  return(loss_tbl_noNA)
}