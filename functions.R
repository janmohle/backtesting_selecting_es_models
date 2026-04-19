####################################################################
################### Function calculates VaR and ES #################
####################################################################
VaR_ES_forecast <- function(data_zoo,
                            var_spec,
                            mean_spec,
                            dist_spec,
                            tolerance_lvl,
                            estimate=TRUE,
                            par_corr=FALSE,
                            fixed_pars=NULL,
                            empirical=FALSE,
                            seed_opt=NULL,
                            current_loop=NA_integer_){
  
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
                                        rseed = if (!is.null(seed_opt)) seed_opt + 111 else NULL, # gosolnp: Seed for inital parameter simulation
                                        trace = 0))
      }, error = function(e){
        write(paste0('\n1: loop ', current_loop),
              file = 'run_logs/solver_problems.txt', append = TRUE)
        return(NULL)
      }
    )
    
    # Assign NULL to garchfit if solver did not converge
    if(!is.null(garchfit) && garchfit@fit$convergence == 1){
      write(paste0('\n2: loop ', current_loop),
            file = 'run_logs/solver_problems.txt', append = TRUE)
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
                                          ranseed = if (!is.null(seed_opt)) seed_opt + 111 else NULL, # Seed for random processes
                                          print_level = 0))
        }, error = function(e){
          write(paste0('\n3: loop ', current_loop),
                file = 'run_logs/solver_problems.txt', append = TRUE)
          return(NULL)
        }
      )
    }
    
    # Assign NULL to garchfit if solver did not converge
    if(!is.null(garchfit) && garchfit@fit$convergence == 1){
      write(paste0('\n4: loop ', current_loop),
            file = 'run_logs/solver_problems.txt', append = TRUE)
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
          write(paste0('\n5: loop ', current_loop),
                file = 'run_logs/solver_problems.txt', append = TRUE)
          return(NULL)
        }
      )
    }
    
    # Assign NULL to garchfit if solver did not converge
    if(!is.null(garchfit) && garchfit@fit$convergence == 1){
      write(paste0('\n6: loop ', current_loop),
            file = 'run_logs/solver_problems.txt', append = TRUE)
      garchfit <- NULL
    }
    
    # Assign NULL to garchfit if coef contains NAs
    if(!is.null(garchfit) && anyNA(rugarch::coef(garchfit))){
      write(paste0('\n7: loop ', current_loop),
            file = 'run_logs/solver_problems.txt', append = TRUE)
      garchfit <- NULL
    }
    
    # Assign NULL to garchfit if coef contains infinite numbers
    if(!is.null(garchfit) && any(is.infinite(rugarch::coef(garchfit)))){
      write(paste0('\n8: loop ', current_loop),
            file = 'run_logs/solver_problems.txt', append = TRUE)
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
    
    # Extracting last standard deviation from fit
    st_dev_t_min_1 <- tail(rugarch::sigma(garchfit), 1)

    # Extracting last residual from fit
    residual_t_min_1 <- tail(rugarch::residuals(garchfit), 1)
    
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
    
    # Extract covariance matrix of parameters and calculate gradients of mu and sigma evaluated at coef_fit
    if(par_corr & estimate){
      vcov_matrix_par <- matrix_to_string(rugarch::vcov(garchfit))

      gradients <- numDeriv::jacobian(func = function(x) {
        
        x_named <- setNames(as.list(x), names(coef_fit))
        
        spec_tmp <- ugarchspec(variance.model = var_spec,
                               mean.model = mean_spec,
                               distribution.model = dist_spec,
                               fixed.pars = x_named)
        
        forecast_tmp <- ugarchforecast(fitORspec = spec_tmp,
                                       data = data_zoo,
                                       n.ahead = 1)
        
        return(c(mu = as.double(rugarch::fitted(forecast_tmp)),
                 sigma = as.double(rugarch::sigma(forecast_tmp))))
      },
      x = coef_fit,
      method = 'Richardson')
      
      mu_grad <- vector_to_string(gradients[1,])
      sigma_grad <- vector_to_string(gradients[2,])
    }
    
    # Extracting mu and sigma
    mu <- rugarch::fitted(garchforecast)
    sigma <- rugarch::sigma(garchforecast)
    
    if(empirical){
      
      # Empirical distribution
      emp_dist <- as.vector(rugarch::residuals(garchfit, standardize = TRUE))
      emp_dist_vec <- vector_to_string(emp_dist)

      # pth quantile of empirical distribution
      q_emp_dist <- stats::quantile(emp_dist,
                                    probs = tolerance_lvl,
                                    na.rm = TRUE,
                                    names = FALSE,
                                    type = 1)

      # VaR and ES calculation
      VaR <- mu + sigma * q_emp_dist
      ES <- mu + sigma * mean(emp_dist[emp_dist <= q_emp_dist])
      
    } else{
      
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
    }
    
    # Return results
    results <- tibble(last_date = last_date,
                      mu = as.double(mu),
                      sigma = as.double(sigma),
                      variance = as.double(sigma)^2,
                      skew = as.double(skew),
                      shape = as.double(shape),
                      lambda = as.double(lambda),
                      dist_spec = dist_spec,
                      st_dev_t_min_1 = as.double(st_dev_t_min_1),
                      variance_t_min_1 = as.double(st_dev_t_min_1)^2,
                      residual_t_min_1 = as.double(residual_t_min_1),
                      residual_t_min_1_quadr = as.double(residual_t_min_1)^2,
                      VaR = as.double(VaR),
                      ES = as.double(ES),
                      vcov_matrix_par = if(par_corr) vcov_matrix_par else NA,
                      mu_grad = if(par_corr) mu_grad else NA,
                      sigma_grad = if(par_corr) sigma_grad else NA,
                      emp_dist_vec = if(empirical) emp_dist_vec else NA)
  } else {
    
    # Log complete fail of estimation
    write(paste0('\n\nComplete fail: loop ', current_loop, '\n\n'),
          file = 'run_logs/solver_problems.txt', append = TRUE)
    
    
    results <- tibble(last_date = last_date,
                      mu = NA_real_,
                      sigma = NA_real_,
                      variance = NA_real_,
                      skew = NA_real_,
                      shape = NA_real_,
                      lambda = NA_real_,
                      dist_spec = NA_character_,
                      st_dev_t_min_1 = NA_real_,
                      variance_t_min_1 = NA_real_,
                      residual_t_min_1 = NA_real_,
                      residual_t_min_1_quadr = NA_real_,
                      VaR = NA_real_,
                      ES = NA_real_,
                      vcov_matrix_par = NA,
                      mu_grad = NA,
                      sigma_grad = NA,
                      emp_dist_vec = NA)
  }
  
  results_zoo <- zoo(x = dplyr::select(results, -'last_date'),
                     order.by = results[['last_date']])
  return(results_zoo)
}


#############################################################################
###################     Mincer Regression    ################################
#############################################################################
mincer_regression <- function(formula,
                              shortfall_tbl,
                              h0,
                              white_adjust){
  mincer_reg <- lm(formula = formula,
                   data = shortfall_tbl)
  mincer_reg_result <- car::linearHypothesis(model = mincer_reg,
                                             hypothesis.matrix = h0,
                                             test = 'F',
                                             white.adjust = white_adjust)
  return(mincer_reg_result$'Pr(>F)'[2])
}


##########################################################################
############# Unconditional coverage test for ES #########################
##########################################################################
ES_uc_backtest <- function(CumVio,
                           tolerance_lvl,
                           par_corr,
                           est_window,
                           Return,
                           mu,
                           sigma,
                           dist_spec,
                           skew,
                           shape,
                           lambda,
                           vcov_matrix_par_code,
                           mu_grad_code,
                           sigma_grad_code,
                           empirical,
                           emp_dist_code){
  
  # Extracting non-NA cumulative violations
  H_hut <- CumVio[!is.na(CumVio)]
  
  # Mean and length of H_hut
  H_mean <- mean(H_hut)
  n <- length(H_hut)
  
  # Parameter uncertainty correction
  if(par_corr){
    
    if(empirical){
      
      # Create tibble with required data
      tibble_tmp <- tibble(Return = Return,
                           mu = mu,
                           sigma = sigma,
                           vcov_matrix_par_code = vcov_matrix_par_code,
                           mu_grad_code = mu_grad_code,
                           sigma_grad_code = sigma_grad_code,
                           emp_dist_code = emp_dist_code) %>%
        mutate(eps = (Return - mu) / sigma,
               vcov_matrix_par = map(vcov_matrix_par_code, eval_string_code),
               mu_grad         = map(mu_grad_code, eval_string_code),
               sigma_grad      = map(sigma_grad_code, eval_string_code),
               emp_dist        = map(emp_dist_code, eval_string_code))
      
      # Calculate R_ES
      tibble_tmp <- tibble_tmp %>%
        mutate(R_ES_t = pmap(
          list(eps, mu_grad, sigma_grad, sigma, emp_dist),
          function(eps, mu_grad, sigma_grad, sigma, emp_dist){
            
            # Check for any NA values
            if (anyNA(c(eps, mu_grad, sigma_grad, sigma))) {return(NA_real_)}
            
            # Compute the density of eps
            kernel_dens_eps <- density(emp_dist,
                                       na.rm = TRUE)
            density_eps <- approx(x = kernel_dens_eps[['x']],
                                  y = kernel_dens_eps[['y']],
                                  xout = eps)[['y']]
            
            # Compute indicator: 1 if eps <= quantile of tolerance_lvl, else 0
            indicator <- ifelse(
              eps <= stats::quantile(emp_dist,
                                     probs = tolerance_lvl,
                                     na.rm = TRUE,
                                     names = FALSE,
                                     type = 8), # Approximately median unbiased regardless of distribution
              1, 0)
            
            # Compute gradient expression
            grad_term <- (mu_grad + eps * sigma_grad) / sigma
            
            # Final result
            return(density_eps * indicator * grad_term)
          }
        ))
      
    } else {
      
      # Create tibble with required data
      tibble_tmp <- tibble(Return = Return,
                           mu = mu,
                           sigma = sigma,
                           dist_spec = dist_spec,
                           skew = skew,
                           shape = shape,
                           lambda = lambda,
                           vcov_matrix_par_code = vcov_matrix_par_code,
                           mu_grad_code = mu_grad_code,
                           sigma_grad_code = sigma_grad_code) %>%
        mutate(eps = (Return - mu) / sigma,
               vcov_matrix_par = map(vcov_matrix_par_code, eval_string_code),
               mu_grad         = map(mu_grad_code, eval_string_code),
               sigma_grad      = map(sigma_grad_code, eval_string_code))
      
      # Calculate R_ES
      tibble_tmp <- tibble_tmp %>%
        mutate(R_ES_t = pmap(
          list(eps, dist_spec, skew, shape, lambda, mu_grad, sigma_grad, sigma),
          function(eps, dist_spec, skew, shape, lambda, mu_grad, sigma_grad, sigma){
            
            # Check for any NA values
            if (anyNA(c(eps, dist_spec, mu_grad, sigma_grad, sigma))) {return(NA_real_)}
            
            # Compute the density of eps
            density_eps <- rugarch::ddist(distribution = dist_spec,
                                          y = eps,
                                          skew = skew,
                                          shape = shape,
                                          lambda = lambda)
            
            # Compute indicator: 1 if eps <= quantile of tolerance_lvl, else 0
            indicator <- ifelse(
              eps <= rugarch::qdist(distribution = dist_spec,
                                    p = tolerance_lvl,
                                    skew = skew,
                                    shape = shape,
                                    lambda = lambda),
              1, 0)
            
            # Compute gradient expression
            grad_term <- (mu_grad + eps * sigma_grad) / sigma
            
            # Final result
            return(density_eps * indicator * grad_term)
          }
        ))
    }
    
    # Sum R_ES_t vectors across all rows (ignoring NA entries)
    R_ES_sum <- tibble_tmp %>%
      dplyr::filter(!map_lgl(R_ES_t, ~ anyNA(.))) %>%
      pull(R_ES_t) %>%
      purrr::reduce(`+`)
    
    # Number of R_ES_t without NA
    n_no_NA <- tibble_tmp %>%
      dplyr::filter(!map_lgl(R_ES_t, ~ anyNA(.))) %>%
      nrow()
    
    # Calculate final R_ES
    R_ES <- (1 / (tolerance_lvl * n_no_NA)) * R_ES_sum
    
    # Calculate mean of vcov_matrix_par
    valid_vcov_mat <- tibble_tmp %>%
      dplyr::filter(!map_lgl(vcov_matrix_par, ~ anyNA(.) || is.null(.))) %>%
      pull(vcov_matrix_par)
    W_mean <- purrr::reduce(valid_vcov_mat, `+`) / length(valid_vcov_mat)
    
    # Final calculation of correction factor
    corr_fac <- (n/est_window) * t(as.vector(R_ES)) %*% as.matrix(W_mean) %*% as.vector(R_ES)
    
    ##TEST
    #tibble_uc <<- tibble_tmp
    #W_mean_save_uc <<- W_mean
    #R_ES_save <<- R_ES
    
  } else {
    
    # No correction if estimate = FALSE
    corr_fac <- 0
  }
  
  ##TEST
  #corr_fac_save <<- corr_fac
  
  # Test statistic calculation
  U <- (sqrt(n) * (H_mean - (tolerance_lvl / 2))) / sqrt(tolerance_lvl * ((1 / 3) - (tolerance_lvl / 4)) + corr_fac)
  
  # p-value calculation
  p <- 2 * pnorm(q = abs(U),
                 lower.tail = FALSE)
  
  # Return p-value
  return(as.numeric(p))
}


##########################################################################
############# Conditional coverage test for ES ###########################
##########################################################################
ES_cc_backtest <- function(CumVio,
                           tolerance_lvl,
                           lags,                           
                           par_corr,
                           est_window,
                           Return,
                           mu,
                           sigma,
                           dist_spec,
                           skew,
                           shape,
                           lambda,
                           vcov_matrix_par_code,
                           mu_grad_code,
                           sigma_grad_code,
                           empirical,
                           emp_dist_code){
  
  # Extracting non-NA cumulative violations
  H_hut <- CumVio[!is.na(CumVio)]
  
  # Length of H_hut
  n <- length(H_hut)
  
  # Variance of H_huts
  gamma_n0 <- 1 / n * sum((H_hut - (tolerance_lvl / 2)) * (H_hut - (tolerance_lvl / 2)))
  
  # Vector with covariance between H_hut and j-laged H_hut
  gamma_nj <- vector(length = lags)
  
  for(j in 1:lags){
    gamma_nj[j] <- 1 / (n - j) * sum((H_hut[(j+1):n] - (tolerance_lvl / 2)) * (H_hut[1:(n-j)] - (tolerance_lvl / 2)))
  }
  
  # Vector with correlations between H_hut and j-laged H_hut
  rho_nj <- as.vector(gamma_nj / gamma_n0)
  
  # Parameter uncertainty correction
  if(par_corr){
    
    if(empirical){
      
      # Create tibble with required data
      tibble_tmp <- tibble(Return = Return,
                           mu = mu,
                           sigma = sigma,
                           vcov_matrix_par_code = vcov_matrix_par_code,
                           mu_grad_code = mu_grad_code,
                           sigma_grad_code = sigma_grad_code,
                           emp_dist_code = emp_dist_code,
                           CumVio = CumVio) %>%
        mutate(eps = (Return - mu) / sigma,
               vcov_matrix_par = map(vcov_matrix_par_code, eval_string_code),
               mu_grad         = map(mu_grad_code, eval_string_code),
               sigma_grad      = map(sigma_grad_code, eval_string_code),
               emp_dist        = map(emp_dist_code, eval_string_code))
      
      # Compute R_j columns
      R_j_list <- list()
      for(j in 1:lags){
        
        # Create lagged CumVio column
        CumVio_lag_j <- paste0('CumVio_lag_', j)
        tibble_tmp <- tibble_tmp %>%
          mutate(!!CumVio_lag_j := dplyr::lag(x = CumVio, n = j))
        
        # Create R_j_t column
        R_j_t <- paste0('R_', j, '_t')
        tibble_tmp <- tibble_tmp %>%
          mutate(!!R_j_t := pmap(
            list(eps, mu_grad, sigma_grad, sigma, emp_dist, .data[[CumVio_lag_j]]),
            function(eps, mu_grad, sigma_grad, sigma, emp_dist, CumVio_lag_j) {
              
              # Check for any NA values
              if (anyNA(c(eps, mu_grad, sigma_grad, sigma, CumVio_lag_j))) {return(NA_real_)}
              
              # Compute the empirical density of eps
              kernel_dens_eps <- density(emp_dist,
                                         na.rm = TRUE)
              density_eps <- approx(x = kernel_dens_eps[['x']],
                                    y = kernel_dens_eps[['y']],
                                    xout = eps)[['y']]
              
              # Compute indicator: 1 if eps <= quantile of tolerance_lvl, else 0
              indicator <- ifelse(
                eps <= stats::quantile(emp_dist,
                                       probs = tolerance_lvl,
                                       na.rm = TRUE,
                                       names = FALSE,
                                       type = 8), # Approximately median unbiased regardless of distribution
                1, 0)
              
              # Compute gradient expression
              grad_term <- (mu_grad + eps * sigma_grad) / sigma
              
              # Final result
              return((CumVio_lag_j - (tolerance_lvl / 2)) * density_eps * indicator * grad_term)
            }
          ))
        
        # Sum R_j_t vectors across all rows (ignoring NA entries)
        R_j_t_sum <- tibble_tmp %>%
          dplyr::filter(!map_lgl(.data[[R_j_t]], ~ anyNA(.))) %>%
          pull(.data[[R_j_t]]) %>%
          purrr::reduce(`+`)
        
        # Number of R_j_t without NA
        n_no_NA <- tibble_tmp %>%
          dplyr::filter(!map_lgl(.data[[R_j_t]], ~ anyNA(.))) %>%
          nrow()
        
        # Calculate final R_j
        R_j_list[[j]] <- (1/(tolerance_lvl * (1/3 - tolerance_lvl/4))) * (1 / n_no_NA) * R_j_t_sum
      }
      
    } else {
      
      # Create tibble with required data
      tibble_tmp <- tibble(Return = Return,
                           mu = mu,
                           sigma = sigma,
                           dist_spec = dist_spec,
                           skew = skew,
                           shape = shape,
                           lambda = lambda,
                           vcov_matrix_par_code = vcov_matrix_par_code,
                           mu_grad_code = mu_grad_code,
                           sigma_grad_code = sigma_grad_code,
                           CumVio = CumVio) %>%
        mutate(eps = (Return - mu) / sigma,
               vcov_matrix_par = map(vcov_matrix_par_code, eval_string_code),
               mu_grad         = map(mu_grad_code, eval_string_code),
               sigma_grad      = map(sigma_grad_code, eval_string_code))
      
      # Compute R_j columns
      R_j_list <- list()
      for(j in 1:lags){
        
        # Create lagged CumVio column
        CumVio_lag_j <- paste0('CumVio_lag_', j)
        tibble_tmp <- tibble_tmp %>%
          mutate(!!CumVio_lag_j := dplyr::lag(x = CumVio, n = j))
        
        # Create R_j_t column
        R_j_t <- paste0('R_', j, '_t')
        tibble_tmp <- tibble_tmp %>%
          mutate(!!R_j_t := pmap(
            list(eps, dist_spec, skew, shape, lambda, mu_grad, sigma_grad, sigma, .data[[CumVio_lag_j]]),
            function(eps, dist_spec, skew, shape, lambda, mu_grad, sigma_grad, sigma, CumVio_lag_j) {
              
              # Check for any NA values
              if (anyNA(c(eps, dist_spec, mu_grad, sigma_grad, sigma, CumVio_lag_j))) {return(NA_real_)}
              
              # Compute the density of eps
              density_eps <- rugarch::ddist(distribution = dist_spec,
                                            y = eps,
                                            skew = skew,
                                            shape = shape,
                                            lambda = lambda)
              
              # Compute indicator: 1 if eps <= quantile of tolerance_lvl, else 0
              indicator <- ifelse(
                eps <= rugarch::qdist(distribution = dist_spec,
                                      p = tolerance_lvl,
                                      skew = skew,
                                      shape = shape,
                                      lambda = lambda),
                1, 0)
              
              # Compute gradient expression
              grad_term <- (mu_grad + eps * sigma_grad) / sigma
              
              # Final result
              (CumVio_lag_j - (tolerance_lvl / 2)) * density_eps * indicator * grad_term
            }
          ))
        
        # Sum R_j_t vectors across all rows (ignoring NA entries)
        R_j_t_sum <- tibble_tmp %>%
          dplyr::filter(!map_lgl(.data[[R_j_t]], ~ anyNA(.))) %>%
          pull(.data[[R_j_t]]) %>%
          purrr::reduce(`+`)
        
        # Number of R_j_t without NA
        n_no_NA <- tibble_tmp %>%
          dplyr::filter(!map_lgl(.data[[R_j_t]], ~ anyNA(.))) %>%
          nrow()
        
        # Calculate final R_j
        R_j_list[[j]] <- (1/(tolerance_lvl * (1/3 - tolerance_lvl/4))) * (1 / n_no_NA) * R_j_t_sum
      }
    }
    
    # Calculate mean of vcov_matrix_par
    valid_vcov_mat <- tibble_tmp %>%
      dplyr::filter(!map_lgl(vcov_matrix_par, ~ anyNA(.) || is.null(.))) %>%
      pull(vcov_matrix_par)
    W_mean <- purrr::reduce(valid_vcov_mat, `+`) / length(valid_vcov_mat)
    
    # Create sigma matrix
    sigma_mat <- matrix(nrow = lags, ncol = lags)
    for(i in 1:lags){
      for(j in 1:lags){
        sigma_mat[i,j] <- ifelse(i==j,1,0) + (n/est_window) * (t(as.vector(R_j_list[[i]])) %*% as.matrix(W_mean) %*% as.vector(R_j_list[[j]]))
      }
    }
    
    # Compute final correction matrix
    corr_matrix <- solve(sigma_mat)
    
    ##TEST
    #tibble_cc <<- tibble_tmp
    #W_mean_save_cc <<- W_mean
    #R_j_list_save <<- R_j_list
    
  } else {
    
    # No correction if estimate = FALSE
    corr_matrix <- diag(x = 1, nrow = lags)
  }
  
  ##TEST
  #corr_matrix_save <<- corr_matrix
  
  # Test statistic calculation
  C <- n * t(rho_nj) %*% corr_matrix %*% rho_nj
  
  # p-value calculation
  p <- pchisq(q = C,
              df = lags,
              lower.tail = FALSE)
  
  # Return p-value
  return(as.numeric(p))
}


########################################################################################
################### Functions convert matrix or vector into string and read it #########
########################################################################################
matrix_to_string <- function(m) {
  res_str <- paste0('matrix(nrow = ', nrow(m), ', ncol = ', ncol(m), ', byrow = FALSE, data = c(')
  elements <- as.vector(m)
  res_str <- paste0(res_str, paste(elements, collapse = ", "), '))')
  return(res_str)
}

vector_to_string <- function(v){
  res_vec <- paste0('c(', paste(v, collapse = ", "), ')')
  return(res_vec)
}

eval_string_code <- function(str_code){
  str_res <- eval(parse(text = str_code))
  return(str_res)
}


####################################################################################
################### For parallel loop: Create vector with names of test ############
####################################################################################
names_of_add_tests <- function(n_boot_de,
                               n_boot_esr,
                               estimate,
                               par_corr){
  
  add_tst <- c('UC', 'CC', 'Auxiliary_ESR', 'Strict_ESR')
  if(!is.null(n_boot_de)){
    add_tst <- c(add_tst, 'UC_boot', 'CC_boot')
  }
  
  if(!is.null(n_boot_esr)){
    add_tst <- c(add_tst, 'Auxiliary_ESR_boot', 'Strict_ESR_boot')
  }
  
  if(estimate & par_corr){
    add_tst <- c(add_tst, 'UC_par_corr', 'CC_par_corr')
  }
  
  add_tst <- c(sort(add_tst[grepl('^UC', add_tst)]),
               sort(add_tst[grepl('^CC', add_tst)]),
               sort(add_tst[grepl('^Auxiliary_ESR', add_tst)]),
               sort(add_tst[grepl('^Strict_ESR', add_tst)]))
  return(add_tst)
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


#############################################################################
################### Parallel estimation Loop ################################
#############################################################################
estimation_loop_par <- function(n_loop,
                                est_window,
                                oos_window,
                                tolerance_lvl,
                                var_spec_sim,
                                mean_spec_sim,
                                dist_spec_sim,
                                fixed_pars_sim,
                                estimate=TRUE,
                                par_corr=FALSE,
                                var_spec_est,
                                mean_spec_est,
                                dist_spec_est,
                                fixed_pars_est=NULL,
                                cores=1,
                                white_adjust='hc3',
                                seed=NULL,
                                mincer_spec,
                                execute_additional_tsts=TRUE,
                                lags_es_cc=NULL,
                                max_lag_es_cc=NULL,
                                n_boot_de=NULL,
                                n_boot_esr=NULL,
                                empirical=FALSE){
  
  # Clear started_loops.txt and finished_loops.txt
  write(x = ' ',
        file = 'run_logs/started_loops.txt',
        append = FALSE)
  write(x = ' ',
        file = 'run_logs/finished_loops.txt',
        append = FALSE)
  
  # Register cores
  registerDoParallel(cores = cores)
  
  # Estimation loop
  result_foreach <- foreach(i=1:n_loop) %dopar% {
    
    # Register started loops
    write(x = paste0(i),
          file = 'run_logs/started_loops.txt',
          append = TRUE)
    
    # Set seed for reproducibility of results
    if(!is.null(seed)){
      current_seed <- (seed+i^2-i) %% 2147483647
      set.seed(current_seed)
    } else {
      current_seed <- NULL
    }
    
    # Progress messages
    if(i %% cores == 0){
      cat(paste0('Loop ', i, ' started at ', Sys.time(), '\n'))
      write(x = paste0('Loop ', i, ' started at ', Sys.time(), '\n'),
            file = 'run_logs/progress_information.txt',
            append = TRUE)
    }
    
    # Simulate garch process
    garchspec_sim <- ugarchspec(variance.model = var_spec_sim,
                                mean.model = mean_spec_sim,
                                distribution.model = dist_spec_sim,
                                fixed.pars = fixed_pars_sim)
    
    garchsimulation <- ugarchpath(spec = garchspec_sim,
                                  n.sim = est_window + oos_window,
                                  n.start = 1000,
                                  m.sim = 1)
    
    # Create tibble and zoo object of simulation
    sim_matrix <- garchsimulation@path$seriesSim
    colnames(sim_matrix) <- 'Return'
    garchsimulation_tbl <- as_tibble(sim_matrix) %>%
      mutate(Date = base::as.Date('2000-01-01') + 1:nrow(sim_matrix))
    garchsimulation_zoo <- zoo(x = garchsimulation_tbl[['Return']],
                               order.by = garchsimulation_tbl[['Date']])
    
    ##TEST
    # garchsimulation_zoo_s <<- garchsimulation_zoo
    # plot_s <<- plot(garchsimulation_zoo)
    # Sys.sleep(time = 2)
    # plot_s2 <<- garchsimulation_zoo %>% as_tibble() %>%
    #   mutate(Price = 100 * exp(cumsum(value))) %>%
    #   select(Price) %>%
    #   zoo(order.by = index(garchsimulation_zoo)) %>%
    #   plot()
    
    ##TEST
    #garchsimulation_zoo_s <<- garchsimulation_zoo
    
    # Calculate GARCH model
    rolling_VaR_ES <- rollapply(data = garchsimulation_zoo,
                                width = est_window,
                                FUN = function(x) VaR_ES_forecast(data_zoo = x,
                                                                  var_spec = var_spec_est,
                                                                  mean_spec = mean_spec_est,
                                                                  dist_spec = dist_spec_est,
                                                                  tolerance_lvl = tolerance_lvl,
                                                                  estimate = estimate,
                                                                  par_corr = par_corr,
                                                                  fixed_pars = fixed_pars_est,
                                                                  empirical = empirical,
                                                                  seed_opt = current_seed,
                                                                  current_loop = i),
                                align = 'right',
                                coredata = FALSE)
    
    rolling_VaR_ES_lead <- stats::lag(rolling_VaR_ES,
                                      k = -1)
    rolling_VaR_ES_tbl <- as_tibble(rolling_VaR_ES_lead) %>%
      mutate(across(-c(dist_spec, vcov_matrix_par, mu_grad, sigma_grad, emp_dist_vec), as.numeric),
             Date = index(rolling_VaR_ES_lead))
    VaR_ES_results_tbl <- inner_join(garchsimulation_tbl, rolling_VaR_ES_tbl, by = 'Date')
    
    ##TEST
    #VaR_ES_results_tbl_save <<- VaR_ES_results_tbl
    
    # Log NAs
    cols_check <- c('Return', 'VaR', 'ES')
    na_counts <- VaR_ES_results_tbl %>%
      dplyr::summarise(dplyr::across(all_of(cols_check), ~sum(is.na(.))))
    
    if(any(na_counts > 0)){
      write(paste0('\nNumber of NAs in loop ', i, ': ', paste(names(na_counts), na_counts, collapse = ', '), '\n'),
            file = 'run_logs/NA_information.txt',
            append = TRUE)
    }
    
    # Filter for shortfalls (NAs are removed automatically)
    shortfall_tbl <- VaR_ES_results_tbl %>%
      dplyr::filter(Return <= VaR) %>%
      mutate(shortfall = Return - ES,
             residual_t_min_1_quadr_lower_0 = ifelse(residual_t_min_1 < 0, residual_t_min_1_quadr, 0),
             indicator_residual_lower_0 = ifelse(residual_t_min_1 < 0, 1, 0))
    
    # Save number of observations that go into Mincer regression
    n_obs_mincer <- nrow(shortfall_tbl)
    
    # Loop over specified white adjustments
    for(i_wa in 1:length(white_adjust)){
      wa <- white_adjust[i_wa]
      
      # Mincer regression execution
      p_values_vec <- vector()
      for(j in 1:length(mincer_spec)){
        p_values_vec[j] <- tryCatch(
          {
            mincer_regression(formula = mincer_spec[[j]][['formula']],
                              shortfall_tbl = shortfall_tbl,
                              h0 = mincer_spec[[j]][['h0']],
                              white_adjust = wa)
          }, error = function(e1){
            
            # Log fail of first try
            write(paste(tolerance_lvl, var_spec_sim[['model']], paste(mean_spec_sim, collapse = ' '), dist_spec_sim, 'loop', i, '\n', paste(mincer_spec[[j]][['h0']], collapse = ' ')),
                  file = 'run_logs/mincer_fail_1.txt', append = TRUE)
            
            tryCatch(
              {
                # Try again without white adjustment
                mincer_regression(formula = mincer_spec[[j]][['formula']],
                                  shortfall_tbl = shortfall_tbl,
                                  h0 = mincer_spec[[j]][['h0']],
                                  white_adjust = FALSE)
              }, error = function(e2){
                
                # Log fail of second try
                write(paste(tolerance_lvl, var_spec_sim[['model']], paste(mean_spec_sim, collapse = ' '), dist_spec_sim, 'loop', i, '\n', paste(mincer_spec[[j]][['h0']], collapse = ' ')),
                      file = 'run_logs/mincer_fail_2.txt', append = TRUE)
                # Return NA if error still occures
                NA_real_
              }
            )
          }
        )
      }
      
      # Save results
      if(i_wa == 1){
        result_lst <- list()
        result_lst_names <- names(mincer_spec)
      }
      
      for(j in 1:length(mincer_spec)){
        result_lst[[result_lst_names[j]]][[paste0('p_', wa)]] <- p_values_vec[j]
        result_lst[[result_lst_names[j]]][[paste0('p0_01_', wa)]] <- ifelse(p_values_vec[j] <= 0.01, 1, 0)
        result_lst[[result_lst_names[j]]][[paste0('p0_05_', wa)]] <- ifelse(p_values_vec[j] <= 0.05, 1, 0)
        result_lst[[result_lst_names[j]]][[paste0('p0_1_', wa)]] <- ifelse(p_values_vec[j] <= 0.1, 1, 0)
      }
      
      # Compute additional tests if execute_additional_tsts = TRUE
      if(execute_additional_tsts){
        
        # Compute cdf value at return
        if(empirical){
          
          # With empirical distribution
          u <- VaR_ES_results_tbl %>%
            mutate(standardized_Return = (Return - mu) / sigma,
                   emp_dist = map(emp_dist_vec, eval_string_code)) %>%
            mutate(u = map2_dbl(.x = emp_dist,
                                .y = standardized_Return,
                                .f = ~ mean(.x <= .y))) %>%
            pull(u)
          
        } else {
          
          # With parametric distribution assumption
          u <- rugarch::pdist(distribution = dist_spec_est,
                              q = VaR_ES_results_tbl[['Return']],
                              mu = VaR_ES_results_tbl[['mu']],
                              sigma = VaR_ES_results_tbl[['sigma']],
                              skew = VaR_ES_results_tbl[['skew']],
                              shape = VaR_ES_results_tbl[['shape']],
                              lambda = VaR_ES_results_tbl[['lambda']])
        }
        
        ##TEST
        #u_save <<- u
        
        # Compute cumulative violations
        CumVio <- (1 / tolerance_lvl) * (tolerance_lvl - u) * ifelse(u <= tolerance_lvl, 1, 0)
        
        # Log NAs in u and CumVio
        NAs_u <- sum(is.na(u))
        NAs_CumVio <- sum(is.na(CumVio))
        
        if((NAs_u > 0) || (NAs_CumVio > 0)){
          write(paste0('\nNumber of NAs in loop ', i, ': u ', NAs_u, ' , CumVio ', NAs_CumVio, '\n'),
                file = 'run_logs/NA_information.txt',
                append = TRUE)
        }
        
        ##TEST
        # u_save <<- u
        # CumVio_save <<- CumVio
        # print(if(!is.null(lags_es_cc)) lags_es_cc else auto_lag(y = CumVio,
        #                                                   max_lag = max_lag_es_cc))
        # plot(acf(CumVio, lag.max = 15))
        
        # Create vector with names of additional tests
        add_tests <- names_of_add_tests(n_boot_de = n_boot_de,
                                        n_boot_esr = n_boot_esr,
                                        estimate = estimate,
                                        par_corr = par_corr)
        p_add <- vector(length = length(add_tests))
        names(p_add) <- add_tests
        
        # Execute unconditional and conditional coverage backtest for ES
        p_add['UC'] <- ES_uc_backtest(CumVio = CumVio,
                                      tolerance_lvl = tolerance_lvl,
                                      par_corr = FALSE)
        
        p_add['CC'] <- ES_cc_backtest(CumVio = CumVio,
                                      tolerance_lvl = tolerance_lvl,
                                      lags = if(!is.null(lags_es_cc)) lags_es_cc else auto_lag(y = CumVio,
                                                                                               max_lag = max_lag_es_cc),
                                      par_corr = FALSE)
        
        # Execute unconditional and conditional coverage backtest with bootstrap for ES
        if(!is.null(n_boot_de)){
          if(!is.null(seed)){set.seed(current_seed)}
          shortfall_de_test_boot_res <- tstests::shortfall_de_test(x = u[!is.na(u)],
                                                                   alpha = tolerance_lvl,
                                                                   lags = if(!is.null(lags_es_cc)) lags_es_cc else auto_lag(y = CumVio,
                                                                                                                            max_lag = max_lag_es_cc),
                                                                   boot = TRUE,
                                                                   n_boot = n_boot_de)
          
          p_add['UC_boot'] <- shortfall_de_test_boot_res[['p_value']][1]
          p_add['CC_boot'] <- shortfall_de_test_boot_res[['p_value']][2]
        }
        
        # Execute unconditional and conditional coverage backtest with parameter uncertainty correction for ES
        if(estimate & par_corr){
          p_add['UC_par_corr'] <- ES_uc_backtest(CumVio = CumVio,
                                                 tolerance_lvl = tolerance_lvl,
                                                 par_corr = TRUE,
                                                 est_window = est_window,
                                                 Return = VaR_ES_results_tbl[['Return']],
                                                 mu = VaR_ES_results_tbl[['mu']],
                                                 sigma = VaR_ES_results_tbl[['sigma']],
                                                 dist_spec = VaR_ES_results_tbl[['dist_spec']],
                                                 skew = VaR_ES_results_tbl[['skew']],
                                                 shape = VaR_ES_results_tbl[['shape']],
                                                 lambda = VaR_ES_results_tbl[['lambda']],
                                                 vcov_matrix_par_code = VaR_ES_results_tbl[['vcov_matrix_par']],
                                                 mu_grad_code = VaR_ES_results_tbl[['mu_grad']],
                                                 sigma_grad_code = VaR_ES_results_tbl[['sigma_grad']],
                                                 empirical = empirical,
                                                 emp_dist_code = VaR_ES_results_tbl[['emp_dist_vec']])
          
          p_add['CC_par_corr'] <- ES_cc_backtest(CumVio = CumVio,
                                                 tolerance_lvl = tolerance_lvl,
                                                 lags = if(!is.null(lags_es_cc)) lags_es_cc else auto_lag(y = CumVio,
                                                                                                          max_lag = max_lag_es_cc),
                                                 par_corr = TRUE,
                                                 est_window = est_window,
                                                 Return = VaR_ES_results_tbl[['Return']],
                                                 mu = VaR_ES_results_tbl[['mu']],
                                                 sigma = VaR_ES_results_tbl[['sigma']],
                                                 dist_spec = VaR_ES_results_tbl[['dist_spec']],
                                                 skew = VaR_ES_results_tbl[['skew']],
                                                 shape = VaR_ES_results_tbl[['shape']],
                                                 lambda = VaR_ES_results_tbl[['lambda']],
                                                 vcov_matrix_par_code = VaR_ES_results_tbl[['vcov_matrix_par']],
                                                 mu_grad_code = VaR_ES_results_tbl[['mu_grad']],
                                                 sigma_grad_code = VaR_ES_results_tbl[['sigma_grad']],
                                                 empirical = empirical,
                                                 emp_dist_code = VaR_ES_results_tbl[['emp_dist_vec']])
        }
        
        # Create VaR_ES_results_tbl without NAs in Return, VaR and ES
        #if(!exists('VaR_ES_results_tbl_noNA')){#WHY do I have this? Test if I can remove it
          VaR_ES_results_tbl_noNA <- VaR_ES_results_tbl %>%
            tidyr::drop_na(Return, VaR, ES)
        #}
        
        # Execute Auxiliary ESR backtest
        if(!is.null(seed)){set.seed(current_seed)}
        p_add['Auxiliary_ESR'] <- esback::esr_backtest(r = VaR_ES_results_tbl_noNA[['Return']],
                                                       q = VaR_ES_results_tbl_noNA[['VaR']],
                                                       e = VaR_ES_results_tbl_noNA[['ES']],
                                                       alpha = tolerance_lvl,
                                                       version = 2)[['pvalue_twosided_asymptotic']]
        
        # Execute Strict ESR backtest
        if(!is.null(seed)){set.seed(current_seed)}
        p_add['Strict_ESR'] <- esback::esr_backtest(r = VaR_ES_results_tbl_noNA[['Return']],
                                                    e = VaR_ES_results_tbl_noNA[['ES']],
                                                    alpha = tolerance_lvl,
                                                    version = 1)[['pvalue_twosided_asymptotic']]
        
        # Execute ESR tests with bootstrapping if n_boot_esr is not NULL
        if(!is.null(n_boot_esr)){
          
          # Execute Auxiliary ESR backtest with bootstrap
          if(!is.null(seed)){set.seed(current_seed)}
          p_add['Auxiliary_ESR_boot'] <- tryCatch(
            {
              esback::esr_backtest(r = VaR_ES_results_tbl_noNA[['Return']],
                                   q = VaR_ES_results_tbl_noNA[['VaR']],
                                   e = VaR_ES_results_tbl_noNA[['ES']],
                                   alpha = tolerance_lvl,
                                   version = 2,
                                   B = n_boot_esr)[['pvalue_twosided_bootstrap']]
            }, error = function(e1){
              
              # Retry with different seed
              tryCatch(
                {
                  if(!is.null(seed)){set.seed(current_seed+11)}
                  esback::esr_backtest(r = VaR_ES_results_tbl_noNA[['Return']],
                                       q = VaR_ES_results_tbl_noNA[['VaR']],
                                       e = VaR_ES_results_tbl_noNA[['ES']],
                                       alpha = tolerance_lvl,
                                       version = 2,
                                       B = n_boot_esr)[['pvalue_twosided_bootstrap']]
                }, error = function(e2){
                  # Log fail of second try
                  write(paste(tolerance_lvl, var_spec_sim[['model']], paste(mean_spec_sim, collapse = ' '), dist_spec_sim, 'loop', i, '\n'),
                        file = 'run_logs/auxiliary_esr_boot_fail.txt', append = TRUE)
                  
                  # Return NA if error appears again
                  NA_real_
                }
              )
            }
          )
          
          # Execute Strict ESR backtest with bootstrap
          if(!is.null(seed)){set.seed(current_seed)}
          p_add['Strict_ESR_boot'] <- tryCatch(
            {
              esback::esr_backtest(r = VaR_ES_results_tbl_noNA[['Return']],
                                   e = VaR_ES_results_tbl_noNA[['ES']],
                                   alpha = tolerance_lvl,
                                   version = 1,
                                   B = n_boot_esr)[['pvalue_twosided_bootstrap']]
            }, error = function(e1){
              
              # Retry with different seed
              tryCatch(
                {
                  if(!is.null(seed)){set.seed(current_seed+11)}
                  esback::esr_backtest(r = VaR_ES_results_tbl_noNA[['Return']],
                                       e = VaR_ES_results_tbl_noNA[['ES']],
                                       alpha = tolerance_lvl,
                                       version = 1,
                                       B = n_boot_esr)[['pvalue_twosided_bootstrap']]
                }, error = function(e2){
                  # Log fail of second try
                  write(paste(tolerance_lvl, var_spec_sim[['model']], paste(mean_spec_sim, collapse = ' '), dist_spec_sim, 'loop', i, '\n'),
                        file = 'run_logs/strict_esr_boot_fail.txt', append = TRUE)
                  
                  # Return NA if error appears again
                  NA_real_
                }
              )
            }
          )
        }
        
        # Store result
        for(add_tst in add_tests){
          result_lst[[add_tst]][[paste0('p_', wa)]] <- p_add[add_tst]
          result_lst[[add_tst]][[paste0('p0_01_', wa)]] <- ifelse(p_add[add_tst] <= 0.01, 1, 0)
          result_lst[[add_tst]][[paste0('p0_05_', wa)]] <- ifelse(p_add[add_tst] <= 0.05, 1, 0)
          result_lst[[add_tst]][[paste0('p0_1_', wa)]] <- ifelse(p_add[add_tst] <= 0.1, 1, 0)
        }
      }
    }
    
    # Save seed number and number of observations that entered Mincer regressions
    result_lst[['n_obs_mincer']] <- n_obs_mincer
    result_lst[['seed']] <- if(!is.null(seed)) current_seed else NA_integer_
    
    # Progress messages
    if(i %% cores == 0){
      cat(paste0('Loop ', i, ' finished at ', Sys.time(), '\n'))
      write(x = paste0('Loop ', i, ' finished at ', Sys.time(), '\n'),
            file = 'run_logs/progress_information.txt',
            append = TRUE)
    }
    
    # Register finished loops
    write(x = paste0(i),
          file = 'run_logs/finished_loops.txt',
          append = TRUE)
    
    # Result list
    result_lst
  }
  
  ##TEST
  # saveRDS(object = result_foreach,
  #         file = 'result_foreach.RData')
  
  # Organize results
  result <- list()
  
  # Create vector with names of executed tests
  if(execute_additional_tsts){
    add_tests <- names_of_add_tests(n_boot_de = n_boot_de,
                                    n_boot_esr = n_boot_esr,
                                    estimate = estimate,
                                    par_corr = par_corr)
    test_names <- c(names(mincer_spec), add_tests)
  } else {
    test_names <- names(mincer_spec)
  }
  
  for(test_name in test_names){
    result[[test_name]] <- list()
    
    for(wa in white_adjust){
      keys <- c(paste0('p_', wa),
                paste0('p0_01_', wa),
                paste0('p0_05_', wa),
                paste0('p0_1_', wa))
      
      # Assign empty vectors to dynamically named list elements
      result[[test_name]][keys] <- replicate(length(keys), vector(), simplify = FALSE)
    }
  }
  result[['n_obs_mincer']] <- vector()
  result[['seed']] <- vector()
  
  for(i in 1:n_loop){
    for(method in names(result_foreach[[i]])){
      if(method %in% c('n_obs_mincer', 'seed')){
        result[[method]] <- c(result[[method]], result_foreach[[i]][[method]])
      } else {
        for(metrix in names(result_foreach[[i]][[method]])){
          result[[method]][[metrix]] <- c(result[[method]][[metrix]], result_foreach[[i]][[method]][[metrix]])
        }
      }
    }
  }
  return(result)
}


##############################################################################
###  Kupiec test: Unconditional Coverage VaR  ################################
##############################################################################
VaR_Kupiec_backtest <- function(n_VaR_exceeded,
                                oos_window,
                                tolerance_lvl){
  
  # Calculating number of non-exceedences, exceedences and proportion of exceedences
  n1 <- n_VaR_exceeded
  n0 <- oos_window - n1
  prop_exceeded <- n1 / (n1 + n0)
  
  # Likelihood Ratio Test
  # Likelihood of unrestricted (L_ur) and restricted (L_r) model
  L_ur <- prop_exceeded ^ n1 * (1 - prop_exceeded) ^ n0
  L_r <- tolerance_lvl ^ n1 * (1 - tolerance_lvl) ^ n0
  
  # Test statistic and assign name speci_dist
  LR <- -2 * log(L_r / L_ur)
  
  # p value (LR ~ X^2(1)) and assign name speci_dist
  p <- pchisq(q = LR,
              df = 1,
              lower.tail = FALSE)
  
  #Return results
  results <- list(LR = LR,
                  p = p)
  return(results)
}