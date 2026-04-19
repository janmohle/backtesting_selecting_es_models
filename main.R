# Before first run:
# renv::restore()

# Initial cleaning
rm(list = ls())
if (dev.cur() != 1) {
  dev.off()
}
cat('\14')

# Check dependencies
required_pkgs <- c(
  'tibble',
  'dplyr',
  'purrr',
  'rugarch',
  'zoo',
  'foreach',
  'doParallel',
  'tidyr',
  'tstests',
  'esback',
  'numDeriv')

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if(length(missing_pkgs) > 0){
  stop('Missing packages: ', paste(missing_pkgs, collapse=', '),
       '\nRun renv::restore()', call.=FALSE)
} else {
  cat('\nDependencies are loaded\n')
}

# Library packages and load functions
library(tibble)
library(dplyr)
library(purrr)
library(rugarch)
library(zoo)
library(foreach)
library(doParallel)
source('functions.R')


#############################################
########### Global parameters ###############
#############################################
tolerance_lvls = c(0.05, 0.025, 0.01)
n_loop = 1000
oos_window = 5000 #2500
est_window = 750
cores = 100
seed = 78
size_or_power = 'both' # size, power or both
sim_size = 1:8 #c(1,2,5,6,7,8)
sim_power = 2:8
est = 1
white_adjust = 'hc3' #c('hc3', FALSE)
execute_additional_tsts = TRUE
lags_es_cc = NULL #if NULL, then automatic lag selection
max_lag_es_cc = 10 #NULL #only effective in lags_es_cc=NULL, if max_lag_es_cc=NULL then max_lag_es_cc=floor(sqrt(n))
n_boot_de = 10000 #NULL
n_boot_esr = NULL #NULL
estimation = TRUE
par_corr = FALSE
empirical = FALSE


########################################################
########### Model Specifications #######################
########################################################

# Spec 1
var_spec_1 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_1 = list(armaOrder = c(0,0))
dist_spec_1 = 'norm'
fixed_pars_1 = list(mu = 0,
                    omega = 0.00001,
                    alpha1 = 0.1,
                    beta1 = 0.85)

# Spec 2
var_spec_2 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_2 = list(armaOrder = c(0,0))
dist_spec_2 = 'std'
fixed_pars_2 = list(mu = 0,
                    omega = 0.00001,
                    alpha1 = 0.1,
                    beta1 = 0.85,
                    shape = 5)

# Spec 3
var_spec_3 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_3 = list(armaOrder = c(0,0))
dist_spec_3 = 'snorm'
fixed_pars_3 = list(mu = 0,
                    omega = 0.00001,
                    alpha1 = 0.1,
                    beta1 = 0.85,
                    skew = 0.85)

# Spec 4
var_spec_4 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_4 = list(armaOrder = c(1,0))
dist_spec_4 = 'norm'
fixed_pars_4 = list(mu = 0,
                    omega = 0.00001,
                    alpha1 = 0.1,
                    beta1 = 0.85,
                    ar1 = 0.1)

# Spec 5
var_spec_5 = list(model = 'gjrGARCH',garchOrder = c(1,1))
mean_spec_5 = list(armaOrder = c(0,0))
dist_spec_5 = 'norm'
fixed_pars_5 = list(mu = 0,
                    omega = 0.00001,
                    alpha1 = 0.04,
                    beta1 = 0.85,
                    gamma1 = 0.12)

# Spec 6
var_spec_6 = list(model = 'gjrGARCH',garchOrder = c(1,1))
mean_spec_6 = list(armaOrder = c(0,0))
dist_spec_6 = 'std'
fixed_pars_6 = list(mu = 0,
                    omega = 0.00001,
                    alpha1 = 0.04,
                    beta1 = 0.85,
                    gamma1 = 0.12,
                    shape = 5)

# Spec 7
var_spec_7 = list(model = 'eGARCH',garchOrder = c(1,1))
mean_spec_7 = list(armaOrder = c(0,0))
dist_spec_7 = 'norm'
fixed_pars_7 = list(mu = 0,
                    omega = -0.447,#-0.305
                    alpha1 = -0.077,#-0.182
                    beta1 = 0.948,#0.969
                    gamma1 = 0.191)#0.187

# Spec 8
var_spec_8 = list(model = 'eGARCH',garchOrder = c(1,1))
mean_spec_8 = list(armaOrder = c(0,0))
dist_spec_8 = 'std'
fixed_pars_8 = list(mu = 0,
                    omega = -0.447,#-0.305
                    alpha1 = -0.077,#-0.182
                    beta1 = 0.948,#0.969
                    gamma1 = 0.191,#0.187
                    shape = 5)

###########################################################
########### Mincer Regression Specifications ##############
###########################################################

mincer_spec <- list(simple_shortfall = list(formula = shortfall ~ 1,# initially considered in paper (but not anymore)
                                            h0 = c('(Intercept) = 0')),
                    simple_return = list(formula = Return ~ ES,
                                         h0 = c('(Intercept) = 0', 'ES = 1')),
                    variance_shortfall = list(formula = shortfall ~ variance_t_min_1,
                                              h0 = c('(Intercept) = 0', 'variance_t_min_1 = 0')),
                    variance_return = list(formula = Return ~ variance_t_min_1 + ES,
                                           h0 = c('(Intercept) = 0', 'variance_t_min_1 = 0', 'ES = 1')),
                    residual_sqrt_shortfall = list(formula = shortfall ~ residual_t_min_1_quadr,
                                                   h0 = c('(Intercept) = 0', 'residual_t_min_1_quadr = 0')),
                    residual_sqrt_return = list(formula = Return ~ residual_t_min_1_quadr + ES,
                                                h0 = c('(Intercept) = 0', 'residual_t_min_1_quadr = 0', 'ES = 1')),
                    indicator_lower_0_shortfall = list(formula = shortfall ~ indicator_residual_lower_0,# initially considered in paper (but not anymore)
                                                       h0 = c('(Intercept) = 0', 'indicator_residual_lower_0 = 0')),
                    indicator_lower_0_return = list(formula = Return ~ indicator_residual_lower_0 + ES,
                                                    h0 = c('(Intercept) = 0', 'indicator_residual_lower_0 = 0', 'ES = 1')),
                    residual_sqrt0_shortfall = list(formula = shortfall ~ residual_t_min_1_quadr_lower_0,
                                                    h0 = c('(Intercept) = 0', 'residual_t_min_1_quadr_lower_0 = 0')),
                    residual_sqrt0_return = list(formula = Return ~ residual_t_min_1_quadr_lower_0 + ES,
                                                 h0 = c('(Intercept) = 0', 'residual_t_min_1_quadr_lower_0 = 0', 'ES = 1'))#,
                    #full_shortfall = list(formula = shortfall ~ variance_t_min_1 + residual_t_min_1_quadr + residual_t_min_1_quadr_lower_0,
                    #                      h0 = c('(Intercept) = 0', 'variance_t_min_1 = 0', 'residual_t_min_1_quadr = 0', 'residual_t_min_1_quadr_lower_0 = 0')),
                    #full_return = list(formula = Return ~ variance_t_min_1 + residual_t_min_1_quadr + residual_t_min_1_quadr_lower_0 + ES,
                    #                   h0 = c('(Intercept) = 0', 'variance_t_min_1 = 0', 'residual_t_min_1_quadr = 0', 'residual_t_min_1_quadr_lower_0 = 0', 'ES = 1'))
)

# Create required folders
for(dir_name in c('run_logs', 'results')){
  dir.create(dir_name, showWarnings = FALSE)
}
rm(dir_name)

if((size_or_power == 'size') || (size_or_power == 'both')){
  
  
  ##########################################
  ########### Size loops execution #########
  ##########################################
  
  # Initiate progress_information.txt
  write(paste0('Size execution started at ', Sys.time(), '\n\n'),
        file = 'run_logs/progress_information.txt',
        append = FALSE)
  cat(paste0('\nSize execution started at ', Sys.time(), '\n\n'))
  
  # Initiate solver_problems.txt
  write('Size solver problems (Check functions.R for location of problems encoded in numbers):\n\n',
        file = 'run_logs/solver_problems.txt',
        append = FALSE)
  
  # Initiate NA_information.txt
  write('Size number of NAs:\n\n',
        file = 'run_logs/NA_information.txt',
        append = FALSE)
  
  for(sim_i in sim_size){
    
    # Log current simulation
    for(file in c('progress_information', 'solver_problems', 'NA_information')){
      write(paste0('\nSimulation ', sim_i, ':\n'),
            file = paste0('run_logs/', file,'.txt'),
            append = TRUE)
    }
    cat(paste0('\nSimulation ', sim_i,':\n'))
    
    for(tolerance_lvl in tolerance_lvls){
      
      # Log current tolerance level
      for(file in c('progress_information', 'solver_problems', 'NA_information')){
        write(paste0('\nTolerance level ', tolerance_lvl, ':\n'),
              file = paste0('run_logs/', file,'.txt'),
              append = TRUE)
      }
      cat(paste0('\nTolerance level ', tolerance_lvl, ' started at ', Sys.time(),'\n'))
      
      result_lst <- estimation_loop_par(n_loop = n_loop,
                                        est_window = est_window,
                                        oos_window = oos_window,
                                        tolerance_lvl = tolerance_lvl,
                                        var_spec_sim = get(paste0('var_spec_', sim_i)),
                                        mean_spec_sim = get(paste0('mean_spec_', sim_i)),
                                        dist_spec_sim = get(paste0('dist_spec_', sim_i)),
                                        fixed_pars_sim = get(paste0('fixed_pars_', sim_i)),
                                        estimate = estimation,
                                        par_corr = par_corr,
                                        var_spec_est = get(paste0('var_spec_', sim_i)),
                                        mean_spec_est = get(paste0('mean_spec_', sim_i)),
                                        dist_spec_est = get(paste0('dist_spec_', sim_i)),
                                        fixed_pars_est = get(paste0('fixed_pars_', sim_i)),
                                        cores = cores,
                                        white_adjust = white_adjust,
                                        seed = seed,
                                        mincer_spec = mincer_spec,
                                        execute_additional_tsts = execute_additional_tsts,
                                        lags_es_cc = lags_es_cc,
                                        max_lag_es_cc = max_lag_es_cc,
                                        n_boot_de = n_boot_de,
                                        n_boot_esr = n_boot_esr,
                                        empirical = empirical)
      
      # Number of NAs in UC_par_corr and CC_par_corr printed to NA_information.txt
      if(estimation && par_corr){
        write(paste0('Number of NAs in UC_par_corr: ', sum(is.na(result_lst[['UC_par_corr']][[white_adjust[1]]])), '\n'),
              file = 'run_logs/NA_information.txt',
              append = TRUE)
        write(paste0('Number of NAs in CC_par_corr: ', sum(is.na(result_lst[['CC_par_corr']][[white_adjust[1]]])), '\n'),
              file = 'run_logs/NA_information.txt',
              append = TRUE)
      }
      
      saveRDS(result_lst,
              file = paste0('results/size_', gsub(pattern = '\\.', replacement = '_', paste0(tolerance_lvl)), '_sim_', sim_i, '_est_', sim_i, ifelse(estimation, '_est', '_fix'), '.rds'))
    }
  }
  
  write(paste0('\n\nSize execution finished at ', Sys.time(), '\n\n'),
        file = 'run_logs/progress_information.txt',
        append = TRUE)
  cat(paste0('\n\nSize execution finished at ', Sys.time(), '\n\n'))
}


if((size_or_power == 'power') || (size_or_power == 'both')){
  
  
  ##########################################
  ########### Power loops execution ########
  ##########################################
  
  # Initiate progress_information.txt
  write(paste0('Power execution started at ', Sys.time(), '\n\n'),
        file = 'run_logs/progress_information.txt',
        append = (size_or_power == 'both'))
  cat(paste0('\nPower execution started at ', Sys.time(), '\n\n'))
  
  # Initiate solver_problems.txt
  write('Power solver problems (Check functions.R for location of problems encoded in numbers):\n\n',
        file = 'run_logs/solver_problems.txt',
        append = (size_or_power == 'both'))
  
  # Initiate NA_information.txt
  write('Power number of NAs:\n\n',
        file = 'run_logs/NA_information.txt',
        append = (size_or_power == 'both'))
  
  for(sim_i in sim_power){
    
    # Log current simulation
    for(file in c('progress_information', 'solver_problems', 'NA_information')){
      write(paste0('\nSimulation ', sim_i, ':\n'),
            file = paste0('run_logs/', file,'.txt'),
            append = TRUE)
    }
    cat(paste0('\nSimulation ', sim_i,':\n'))
    
    for(tolerance_lvl in tolerance_lvls){
      
      # Log current tolerance level
      for(file in c('progress_information', 'solver_problems', 'NA_information')){
        write(paste0('\nTolerance level ', tolerance_lvl, ':\n'),
              file = paste0('run_logs/', file,'.txt'),
              append = TRUE)
      }
      cat(paste0('\nTolerance level ', tolerance_lvl, ' started at ', Sys.time(),'\n'))
      
      result_lst <- estimation_loop_par(n_loop = n_loop,
                                        est_window = est_window,
                                        oos_window = oos_window,
                                        tolerance_lvl = tolerance_lvl,
                                        var_spec_sim = get(paste0('var_spec_', sim_i)),
                                        mean_spec_sim = get(paste0('mean_spec_', sim_i)),
                                        dist_spec_sim = get(paste0('dist_spec_', sim_i)),
                                        fixed_pars_sim = get(paste0('fixed_pars_', sim_i)),
                                        estimate = estimation,
                                        par_corr = par_corr,
                                        var_spec_est = get(paste0('var_spec_', est)),
                                        mean_spec_est = get(paste0('mean_spec_', est)),
                                        dist_spec_est = get(paste0('dist_spec_', est)),
                                        fixed_pars_est = get(paste0('fixed_pars_', est)),
                                        cores = cores,
                                        white_adjust = white_adjust,
                                        seed = seed,
                                        mincer_spec = mincer_spec,
                                        execute_additional_tsts = execute_additional_tsts,
                                        lags_es_cc = lags_es_cc,
                                        max_lag_es_cc = max_lag_es_cc,
                                        n_boot_de = n_boot_de,
                                        n_boot_esr = n_boot_esr,
                                        empirical = empirical)
      
      # Number of NAs in UC_par_corr and CC_par_corr printed to NA_information.txt
      if(estimation && par_corr){
        write(paste0('Number of NAs in UC_par_corr: ', sum(is.na(result_lst[['UC_par_corr']][[white_adjust[1]]])), '\n'),
              file = 'run_logs/NA_information.txt',
              append = TRUE)
        write(paste0('Number of NAs in CC_par_corr: ', sum(is.na(result_lst[['CC_par_corr']][[white_adjust[1]]])), '\n'),
              file = 'run_logs/NA_information.txt',
              append = TRUE)
      }
      
      saveRDS(result_lst,
              file = paste0('results/power_', gsub(pattern = '\\.', replacement = '_', paste0(tolerance_lvl)), '_sim_', sim_i, '_est_', est, ifelse(estimation, '_est', '_fix'), '.rds'))
    }
  }
  
  write(paste0('\n\nPower execution finished at ', Sys.time()),
        file = 'run_logs/progress_information.txt',
        append = TRUE)
  cat(paste0('\n\nPower execution finished at ', Sys.time(), '\n\n'))
}
rm(sim_i, tolerance_lvl, result_lst)