# TO-DO:
# still to decide which backtests to use
# decide which specifications to use
# decide process parameter
# github update

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
  'MCS')

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
source('functions_decision_process.R')


#############################################
########### Global parameters ###############
#############################################

fresh_start = TRUE
seed = 78+22 # NULL if no seed #79
cores = 1#00
n_loop = 1#00    #1000
oos_window = 500  #1750
est_window = 750
tolerance_lvl = 0.025
sign_lvl_backtests = 0.05
estimate = TRUE
spec_sim_all = 1:8
spec_est_all = 1#:8
backtests = c('CC')#c('A_ESR', , 'CC') #DECIDE # maybe UC, CC, A_ESR?
n_boot_de = 2000 # NULL if no bootstrap for UC and CC backtest
lags_es_cc = NULL#5  # if NULL then automatic lag selection
max_lag_es_cc = NULL # if NULL then max_lag=floor(sqrt(n))
MCS_alpha = 0.3#DECIDE
MCS_statistic = 'Tmax' # Tmax is primary statistic in MCS paper
n_boot_MCS = 50000

#########################################################################
################# Process specifications ################################
#########################################################################

var_spec_1 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_1 = list(armaOrder = c(0,0))
dist_spec_1 = 'norm'
fixed_pars_1 = list(mu = 0,
                    omega = 0.00001,
                    alpha1 = 0.1,
                    beta1 = 0.85)

var_spec_2 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_2 = list(armaOrder = c(0,0))
dist_spec_2 = 'std'
fixed_pars_2 = list(mu = 0,
                    omega = 0.00001,
                    alpha1 = 0.1,
                    beta1 = 0.85,
                    shape = 5)

var_spec_3 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_3 = list(armaOrder = c(0,0))
dist_spec_3 = 'snorm'
fixed_pars_3 = list(mu = 0,
                    omega = 0.00001,
                    alpha1 = 0.1,
                    beta1 = 0.85,
                    skew = 0.85)

var_spec_4 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_4 = list(armaOrder = c(1,0))
dist_spec_4 = 'norm'
fixed_pars_4 = list(mu = 0,
                    omega = 0.00001,
                    alpha1 = 0.1,
                    beta1 = 0.85,
                    ar1 = 0.1)#CHANGED FROM 0.05!!!

var_spec_5 = list(model = 'gjrGARCH',garchOrder = c(1,1))
mean_spec_5 = list(armaOrder = c(0,0))
dist_spec_5 = 'norm'
fixed_pars_5 = list(mu = 0,
                    omega = 0.00001,
                    alpha1 = 0.04,
                    beta1 = 0.85,
                    gamma1 = 0.12)

var_spec_6 = list(model = 'gjrGARCH',garchOrder = c(1,1))
mean_spec_6 = list(armaOrder = c(0,0))
dist_spec_6 = 'std'
fixed_pars_6 = list(mu = 0,
                    omega = 0.00001,
                    alpha1 = 0.04,
                    beta1 = 0.85,
                    gamma1 = 0.12,
                    shape = 5)

var_spec_7 = list(model = 'eGARCH',garchOrder = c(1,1))
mean_spec_7 = list(armaOrder = c(0,0))
dist_spec_7 = 'norm'
fixed_pars_7 = list(mu = 0,
                    omega = -0.447,
                    alpha1 = -0.077,
                    beta1 = 0.948,
                    gamma1 = 0.191)

var_spec_8 = list(model = 'eGARCH',garchOrder = c(1,1))
mean_spec_8 = list(armaOrder = c(0,0))
dist_spec_8 = 'std'
fixed_pars_8 = list(mu = 0,
                    omega = -0.447,
                    alpha1 = -0.077,
                    beta1 = 0.948,
                    gamma1 = 0.191,
                    shape = 5)

# Combine specifications into lists
spec_sim_lst <- list()
for(i in spec_sim_all){
  spec_sim_lst[[paste0('spec_', i)]] <- list(var_spec = get(paste0('var_spec_', i)),
                                             mean_spec = get(paste0('mean_spec_', i)),
                                             dist_spec = get(paste0('dist_spec_', i)),
                                             fixed_pars = get(paste0('fixed_pars_', i)))
}
rm(i)

spec_est_lst <- list()
for(i in spec_est_all){
  spec_est_lst[[paste0('spec_', i)]] <- list(var_spec = get(paste0('var_spec_', i)),
                                             mean_spec = get(paste0('mean_spec_', i)),
                                             dist_spec = get(paste0('dist_spec_', i)),
                                             fixed_pars = get(paste0('fixed_pars_', i)))
}
rm(i)


#############################################################
################# Process execution #########################
#############################################################

# Create required folders
for(dir_name in c('run_logs', 'results', 'checkpoints')){
  dir.create(dir_name, showWarnings = FALSE)
}
rm(dir_name)

# Empty checkpoints folder
if(fresh_start){
  unlink(list.files('checkpoints', full.names = TRUE))
}

# Initialize logs
write(paste0('Execution started at ', Sys.time(), '\n'),
      file = 'run_logs/progress_information.txt',
      append = FALSE)

write('Solver problems (Check estimation code for location of problems encoded in numbers):\n',
      file = 'run_logs/solver_problems.txt',
      append = FALSE)

write('Number of NAs:\n',
      file = 'run_logs/NA_information.txt',
      append = FALSE)

# Initialize list for results
results_collapsed <- list()

# Register cores
registerDoParallel(cores = cores)

# Loop over simulation specifications
for(spec_sim_name in names(spec_sim_lst)){
  
  # Current simulation specification
  cat(paste0('\nCurrent simulation: ', spec_sim_name, ' ', Sys.time(), '\n\n'))
  
  write(paste0('\nCurrent simulation: ', spec_sim_name, ' ', Sys.time(), '\n\n'),
        file = 'run_logs/progress_information.txt',
        append = TRUE)
  
  write(paste0('\nCurrent simulation: ', spec_sim_name, '\n'),
        file = 'run_logs/solver_problems.txt',
        append = TRUE)
  
  write(paste0('\nCurrent simulation: ', spec_sim_name, '\n'),
        file = 'run_logs/NA_information.txt',
        append = TRUE)
  
  # Clear started_loops.txt and finished_loops.txt
  write(x = ' ',
        file = 'run_logs/started_loops.txt',
        append = FALSE)
  write(x = ' ',
        file = 'run_logs/finished_loops.txt',
        append = FALSE)
  
  # Extract current specification
  spec_sim_i <- spec_sim_lst[[spec_sim_name]]
  
  # Parallel estimation loop
  result_foreach <- foreach(i=1:n_loop, .errorhandling = 'pass') %dopar% {
    
    # Log started loops
    write(x = paste0(i),
          file = 'run_logs/started_loops.txt',
          append = TRUE)
    
    if(i %% cores == 0){
      write(x = paste0('Loop ', i, ' started at ', Sys.time(), '\n'),
            file = 'run_logs/progress_information.txt',
            append = TRUE)
    }
    
    # Set seed for reproducibility of results
    if(!is.null(seed)){
      current_seed <- (seed+i^2-i) %% 2147483647
      set.seed(current_seed)
    } else {
      current_seed <- NULL
    }
    
    ######## Simulation ########
    garchspec_sim <- ugarchspec(variance.model = spec_sim_i[['var_spec']],
                                mean.model = spec_sim_i[['mean_spec']],
                                distribution.model = spec_sim_i[['dist_spec']],
                                fixed.pars = spec_sim_i[['fixed_pars']])
    
    garchsimulation <- ugarchpath(spec = garchspec_sim,
                                  n.sim = est_window + oos_window,
                                  n.start = 1000,
                                  m.sim = 1)
    
    # Create tibble of simulation
    sim_matrix <- garchsimulation@path$seriesSim
    colnames(sim_matrix) <- 'Return'
    garchsimulation_tbl <- as_tibble(sim_matrix) %>%
      mutate(Date = base::as.Date('2000-01-01') + 1:nrow(sim_matrix))
    
    ######## Forecasting, backtesting and loss function calculation ########
    loss_tbl_noNA <- evaluate_var_es_backtests_and_loss(return_data_tbl = garchsimulation_tbl,
                                                        spec_est_lst = spec_est_lst,
                                                        tolerance_lvl = tolerance_lvl,
                                                        sign_lvl_backtests = sign_lvl_backtests,
                                                        est_window = est_window,
                                                        estimate = estimate,
                                                        lags_es_cc = lags_es_cc,
                                                        max_lag_es_cc = max_lag_es_cc,
                                                        n_boot_de = n_boot_de,
                                                        current_seed = current_seed,
                                                        current_loop = i,
                                                        run_log = TRUE,
                                                        backtests = backtests)
    
    # Test which specifications passed backtesting
    passed_backtesting <- setNames(
      as.integer(names(spec_est_lst) %in% names(loss_tbl_noNA)),
      names(spec_est_lst)
    )
    
    # Test which specification has rank 1
    rank1 <- setNames(integer(length(spec_est_lst)), names(spec_est_lst))
    
    mean_loss <- colMeans(loss_tbl_noNA)
    ranks <- rank(mean_loss,
                  ties.method = 'min')
    
    rank1[names(ranks)] <- as.integer(ranks == 1) # if no model passes backtesting all models get rank 0 assigned
    
    # Test which specifications are in MCS
    # Initiate binary vector with 0s which signals MCS inclusion with 1s
    MCS_included <- setNames(integer(length(spec_est_lst)), names(spec_est_lst))
    
    if(ncol(loss_tbl_noNA) >= 2){
      
      # MCS procedure
      if(!is.null(seed)){set.seed(current_seed+1)}
      
      MCS_output <- MCS::MCSprocedure(Loss = loss_tbl_noNA,
                                      alpha = MCS_alpha,
                                      B = n_boot_MCS,
                                      statistic = MCS_statistic,
                                      verbose = FALSE)
      
      # Extract MCS
      MCS_vec <- MCS_output@Info$model.names
      
      # Assign 1 to specifications which are in the MCS
      MCS_included[MCS_vec] <- 1
      
    } else if (ncol(loss_tbl_noNA) == 1){
      
      # Assign 1 to single model in loss_tbl_noNA
      MCS_included[colnames(loss_tbl_noNA)] <- 1
    }
    # if loss_tbl_noNA empty then MCS_included contains only 0s
    
    # Log finished loop
    write(x = paste0(i),
          file = 'run_logs/finished_loops.txt',
          append = TRUE)
    
    if(i %% cores == 0){
      write(x = paste0('Loop ', i, ' finished at ', Sys.time(), '\n'),
            file = 'run_logs/progress_information.txt',
            append = TRUE)
    }
    
    # Result list
    result_lst <- list(passed_backtesting = passed_backtesting,
                       rank1 = rank1,
                       MCS_included = MCS_included)
    
    # Checkpoint save
    saveRDS(result_lst, sprintf('checkpoints/results_sim_%s_iter_%04d.rds', spec_sim_name, i))
    
    # Return list
    result_lst
  }
  
  # Checking for errors in iterations
  is_error <- vapply(result_foreach, inherits, logical(1), what = 'error')
  
  if(sum(is_error) > 0){
    
    # Save result_foreach with error iterations
    saveRDS(result_foreach,
            file = paste0('results/error_result_foreach_', spec_sim_name, '.rds'))
    
    # Write error messages to progress_information.txt and console
    error_loops <- which(is_error)
    error_messages <- lapply(result_foreach[error_loops], conditionMessage)
    
    log_message <- paste0('\nProportion of loops with error: ', sum(is_error), '/', length(result_foreach), ' = ', sum(is_error)/length(result_foreach), ',\n',
                          paste('Loop', error_loops, ':', error_messages, collapse = '\n'), '\n')
    
    cat(log_message)
    write(log_message,
          file = 'run_logs/progress_information.txt',
          append = TRUE)
    
    # Remove iterations with error for further processing
    result_foreach <- result_foreach[!is_error]
  }
  
  # Collapse results from parallel loop and divide by number of successful iterations
  results_collapsed[[spec_sim_name]] <- list(
    passed_backtesting = Reduce(`+`, lapply(result_foreach, `[[`, 'passed_backtesting')) / length(result_foreach),
    rank1              = Reduce(`+`, lapply(result_foreach, `[[`, 'rank1')) / length(result_foreach),
    MCS_included       = Reduce(`+`, lapply(result_foreach, `[[`, 'MCS_included')) / length(result_foreach)
  )
  
  # Checkpoint save
  saveRDS(results_collapsed,
          file = 'checkpoints/decision_process_sim_checkpoint.rds')
  
  # Remove checkpoint files from parallel loop
  files_to_remove <- list.files(
    path = 'checkpoints',
    pattern = paste0('results_sim_', spec_sim_name, '_iter_'),
    full.names = TRUE
  )
  file.remove(files_to_remove)
}
rm(spec_sim_name, spec_sim_i, files_to_remove)

# Create result matrix
passed_backtesting_results <- t(sapply(results_collapsed, `[[`, 'passed_backtesting'))
rank1_results <- t(sapply(results_collapsed, `[[`, 'rank1'))
MCS_included_results <- t(sapply(results_collapsed, `[[`, 'MCS_included'))

# Save results
save(passed_backtesting_results, rank1_results, MCS_included_results,
     file = 'results/decision_process_sim.RData')

# Remove checkpoint file
#file.remove('checkpoints/decision_process_sim_checkpoint.rds')

# Log finishing of execution 
write(paste0('\n\nExecution finished at ', Sys.time()),
      file = 'run_logs/progress_information.txt',
      append = TRUE)

# passed backtesting?
# in MCS?
# rank 1?
# -> answer these with binary vectors

# Finished loops:
# ls checkpoints > checkpoint_files.txt