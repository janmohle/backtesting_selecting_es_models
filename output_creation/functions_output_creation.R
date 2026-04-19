####################################################################
################### Generate matrix with results ###################
####################################################################
create_result_matrix <- function(folder_path,
                                 tolerance_lvls,
                                 spec_idx,
                                 est_idx,
                                 tests,
                                 estimation,
                                 alpha,
                                 white_adjust='hc3',
                                 size_or_power,
                                 size_correct_power,
                                 na_rm=TRUE){
  
  # Initiate result matrix
  results <- matrix(nrow = length(tolerance_lvls) * length(spec_idx), ncol = length(tests))
  
  # Initiate matrix with result names and used alpha levels
  results_names <- matrix(nrow = length(tolerance_lvls) * length(spec_idx), ncol = length(tests))
  alpha_lst <- matrix(nrow = length(tolerance_lvls) * length(spec_idx), ncol = length(tests))
  
  ##TEST
  #size_corr_names <- matrix(nrow = length(tolerance_lvls) * length(spec_idx), ncol = length(tests))
  
  row_idx <- 0
  
  for(tolerance_lvl in tolerance_lvls){
    
    for(i in spec_idx){
      
      # Read in size result list for size correction
      if(size_or_power == 'power' & size_correct_power){
        size_lst_name <- paste0(folder_path, '/size_', gsub(pattern = '\\.', replacement = '_', paste0(tolerance_lvl)), '_sim_', est_idx, '_est_', est_idx, ifelse(estimation, '_est', '_fix'), '.rds')
        size_lst <- readRDS(size_lst_name)
      }
      
      # Read in result list
      lst_name <- paste0(folder_path, '/', size_or_power, '_', gsub(pattern = '\\.', replacement = '_', paste0(tolerance_lvl)), '_sim_', i, '_est_', ifelse(size_or_power == 'power', est_idx, i), ifelse(estimation, '_est', '_fix'), '.rds')
      lst <- readRDS(lst_name)
      
      row_idx <- row_idx + 1
      
      for(test_idx in 1:length(tests)){
        test <- tests[test_idx]
        
        # Store rejection proportion in result matrix
        if(size_or_power == 'power' & size_correct_power){
          
          # Modify alpha for size adjusted power
          alpha_mod <- quantile(x = size_lst[[test]][[paste0('p_', white_adjust)]],
                                probs = alpha,
                                type = 8,
                                na.rm = na_rm)
          
          results[row_idx, test_idx] <- mean(lst[[test]][[paste0('p_', white_adjust)]] <= alpha_mod,
                                             na.rm = na_rm)
          
        } else {
          
          results[row_idx, test_idx] <- mean(lst[[test]][[paste0('p_', white_adjust)]] <= alpha,
                                             na.rm = na_rm)
        }
        
        # Store result names and used alpha levels
        results_names[row_idx, test_idx] <- paste0(lst_name, ' | ', tests[test_idx])
        alpha_lst[row_idx, test_idx] <- if((size_or_power == 'power') & size_correct_power) alpha_mod else alpha
        
        ##TEST
        #size_corr_names[row_idx, test_idx] <- paste0(size_lst_name, ' | ', tests[test_idx])
      }
    }
  }
  
  ##TEST
  #size_corr_names_s<<-size_corr_names
  
  # Combine results
  result_lst <- list(results = results,
                     names = results_names,
                     alphas = alpha_lst)
  
  return(result_lst)
}


###########################################################################
################### Generate latex code for size tables ###################
###########################################################################
generate_latex_size_table <- function(results,
                                      tolerance_lvls,
                                      garch_specs,  # c('GARCH', 'GJR-GARCH', 'E-GARCH')
                                      dist_specs){  # c('Normal', 'Student's $t$')
  
  # Checks
  stopifnot(nrow(results) == length(tolerance_lvls) * 6L)
  stopifnot(length(garch_specs) == 3L)
  stopifnot(length(dist_specs) == 2L)
  
  lines <- c()
  
  # Table header
  lines <- c(
    lines,
    '\\begin{tabular}{c c l l',
    '                S[table-format=1.3]',
    '                S[table-format=1.3]',
    '                S[table-format=1.3]',
    '                S[table-format=1.3]}',
    '\\toprule',
    '$\\alpha$ & P & Specification & Distribution & {UC} & {CC} & {A-ESR} & {S-ESR} \\\\',
    '\\midrule'
  )
  
  for (a_idx in seq_along(tolerance_lvls)) {
    alpha <- tolerance_lvls[a_idx]
    
    # starting row index in results for current alpha
    base <- (a_idx - 1L) * 6L
    
    # multirow for alpha
    lines <- c(lines, sprintf('\\multirow{6}{*}{%.3f}', alpha))
    
    # (P1): GARCH-Normal
    vals <- paste(sprintf('%.3f', results[base + 1L, ]), collapse = ' & ')
    lines <- c(
      lines,
      sprintf('& 1 & \\multirow{2}{*}{%s} & %s & %s \\\\',
              garch_specs[1], dist_specs[1], vals)
    )
    
    # (P2): GARCH-Student's t
    vals <- paste(sprintf('%.3f', results[base + 2L, ]), collapse = ' & ')
    lines <- c(
      lines,
      sprintf('& 2 & & %s & %s \\\\',
              dist_specs[2], vals)
    )
    
    # (P5): GJR-GARCH-Normal
    vals <- paste(sprintf('%.3f', results[base + 3L, ]), collapse = ' & ')
    lines <- c(
      lines,
      sprintf('& 5 & \\multirow{2}{*}{%s} & %s & %s \\\\',
              garch_specs[2], dist_specs[1], vals)
    )
    
    # (P6): GJR-GARCH-Student's t
    vals <- paste(sprintf('%.3f', results[base + 4L, ]), collapse = ' & ')
    lines <- c(
      lines,
      sprintf('& 6 & & %s & %s \\\\',
              dist_specs[2], vals)
    )
    
    # (P7): E-GARCH-Normal
    vals <- paste(sprintf('%.3f', results[base + 5L, ]), collapse = ' & ')
    lines <- c(
      lines,
      sprintf('& 7 & \\multirow{2}{*}{%s} & %s & %s \\\\',
              garch_specs[3], dist_specs[1], vals)
    )
    
    # (P8): E-GARCH-Student's t
    vals <- paste(sprintf('%.3f', results[base + 6L, ]), collapse = ' & ')
    lines <- c(
      lines,
      sprintf('& 8 & & %s & %s \\\\',
              dist_specs[2], vals)
    )
    
    # midrule between alpha blocks
    if (a_idx != length(tolerance_lvls)) {
      lines <- c(lines, '', '\\midrule', '')
    }
  }
  
  # bottonrule at the end
  lines <- c(lines, '\\bottomrule', '\\end{tabular}')
  paste(lines, collapse = '\n')
}


############################################################################
################### Generate latex code for power tables ###################
############################################################################
generate_latex_power_table <- function(results,
                                       tolerance_lvls,
                                       garch_specs,  # c('GARCH', 'AR(1)-GARCH', 'GJR-GARCH', 'E-GARCH')
                                       dist_specs){ # c('Normal', 'Student's $t$', 'Skewed Normal')
  
  # Checks
  stopifnot(nrow(results) == length(tolerance_lvls) * 7L)
  stopifnot(length(garch_specs) == 4L)
  stopifnot(length(dist_specs) == 3L)
  
  lines <- c()
  
  # Table header
  lines <- c(
    lines,
    '\\begin{tabular}{c c l l',
    '                S[table-format=1.3]',
    '                S[table-format=1.3]',
    '                S[table-format=1.3]',
    '                S[table-format=1.3]}',
    '\\toprule',
    '$\\alpha$ & P & Specification & Distribution & {UC} & {CC} & {A-ESR} & {S-ESR} \\\\',
    '\\midrule'
  )
  
  for (a_idx in seq_along(tolerance_lvls)) {
    alpha <- tolerance_lvls[a_idx]
    
    # starting row index in results for current alpha
    base  <- (a_idx - 1L) * 7L
    
    # multirow for alpha
    lines <- c(lines, sprintf('\\multirow{7}{*}{%.3f}', alpha))
    
    # (P2): GARCH-Student's t
    vals <- paste(sprintf('%.3f', results[base + 1L, ]), collapse = ' & ')
    lines <- c(
      lines,
      sprintf('& 2 & \\multirow{2}{*}{%s} & %s & %s \\\\',
              garch_specs[1], dist_specs[2], vals)
    )
    
    # (P3): GARCH-Skewed Normal
    vals <- paste(sprintf('%.3f', results[base + 2L, ]), collapse = ' & ')
    lines <- c(
      lines,
      sprintf('& 3 & & %s & %s \\\\',
              dist_specs[3], vals)
    )
    
    # (P4): AR(1)-GARCH-Normal
    vals <- paste(sprintf('%.3f', results[base + 3L, ]), collapse = ' & ')
    lines <- c(
      lines,
      sprintf('& 4 & %s & %s & %s \\\\',
              garch_specs[2], dist_specs[1], vals)
    )
    
    # (P5): GJR-GARCH-Normal
    vals <- paste(sprintf('%.3f', results[base + 4L, ]), collapse = ' & ')
    lines <- c(
      lines,
      sprintf('& 5 & \\multirow{2}{*}{%s} & %s & %s \\\\',
              garch_specs[3], dist_specs[1], vals)
    )
    
    # (P6) GJR-GARCH-Student's t
    vals <- paste(sprintf('%.3f', results[base + 5L, ]), collapse = ' & ')
    lines <- c(
      lines,
      sprintf('& 6 & & %s & %s \\\\',
              dist_specs[2], vals)
    )
    
    # (P7): E-GARCH-Normal
    vals <- paste(sprintf('%.3f', results[base + 6L, ]), collapse = ' & ')
    lines <- c(
      lines,
      sprintf('& 7 & \\multirow{2}{*}{%s} & %s & %s \\\\',
              garch_specs[4], dist_specs[1], vals)
    )
    
    # (P8): E-GARCH-Student's t
    vals <- paste(sprintf('%.3f', results[base + 7L, ]), collapse = ' & ')
    lines <- c(
      lines,
      sprintf('& 8 & & %s & %s \\\\',
              dist_specs[2], vals)
    )
    
    # midrule between alpha blocks
    if (a_idx != length(tolerance_lvls)) {
      lines <- c(lines, '', '\\midrule', '')
    }
  }
  
  # bottonrule at the end
  lines <- c(lines, '\\bottomrule', '\\end{tabular}')
  paste(lines, collapse = '\n')
}


#############################################################################
################### Plot ROC of size and power ##############################
#############################################################################

plot_emp_size_power <- function(plot_data_tbl, # table with columns sign_lvl, UC, CC, A_ESR, S_ESR
                                ylim=NULL,
                                x_lab=NULL,
                                y_lab=NULL,
                                title=NULL,
                                print_legend=TRUE,
                                linewidth=1,
                                text_size_basis=6){
  
  plot_data_long <- plot_data_tbl %>%
    tidyr::pivot_longer(cols = c(UC, CC, A_ESR, S_ESR),#(parameter)
                        names_to = 'Backtests',
                        values_to = 'Rejection_prop')
  
  points_data <- plot_data_long %>%
    group_by(Backtests) %>%
    slice(round(seq(1, n(), length.out = 10))) %>%
    ungroup()
  
  plot_res <- ggplot(plot_data_long,
                     aes(x = sign_lvl,
                         y = Rejection_prop,
                         color = Backtests,
                         shape = Backtests)) +
    geom_line(linewidth = linewidth) +
    geom_point(data = points_data, size = 1.8) +
    geom_abline(slope = 1, linetype = 'dashed', linewidth = linewidth) +
    scale_color_manual(name = 'Backtests',
                       values = c(UC = 'blue',#(parameter)
                                  CC = 'red',#(parameter)
                                  A_ESR = 'darkorange',#(parameter)
                                  S_ESR = 'purple'),#(parameter)
                       breaks = c('UC', 'CC', 'A_ESR', 'S_ESR'),#(parameter)
                       labels = c('UC', 'CC', 'A-ESR', 'S-ESR')) +#(parameter)
    scale_shape_manual(name = 'Backtests',
                       values = c(UC = 16,#(parameter)
                                  CC = 17,#(parameter)
                                  A_ESR = 15,#(parameter)
                                  S_ESR = 18),#(parameter)
                       breaks = c('UC', 'CC', 'A_ESR', 'S_ESR'),#(parameter)
                       labels = c('UC', 'CC', 'A-ESR', 'S-ESR')) +#(parameter)
    theme_minimal() + 
    theme(axis.text = element_text(size = text_size_basis),
          legend.text = element_text(size = text_size_basis+1),
          legend.title = element_text(size = text_size_basis+2))
  
  if(!is.null(x_lab)){
    plot_res <- plot_res + 
      labs(x = x_lab) +
      theme(axis.title.x = element_text(size = text_size_basis))
  }
  
  if(!is.null(y_lab)){
    plot_res <- plot_res + 
      labs(y = y_lab) +
      theme(axis.title.y = element_text(size = text_size_basis))
  }
  
  if(!is.null(title)){
    plot_res <- plot_res + 
      labs(title = title) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size = text_size_basis+2))
  }
  
  if(!is.null(ylim)){
    plot_res <- plot_res + 
      scale_y_continuous(limits = ylim)
  }
  
  if(!print_legend){
    plot_res <- plot_res +
      theme(legend.position = 'none')
  }
  
  return(plot_res)
}


#############################################################################
################### Result plot creation ####################################
#############################################################################

result_plots_creation <- function(folder_path,
                                  tolerance_lvl,
                                  sign_lvls,
                                  spec_idx,
                                  size_or_power,
                                  size_correct_power,
                                  ylim=NULL,
                                  x_lab=NULL,
                                  y_lab=NULL,
                                  titles_vec=NULL,
                                  print_legend=TRUE,
                                  estimation=TRUE,
                                  na_rm=TRUE,
                                  linewidth=1,
                                  text_size_basis=6){
  
  plot_data_tbl_lst <- list()
  for(i in seq_along(spec_idx)){
    plot_data_tbl_lst[[paste0('P', spec_idx[i])]] <- tibble(sign_lvl = sign_lvls,
                                                            UC = double(length(sign_lvls)),#(parameter)
                                                            CC = double(length(sign_lvls)),#(parameter)
                                                            A_ESR = double(length(sign_lvls)),#(parameter)
                                                            S_ESR = double(length(sign_lvls)))#(parameter)
  }
  
  for(i_sign_lvl in seq_along(sign_lvls)){
    
    sign_lvl <- sign_lvls[i_sign_lvl]
    
    print(sign_lvl)
    
    results <- create_result_matrix(folder_path = folder_path,
                                    tolerance_lvls = tolerance_lvl,
                                    spec_idx = spec_idx,
                                    est_idx = 1,
                                    tests = c('UC',#(parameter) ##TEST
                                              'CC_boot',#(parameter) ##TEST
                                              'Auxiliary_ESR',#(parameter) ##TEST
                                              'Strict_ESR'),#(parameter) ##TEST
                                    estimation = estimation,
                                    alpha = sign_lvl,
                                    white_adjust = 'hc3',
                                    size_or_power = size_or_power,
                                    size_correct_power = size_correct_power,
                                    na_rm = na_rm)$results
    
    colnames(results) <- c('UC', 'CC', 'A_ESR', 'S_ESR')#(parameter)
    
    for(i_spec in seq_along(spec_idx)){
      for(test in c('UC', 'CC', 'A_ESR', 'S_ESR')){#(parameter)
        plot_data_tbl_lst[[paste0('P', spec_idx[i_spec])]][[test]][i_sign_lvl] <- results[i_spec,test]
      }
    }
  }
  
  ##TEST
  #plot_data_tbl_lst_save <<- plot_data_tbl_lst
  
  plots_return <- list()
  for(i_spec in seq_along(spec_idx)){
    i <- spec_idx[i_spec]
    plots_return[[paste0('P', i)]] <- plot_emp_size_power(plot_data_tbl = plot_data_tbl_lst[[paste0('P', i)]],
                                                          ylim=ylim,
                                                          x_lab=x_lab,
                                                          y_lab=y_lab,
                                                          title= if(!is.null(titles_vec)) titles_vec[i_spec] else NULL,
                                                          print_legend=print_legend,
                                                          linewidth = linewidth,
                                                          text_size_basis = text_size_basis)
  }
  return(plots_return)
}