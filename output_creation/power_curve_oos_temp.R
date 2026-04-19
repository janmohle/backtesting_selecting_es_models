#############################################################################
################### Plot ROC of size and power ##############################
#############################################################################

plot_emp_size_power <- function(plot_data_tbl, # table with columns oos_window, UC, CC, A_ESR, S_ESR
                                sign_lvl,
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
  
  ##TEST
  #print(plot_data_long)
  
  points_data <- plot_data_long %>%
    group_by(Backtests) %>%
    slice(round(seq(1, n(), length.out = nrow(plot_data_tbl)))) %>%##TEST
    ungroup()
  
  plot_res <- ggplot(plot_data_long,
                     aes(x = oos_window,
                         y = Rejection_prop,
                         color = Backtests,
                         shape = Backtests)) +
    geom_line(linewidth = linewidth) +
    geom_point(data = points_data, size = 1.8) +
    geom_abline(intercept = sign_lvl, slope = 0, linetype = 'dashed', linewidth = linewidth) +
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
                                  oos_windows,
                                  tolerance_lvl,
                                  sign_lvl,
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
    plot_data_tbl_lst[[paste0('P', spec_idx[i])]] <- tibble(oos_window = oos_windows,
                                                            UC = double(length(oos_windows)),#(parameter)
                                                            CC = double(length(oos_windows)),#(parameter)
                                                            A_ESR = double(length(oos_windows)),#(parameter)
                                                            S_ESR = double(length(oos_windows)))#(parameter)
  }
  
  for(i_oos_windows in seq_along(oos_windows)){
    
    oos_window <- oos_windows[i_oos_windows]
    
    print(oos_window)
    
    results <- create_result_matrix(folder_path = paste0(folder_path, oos_window),
                                    tolerance_lvls = tolerance_lvl,
                                    spec_idx = spec_idx,
                                    est_idx = 1,
                                    tests = c('UC_boot',#(parameter) ##TEST
                                              'CC_boot',#(parameter) ##TEST
                                              'Auxiliary_ESR',#(parameter)
                                              'Strict_ESR'),#(parameter)
                                    estimation = estimation,
                                    alpha = sign_lvl,
                                    white_adjust = 'hc3',
                                    size_or_power = size_or_power,
                                    size_correct_power = size_correct_power,
                                    na_rm = na_rm)$results
    
    colnames(results) <- c('UC', 'CC', 'A_ESR', 'S_ESR')#(parameter)
    
    for(i_spec in seq_along(spec_idx)){
      for(test in c('UC', 'CC', 'A_ESR', 'S_ESR')){#(parameter)
        plot_data_tbl_lst[[paste0('P', spec_idx[i_spec])]][[test]][i_oos_windows] <- results[i_spec,test]
      }
    }
  }
  
  ##TEST
  #print(plot_data_tbl_lst)
  
  plots_return <- list()
  for(i_spec in seq_along(spec_idx)){
    i <- spec_idx[i_spec]
    plots_return[[paste0('P', i)]] <- plot_emp_size_power(plot_data_tbl = plot_data_tbl_lst[[paste0('P', i)]],
                                                          sign_lvl = sign_lvl,
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



#############################################
########### Global parameters ###############
#############################################
folder_path = 'results_final/oos_'
size_correct_power = TRUE
estimation = TRUE
na_rm = TRUE

####### plots only ##########################
tolerance_lvl = 0.01
create_plots = TRUE
plots_to_latex = FALSE
oos_windows = c(500, 1000, 5000)
sign_lvl = 0.05
linewidth = 1
text_size_basis = 6

size_result_plots <- result_plots_creation(folder_path = folder_path,
                                           oos_windows =oos_windows,
                                           tolerance_lvl = tolerance_lvl,
                                           sign_lvl = sign_lvl,
                                           spec_idx = c(1,2,5,6,7,8),
                                           size_or_power = 'size',
                                           size_correct_power = size_correct_power,
                                           ylim = c(0,0.19),
                                           x_lab = 'Backtesting sample length',
                                           y_lab = 'Empirical size',
                                           titles_vec = c('(P1): GARCH-Normal', '(P2): GARCH-Student\'s t', '(P5): GJR-GARCH-Normal', '(P6): GJR-GARCH-Student\'s t', '(P7): E-GARCH-Normal', '(P8): E-GARCH-Student\'s t'),
                                           print_legend = TRUE,
                                           estimation = estimation,
                                           na_rm = na_rm,
                                           linewidth = linewidth,
                                           text_size_basis = text_size_basis)


size_plots_combined <- (size_result_plots$P1 | size_result_plots$P2) /
  (size_result_plots$P5 | size_result_plots$P6) /
  (size_result_plots$P7 | size_result_plots$P8) +
  plot_layout(guides = 'collect')
plot(size_plots_combined)



power_result_plots <- result_plots_creation(folder_path = folder_path,
                                            oos_windows =oos_windows,
                                            tolerance_lvl = tolerance_lvl,
                                            sign_lvl = sign_lvl,
                                            spec_idx = 2:8,
                                            size_or_power = 'power',
                                            size_correct_power = size_correct_power,
                                            ylim = c(0,1),
                                            x_lab = 'Backtesting sample length',
                                            y_lab = 'Empirical size-adjusted power',
                                            titles_vec = c('(P2): GARCH-Student\'s t', '(P3): GARCH-skewed Normal', '(P4): AR(1)-GARCH-Normal', '(P5): GJR-GARCH-Normal', '(P6): GJR-GARCH-Student\'s t', '(P7): E-GARCH-Normal', '(P8): E-GARCH-Student\'s t'),
                                            print_legend = TRUE,
                                            estimation = estimation,
                                            na_rm = na_rm,
                                            linewidth = linewidth,
                                            text_size_basis = text_size_basis)

power_plots_combined <- wrap_plots(
  power_result_plots$P2, power_result_plots$P3,
  power_result_plots$P4, guide_area(),
  power_result_plots$P5, power_result_plots$P6,
  power_result_plots$P7, power_result_plots$P8,
  ncol = 2,
  guides = 'collect'
) &
  theme(legend.position = 'right')
plot(power_plots_combined)



# Set plot options
textwidth = 15.5/2.54  # LaTeX textwidth from "cm" into "inch"
options(tikzMetricsDictionary = 'tikzMetricsDictionary') # create MetricsDictionary in file

# Export plots
tikz(file=paste0('../final_paper/plots/size_oos_plots_combined_', gsub(pattern = '\\.', replacement = '_', paste0(tolerance_lvl)), '.tex'), width=textwidth, height=1.2*textwidth)
plot(size_plots_combined)
dev.off()

tikz(file=paste0('../final_paper/plots/power_oos_plots_combined_', gsub(pattern = '\\.', replacement = '_', paste0(tolerance_lvl)), '.tex'), width=textwidth, height=1.2*textwidth)
plot(power_plots_combined)
dev.off()