# Initial cleaning
#rm(list = ls())
# if (dev.cur() != 1) {
#   dev.off()
# }
cat('\14')

# Load packages and functions
library(ggplot2)
library(dplyr)
library(patchwork)
library(tikzDevice)
source('output_creation/functions_output_creation.R')


#############################################
########### Global parameters ###############
#############################################
folder_path = 'results_paper/1000'
size_correct_power = TRUE
estimation = TRUE
na_rm = TRUE

####### plots only ##########################
create_plots = TRUE
plots_to_latex = FALSE
tolerance_lvl = 0.05
sign_lvls = c((1:100)/1000)
linewidth = 0.7#1
text_size_basis = 8#6

####### tables only #########################
create_tables = TRUE
tolerance_lvls = c(0.05, 0.025, 0.01)
sign_lvl = 0.05
tables_to_latex = TRUE



##### Plot creation ##########################
if(create_plots){
  size_result_plots <- result_plots_creation(folder_path = folder_path,
                                             tolerance_lvl = tolerance_lvl,
                                             sign_lvls = sign_lvls,
                                             spec_idx = c(1,2,5,6,7,8),#1:8
                                             size_or_power = 'size',
                                             size_correct_power = size_correct_power,
                                             ylim = c(0,0.19),
                                             x_lab = 'nominal',
                                             y_lab = 'empirical',
                                             title = c('(P1): GARCH-Normal', '(P2): GARCH-Student\'s t', '(P5): GJR-GARCH-Normal', '(P6): GJR-GARCH-Student\'s t', '(P7): E-GARCH-Normal', '(P8): E-GARCH-Student\'s t'),
                                             print_legend = TRUE,
                                             estimation = estimation,
                                             na_rm = na_rm,
                                             linewidth = linewidth,
                                             text_size_basis = text_size_basis)
  
  size_plots_combined <- 
    (size_result_plots$P1 | size_result_plots$P2) /
    #(size_result_plots$P3 | size_result_plots$P4) /
    (size_result_plots$P5 | size_result_plots$P6) /
    (size_result_plots$P7 | size_result_plots$P8) +
    plot_layout(guides = 'collect')
  
  ##TEST
  # size_plots_combined <- wrap_plots(
  #   size_result_plots$P1, guide_area(),
  #   size_result_plots$P5, size_result_plots$P6,
  #   ncol = 2,
  #   guides = 'collect'
  # ) &
  #   theme(legend.position = 'right')
  
  plot(size_plots_combined)
  
  power_result_plots <- result_plots_creation(folder_path = folder_path,
                                              tolerance_lvl = tolerance_lvl,
                                              sign_lvls = sign_lvls,
                                              spec_idx = 2:8,#c(2,5,6),#,#
                                              size_or_power = 'power',
                                              size_correct_power = size_correct_power,
                                              ylim = c(0,1),
                                              x_lab = 'Target size',
                                              y_lab = 'Empirical size-adjusted power',
                                              title = c('(P2): GARCH-Student\'s t', '(P3): GARCH-skewed Normal', '(P4): AR(1)-GARCH-Normal', '(P5): GJR-GARCH-Normal', '(P6): GJR-GARCH-Student\'s t', '(P7): E-GARCH-Normal', '(P8): E-GARCH-Student\'s t'),
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
  
  ##TEST
  # power_plots_combined <- wrap_plots(
  #   power_result_plots$P2, guide_area(),
  #   power_result_plots$P5, power_result_plots$P6,
  #   ncol = 2,
  #   guides = 'collect'
  # ) &
  #   theme(legend.position = 'right')
  
  plot(power_plots_combined)
  
  if(plots_to_latex){
    # Set plot options
    textwidth = 15.5/2.54  # LaTeX textwidth from "cm" into "inch"
    options(tikzMetricsDictionary = 'tikzMetricsDictionary') # create MetricsDictionary in file
    
    # Export plots
    tikz(file=paste0('../final_paper/plots/size_plots_combined_', gsub(pattern = '\\.', replacement = '_', paste0(tolerance_lvl)), '.tex'), width=textwidth, height=1.2*textwidth)
    plot(size_plots_combined)
    dev.off()
    
    # tikz(file=paste0('../final_paper/plots/power_plots_combined_', gsub(pattern = '\\.', replacement = '_', paste0(tolerance_lvl)), '.tex'), width=textwidth, height=1.2*textwidth)
    # plot(power_plots_combined)
    # dev.off()
  }
}


#################### Result tables creation #####################################

if(create_tables){
  # Size results
  results_size_est <- create_result_matrix(folder_path = folder_path,
                                           tolerance_lvls = tolerance_lvls,
                                           spec_idx = 1:8,#c(1,2,5,6,7,8),
                                           est_idx = 1,
                                           tests = c('UC', ##TEST
                                                     'CC_boot', ##TEST
                                                     'Auxiliary_ESR',
                                                     'Strict_ESR'),
                                           estimation = estimation,
                                           alpha = sign_lvl,
                                           white_adjust = 'hc3',
                                           size_or_power = 'size',
                                           size_correct_power = size_correct_power,
                                           na_rm = na_rm)
  
  if(tables_to_latex){
    latex_tbl_size_est <- generate_latex_size_table(results = results_size_est$results, 
                                                    tolerance_lvls = c(0.05, 0.025, 0.01),
                                                    garch_specs = c('GARCH', 'GJR-GARCH', 'E-GARCH'),
                                                    dist_specs = c('Normal', 'Student\'s $t$'))
  }
  
  # Power results
  results_power_est <- create_result_matrix(folder_path = folder_path,
                                            tolerance_lvls = tolerance_lvls,
                                            spec_idx = 2:8,
                                            est_idx = 1,
                                            tests = c('UC',
                                                      'CC_boot',
                                                      'Auxiliary_ESR',
                                                      'Strict_ESR'),
                                            estimation = estimation,
                                            alpha = sign_lvl,
                                            white_adjust = 'hc3',
                                            size_or_power = 'power',
                                            size_correct_power = size_correct_power,
                                            na_rm = na_rm)
  
  if(tables_to_latex){
    latex_tbl_power_est <- generate_latex_power_table(results = results_power_est$results, 
                                                      tolerance_lvls = c(0.05, 0.025, 0.01),
                                                      garch_specs = c('GARCH', 'AR(1)-GARCH', 'GJR-GARCH', 'E-GARCH'),
                                                      dist_specs = c('Normal', 'Student\'s $t$', 'Skewed Normal'))
  }
}
