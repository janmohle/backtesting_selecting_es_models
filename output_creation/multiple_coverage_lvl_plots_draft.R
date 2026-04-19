# Multiple versions of plots with multiple coverage level in one plot create from ChatGPT as a draft

# Initial cleaning
rm(list = ls())
if (dev.cur() != 1) {
  dev.off()
}
cat('\14')

# Load packages and functions
library(ggplot2)
library(dplyr)
library(patchwork)
library(tikzDevice)
source('output_creation/functions_output_creation.R')   # contains create_result_matrix etc.
# source the file where you put the plot functions, or paste them above:


########## FIRST VERSION ############


#############################################################################
################### Plot ROC of size and power (single tol) #################
#############################################################################

plot_emp_size_power <- function(plot_data_tbl, # table with columns sign_lvl, UC, CC, A_ESR, S_ESR
                                ylim = NULL,
                                x_lab = NULL,
                                y_lab = NULL,
                                title = NULL,
                                print_legend = TRUE,
                                linewidth = 1,
                                text_size_basis = 6) {
  
  plot_data_long <- plot_data_tbl %>%
    tidyr::pivot_longer(cols = c(UC, CC, A_ESR, S_ESR),
                        names_to = 'Backtests',
                        values_to = 'Rejection_prop')
  
  points_data <- plot_data_long %>%
    dplyr::group_by(Backtests) %>%
    dplyr::slice(round(seq(1, dplyr::n(), length.out = 10))) %>%
    dplyr::ungroup()
  
  plot_res <- ggplot(plot_data_long,
                     aes(x = sign_lvl,
                         y = Rejection_prop,
                         color = Backtests,
                         shape = Backtests)) +
    geom_line(linewidth = linewidth) +
    geom_point(data = points_data, size = 1.8) +
    geom_abline(slope = 1, linetype = 'dashed', linewidth = linewidth) +
    scale_color_manual(name = 'Backtests',
                       values = c(UC = 'blue',
                                  CC = 'red',
                                  A_ESR = 'darkorange',
                                  S_ESR = 'purple'),
                       breaks = c('UC', 'CC', 'A_ESR', 'S_ESR'),
                       labels = c('UC', 'CC', 'A-ESR', 'S-ESR')) +
    scale_shape_manual(name = 'Backtests',
                       values = c(UC = 16,
                                  CC = 17,
                                  A_ESR = 15,
                                  S_ESR = 18),
                       breaks = c('UC', 'CC', 'A_ESR', 'S_ESR'),
                       labels = c('UC', 'CC', 'A-ESR', 'S-ESR')) +
    theme_minimal() + 
    theme(axis.text    = element_text(size = text_size_basis),
          legend.text  = element_text(size = text_size_basis + 1),
          legend.title = element_text(size = text_size_basis + 2))
  
  if (!is.null(x_lab)) {
    plot_res <- plot_res + 
      labs(x = x_lab) +
      theme(axis.title.x = element_text(size = text_size_basis))
  }
  
  if (!is.null(y_lab)) {
    plot_res <- plot_res + 
      labs(y = y_lab) +
      theme(axis.title.y = element_text(size = text_size_basis))
  }
  
  if (!is.null(title)) {
    plot_res <- plot_res + 
      labs(title = title) +
      theme(plot.title = element_text(hjust = 0.5,
                                      size  = text_size_basis + 2))
  }
  
  if (!is.null(ylim)) {
    plot_res <- plot_res + 
      scale_y_continuous(limits = ylim)
  }
  
  if (!print_legend) {
    plot_res <- plot_res +
      theme(legend.position = 'none')
  }
  
  return(plot_res)
}


#############################################################################
################### Plot ROC with multiple ES tolerances ####################
#############################################################################

plot_emp_size_power_tols <- function(plot_data_tbl,  # columns: sign_lvl, tolerance_lvl, UC, CC, A_ESR, S_ESR
                                     ylim = NULL,
                                     x_lab = NULL,
                                     y_lab = NULL,
                                     title = NULL,
                                     print_legend = TRUE,
                                     linewidth = 0.7,
                                     text_size_basis = 6) {
  
  # Long format + tolerance levels
  plot_data_long <- plot_data_tbl %>%
    tidyr::pivot_longer(cols = c(UC, CC, A_ESR, S_ESR),
                        names_to = 'Backtests',
                        values_to = 'Rejection_prop')
  
  # Subsample points per Backtests x tolerance level
  points_data <- plot_data_long %>%
    dplyr::group_by(Backtests, tolerance_lvl) %>%
    dplyr::slice(round(seq(1, dplyr::n(), length.out = 10))) %>%
    dplyr::ungroup()
  
  # Define linetypes for the ES tolerance levels
  tol_levels <- sort(unique(plot_data_long$tolerance_lvl))
  tol_labels <- sprintf("ES tol = %.3f", tol_levels)
  names(tol_labels) <- tol_levels
  
  plot_res <- ggplot(plot_data_long,
                     aes(x = sign_lvl,
                         y = Rejection_prop,
                         color = Backtests,
                         shape = Backtests,
                         linetype = factor(tolerance_lvl))) +
    geom_line(linewidth = linewidth) +
    geom_point(data = points_data, size = 1.8) +
    geom_abline(slope = 1, linetype = 'dashed', linewidth = linewidth) +
    scale_color_manual(name = 'Backtests',
                       values = c(UC = 'blue',
                                  CC = 'red',
                                  A_ESR = 'darkorange',
                                  S_ESR = 'purple'),
                       breaks = c('UC', 'CC', 'A_ESR', 'S_ESR'),
                       labels = c('UC', 'CC', 'A-ESR', 'S-ESR')) +
    scale_shape_manual(name = 'Backtests',
                       values = c(UC = 16,
                                  CC = 17,
                                  A_ESR = 15,
                                  S_ESR = 18),
                       breaks = c('UC', 'CC', 'A_ESR', 'S_ESR'),
                       labels = c('UC', 'CC', 'A-ESR', 'S-ESR')) +
    scale_linetype_discrete(name   = "ES tolerance",
                            labels = tol_labels) +
    theme_minimal() +
    theme(axis.text    = element_text(size = text_size_basis),
          legend.text  = element_text(size = text_size_basis + 1),
          legend.title = element_text(size = text_size_basis + 2))
  
  if (!is.null(x_lab)) {
    plot_res <- plot_res +
      labs(x = x_lab) +
      theme(axis.title.x = element_text(size = text_size_basis))
  }
  
  if (!is.null(y_lab)) {
    plot_res <- plot_res +
      labs(y = y_lab) +
      theme(axis.title.y = element_text(size = text_size_basis))
  }
  
  if (!is.null(title)) {
    plot_res <- plot_res +
      labs(title = title) +
      theme(plot.title = element_text(hjust = 0.5,
                                      size  = text_size_basis + 2))
  }
  
  if (!is.null(ylim)) {
    plot_res <- plot_res +
      scale_y_continuous(limits = ylim)
  }
  
  if (!print_legend) {
    plot_res <- plot_res +
      theme(legend.position = 'none')
  }
  
  return(plot_res)
}


#############################################################################
################### Result plot creation (single tol) #######################
#############################################################################

result_plots_creation <- function(folder_path,
                                  tolerance_lvl,
                                  sign_lvls,
                                  spec_idx,
                                  size_or_power,
                                  size_correct_power,
                                  ylim = NULL,
                                  x_lab = NULL,
                                  y_lab = NULL,
                                  titles_vec = NULL,
                                  print_legend = TRUE,
                                  estimation = TRUE,
                                  na_rm = TRUE,
                                  linewidth = 1,
                                  text_size_basis = 6){
  
  plot_data_tbl_lst <- list()
  for (i in seq_along(spec_idx)) {
    plot_data_tbl_lst[[paste0('P', spec_idx[i])]] <- tibble(
      sign_lvl = sign_lvls,
      UC       = double(length(sign_lvls)),
      CC       = double(length(sign_lvls)),
      A_ESR    = double(length(sign_lvls)),
      S_ESR    = double(length(sign_lvls))
    )
  }
  
  for (i_sign_lvl in seq_along(sign_lvls)){
    
    sign_lvl <- sign_lvls[i_sign_lvl]
    print(sign_lvl)
    
    results <- create_result_matrix(
      folder_path        = folder_path,
      tolerance_lvls     = tolerance_lvl,
      spec_idx           = spec_idx,
      est_idx            = 1,
      tests              = c('UC_boot', 'CC_boot', 'Auxiliary_ESR', 'Strict_ESR'),
      estimation         = estimation,
      alpha              = sign_lvl,
      white_adjust       = 'hc3',
      size_or_power      = size_or_power,
      size_correct_power = size_correct_power,
      na_rm              = na_rm
    )$results
    
    colnames(results) <- c('UC', 'CC', 'A_ESR', 'S_ESR')
    
    for (i_spec in seq_along(spec_idx)) {
      for (test in c('UC', 'CC', 'A_ESR', 'S_ESR')) {
        plot_data_tbl_lst[[paste0('P', spec_idx[i_spec])]][[test]][i_sign_lvl] <- results[i_spec, test]
      }
    }
  }
  
  plots_return <- list()
  for (i_spec in seq_along(spec_idx)) {
    i <- spec_idx[i_spec]
    plots_return[[paste0('P', i)]] <- plot_emp_size_power(
      plot_data_tbl   = plot_data_tbl_lst[[paste0('P', i)]],
      ylim            = ylim,
      x_lab           = x_lab,
      y_lab           = y_lab,
      title           = if (!is.null(titles_vec)) titles_vec[i_spec] else NULL,
      print_legend    = print_legend,
      linewidth       = linewidth,
      text_size_basis = text_size_basis
    )
  }
  return(plots_return)
}


#############################################################################
################### Result plot creation (multi ES tol) #####################
#############################################################################

result_plots_creation_multi_tol <- function(folder_path,
                                            tolerance_lvls = c(0.05, 0.025, 0.01),
                                            sign_lvls,
                                            spec_idx,
                                            size_or_power,
                                            size_correct_power,
                                            ylim = NULL,
                                            x_lab = NULL,
                                            y_lab = NULL,
                                            titles_vec = NULL,
                                            print_legend = TRUE,
                                            estimation = TRUE,
                                            na_rm = TRUE,
                                            linewidth = 0.7,
                                            text_size_basis = 6) {
  
  # For each spec, create data frame with rows = sign_lvls * tolerance_lvls
  plot_data_tbl_lst <- list()
  for (i in seq_along(spec_idx)) {
    Pname <- paste0("P", spec_idx[i])
    plot_data_tbl_lst[[Pname]] <- tibble(
      sign_lvl      = rep(sign_lvls, times = length(tolerance_lvls)),
      tolerance_lvl = rep(tolerance_lvls, each = length(sign_lvls)),
      UC            = double(length(sign_lvls) * length(tolerance_lvls)),
      CC            = double(length(sign_lvls) * length(tolerance_lvls)),
      A_ESR         = double(length(sign_lvls) * length(tolerance_lvls)),
      S_ESR         = double(length(sign_lvls) * length(tolerance_lvls))
    )
  }
  
  # Fill in rejection proportions for each sign level and tolerance level
  for (tol_idx in seq_along(tolerance_lvls)) {
    tol <- tolerance_lvls[tol_idx]
    
    for (i_sign_lvl in seq_along(sign_lvls)) {
      sign_lvl <- sign_lvls[i_sign_lvl]
      message("Processing tol = ", tol, ", alpha = ", sign_lvl)
      
      results <- create_result_matrix(
        folder_path        = folder_path,
        tolerance_lvls     = tol,
        spec_idx           = spec_idx,
        est_idx            = 1,
        tests              = c('UC_boot', 'CC_boot', 'Auxiliary_ESR', 'Strict_ESR'),
        estimation         = estimation,
        alpha              = sign_lvl,
        white_adjust       = 'hc3',
        size_or_power      = size_or_power,
        size_correct_power = size_correct_power,
        na_rm              = na_rm
      )$results
      
      colnames(results) <- c('UC', 'CC', 'A_ESR', 'S_ESR')
      
      # Row index in the stacked table for this tolerance & significance level
      row_offset <- (tol_idx - 1) * length(sign_lvls) + i_sign_lvl
      
      for (i_spec in seq_along(spec_idx)) {
        Pname <- paste0("P", spec_idx[i_spec])
        for (test in c('UC', 'CC', 'A_ESR', 'S_ESR')) {
          plot_data_tbl_lst[[Pname]][[test]][row_offset] <- results[i_spec, test]
        }
      }
    }
  }
  
  plots_return <- list()
  for (i_spec in seq_along(spec_idx)) {
    i <- spec_idx[i_spec]
    Pname <- paste0("P", i)
    plots_return[[Pname]] <- plot_emp_size_power_tols(
      plot_data_tbl   = plot_data_tbl_lst[[Pname]],
      ylim            = ylim,
      x_lab           = x_lab,
      y_lab           = y_lab,
      title           = if (!is.null(titles_vec)) titles_vec[i_spec] else NULL,
      print_legend    = print_legend,
      linewidth       = linewidth,
      text_size_basis = text_size_basis
    )
  }
  
  return(plots_return)
}




# source('output_creation/functions_plots.R')

#############################################
########### Global parameters ###############
#############################################
folder_path        <- 'results_final/oos_5000'
size_correct_power <- TRUE
estimation         <- TRUE
na_rm              <- TRUE

####### plots only ##########################
create_plots      <- TRUE
plots_to_latex    <- FALSE
tolerance_lvls_pl <- c(0.05, 0.025)   # three ES tolerance levels in ONE plot
sign_lvls         <- c((1:100)/1000)
linewidth         <- 0.7
text_size_basis   <- 8

####### tables only #########################
create_tables    <- TRUE
tolerance_lvls   <- c(0.05, 0.025, 0.01)
sign_lvl         <- 0.05
tables_to_latex  <- TRUE


##### Plot creation with multiple ES tolerance levels ##########################
if (create_plots) {
  
  ## SIZE PLOTS: one plot per specification, linetype = ES tolerance level
  size_result_plots_multi <- result_plots_creation_multi_tol(
    folder_path        = folder_path,
    tolerance_lvls     = tolerance_lvls_pl,                  # 0.05, 0.025, 0.01
    sign_lvls          = sign_lvls,
    spec_idx           = c(1, 2, 5, 6, 7, 8),
    size_or_power      = 'size',
    size_correct_power = size_correct_power,
    ylim               = c(0, 0.19),
    x_lab              = 'Nominal size',
    y_lab              = 'Empirical size',
    titles_vec         = c('(P1): GARCH-Normal',
                           '(P2): GARCH-Student\'s t',
                           '(P5): GJR-GARCH-Normal',
                           '(P6): GJR-GARCH-Student\'s t',
                           '(P7): E-GARCH-Normal',
                           '(P8): E-GARCH-Student\'s t'),
    print_legend       = TRUE,
    estimation         = estimation,
    na_rm              = na_rm,
    linewidth          = linewidth,
    text_size_basis    = text_size_basis
  )
  
  size_plots_combined_multi <- (size_result_plots_multi$P1 | size_result_plots_multi$P2) /
    (size_result_plots_multi$P5 | size_result_plots_multi$P6) /
    (size_result_plots_multi$P7 | size_result_plots_multi$P8) +
    plot_layout(guides = 'collect')
  
  plot(size_plots_combined_multi)
  
  
  ## POWER PLOTS: one plot per specification, linetype = ES tolerance level
  power_result_plots_multi <- result_plots_creation_multi_tol(
    folder_path        = folder_path,
    tolerance_lvls     = tolerance_lvls_pl,                  # 0.05, 0.025, 0.01
    sign_lvls          = sign_lvls,
    spec_idx           = 2:8,
    size_or_power      = 'power',
    size_correct_power = size_correct_power,
    ylim               = c(0, 1),
    x_lab              = 'Target size',
    y_lab              = 'Empirical size-adjusted power',
    titles_vec         = c('(P2): GARCH-Student\'s t',
                           '(P3): GARCH-skewed Normal',
                           '(P4): AR(1)-GARCH-Normal',
                           '(P5): GJR-GARCH-Normal',
                           '(P6): GJR-GARCH-Student\'s t',
                           '(P7): E-GARCH-Normal',
                           '(P8): E-GARCH-Student\'s t'),
    print_legend       = TRUE,
    estimation         = estimation,
    na_rm              = na_rm,
    linewidth          = linewidth,
    text_size_basis    = text_size_basis
  )
  
  power_plots_combined_multi <- wrap_plots(
    power_result_plots_multi$P2, power_result_plots_multi$P3,
    power_result_plots_multi$P4, guide_area(),
    power_result_plots_multi$P5, power_result_plots_multi$P6,
    power_result_plots_multi$P7, power_result_plots_multi$P8,
    ncol  = 2,
    guides = 'collect'
  ) &
    theme(legend.position = 'right')
  
  plot(power_plots_combined_multi)
  
  
  ####### Export plots to LaTeX (optional) ###################################
  if (plots_to_latex) {
    # Set plot options
    textwidth <- 15.5/2.54  # LaTeX textwidth from "cm" into "inch"
    options(tikzMetricsDictionary = 'tikzMetricsDictionary') # create MetricsDictionary in file
    
    # Note: use some label for tolerance set in filename, e.g. 'multi'
    tikz(file = '../final_paper/plots/size_plots_combined_multi_tol.tex',
         width = textwidth, height = 1.2 * textwidth)
    plot(size_plots_combined_multi)
    dev.off()
    
    tikz(file = '../final_paper/plots/power_plots_combined_multi_tol.tex',
         width = textwidth, height = 1.2 * textwidth)
    plot(power_plots_combined_multi)
    dev.off()
  }
}







############## SECOND VERSION ##############
#############################################
# Multi-tolerance plot with:
# - ES tol = 0.050 -> solid
# - ES tol = 0.025 -> dense dotted
#############################################

plot_emp_size_power_tols <- function(plot_data_tbl,  # sign_lvl, tolerance_lvl, UC, CC, A_ESR, S_ESR
                                     ylim = NULL,
                                     x_lab = NULL,
                                     y_lab = NULL,
                                     title = NULL,
                                     print_legend = TRUE,
                                     linewidth = 0.9,
                                     text_size_basis = 8) {
  
  # Original colors
  base_cols <- c(
    UC    = "blue",
    CC    = "red",
    A_ESR = "darkorange",
    S_ESR = "purple"
  )
  
  # Original shapes per backtest
  base_shapes <- c(
    UC    = 16,
    CC    = 17,
    A_ESR = 15,
    S_ESR = 18
  )
  
  # Long format
  plot_data_long <- plot_data_tbl %>%
    tidyr::pivot_longer(cols = c(UC, CC, A_ESR, S_ESR),
                        names_to = "Backtests",
                        values_to = "Rejection_prop") %>%
    dplyr::mutate(
      tol_str = sprintf("%.3f", as.numeric(tolerance_lvl))  # "0.025" or "0.050"
    )
  
  # Only two tolerance levels
  tol_levels <- c("0.025", "0.050")
  # 0.050 solid, 0.025 dense dotted
  tol_lty <- c(
    "0.025" = "11",    # short dash/short gap → dense dotted look
    "0.050" = "solid"
  )
  
  # Map any unexpected tolerance to 0.050
  plot_data_long <- plot_data_long %>%
    dplyr::mutate(
      tol_str = dplyr::if_else(tol_str %in% tol_levels, tol_str, "0.050")
    )
  
  # Subsample points for clarity (same symbols as your original)
  points_data <- plot_data_long %>%
    dplyr::group_by(Backtests, tol_str) %>%
    dplyr::slice(round(seq(1, dplyr::n(), length.out = 8))) %>%
    dplyr::ungroup()
  
  tol_labels <- c(
    "0.025" = "ES tol = 0.025 (dotted)",
    "0.050" = "ES tol = 0.050 (solid)"
  )
  
  p <- ggplot(plot_data_long,
              aes(x = sign_lvl,
                  y = Rejection_prop,
                  color = Backtests,
                  linetype = tol_str,
                  group = interaction(Backtests, tol_str))) +
    geom_line(linewidth = linewidth) +
    geom_point(data = points_data,
               aes(shape = Backtests),
               size = 1.4,
               stroke = 0.2) +
    geom_abline(slope = 1,
                linetype = "dotted",
                linewidth = 0.7,
                color = "grey40") +
    scale_color_manual(
      name   = "Backtests",
      values = base_cols,
      breaks = c("UC", "CC", "A_ESR", "S_ESR"),
      labels = c("UC", "CC", "A-ESR", "S-ESR")
    ) +
    scale_shape_manual(
      name   = "Backtests",
      values = base_shapes,
      breaks = c("UC", "CC", "A_ESR", "S_ESR"),
      labels = c("UC", "CC", "A-ESR", "S-ESR")
    ) +
    scale_linetype_manual(
      name   = "ES tolerance",
      values = tol_lty,
      breaks = tol_levels,
      labels = tol_labels
    ) +
    theme_minimal() +
    theme(
      axis.text        = element_text(size = text_size_basis),
      axis.title       = element_text(size = text_size_basis),
      legend.text      = element_text(size = text_size_basis),
      legend.title     = element_text(size = text_size_basis + 1),
      legend.box       = "vertical",
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(x_lab)) {
    p <- p + labs(x = x_lab)
  }
  if (!is.null(y_lab)) {
    p <- p + labs(y = y_lab)
  }
  if (!is.null(title)) {
    p <- p + labs(title = title) +
      theme(plot.title = element_text(hjust = 0.5,
                                      size  = text_size_basis + 1))
  }
  if (!is.null(ylim)) {
    p <- p + scale_y_continuous(limits = ylim)
  }
  if (!print_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  p
}








############## OTHER #################
tolerance_lvls_plot <- c(0.05, 0.025)

all_size_figures  <- list()
all_power_figures <- list()

for (tol in tolerance_lvls_plot) {
  message("Creating plots for ES tolerance = ", tol)
  
  # SIZE plots (your original function)
  size_result_plots <- result_plots_creation(
    folder_path        = folder_path,
    tolerance_lvl      = tol,
    sign_lvls          = sign_lvls,
    spec_idx           = c(1,2,5,6,7,8),
    size_or_power      = "size",
    size_correct_power = size_correct_power,
    ylim               = c(0, 0.19),
    x_lab              = "Nominal size",
    y_lab              = "Empirical size",
    titles_vec         = c("(P1): GARCH-Normal",
                           "(P2): GARCH-Student's t",
                           "(P5): GJR-GARCH-Normal",
                           "(P6): GJR-GARCH-Student's t",
                           "(P7): E-GARCH-Normal",
                           "(P8): E-GARCH-Student's t"),
    print_legend       = TRUE,
    estimation         = estimation,
    na_rm              = na_rm,
    linewidth          = linewidth,
    text_size_basis    = text_size_basis
  )
  
  size_fig <- (size_result_plots$P1 | size_result_plots$P2) /
    (size_result_plots$P5 | size_result_plots$P6) /
    (size_result_plots$P7 | size_result_plots$P8) +
    plot_layout(guides = "collect")
  
  all_size_figures[[paste0("tol_", tol)]] <- size_fig
  
  # POWER plots
  power_result_plots <- result_plots_creation(
    folder_path        = folder_path,
    tolerance_lvl      = tol,
    sign_lvls          = sign_lvls,
    spec_idx           = 2:8,
    size_or_power      = "power",
    size_correct_power = size_correct_power,
    ylim               = c(0, 1),
    x_lab              = "Target size",
    y_lab              = "Empirical size-adjusted power",
    titles_vec         = c("(P2): GARCH-Student's t",
                           "(P3): GARCH-skewed Normal",
                           "(P4): AR(1)-GARCH-Normal",
                           "(P5): GJR-GARCH-Normal",
                           "(P6): GJR-GARCH-Student's t",
                           "(P7): E-GARCH-Normal",
                           "(P8): E-GARCH-Student's t"),
    print_legend       = TRUE,
    estimation         = estimation,
    na_rm              = na_rm,
    linewidth          = linewidth,
    text_size_basis    = text_size_basis
  )
  
  power_fig <- wrap_plots(
    power_result_plots$P2, power_result_plots$P3,
    power_result_plots$P4, guide_area(),
    power_result_plots$P5, power_result_plots$P6,
    power_result_plots$P7, power_result_plots$P8,
    ncol  = 2,
    guides = "collect"
  ) &
    theme(legend.position = "right")
  
  all_power_figures[[paste0("tol_", tol)]] <- power_fig
}

# Example: inspect one tolerance at a time
plot(all_size_figures[["tol_0.05"]])
plot(all_power_figures[["tol_0.05"]])

plot(all_size_figures[["tol_0.025"]])
plot(all_power_figures[["tol_0.025"]])
