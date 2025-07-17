#' plot
#'
#' @description This function plots a summary of the results of the model.
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @importFrom tidyr unite
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr with_groups
#' @importFrom magrittr equals
#'
#' @param x A tibble including a cell_group name column | sample name column | read counts column | factor columns | Pvalue column | a significance column
#' @param significance_threshold Numeric value specifying the significance threshold for highlighting differences. Default is 0.025.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' @param significance_statistic Character vector indicating which statistic to highlight. Default is "FDR".
#' @param show_fdr_message Logical. Whether to show the Bayesian FDR interpretation message on the plot. Default is TRUE.
#' @param ... For internal use 
#'
#' @return A `ggplot`
#'
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
#'
#' @export
#'
#' @examples
#'
#' print("cmdstanr is needed to run this example.")
#' # Note: Before running the example, ensure that the 'cmdstanr' package is installed:
#' # install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))
#'
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists()) {
#'     data("counts_obj")
#'
#'     estimate = sccomp_estimate(
#'       counts_obj,
#'       ~ type, ~1, "sample", "cell_group", "count",
#'       cores = 1
#'     )
#'
#'     # estimate |> plot()
#'   }
#' }
#'
plot.sccomp_tbl <- function(x,  significance_threshold = 0.05, test_composition_above_logit_fold_change = .data |> attr("test_composition_above_logit_fold_change"), significance_statistic = c("FDR", "pH0"), show_fdr_message = TRUE, ...) {
  
  # Define the variables as NULL to avoid CRAN NOTES
  parameter <- NULL
  count_data <- x |> attr("count_data")
  v_effect <- NULL
  
  .cell_group = attr(x, ".cell_group")
  .count = attr(x, ".count")
  .sample = attr(x, ".sample")
  
  plots = list()
  
  # Check if test have been done
  if(x |> select(ends_with("FDR")) |> ncol() |> equals(0))
    stop("sccomp says: to produce plots, you need to run the function sccomp_test() on your estimates.")
  
  data_proportion =
    x |>
    
    # Otherwise does not work
    select(-`factor`) |>
    
    pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) |>
    left_join(
      count_data, by = join_by(!!.cell_group)
    ) |> 
    with_groups(!!.sample, ~ mutate(.x, proportion = (!!.count)/sum(!!.count)) )
  
  
  # If I don't have outliers add them
  if(!"outlier" %in% colnames(data_proportion)) data_proportion = data_proportion |> mutate(outlier = FALSE) 
  
  # Select the factors to plot  
  my_factors =
    x |>
    filter(!is.na(`factor`)) |>
    distinct(`factor`) |>
    pull(`factor`) 
  
  if(length(my_factors)==0)
    message("sccomp says: the contrasts you have tested do not represent factors. Therefore, plot of the the posterior predictive check will be omitted.")
  else {
    # Boxplot
    plots$boxplot =
      
      my_factors |>
      map(
        ~ {
          
          # If variable is continuous
          if(data_proportion |> select(all_of(.x)) |> pull(1) |> is("numeric"))
            my_plot = 
              plot_scatterplot(
                .data = x,
                data_proportion = data_proportion,
                factor_of_interest = .x,
                .cell_group = !!.cell_group,
                .sample =  !!.sample,
                my_theme = sccomp_theme(),
                significance_threshold = significance_threshold
              )
          
          # If discrete
          else 
            my_plot = 
              sccomp_boxplot(
                .data = x,
                factor = .x,
                significance_threshold = significance_threshold
              ) 
          
          # Return
          my_plot +
            ggtitle(sprintf("Grouped by %s (for multi-factor models, associations could be hardly observable with unidimensional data stratification)", .x))
        } )
    
  }
  
  # 1D intervals
  plots$credible_intervals_1D = plot_1D_intervals(
    .data = x, 
    significance_threshold = significance_threshold,
    significance_statistic = significance_statistic,
    show_fdr_message = show_fdr_message
  )
  
  # 2D intervals
  if("v_effect" %in% colnames(x) && (x |> filter(!is.na(v_effect)) |> nrow()) > 0)  plots$credible_intervals_2D = plot_2D_intervals(
    .data = x, 
    significance_threshold = significance_threshold,
    significance_statistic = significance_statistic,
    show_fdr_message = show_fdr_message
  )
  
  plots
  
}

#' Plot 1D Intervals for Cell-group Effects
#'
#' This function creates a series of 1D interval plots for cell-group effects, highlighting significant differences based on a given significance threshold.
#'
#' @param .data Data frame containing the main data.
#' @param significance_threshold Numeric value specifying the significance threshold for highlighting differences. 
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' @param show_fdr_message Logical. Whether to show the Bayesian FDR interpretation message on the plot. Default is TRUE.
#' @param significance_statistic Character vector indicating which statistic to highlight. Default is "FDR".
#' @importFrom patchwork wrap_plots
#' @importFrom forcats fct_reorder
#' @importFrom tidyr drop_na
#' 
#' @export
#' 
#' @return A combined plot of 1D interval plots.
#' @examples
#' 
#' print("cmdstanr is needed to run this example.")
#' 
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists()) {
#'     data("counts_obj")
#'
#'     estimate <- sccomp_estimate(
#'       counts_obj,
#'       ~ type,
#'       ~1,
#'       "sample",
#'       "cell_group",
#'       "count",
#'       cores = 1
#'     ) |> 
#'     sccomp_test()
#'     
#'   # Example usage:
#'   my_plot = plot_1D_intervals(estimate)
#'     
#'   }
#' }
#'
#' 
plot_1D_intervals = function(
  .data, 
  significance_threshold = 0.05, 
  test_composition_above_logit_fold_change = .data |> attr("test_composition_above_logit_fold_change"), 
  show_fdr_message = TRUE,
  significance_statistic = c("FDR", "pH0")
) {
  significance_statistic <- match.arg(significance_statistic)
  
  # Define the variables as NULL to avoid CRAN NOTES
  parameter <- NULL
  estimate <- NULL
  value <- NULL
  pH0 <- NULL
  FDR <- NULL
  
  .cell_group = attr(.data, ".cell_group")
  
  # Check if test have been done
  if(.data |> select(ends_with("FDR")) |> ncol() |> equals(0))
    stop("sccomp says: to produce plots, you need to run the function sccomp_test() on your estimates.")
  
  plot_list = 
    .data |>
    filter(parameter != "(Intercept)") |>
    
    # Reshape data
    select(-contains("n_eff"), -contains("R_k_hat"), -contains("_ess"), -contains("_rhat")) |> 
    pivot_longer(c(contains("c_"), contains("v_")), names_sep = "_", names_to = c("which", "estimate")) |>
    pivot_wider(names_from = estimate, values_from = value) |>
    
    # Nest data by parameter and which
    nest(data = -c(parameter, which)) |>
    mutate(plot = pmap(
      list(data, which, parameter),
      ~  {
        plot_data <- ..1
        # Check if there are any statistics to plot
        if(plot_data |> filter(!effect |> is.na()) |> nrow() |> equals(0))
          return(NA)
        
        # Choose color variable and legend
        if (significance_statistic == "FDR") {
          color_var <- plot_data$FDR < significance_threshold
          color_aes <- aes(xmin = lower, xmax = upper, color = FDR < significance_threshold)
          color_scale <- scale_color_manual(values = c("grey40", "red"))
          legend_title <- "FDR < significance_threshold"
        } else {
          color_var <- plot_data$pH0 < significance_threshold
          color_aes <- aes(xmin = lower, xmax = upper, color = pH0 < significance_threshold)
          color_scale <- scale_color_manual(values = c("grey40", "red"))
          legend_title <- "pH0 < significance_threshold"
        }
        
        ggplot(plot_data, aes(x = effect, y = fct_reorder(!!.cell_group, effect))) +
          geom_vline(xintercept = test_composition_above_logit_fold_change, colour = "grey") +
          geom_vline(xintercept = -test_composition_above_logit_fold_change, colour = "grey") +
          geom_errorbar(color_aes) +
          geom_point() +
          color_scale +
          xlab("Credible interval of the slope") +
          ylab("Cell group") +
          ggtitle(sprintf("%s %s", ..2, ..3)) +
          sccomp_theme() +
          theme(legend.position = "bottom") +
          guides(color = guide_legend(title = legend_title))
      }
    )) %>%
    
    # Filter out NA plots
    filter(!plot |> is.na()) |> 
    pull(plot) 
  
  # Combine all individual plots into one plot
  combined_plot <- plot_list |>
    wrap_plots(ncol = plot_list |> length() |> sqrt() |> ceiling())

  # Only show the FDR message if significance_statistic == "FDR" and show_fdr_message is TRUE
  if (significance_statistic == "FDR" && show_fdr_message) {
    combined_plot <- combined_plot + theme(plot.caption = ggplot2::element_text(hjust = 0))
    combined_plot <- combined_plot + patchwork::plot_annotation(
      caption = paste(
        "Bayesian FDR: Stephens' method (doi: 10.1093/biostatistics/kxw041)",
        "\nFDR-significant populations may cross fold change thresholds because Bayesian FDR considers posterior probabilities rather than p-values.",
        "\nThe method sorts null hypothesis probabilities in ascending order and calculates cumulative averages for robust false discovery control.",
        sep = ""
      )
    )
  }
  combined_plot
}


#' Plot 2D Intervals for Mean-Variance Association
#'
#' This function creates a 2D interval plot for mean-variance association, highlighting significant differences based on a given significance threshold.
#'
#' @param .data Data frame containing the main data.
#' @param significance_threshold Numeric value specifying the significance threshold for highlighting differences. Default is 0.025.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' @param show_fdr_message Logical. Whether to show the Bayesian FDR interpretation message on the plot. Default is TRUE.
#' @param significance_statistic Character vector indicating which statistic to highlight. Default is "FDR".
#' 
#' 
#' @importFrom dplyr filter arrange mutate if_else row_number
#' @importFrom ggplot2 ggplot geom_vline geom_hline geom_errorbar geom_point annotate aes facet_wrap
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales trans_new
#' @importFrom stringr str_replace
#' @importFrom stats quantile
#' @importFrom magrittr equals
#' 
#' @export
#' 
#' @return A ggplot object representing the 2D interval plot.
#' 
#' @examples
#' 
#' print("cmdstanr is needed to run this example.")
#' 
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists()) {
#'     data("counts_obj")
#'
#'     estimate <- sccomp_estimate(
#'       counts_obj,
#'       ~ type,
#'       ~type,
#'       "sample",
#'       "cell_group",
#'       "count",
#'       cores = 1
#'     ) |> 
#'     sccomp_test()
#'     
#'   # Example usage:
#'   my_plot = plot_2D_intervals(estimate)
#'     
#'   }
#' }
#' 
plot_2D_intervals = function(
    .data, 
    significance_threshold = 0.05, 
    test_composition_above_logit_fold_change = 
      .data |> attr("test_composition_above_logit_fold_change"),
    show_fdr_message = TRUE,
    significance_statistic = c("FDR", "pH0")
){
  significance_statistic <- match.arg(significance_statistic)
  
  # Define the variables as NULL to avoid CRAN NOTES
  v_effect <- NULL
  parameter <- NULL
  c_effect <- NULL
  c_lower <- NULL
  c_upper <- NULL
  c_FDR <- NULL
  v_lower <- NULL
  v_upper <- NULL
  v_FDR <- NULL
  cell_type_label <- NULL
  pH0 <- NULL
  FDR <- NULL
  
  .cell_group = attr(.data, ".cell_group")
  
  # Check if test have been done
  if(.data |> select(ends_with("FDR")) |> ncol() |> equals(0))
    stop("sccomp says: to produce plots, you need to run the function sccomp_test() on your estimates.")
  
  # Extract prec_coeff parameters from the fitted model for regression line
  fit = attr(.data, "fit")
  prec_coeff_summary = fit$summary("prec_coeff")
  prec_coeff_intercept = prec_coeff_summary$mean[1]
  prec_coeff_slope = prec_coeff_summary$mean[2]
  
  # Add adjusted intercept (unbiased variability) as a new facet
  .data_adjusted <- .data %>%
    filter(parameter == "(Intercept)") %>%
    mutate(
      v_effect = v_effect + prec_coeff_slope * c_effect,
      v_lower = v_lower + prec_coeff_slope * c_effect,
      v_upper = v_upper + prec_coeff_slope * c_effect,
      parameter = "(Intercept, adjusted)"
    )
  
  # Bind adjusted data to original (uncorrected comes first)
  .data_plot <- bind_rows(.data, .data_adjusted)
  
  # Always set parameter factor levels to only those present in the data
  if ("parameter" %in% colnames(.data_plot)) {
    .data_plot$parameter <- factor(.data_plot$parameter, levels = unique(.data_plot$parameter))
  }

  # Use .data_plot instead of .data in the rest of the function
  plot <- .data_plot %>%
    # Filter where variance is inferred
    filter(!is.na(v_effect)) %>%
    
    # Add labels for significant cell groups
    with_groups(
      parameter,
      ~ .x %>%
        arrange(c_FDR) %>%
        mutate(cell_type_label = if_else(row_number() <= 3 & c_FDR < significance_threshold & !parameter %in% c("(Intercept)", "(Intercept, adjusted)"), !!.cell_group, ""))
    ) %>%
    with_groups(
      parameter,
      ~ .x %>%
        arrange(v_FDR) %>%
        mutate(cell_type_label = if_else((row_number() <= 3 & v_FDR < significance_threshold & !parameter %in% c("(Intercept)", "(Intercept, adjusted)") ), !!.cell_group, cell_type_label))
    ) %>%
    {
      .x = (.)
      
      # Choose color variable and legend
      if (significance_statistic == "FDR") {
        color_c_aes <- aes(xmin = c_lower, xmax = c_upper, color = c_FDR < significance_threshold & !parameter %in% c("(Intercept)", "(Intercept, adjusted)"), alpha = c_FDR < significance_threshold & !parameter %in% c("(Intercept)", "(Intercept, adjusted)"))
        color_v_aes <- aes(ymin = v_lower, ymax = v_upper, color = v_FDR < significance_threshold & !parameter %in% c("(Intercept)", "(Intercept, adjusted)"), alpha = v_FDR < significance_threshold & !parameter %in% c("(Intercept)", "(Intercept, adjusted)"))
        color_scale <- scale_color_manual(values = c("#D3D3D3", "#E41A1C"))
        alpha_scale <- scale_alpha_manual(values = c(0.4, 1))
        legend_title <- "FDR < significance_threshold"
      } else {
        color_c_aes <- aes(xmin = c_lower, xmax = c_upper, color = c_pH0 < significance_threshold & !parameter %in% c("(Intercept)", "(Intercept, adjusted)"), alpha = c_pH0 < significance_threshold & !parameter %in% c("(Intercept)", "(Intercept, adjusted)"))
        color_v_aes <- aes(ymin = v_lower, ymax = v_upper, color = v_pH0 < significance_threshold & !parameter %in% c("(Intercept)", "(Intercept, adjusted)"), alpha = v_pH0 < significance_threshold & !parameter %in% c("(Intercept)", "(Intercept, adjusted)"))
        color_scale <- scale_color_manual(values = c("#D3D3D3", "#377EB8"))
        alpha_scale <- scale_alpha_manual(values = c(0.4, 1))
        legend_title <- "pH0 < significance_threshold"
      }
      
      # Calculate range for regression line
      c_range = range(.x$c_effect, na.rm = TRUE)
      c_seq = seq(c_range[1], c_range[2], length.out = 100)
      v_pred = -(prec_coeff_intercept + prec_coeff_slope * c_seq)
      
      # Create regression line data only for (Intercept)
      regression_data = data.frame(
        c_effect = c_seq,
        v_effect = v_pred,
        parameter = "(Intercept)"
      )
      
      # Plot
      p <- ggplot(.x, aes(c_effect, v_effect)) +
        # Add vertical and horizontal lines
        geom_vline(xintercept = c(-test_composition_above_logit_fold_change, test_composition_above_logit_fold_change), colour = "grey", linetype = "dashed", linewidth = 0.3) +
        geom_hline(yintercept = c(-test_composition_above_logit_fold_change, test_composition_above_logit_fold_change), colour = "grey", linetype = "dashed", linewidth = 0.3)
      
       # Add regression line only for (Intercept) facet
      p <- p + geom_line(
        data = regression_data,
        mapping = aes(c_effect, v_effect),
        color = "#0072B2", linewidth = 0.5, alpha = 0.8,
        inherit.aes = FALSE
      )

      # Add horizontal line for (Intercept, adjusted) facet
      if ("(Intercept, adjusted)" %in% unique(.x$parameter)) {
        mean_adjusted <- mean(.x$v_effect[.x$parameter == "(Intercept, adjusted)"])
        c_range_adjusted <- range(.x$c_effect[.x$parameter == "(Intercept, adjusted)"], na.rm = TRUE)
        
        # Create horizontal line data for adjusted intercept only
        horizontal_line_data <- data.frame(
          c_effect = c_range_adjusted,
          v_effect = rep(mean_adjusted, 2),
          parameter = "(Intercept, adjusted)"
        )
        
        p <- p +
          geom_line(
            data = horizontal_line_data,
            mapping = aes(c_effect, v_effect),
            color = "#0072B2", linewidth = 0.5, alpha = 0.8,
            inherit.aes = FALSE
          )
      }
    
      p <- p +
        # Add error bars
        geom_errorbar(color_c_aes, linewidth = 0.2) +
        geom_errorbar(color_v_aes, linewidth = 0.2) +
        # Add points
        geom_point(size = 0.2) +
        # Add annotations
        # annotate("text", x = 0, y = 3.5, label = "Variability", size = 2) + # Disabled temporarily; revisit if annotations are needed for variability.
        # annotate("text", x = 5, y = 0, label = "Abundance", size = 2, angle = 270) + # Disabled temporarily; revisit if annotations are needed for abundance.
        # Add text labels for significant cell groups
        geom_text_repel(aes(c_effect, -v_effect, label = cell_type_label), size = 2.5, data = .x %>% filter(cell_type_label != "")) +
        # Set color and alpha scales
        color_scale +
        alpha_scale +
        # Facet by parameter
        facet_wrap(
          ~ fct_relevel(parameter, c("(Intercept)", "(Intercept, adjusted)")) |>
            
            # I have to understand why this works but throws warning
            suppressWarnings(),
          scales = "free"
        ) +
        xlab("c_effect (Abundance effect)") +
        ylab("v_effect (Variability effect)") +
        # Apply custom theme
        sccomp_theme() +
        theme(legend.position = "bottom") +
        guides(color = guide_legend(title = legend_title), alpha = "none")
      p
    }

  # Only show the FDR message if significance_statistic == "FDR" and show_fdr_message is TRUE
  if (significance_statistic == "FDR" && show_fdr_message) {
    plot <- plot + theme(plot.caption = ggplot2::element_text(hjust = 0))
    plot <- plot + patchwork::plot_annotation(
      caption = paste(
        "Bayesian FDR: Stephens' method (doi: 10.1093/biostatistics/kxw041)",
        "\nFDR-significant populations may cross fold change thresholds because Bayesian FDR considers posterior probabilities rather than p-values.",
        "\nThe method sorts null hypothesis probabilities in ascending order and calculates cumulative averages for robust false discovery control.",
        sep = ""
      )
    )
  }
  plot
}




#' Plot Scatterplot of Cell-group Proportion
#'
#' This function creates a scatterplot of cell-group proportions, optionally highlighting significant differences based on a given significance threshold.
#'
#' @param .data Data frame containing the main data.
#' @param data_proportion Data frame containing proportions of cell groups.
#' @param factor_of_interest A factor indicating the biological condition of interest.
#' @param .cell_group The cell group to be analysed.
#' @param .sample The sample identifier.
#' @param significance_threshold Numeric value specifying the significance threshold for highlighting differences. Default is 0.025.
#' @param my_theme A ggplot2 theme object to be applied to the plot.
#' @importFrom scales trans_new
#' @importFrom stringr str_replace
#' @importFrom stats quantile
#' @importFrom magrittr equals
#' 
#' 
#' @return A ggplot object representing the scatterplot.
#' @examples
#' # Example usage:
#' # plot_scatterplot(.data, data_proportion, "condition", "cell_group", "sample", 0.025, theme_minimal())
plot_scatterplot = function(
    .data, data_proportion, factor_of_interest, .cell_group,
    .sample, significance_threshold = 0.05, my_theme
){
  
  # Define the variables as NULL to avoid CRAN NOTES
  stats_name <- NULL
  parameter <- NULL
  stats_value <- NULL
  count_data <- NULL
  generated_proportions <- NULL
  proportion <- NULL
  name <- NULL
  outlier <- NULL
  
  # Function to remove leading zero from labels
  dropLeadingZero <- function(l){  stringr::str_replace(l, '0(?=.)', '') }
  
  # Define square root transformation and its inverse
  S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
  IS_sqrt <- function(x){x^2*sign(x)}
  S_sqrt_trans <- function() scales::trans_new("S_sqrt",S_sqrt,IS_sqrt)
  
  .cell_group = enquo(.cell_group)
  .sample = enquo(.sample)
  
  # Prepare significance colors
  significance_colors =
    .data %>%
    pivot_longer(
      c(contains("c_"), contains("v_")),
      names_pattern = "([cv])_([a-zA-Z0-9]+)",
      names_to = c("which", "stats_name"),
      values_to = "stats_value"
    ) %>%
    filter(stats_name == "FDR") %>%
    filter(parameter != "(Intercept)") %>%
    filter(stats_value < significance_threshold) %>%
    filter(`factor` == factor_of_interest) 
  
  if(nrow(significance_colors) > 0){
    
    if(.data |> attr("contrasts") |> is.null())
      significance_colors =
        significance_colors %>%
        unite("name", c(which, parameter), remove = FALSE) %>%
        distinct() %>%
        
        # Get clean parameter
        mutate(!!as.symbol(factor_of_interest) := str_replace(parameter, sprintf("^%s", `factor`), "")) %>%
        with_groups(c(!!.cell_group, !!as.symbol(factor_of_interest)), ~ .x %>% summarise(name = paste(name, collapse = ", ")))
    else
      significance_colors =
        significance_colors |>
        mutate(
          factor_values = attr(.data, "count_data") |>
            select(all_of(factor_of_interest)) |>
            distinct() |>
            pull(all_of(factor_of_interest))
        ) |>
        unnest(factor_values) |>
        
        # Filter relevant parameters
        mutate( !!as.symbol(factor_of_interest) := as.character(factor_values) ) |>
        filter(str_detect(parameter, !!as.symbol(factor_of_interest) )) |>
        
        # Rename
        select(!!.cell_group, !!as.symbol(factor_of_interest), name = parameter) |>
        
        # Merge contrasts
        with_groups(c(!!.cell_group, !!as.symbol(factor_of_interest)), ~ .x %>% summarise(name = paste(name, collapse = ", ")))
  }
  
  my_scatterplot = ggplot()
  
  if("fit" %in% names(attributes(.data))){
    
    simulated_proportion =
      .data |>
      sccomp_replicate(number_of_draws = 1000) |>
      left_join(data_proportion %>% distinct(!!as.symbol(factor_of_interest), !!.sample, !!.cell_group))
    
    my_scatterplot = 
      my_scatterplot +
      
      # Add smoothed line for simulated proportions
      geom_smooth(
        aes(!!as.symbol(factor_of_interest), (generated_proportions)),
        lwd=0.2,
        data =
          simulated_proportion %>%
          inner_join(data_proportion %>% distinct(!!as.symbol(factor_of_interest), !!.cell_group, !!.sample)) ,
        color="blue", fill="blue",
        span = 1
      )
  }
  
  if(
    nrow(significance_colors)==0 ||
    
    significance_colors |> 
    pull(!!as.symbol(factor_of_interest)) |> 
    intersect(
      data_proportion |> 
      pull(!!as.symbol(factor_of_interest))
    ) |> 
    length() |> 
    equals(0)
  ) {
    
    my_scatterplot=
      my_scatterplot +
      
      # Add smoothed line without significance colors
      geom_smooth(
        aes(!!as.symbol(factor_of_interest), proportion, fill = NULL),
        data =
          data_proportion ,
        lwd=0.5,
        color = "black",
        span = 1
      )
  } else {
    my_scatterplot=
      my_scatterplot +
      
      # Add smoothed line with significance colors
      geom_smooth(
        aes(!!as.symbol(factor_of_interest), proportion, fill = name),
        data = data_proportion ,
        fatten = 0.5,
        lwd=0.5,
        color = "black",
        span = 1
      )
  }
  
  my_scatterplot +
    
    # Add jittered points for individual data
    geom_point(
      aes(!!as.symbol(factor_of_interest), proportion, shape=outlier, color=outlier),
      data = data_proportion,
      position=position_jitterdodge(jitter.height = 0, jitter.width = 0.2),
      size = 0.5
    ) +
    
    # Facet wrap by cell group
    facet_wrap(
      vars(!!.cell_group),
      scales = "free_y",
      nrow = 4
    ) +
    scale_color_manual(values = c("black", "#e11f28")) +
    scale_y_continuous(trans=S_sqrt_trans(), labels = dropLeadingZero) +
    scale_fill_discrete(na.value = "white") +
    xlab("Biological condition") +
    ylab("Cell-group proportion") +
    guides(color="none", alpha="none", size="none") +
    labs(fill="Significant difference") +
    ggtitle("Note: Be careful judging significance (or outliers) visually for lowly abundant cell groups. \nVisualising proportion hides the uncertainty characteristic of count data, that a count-based statistical model can estimate.") +
    my_theme +
    theme(axis.text.x =  element_text(angle=20, hjust = 1), title = element_text(size = 3))
}
