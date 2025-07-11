#' sccomp_boxplot
#'
#' @description
#' Creates a boxplot visualization of the model results from `sccomp`. This function plots the estimated cell proportions across samples, highlighting significant changes in cell composition according to a specified factor.
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr pivot_longer pivot_wider unite unnest
#' @importFrom dplyr select with_groups mutate left_join rename
#'
#' @param .data A tibble containing the results from `sccomp_estimate` and `sccomp_test`, including the columns: cell_group name, sample name, read counts, factor(s), p-values, and significance indicators.
#' @param factor A character string specifying the factor of interest included in the model for stratifying the boxplot.
#' @param significance_threshold A numeric value indicating the False Discovery Rate (FDR) threshold for labeling significant cell-groups. Defaults to 0.05.
#' @param test_composition_above_logit_fold_change A positive numeric value representing the effect size threshold used in the hypothesis test. A value of 0.2 corresponds to a change in cell proportion of approximately 10% for a cell type with a baseline proportion of 50% (e.g., from 45% to 55%). This threshold is consistent on the logit-unconstrained scale, even when the baseline proportion is close to 0 or 1.
#' @param remove_unwanted_effects A logical value indicating whether to remove unwanted variation from the data before plotting. Defaults to `FALSE`.
#'
#' @return A `ggplot` object representing the boxplot of cell proportions across samples, stratified by the specified factor.
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
#' if (instantiate::stan_cmdstan_exists()) {
#'     data("counts_obj")
#'
#'     estimate <- sccomp_estimate(
#'       counts_obj,
#'       formula_composition = ~ type,
#'       formula_variability = ~ 1,
#'       sample = "sample",
#'       cell_group = "cell_group",
#'       abundance = "count",
#'       cores = 1
#'     ) |>
#'     sccomp_test()
#'
#'     # Plot the boxplot of estimated cell proportions
#'     sccomp_boxplot(
#'         .data = estimate,
#'         factor = "type",
#'         significance_threshold = 0.05
#'     )
#' }
#' }
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
sccomp_boxplot = function(
    .data, 
    factor, 
    significance_threshold = 0.05, 
    test_composition_above_logit_fold_change = .data |> attr("test_composition_above_logit_fold_change"),
    remove_unwanted_effects = FALSE
){
  
  
  .cell_group = attr(.data, ".cell_group")
  .count = attr(.data, ".count")
  .sample = attr(.data, ".sample")
  
  # Check if test have been done
  if(.data |> select(ends_with("FDR")) |> ncol() |> equals(0))
    stop("sccomp says: to produce plots, you need to run the function sccomp_test() on your estimates.")
  
  
  data_proportion =
    .data |>
    
    # Otherwise does not work
    select(-`factor`) |>
    
    pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) |>
    left_join(
      attr(.data, "count_data") ,
      by = quo_name(.cell_group)
    ) |>
    with_groups(!!.sample, ~ mutate(.x, proportion = (!!.count)/sum(!!.count)) ) |> 
    mutate(is_zero = proportion==0) 
  
  if(remove_unwanted_effects){
    .data_adjusted = 
      .data |> 
      sccomp_remove_unwanted_effects(formula_composition_keep = as.formula("~ " |> paste(factor))) |> 
      rename(proportion = adjusted_proportion)
    
    data_proportion = 
      data_proportion |> 
      select(-proportion) |> 
      left_join(.data_adjusted, by = join_by(!!.cell_group, !!.sample))
  }
  else 
    message( "sccomp says: When visualising proportions, especially for complex models, consider setting `remove_unwanted_effects=TRUE`. This will adjust the proportions, preserving only the observed effect.")
  
  # If I don't have outliers add them
  if(!"outlier" %in% colnames(data_proportion)) data_proportion = data_proportion |> mutate(outlier = FALSE) 
  
  plot_boxplot(
    .data,
    data_proportion,
    factor,
    !!.cell_group,
    !!.sample,
    significance_threshold = significance_threshold,
    create_multipanel_theme(),
    remove_unwanted_effects = remove_unwanted_effects
  ) +
    ggtitle(sprintf("Grouped by %s (for multi-factor models, associations could be hardly observable with unidimensional data stratification)", factor))
  
  
}

#' Plot Boxplot of Cell-group Proportion
#'
#' This function creates a boxplot of cell-group proportions, optionally highlighting significant differences based on a given significance threshold.
#'
#' @param .data Data frame containing the main data.
#' @param data_proportion Data frame containing proportions of cell groups.
#' @param factor_of_interest Character string specifying the biological condition of interest.
#' @param .cell_group Character string specifying the cell group to be analyzed.
#' @param .sample Character string specifying the sample identifier.
#' @param significance_threshold Numeric value specifying the significance threshold for highlighting differences. Default is 0.05.
#' @param my_theme A ggplot2 theme object to be applied to the plot.
#' @param remove_unwanted_effects Logical value indicating whether to remove unwanted effects. Default is FALSE.
#' 
#' @importFrom scales trans_new
#' @importFrom stringr str_replace
#' @importFrom stats quantile
#'
#' @return A ggplot object representing the boxplot.
#' 
#' @examples
#' # Example usage:
#' # plot_boxplot(
#' #   .data, data_proportion, factor_of_interest = "condition",
#' #   .cell_group = "cell_group", .sample = "sample",
#' #   significance_threshold = 0.05, my_theme = theme_minimal(),
#' #   remove_unwanted_effects = FALSE
#' # )
#' 
#' @noRd
plot_boxplot = function(
    .data, data_proportion, factor_of_interest, .cell_group,
    .sample, 
    significance_threshold = 0.05, 
    my_theme, 
    remove_unwanted_effects = FALSE
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
  
  # Function to calculate boxplot statistics
  calc_boxplot_stat <- function(x) {
    coef <- 1.5
    n <- sum(!is.na(x))
    
    # Calculate quantiles
    stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
    names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
    iqr <- diff(stats[c(2, 4)])
    
    # Set whiskers
    outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
    if (any(outliers)) {
      stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
    }
    return(stats)
  }
  
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
  
  my_boxplot = ggplot()
  
  if("fit" %in% names(attributes(.data))){
    
    # Remove unwanted effects?
    if(remove_unwanted_effects) formula_composition = as.formula("~ " |> paste(factor_of_interest))
    else formula_composition = NULL
    
    simulated_proportion =
      .data |>
      sccomp_replicate(formula_composition = formula_composition, number_of_draws = 100) |>
      left_join(data_proportion %>% distinct(!!as.symbol(factor_of_interest), !!.sample, !!.cell_group, is_zero))
    
    my_boxplot = my_boxplot +
      
      # Add boxplot for simulated proportions
      stat_summary(
        aes(!!as.symbol(factor_of_interest), (generated_proportions)),
        fun.data = calc_boxplot_stat, geom="boxplot",
        outlier.shape = NA, outlier.color = NA, outlier.size = 0,
        fatten = 0.5, lwd=0.2,
        data =
          simulated_proportion %>%
          inner_join(data_proportion %>% distinct(!!as.symbol(factor_of_interest), !!.cell_group)),
        color="blue"
      )
  }
  
  if(nrow(significance_colors) == 0 ||
     length(intersect(
       significance_colors |> pull(!!as.symbol(factor_of_interest)),
       data_proportion |> pull(!!as.symbol(factor_of_interest))
     )) == 0){
    
    my_boxplot=
      my_boxplot +
      
      # Add boxplot without significance colors
      geom_boxplot(
        aes(!!as.symbol(factor_of_interest), proportion,  group=!!as.symbol(factor_of_interest), fill = NULL),
        outlier.shape = NA, outlier.color = NA,outlier.size = 0,
        data =
          data_proportion |>
          mutate(!!as.symbol(factor_of_interest) := as.character(!!as.symbol(factor_of_interest))) ,
        fatten = 0.5,
        lwd=0.5,
      )
  } else {
    my_boxplot=
      my_boxplot +
      
      # Add boxplot with significance colors
      geom_boxplot(
        aes(!!as.symbol(factor_of_interest), proportion,  group=!!as.symbol(factor_of_interest), fill = name),
        outlier.shape = NA, outlier.color = NA,outlier.size = 0,
        data =
          data_proportion |>
          mutate(!!as.symbol(factor_of_interest) := as.character(!!as.symbol(factor_of_interest))) %>%
          left_join(significance_colors, by = c(quo_name(.cell_group), factor_of_interest)),
        fatten = 0.5,
        lwd=0.5,
      )
  }
  
  my_boxplot +
    
    # Add jittered points for individual data
    geom_jitter(
      aes(!!as.symbol(factor_of_interest), proportion, shape=is_zero, color=outlier,  group=!!as.symbol(factor_of_interest)),
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
    scale_shape_manual(values = c(16, 21)) +
    scale_y_continuous(trans=S_sqrt_trans(), labels = dropLeadingZero) +
    scale_fill_discrete(na.value = "white") +
    xlab("Biological condition") +
    ylab("Cell-group proportion") +
    guides( alpha="none", size="none") +
    labs(fill="Significant difference") +
    my_theme +
    theme(axis.text.x =  element_text(angle=20, hjust = 1), title = element_text(size = 3))
}