

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
#'       .sample = sample,
#'       .cell_group = cell_group,
#'       .count = count,
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
    unnest(count_data) |>
    with_groups(!!.sample, ~ mutate(.x, proportion = (!!.count)/sum(!!.count)) ) |> 
    mutate(is_zero = proportion==0) 
  
  if(remove_unwanted_effects){
    .data_adjusted = 
      .data |> 
      sccomp_remove_unwanted_variation(formula_composition_keep = as.formula("~ " |> paste(factor))) |> 
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
    multipanel_theme,
    remove_unwanted_effects = remove_unwanted_effects
  ) +
    ggtitle(sprintf("Grouped by %s (for multi-factor models, associations could be hardly observable with unidimensional data stratification)", factor))
  
  
}