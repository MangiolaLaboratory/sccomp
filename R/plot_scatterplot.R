#' sccomp_scatterplot
#'
#' @description
#' Creates a scatterplot of cell-group proportions against a continuous covariate,
#' optionally overlaying the posterior predictive trend from the fitted model.
#' Use this visualisation when the factor of interest is numeric; for discrete
#' factors prefer [sccomp_boxplot()].
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer pivot_wider unite unnest
#' @importFrom dplyr select with_groups mutate left_join rename
#'
#' @param .data A tibble containing the results from `sccomp_estimate` and
#'   `sccomp_test`, including the columns: cell_group name, sample name, read
#'   counts, factor(s), p-values, and significance indicators.
#' @param factor A character string specifying the continuous covariate of
#'   interest included in the model.
#' @param significance_threshold A numeric value indicating the threshold for
#'   highlighting significant cell-groups. Defaults to 0.05.
#' @param remove_unwanted_effects A logical value indicating whether to remove
#'   unwanted variation from the data before plotting. Defaults to `FALSE`.
#'
#' @return A `ggplot` object representing the scatterplot of cell proportions
#'   against the continuous covariate, faceted by cell group.
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
#'     data("seurat_obj")
#'
#'     estimate <- sccomp_estimate(
#'       seurat_obj,
#'       formula_composition = ~ continuous_covariate,
#'       formula_variability = ~ 1,
#'       sample = "sample",
#'       cell_group = "cell_group",
#'       cores = 1
#'     ) |>
#'     sccomp_test()
#'
#'     # Plot proportions against the continuous covariate
#'     sccomp_scatterplot(
#'         .data = estimate,
#'         factor = "continuous_covariate",
#'         significance_threshold = 0.05
#'     )
#' }
#' }
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
sccomp_scatterplot = function(
    .data,
    factor,
    significance_threshold = 0.05,
    remove_unwanted_effects = FALSE
){

  .cell_group = attr(.data, ".cell_group")
  .count = attr(.data, ".count")
  .sample = attr(.data, ".sample")

  # Check if test have been done
  if(.data |> select(ends_with("FDR")) |> ncol() |> equals(0))
    stop("sccomp says: to produce plots, you need to run the function sccomp_test() on your estimates.")

  .data_filtered =
    subset_results_by_factor(.data, factor, keep_intercept = FALSE)

  data_proportion =
    .data_filtered |>

    # Otherwise does not work
    select(-`factor`) |>

    pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) |>
    left_join(
      attr(.data, "count_data") ,
      by = quo_name(.cell_group)
    ) |>
    with_groups(!!.sample, ~ mutate(.x, proportion = (!!.count)/sum(!!.count)) )

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
  if(.data |> attr("outliers") |> is.null() |> not())
    data_proportion =
    data_proportion |>
    left_join(.data |> attr("outliers"), by = c(quo_name(.sample), quo_name(.cell_group)))
 else
  data_proportion = data_proportion |> mutate(outlier = FALSE)

  plot_scatterplot(
    .data = .data,
    data_proportion = data_proportion,
    factor_of_interest = factor,
    significance_threshold = significance_threshold,
    my_theme = sccomp_theme()
  ) +
    ggtitle(sprintf("Grouped by %s (for multi-factor models, associations could be hardly observable with unidimensional data stratification)", factor))

}

#' Plot Scatterplot of Cell-group Proportion
#'
#' This function creates a scatterplot of cell-group proportions, optionally
#' highlighting significant differences based on a given significance threshold.
#'
#' @param .data Data frame containing the main data.
#' @param data_proportion Data frame containing proportions of cell groups.
#' @param factor_of_interest A factor indicating the biological condition of interest.
#' @param significance_threshold Numeric value specifying the significance threshold
#'   for highlighting differences. Default is 0.05.
#' @param my_theme A ggplot2 theme object to be applied to the plot.
#' @importFrom scales trans_new
#' @importFrom stringr str_replace str_detect
#'
#' @return A ggplot object representing the scatterplot.
#'
#' @noRd
plot_scatterplot = function(
    .data, data_proportion, factor_of_interest,
    significance_threshold = 0.05, my_theme
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
  
  .cell_group = attr(.data, ".cell_group")
  .count = attr(.data, ".count")
  .sample = attr(.data, ".sample")
  
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
        linewidth = 0.5,
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
