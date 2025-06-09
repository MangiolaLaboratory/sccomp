

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
#' @param ... For internal use 
#'
#' @return A `ggplot`
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
#'       ~ type, ~1, sample, cell_group, count,
#'       cores = 1
#'     )
#'
#'     # estimate |> plot()
#'   }
#' }
#'
plot.sccomp_tbl <- function(x,  significance_threshold = 0.05, test_composition_above_logit_fold_change = .data |> attr("test_composition_above_logit_fold_change"), ...) {
  
  # Define the variables as NULL to avoid CRAN NOTES
  parameter <- NULL
  count_data <- NULL
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
    unnest(count_data) |>
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
                my_theme = multipanel_theme,
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
  plots$credible_intervals_1D = plot_1D_intervals(.data = x, significance_threshold = significance_threshold)
  
  # 2D intervals
  if("v_effect" %in% colnames(x) && (x |> filter(!is.na(v_effect)) |> nrow()) > 0)  plots$credible_intervals_2D = plot_2D_intervals(.data = x, significance_threshold = significance_threshold)
  
  plots
  
}
