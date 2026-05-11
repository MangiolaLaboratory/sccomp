#' plot
#'
#' @description This function plots a summary of the results of the model.
#'
#' @importFrom dplyr filter select
#' @importFrom magrittr equals
#'
#' @param x A tibble including a cell_group name column | sample name column | read counts column | factor columns | Pvalue column | a significance column
#' @param significance_threshold Numeric value specifying the significance threshold for highlighting differences. Default is 0.05.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' @param significance_statistic Character vector indicating which statistic to highlight. Default is "pH0".
#' @param show_fdr_message Logical. Whether to show the Bayesian FDR interpretation message on the plot. Default is TRUE.
#' @param add_marginal_density Logical. Whether to add marginal density plots on adjusted panels in 2D intervals. Default is TRUE.
#' @param omit_ci Logical. Whether to omit credible interval error bars from 2D interval plots. Default is FALSE.
#' @param ... For internal use
#'
#' @return A list containing ggplot objects
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
#'     ) |>
#'     sccomp_test()
#'
#'     plots = estimate |> plot()
#'   }
#' }
#'
plot.sccomp_tbl <- function(
    x,
    significance_threshold = 0.05,
    test_composition_above_logit_fold_change = x |> attr("test_composition_above_logit_fold_change"),
    significance_statistic = c("pH0", "FDR"),
    show_fdr_message = TRUE,
    add_marginal_density = TRUE,
    omit_ci = FALSE,
    sort_by = c("none", "effect", "significance", "alphabetical"),
    ...
) {

  significance_statistic <- match.arg(significance_statistic)
  sort_by <- match.arg(sort_by)

  # Quosures from estimate (same names as legacy `.cell_group` args; stored on tbl attributes)
  .sample <- attr(x, ".sample")
  .cell_group <- attr(x, ".cell_group")
  .count <- attr(x, ".count")

  # Define the variables as NULL to avoid CRAN NOTES
  v_effect <- NULL

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
      attr(x, "count_data"),
      by = join_by(!!.cell_group)
    ) |> 
    with_groups(!!.sample, ~ mutate(.x, proportion = (!!.count)/sum(!!.count)) )
  
  
  # If I don't have outliers add them
  if(x |> attr("outliers") |> is.null() |> not())
    data_proportion = 
    data_proportion |> 
    left_join(x |> attr("outliers"), by = join_by(!!.sample, !!.cell_group))
 else 
  data_proportion = data_proportion |> mutate(outlier = FALSE) 
  
 

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
                significance_threshold = significance_threshold,
                my_theme = sccomp_theme()
              )
          
          # If discrete
          else 
            my_plot = 
              sccomp_boxplot(
                .data = x,
                factor = .x,
                significance_threshold = significance_threshold,
                significance_statistic = significance_statistic
              ) 
          
          # Return
          my_plot +
            ggtitle(sprintf("Grouped by %s (for multi-factor models, associations could be hardly observable with unidimensional data stratification)", .x))
        } )
    
  }
  

  # 1D intervals
  plots$credible_intervals_1D = sccomp_plot_intervals_1D(
    .data = x,
    significance_threshold = significance_threshold,
    test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
    significance_statistic = significance_statistic,
    show_fdr_message = show_fdr_message,
    sort_by = sort_by
  )

  # 2D intervals (only if variance effects exist)
  if("v_effect" %in% colnames(x) && (x |> filter(!is.na(v_effect)) |> nrow()) > 0) {
    plots$credible_intervals_2D = sccomp_plot_intervals_2D(
      .data = x,
      significance_threshold = significance_threshold,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      significance_statistic = significance_statistic,
      show_fdr_message = show_fdr_message,
      add_marginal_density = add_marginal_density,
      omit_ci = omit_ci
    )
  }

  plots

}
