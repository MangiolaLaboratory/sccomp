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
      add_marginal_density = add_marginal_density
    )
  }

  plots

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

