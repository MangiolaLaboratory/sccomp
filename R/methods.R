
#' sccomp_glm main
#'
#' @description This function runs the data modelling and statistical test for the hypothesis that a cell_type includes outlier biological replicate.
#'
#' @importFrom tibble as_tibble
#' @import dplyr
#' @importFrom tidyr spread
#' @import tidybayes
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom purrr map2
#' @importFrom purrr map_int
#' @importFrom tidybulk scale_abundance
#' @importFrom benchmarkme get_ram
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom purrr map
#' @importFrom tibble rowid_to_column
#' @importFrom furrr future_map
#' @importFrom purrr map_lgl
#'
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | covariate columns | Pvaue column | a significance column
#' @param formula A formula. The sample formula used to perform the differential cell_type abundance analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_type identifier
#' @param .count A column name as symbol. The cell_type abundance (read count)
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param seed An integer. Used for development and testing purposes
#'
#' @return A nested tibble `tbl` with cell_type-wise information: `sample wise data` | plot | `ppc samples failed` | `exposure deleterious outliers`
#'
#' @export
#'
sccomp_glm = function(.data,
                      formula = ~ 1,
                      .sample,
                      .cell_group,

                      # Main arguments
                      check_outliers = TRUE,
                      approximate_posterior_inference = TRUE,
                      verbose = FALSE,
                      noise_model = "multi_beta_binomial",
                      seed = sample(1:99999, size = 1)
) {

  cores = 4 #detect_cores()

  # See https://community.rstudio.com/t/how-to-make-complete-nesting-work-with-quosures-and-tidyeval/16473
  # See https://github.com/tidyverse/tidyr/issues/506
  .sample_for_tidyr = enexpr(.sample)
  .cell_group_for_tidyr = enexpr(.cell_group)

  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)



  # Check if columns exist
  check_columns_exist(.data, !!.sample, !!.cell_group,  parse_formula(formula))

  # Check if any column is NA or null
  check_if_any_NA(.data, !!.sample, !!.cell_group, parse_formula(formula))

  # Make counts
  .data_count =
    .data %>%
    count(!!.sample, !!.cell_group, !!as.symbol(parse_formula(formula)), name = "count") %>%
    complete(
      nesting( !!.sample_for_tidyr,  !!as.symbol(parse_formula(formula))),
      !!.cell_group_for_tidyr,
      fill = list(count = 0)
    ) %>%
    mutate(count = as.integer(count))

  # Choose linear model
  my_glm =
    noise_model %>%
    when(
      equals(., "multi_beta_binomial") ~ multi_beta_binomial_glm,
      #equals(., "multi_beta") ~ multi_beta_glm,
      equals(., "dirichlet_multinomial") ~ dirichlet_multinomial_glm
    )

  # Return
  .data_count %>%
    my_glm(
      formula = formula,
      .sample = !!.sample,
      .cell_type = !!.cell_group,
      .count = count,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      cores = cores,
      verbose = verbose,
      seed = seed
    )


}
