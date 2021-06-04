

#' sccomp_glm main
#'
#' @description This function runs the data modelling and statistical test for the hypothesis that a cell_type includes outlier biological replicate.
#'
#' @import dplyr
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom rlang quo_is_null
#' @importFrom SingleCellExperiment colData
#' @importFrom parallel detectCores
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
#' @return A nested tibble `tbl` with cell_group-wise statistics
#'
#' @export
#'
#'
sccomp_glm <- function(.data,
                       formula ,
                       .sample,
                       .cell_group,
                       .count = NULL,
                       # Secondary arguments
                       check_outliers = TRUE,
                       approximate_posterior_inference = FALSE,
                       verbose = FALSE,
                       noise_model = "multi_beta_binomial",
                       cores = detectCores(),
                       seed = 42) {
  UseMethod("sccomp_glm", .data)
}

#' @export
sccomp_glm.Seurat = function(.data,
                             formula ,
                             .sample,
                             .cell_group,
                             .count = NULL,
                             # Secondary arguments
                             check_outliers = TRUE,
                             approximate_posterior_inference = FALSE,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             cores = detectCores(),
                             seed = 42) {

  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")

  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)

  .data[[]] %>%
    sccomp_glm(
      formula = formula,!!.sample,!!.cell_group,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      verbose = verbose,
      noise_model = noise_model,
      cores = cores,
      seed = seed
    )


}

#' @export
sccomp_glm.SingleCellExperiment = function(.data,
                                           formula ,
                                           .sample,
                                           .cell_group,
                                           .count = NULL,
                                           # Secondary arguments
                                           check_outliers = TRUE,
                                           approximate_posterior_inference = FALSE,
                                           verbose = FALSE,
                                           noise_model = "multi_beta_binomial",
                                           cores = detectCores(),
                                           seed = 42) {

  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")


  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)

  .data %>%
    colData() %>%
    sccomp_glm(
      formula = formula,!!.sample,!!.cell_group,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      verbose = verbose,
      noise_model = noise_model,
      cores = cores,
      seed = seed
    )


}

#' @export
sccomp_glm.DFrame = function(.data,
                             formula ,
                             .sample,
                             .cell_group,
                             .count = NULL,

                             # Secondary arguments
                             check_outliers = TRUE,
                             approximate_posterior_inference = FALSE,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             cores = detectCores(),
                             seed = 42) {

  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")


  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)

  .data %>%
    as.data.frame %>%
    sccomp_glm(
      formula = formula,!!.sample,!!.cell_group,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      verbose = verbose,
      noise_model = noise_model,
      cores = cores,
      seed = seed
    )
}

#' @importFrom purrr when
#' @export
sccomp_glm.data.frame = function(.data,
                                 formula ,
                                 .sample,
                                 .cell_group,
                                 .count = NULL,

                                 # Secondary arguments
                                 check_outliers = TRUE,
                                 approximate_posterior_inference = FALSE,
                                 verbose = FALSE,
                                 noise_model = "multi_beta_binomial",
                                 cores = detectCores(),
                                 seed = 42) {

  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)

  # Choose linear model
  my_glm_model =
    noise_model %>%
    when(
      equals(., "multi_beta_binomial") ~ multi_beta_binomial_glm,
      equals(., "dirichlet_multinomial") ~ dirichlet_multinomial_glm
    )

  .count %>%
    when(
      # If the dataframe does not include counts, but is metadata
      quo_is_null(.) ~ sccomp_glm_data_frame_raw(
        .data,
        formula = formula,
        !!.sample,
        !!.cell_group,
        check_outliers = check_outliers,
        approximate_posterior_inference = approximate_posterior_inference,
        verbose = verbose,
        my_glm_model = my_glm_model,
        cores = cores,
        seed = seed
      ),

      # If the dataframe does includes counts
      ~ sccomp_glm_data_frame_counts(
        .data,
        formula = formula,
        !!.sample,
        !!.cell_group,
        !!.count,
        check_outliers = check_outliers,
        approximate_posterior_inference = approximate_posterior_inference,
        verbose = verbose,
        my_glm_model = my_glm_model,
        cores = cores,
        seed = seed
      )
    )
}

#' @importFrom tidyr complete
#' @importFrom tidyr nesting
sccomp_glm_data_frame_raw = function(.data,
                                     formula,
                                     .sample,
                                     .cell_group,
                                     my_glm_model,
                                     # Secondary arguments
                                     check_outliers = TRUE,
                                     approximate_posterior_inference = FALSE,
                                     verbose = FALSE,
                                     noise_model = "multi_beta_binomial",
                                     cores = 4,
                                     seed = sample(1:99999, size = 1)) {

  # See https://community.rstudio.com/t/how-to-make-complete-nesting-work-with-quosures-and-tidyeval/16473
  # See https://github.com/tidyverse/tidyr/issues/506
  .sample_for_tidyr = ensym(.sample)
  .cell_group_for_tidyr = ensym(.cell_group)

  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)

  # Check if columns exist
  check_columns_exist(.data, c(
    quo_name(.sample),
    quo_name(.cell_group),
    parse_formula(formula)
  ))

  # Check if any column is NA or null
  check_if_any_NA(.data, c(
    quo_name(.sample),
    quo_name(.cell_group),
    parse_formula(formula)
  ))

  # Make counts
  .data %>%
    count(!!.sample,
          !!.cell_group,
          !!as.symbol(parse_formula(formula)),
          name = "count") %>%
    complete(
      nesting(!!.sample_for_tidyr,!!as.symbol(parse_formula(formula))),!!.cell_group_for_tidyr,
      fill = list(count = 0)
    ) %>%
    mutate(count = as.integer(count)) %>%

    # Return
    sccomp_glm_data_frame_counts(
      formula = formula,
      .sample = !!.sample,
      .cell_group = !!.cell_group,
      .count = count,
      my_glm_model = my_glm_model,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      cores = cores,
      verbose = verbose,
      seed = seed
    )
}

sccomp_glm_data_frame_counts = function(.data,
                                        formula,
                                        .sample,
                                        .cell_group,
                                        .count,
                                        my_glm_model,
                                        # Secondary arguments
                                        check_outliers = TRUE,
                                        approximate_posterior_inference = FALSE,
                                        verbose = FALSE,
                                        noise_model = "multi_beta_binomial",
                                        cores = 4,
                                        seed = sample(1:99999, size = 1)) {

  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)

  # Check if columns exist
  check_columns_exist(.data, c(
    quo_name(.sample),
    quo_name(.cell_group),
    quo_name(.count),
    parse_formula(formula)
  ))

  # Check if any column is NA or null
  check_if_any_NA(.data, c(
    quo_name(.sample),
    quo_name(.cell_group),
    quo_name(.count),
    parse_formula(formula)
  ))

  # Return
  .data %>%
    my_glm_model(
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
