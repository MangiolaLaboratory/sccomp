

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
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | covariate columns | Pvalue column | a significance column
#' @param formula A formula. The sample formula used to perform the differential cell_type abundance analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_type identifier
#' @param .count A column name as symbol. The cell_type abundance (read count). Used only for data frame count output.
#' @param prior_mean_variable_association A list of the form list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(0.4527292, 0.3318759)). Where for each parameter, we specify mean and standard deviation. This is used to incorporate prior knowledge about the mean/variability association of cell-type proportions.
#' @param percent_false_positive A real between 0 and 100. It is the aimed percent of cell types being a false positive. For example, percent_false_positive_genes = 1 provide 1 percent of the calls for significant changes that are actually not significant.
#' @param check_outliers A boolean. Whether to check for outliers before the fit.
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param verbose A boolean. Prints progression.
#' @param noise_model A character string. The two noise models available are multi_beta_binomial (default) and dirichlet_multinomial.
#' @param variance_association A boolean. Whether the variance should depend on the factor of interest.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param seed An integer. Used for development and testing purposes
#'
#' @return A nested tibble `tbl` with cell_group-wise statistics
#'
#' @examples
#'
#' data("counts_obj")
#'
#' estimate =
#'   sccomp_glm(
#'   counts_obj ,
#'    ~ type,  sample, cell_group, count,
#'     approximate_posterior_inference = FALSE
#'   )
#'
#' @export
#'
#'
sccomp_glm <- function(.data,
                       formula = ~ 1 ,
                       .sample,
                       .cell_group,
                       .count = NULL,
                       # Secondary arguments
                       prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(0.4527292, 0.3318759)),
                       percent_false_positive = 5,
                       check_outliers = TRUE,
                       approximate_posterior_inference = FALSE,
                       verbose = FALSE,
                       noise_model = "multi_beta_binomial",
                       variance_association = FALSE,
                       cores = detectCores(),
                       seed = 42) {
  UseMethod("sccomp_glm", .data)
}

#' @export
sccomp_glm.Seurat = function(.data,
                             formula = ~ 1 ,
                             .sample,
                             .cell_group,
                             .count = NULL,
                             # Secondary arguments
                             prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(0.4527292, 0.3318759)),
                             percent_false_positive = 5,
                             check_outliers = TRUE,
                             approximate_posterior_inference = FALSE,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             variance_association = FALSE,
                             cores = detectCores(),
                             seed = 42) {

  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")

  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)

  .data[[]] %>%
    sccomp_glm(
      formula = formula,!!.sample,!!.cell_group,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      verbose = verbose,
      noise_model = noise_model,
      variance_association = variance_association,
      cores = cores,
      seed = seed
    )


}

#' @export
sccomp_glm.SingleCellExperiment = function(.data,
                                           formula = ~ 1 ,
                                           .sample,
                                           .cell_group,
                                           .count = NULL,

                                           # Secondary arguments
                                           prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(0.4527292, 0.3318759)),
                                           percent_false_positive = 5,
                                           check_outliers = TRUE,
                                           approximate_posterior_inference = FALSE,
                                           verbose = FALSE,
                                           noise_model = "multi_beta_binomial",
                                           variance_association = FALSE,
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
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      approximate_posterior_inference = approximate_posterior_inference,
      verbose = verbose,
      noise_model = noise_model,
      variance_association = variance_association,
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
                             prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(0.4527292, 0.3318759)),
                             percent_false_positive = 5,
                             check_outliers = TRUE,
                             approximate_posterior_inference = FALSE,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             variance_association = FALSE,
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
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      verbose = verbose,
      noise_model = noise_model,
      variance_association = variance_association,
      cores = cores,
      seed = seed
    )
}

#' @importFrom purrr when
#' @export
sccomp_glm.data.frame = function(.data,
                                 formula = ~ 1 ,
                                 .sample,
                                 .cell_group,
                                 .count = NULL,

                                 # Secondary arguments
                                 prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(0.4527292, 0.3318759)),
                                 percent_false_positive =  5,
                                 check_outliers = TRUE,
                                 approximate_posterior_inference = FALSE,
                                 verbose = FALSE,
                                 noise_model = "multi_beta_binomial",
                                 variance_association = FALSE,
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
        prior_mean_variable_association = prior_mean_variable_association,
        percent_false_positive = percent_false_positive ,
        check_outliers = check_outliers,
        approximate_posterior_inference = approximate_posterior_inference,
        verbose = verbose,
        my_glm_model = my_glm_model,
        variance_association = variance_association,
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
        prior_mean_variable_association = prior_mean_variable_association,
        percent_false_positive = percent_false_positive ,
        check_outliers = check_outliers,
        approximate_posterior_inference = approximate_posterior_inference,
        verbose = verbose,
        my_glm_model = my_glm_model,
        variance_association = variance_association,
        cores = cores,
        seed = seed
      )
    )
}

#' @importFrom tidyr complete
#' @importFrom tidyr nesting
sccomp_glm_data_frame_raw = function(.data,
                                     formula = ~ 1,
                                     .sample,
                                     .cell_group,
                                     my_glm_model,

                                     # Secondary arguments
                                     prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(0.4527292, 0.3318759)),
                                     percent_false_positive =  5,
                                     check_outliers = TRUE,
                                     approximate_posterior_inference = FALSE,
                                     verbose = FALSE,
                                     noise_model = "multi_beta_binomial",
                                     variance_association = FALSE,
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
          name = "count") %>%

    complete(
      nesting(!!.sample_for_tidyr),!!.cell_group_for_tidyr,
      fill = list(count = 0)
    ) %>%
    mutate(count = as.integer(count)) %>%

    # Add formula information
    when(
      length(parse_formula(formula))>0 ~ left_join(., .data %>% distinct(!!.sample,!!as.symbol(parse_formula(formula)) )),
      ~ (.)
    ) %>%

    # Return
    sccomp_glm_data_frame_counts(
      formula = formula,
      .sample = !!.sample,
      .cell_group = !!.cell_group,
      .count = count,
      my_glm_model = my_glm_model,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive =  percent_false_positive,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      variance_association = variance_association,
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
                                        prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(0.4527292, 0.3318759)),
                                        percent_false_positive = 5,
                                        check_outliers = TRUE,
                                        approximate_posterior_inference = FALSE,
                                        verbose = FALSE,
                                        noise_model = "multi_beta_binomial",
                                        variance_association = FALSE,
                                        cores = 4,
                                        seed = sample(1:99999, size = 1)) {

  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)

  #Check column class
  check_if_columns_right_class(.data, !!.sample, !!.cell_group)

  # Check that count is integer
  if(.data %>% pull(!!.count) %>% is("integer") %>% not())
    stop(sprintf("sccomp: %s column must be an integer", quo_name(.count)))

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
      .count = !!.count,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      variance_association = variance_association,
      cores = cores,
      verbose = verbose,
      seed = seed
    )
}


#' simulate_data
#'
#' @description This function simulates counts from a linear model.
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
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | covariate columns | Pvalue column | a significance column
#' @param formula A formula. The sample formula used to perform the differential cell_type abundance analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_type identifier
#' @param .sample_cell_count A integer vector. The total number of cells for each sample.
#' @param .coefficients A matrix of coefficients.
#' @param mean_variable_association A numeric vector of size 3. The intercept, slope and standard deviation of the proportion mean/variability association.
#' @param percent_false_positive A real between 0 and 100. It is the aimed percent of cell types being a false positive. For example, percent_false_positive_genes = 1 provide 1 percent of the calls for significant changes that are actually not significant.
#' @param check_outliers A boolean. Whether to check for outliers before the fit.
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param verbose A boolean. Prints progression.
#' @param noise_model A character string. The two noise models available are multi_beta_binomial (default) and dirichlet_multinomial.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param seed An integer. Used for development and testing purposes
#'
#' @return A nested tibble `tbl` with cell_group-wise statistics
#'
#' @export
#'
#' @examples
#'
#' data("counts_obj")
#'
#' estimate =
#'   sccomp_glm(
#'   counts_obj ,
#'    ~ type,  sample, cell_group, count,
#'     approximate_posterior_inference = FALSE
#'   )
#'
simulate_data <- function(.data,
                       formula = ~ 1 ,
                       .sample,
                       .cell_group,
                       .sample_cell_count,
                       .coefficients,
                       # Secondary arguments
                       mean_variable_association = c( 5.6260004, -0.6940178, 0.816423129),
                       percent_false_positive = 5,
                       check_outliers = TRUE,
                       approximate_posterior_inference = FALSE,
                       verbose = FALSE,
                       noise_model = "multi_beta_binomial",
                       cores = detectCores(),
                       seed = 42) {
  UseMethod("simulate_data", .data)
}

#' @export
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom SingleCellExperiment counts
#' @importFrom purrr map_dbl
#' @importFrom purrr pmap
#'
simulate_data.data.frame = function(.data,
                                    formula = ~ 1 ,
                                    .sample,
                                    .cell_group,
                                    .sample_cell_count,
                                    .coefficients,
                                    # Secondary arguments
                                    mean_variable_association = c( 5.6260004, -0.6940178, 0.816423129),
                                    percent_false_positive = 5,
                                    check_outliers = TRUE,
                                    approximate_posterior_inference = FALSE,
                                    verbose = FALSE,
                                    noise_model = "multi_beta_binomial",
                                    cores = detectCores(),
                                    seed = 42){


  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .sample_cell_count = enquo(.sample_cell_count)
  .coefficients = enquo(.coefficients)

  #Check column class
  check_if_columns_right_class(.data, !!.sample, !!.cell_group)

  model_data =
    .data %>%
    data_simulation_to_model_input(formula, !!.sample, !!.cell_group, !!.sample_cell_count, !!.coefficients )

  model_data$prec_coeff = mean_variable_association[1:2]
  model_data$prec_sd  = mean_variable_association[3]

    # [1]  5.6260004 -0.6940178
    # prec_sd  = 0.816423129

  fit =
    sampling(
    stanmodels$glm_multi_beta_binomial_simulate_data,
    #stan_model("inst/stan/glm_multi_beta_binomial_simulate_data.stan"),
    data = model_data,
    chains = 1,
    cores = 1,
    iter = 151,
    warmup = 150,
    #refresh = ifelse(verbose, 1000, 0),
    seed = seed,
    #pars = pars,
    save_warmup = FALSE
  )

  .data %>%
    arrange(!!.sample, !!.cell_group) %>%
    bind_cols(
      fit %>%
        draws_to_tibble_x_y("counts", "N", "M") %>%
        ungroup() %>%
        arrange(N, M) %>%
        select(.value)
    ) %>%

    # Scale counts for exposure
    nest(data = -c(sample, tot_count )) %>%
    mutate(tot_count_simulated = map_dbl(data, ~ sum(.x$.value))) %>%
    mutate(data = pmap(list(data, tot_count, tot_count_simulated),  ~ ..1 %>%  mutate(.value = .value %>% divide_by(..3) %>% multiply_by(..2) %>% round %>% as.integer ))) %>%
    unnest(data)


}
