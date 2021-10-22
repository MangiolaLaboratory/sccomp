

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
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' @param verbose A boolean. Prints progression.
#' @param noise_model A character string. The two noise models available are multi_beta_binomial (default) and dirichlet_multinomial.
#' @param variance_association A boolean. Whether the variance should depend on the factor of interest.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
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
#'     approximate_posterior_inference = TRUE,
#'     check_outliers = FALSE,
#'     cores = 1
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
                       test_composition_above_logit_fold_change = 0.2,
                       verbose = FALSE,
                       noise_model = "multi_beta_binomial",
                       variance_association = FALSE,
                       cores = detectCores(),
                       mcmc_seed = seq_len(99999) |> sample(size = 1)) {
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
                             test_composition_above_logit_fold_change = 0.2,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             variance_association = FALSE,
                             cores = detectCores(),
                             mcmc_seed = seq_len(99999) |> sample(size = 1)) {

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
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      noise_model = noise_model,
      variance_association = variance_association,
      cores = cores,
      mcmc_seed = mcmc_seed
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
                                           test_composition_above_logit_fold_change = 0.2,
                                           verbose = FALSE,
                                           noise_model = "multi_beta_binomial",
                                           variance_association = FALSE,
                                           cores = detectCores(),
                                           mcmc_seed = seq_len(99999) |> sample(size = 1)) {

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
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      noise_model = noise_model,
      variance_association = variance_association,
      cores = cores,
      mcmc_seed = mcmc_seed
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
                             test_composition_above_logit_fold_change = 0.2,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             variance_association = FALSE,
                             cores = detectCores(),
                             mcmc_seed = seq_len(99999) |> sample(size = 1)) {

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
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      noise_model = noise_model,
      variance_association = variance_association,
      cores = cores,
      mcmc_seed = mcmc_seed
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
                                 test_composition_above_logit_fold_change = 0.2,
                                 verbose = FALSE,
                                 noise_model = "multi_beta_binomial",
                                 variance_association = FALSE,
                                 cores = detectCores(),
                                 mcmc_seed = seq_len(99999) |> sample(size = 1)) {

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
        test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
        verbose = verbose,
        my_glm_model = my_glm_model,
        variance_association = variance_association,
        cores = cores,
        mcmc_seed = mcmc_seed
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
        test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
        verbose = verbose,
        my_glm_model = my_glm_model,
        variance_association = variance_association,
        cores = cores,
        mcmc_seed = mcmc_seed
      )
    ) %>%

    # Track input parameters
    add_attr(noise_model, "noise_model") %>%
    add_attr(.sample, ".sample") %>%
    add_attr(.cell_group, ".cell_group")
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
                                     test_composition_above_logit_fold_change = 0.2,
                                     verbose = FALSE,
                                     noise_model = "multi_beta_binomial",
                                     variance_association = FALSE,
                                     cores = 4,
                                     mcmc_seed = seq_len(99999) |> sample(size = 1) ) {

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
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      mcmc_seed = mcmc_seed
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
                                        test_composition_above_logit_fold_change = 0.2,
                                        verbose = FALSE,
                                        noise_model = "multi_beta_binomial",
                                        variance_association = FALSE,
                                        cores = 4,
                                        mcmc_seed = seq_len(99999) |> sample(size = 1)) {

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
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      seed = mcmc_seed
    )
}

#' replicate_data
#'
#' @description This function replicates counts from a real-world dataset.
#'
#'
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | covariate columns | Pvalue column | a significance column
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#'
#' @return A nested tibble `tbl` with cell_group-wise statistics
#'
#' @export
#'
#' @examples
#'
#' data("counts_obj")
#'
#'   sccomp_glm(
#'   counts_obj ,
#'    ~ type,  sample, cell_group, count,
#'     approximate_posterior_inference = TRUE,
#'     check_outliers = FALSE,
#'     cores = 1
#'   ) |>
#'
#'   replicate_data()
#'
replicate_data <- function(.data,
                           number_of_draws = 1,
                           mcmc_seed = seq_len(99999) |> sample(size = 1)) {
  UseMethod("replicate_data", .data)
}

#' @export
#'
replicate_data.data.frame = function(.data,
                                     number_of_draws = 1,
                                     mcmc_seed = seq_len(99999) |> sample(size = 1)){


  # Select model based on noise model
  my_model = attr(.data, "noise_model") %>% when(
    (.) == "multi_beta_binomial" ~ stanmodels$glm_multi_beta_binomial_generate_date,
    (.) == "dirichlet_multinomial" ~ get_model_from_data("model_glm_dirichlet_multinomial_generate_quantities.rds", glm_dirichlet_multinomial_generate_quantities)
  )

  model_input = attr(.data, "model_input")
  .sample = attr(.data, ".sample")
  .cell_group = attr(.data, ".cell_group")

  # Generate quantities
  rstan::gqs(
    my_model,
    draws =  as.matrix(attr(.data, "fit") ),
    data = model_input,
    seed = mcmc_seed
  ) %>%

    # Parse
    parse_generated_quantities(number_of_draws = number_of_draws) %>%

    # Get sample name
    nest(data = -N) %>%
    arrange(N) %>%
    mutate(!!.sample := rownames(model_input$y)) %>%
    unnest(data) %>%

    # get cell type name
    nest(data = -M) %>%
    mutate(!!.cell_group := colnames(model_input$y)) %>%
    unnest(data) %>%

    select(-N, -M)

  # %>%
  #   nest(generated_data = -c(!!.sample, !!.cell_group))

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
#' @param .estimate_object The result of sccomp_glm execution. This is used for sampling from real-data properies.
#' @param formula A formula. The sample formula used to perform the differential cell_type abundance analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_type identifier
#' @param .coefficients Th column names for coefficients, for example, c(b_0, b_1)
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#'
#' @return A nested tibble `tbl` with cell_group-wise statistics
#'
#' @export
#'
#' @examples
#'
#' data("counts_obj")
#' library(dplyr)
#'
#' estimate =
#'   sccomp_glm(
#'   counts_obj ,
#'    ~ type,  sample, cell_group, count,
#'     approximate_posterior_inference = TRUE,
#'     check_outliers = FALSE,
#'     cores = 1
#'   )
#'
#' # Set coefficients for cell_types. In this case all coefficients are 0 for simplicity.
#' counts_obj = counts_obj |> mutate(b_0 = 1, b_1 = 1)
#'
#' # Simulate data
#' simulate_data(counts_obj, estimate, ~type, sample, cell_group, c(b_0, b_1))
#'
simulate_data <- function(.data,
                          .estimate_object,
                          formula,
                       .sample = NULL,
                       .cell_group = NULL,
                       .coefficients = NULL,
                       number_of_draws = 1,
                       mcmc_seed = seq_len(99999) |> sample(size = 1)) {
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
                                    .estimate_object,
                                    formula,
                                    .sample = NULL,
                                    .cell_group = NULL,
                                    .coefficients = NULL,
                                    number_of_draws = 1,
                                    mcmc_seed = seq_len(99999) |> sample(size = 1)){


  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .coefficients = enquo(.coefficients)

  #Check column class
  check_if_columns_right_class(.data, !!.sample, !!.cell_group)

  model_data = attr(.estimate_object, "model_input")

  # Select model based on noise model
  my_model = attr(.estimate_object, "noise_model") %>% when(
    (.) == "multi_beta_binomial" ~ stanmodels$glm_multi_beta_binomial_simulate_data,
    (.) == "dirichlet_multinomial" ~ get_model_from_data("model_glm_dirichlet_multinomial_generate_quantities.rds", glm_dirichlet_multinomial_generate_quantities)
  )


  model_input =
    .data %>%
    nest(data___ = -!!.sample) %>%
    mutate(.exposure = sample(model_data$exposure, size = n(), replace = TRUE )) %>%
    unnest(data___) %>%
    data_simulation_to_model_input(formula, !!.sample, !!.cell_group, .exposure, !!.coefficients )



    # [1]  5.6260004 -0.6940178
    # prec_sd  = 0.816423129

  fit =
    rstan::gqs(
    my_model,
    draws =  as.matrix(attr(.estimate_object, "fit") ),
    data = model_input,
    seed = mcmc_seed
  )

  parsed_fit =
    parse_generated_quantities(fit, number_of_draws = number_of_draws) %>%

    # Get sample name
    nest(data = -N) %>%
    arrange(N) %>%
    mutate(!!.sample := rownames(model_input$X)) %>%
    unnest(data) %>%

    # get cell type name
    nest(data = -M) %>%
    mutate(!!.cell_group := colnames(model_input$beta)) %>%
    unnest(data) %>%

    select(-N, -M)

  .data %>%
    left_join(
      parsed_fit,
      by = c(quo_name(.sample), quo_name(.cell_group))
    )


}


