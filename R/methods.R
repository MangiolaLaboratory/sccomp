

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
#' @param prior_mean_variable_association A list of the form list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)). Where for each parameter, we specify mean and standard deviation. This is used to incorporate prior knowledge about the mean/variability association of cell-type proportions.
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
#'     approximate_posterior_inference = "all",
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
                       formula_variability = ~ 1,
                       # Secondary arguments
                       prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                       percent_false_positive = 5,
                       check_outliers = TRUE,
                       approximate_posterior_inference = "outlier_detection",
                       test_composition_above_logit_fold_change = 0.2,
                       verbose = FALSE,
                       noise_model = "multi_beta_binomial",
                       variance_association = FALSE,
                       exclude_priors = FALSE,
                       cores = detectCores(),
                       mcmc_seed = sample(1e5, 1),
                       max_sampling_iterations = 20000,
                       pass_fit = TRUE) {
  UseMethod("sccomp_glm", .data)
}

#' @export
sccomp_glm.Seurat = function(.data,
                             formula = ~ 1 ,
                             .sample,
                             .cell_group,
                             .count = NULL,
                             formula_variability = ~ 1,
                             # Secondary arguments
                             prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                             percent_false_positive = 5,
                             check_outliers = TRUE,
                             approximate_posterior_inference = "outlier_detection",
                             test_composition_above_logit_fold_change = 0.2,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             variance_association = FALSE,
                             exclude_priors = FALSE,
                             cores = detectCores(),
                             mcmc_seed = sample(1e5, 1),
                             max_sampling_iterations = 20000,
                             pass_fit = TRUE) {

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
      exclude_priors = exclude_priors,
      cores = cores,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )


}

#' @export
sccomp_glm.SingleCellExperiment = function(.data,
                                           formula = ~ 1 ,
                                           .sample,
                                           .cell_group,
                                           .count = NULL,
                                           formula_variability = ~ 1,

                                           # Secondary arguments
                                           prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                                           percent_false_positive = 5,
                                           check_outliers = TRUE,
                                           approximate_posterior_inference = "outlier_detection",
                                           test_composition_above_logit_fold_change = 0.2,
                                           verbose = FALSE,
                                           noise_model = "multi_beta_binomial",
                                           variance_association = FALSE,
                                           exclude_priors = FALSE,
                                           cores = detectCores(),
                                           mcmc_seed = sample(1e5, 1),
                                           max_sampling_iterations = 20000,
                                           pass_fit = TRUE) {

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
      exclude_priors = exclude_priors,
      cores = cores,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )


}

#' @export
sccomp_glm.DFrame = function(.data,
                             formula ,
                             .sample,
                             .cell_group,
                             .count = NULL,
                             formula_variability = ~ 1,

                             # Secondary arguments
                             prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                             percent_false_positive = 5,
                             check_outliers = TRUE,
                             approximate_posterior_inference = "outlier_detection",
                             test_composition_above_logit_fold_change = 0.2,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             variance_association = FALSE,
                             exclude_priors = FALSE,
                             cores = detectCores(),
                             mcmc_seed = sample(1e5, 1),
                             max_sampling_iterations = 20000,
                             pass_fit = TRUE) {

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
      exclude_priors = exclude_priors,
      cores = cores,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )
}

#' @importFrom purrr when
#' @export
sccomp_glm.data.frame = function(.data,
                                 formula = ~ 1 ,
                                 .sample,
                                 .cell_group,
                                 .count = NULL,
                                 formula_variability = ~ 1,

                                 # Secondary arguments
                                 prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                                 percent_false_positive =  5,
                                 check_outliers = TRUE,
                                 approximate_posterior_inference = "outlier_detection",
                                 test_composition_above_logit_fold_change = 0.2,
                                 verbose = FALSE,
                                 noise_model = "multi_beta_binomial",
                                 variance_association = FALSE,
                                 exclude_priors = FALSE,
                                 cores = detectCores(),
                                 mcmc_seed = sample(1e5, 1),
                                 max_sampling_iterations = 20000,
                                 pass_fit = TRUE) {

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
        formula_variability = formula_variability,
        prior_mean_variable_association = prior_mean_variable_association,
        percent_false_positive = percent_false_positive ,
        check_outliers = check_outliers,
        approximate_posterior_inference = approximate_posterior_inference,
        test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
        verbose = verbose,
        my_glm_model = my_glm_model,
        variance_association = variance_association,
        exclude_priors = exclude_priors,
        cores = cores,
        mcmc_seed = mcmc_seed,
        max_sampling_iterations = max_sampling_iterations,
        pass_fit = pass_fit
      ),

      # If the dataframe does includes counts
      ~ sccomp_glm_data_frame_counts(
        .data,
        formula = formula,
        !!.sample,
        !!.cell_group,
        !!.count,
        formula_variability = formula_variability,
        prior_mean_variable_association = prior_mean_variable_association,
        percent_false_positive = percent_false_positive ,
        check_outliers = check_outliers,
        approximate_posterior_inference = approximate_posterior_inference,
        test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
        verbose = verbose,
        my_glm_model = my_glm_model,
        variance_association = variance_association,
        exclude_priors = exclude_priors,
        cores = cores,
        mcmc_seed = mcmc_seed,
        max_sampling_iterations = max_sampling_iterations,
        pass_fit = pass_fit
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
                                     formula_variability = ~ 1,
                                     my_glm_model,

                                     # Secondary arguments
                                     prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                                     percent_false_positive =  5,
                                     check_outliers = TRUE,
                                     approximate_posterior_inference = "outlier_detection",
                                     test_composition_above_logit_fold_change = 0.2,
                                     verbose = FALSE,
                                     noise_model = "multi_beta_binomial",
                                     variance_association = FALSE,
                                     exclude_priors = FALSE,
                                     cores = 4,
                                     mcmc_seed = sample(1e5, 1),
                                     max_sampling_iterations = 20000,
                                     pass_fit = TRUE ) {

  # See https://community.rstudio.com/t/how-to-make-complete-nesting-work-with-quosures-and-tidyeval/16473
  # See https://github.com/tidyverse/tidyr/issues/506


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
    class_list_to_counts(!!.sample, !!.cell_group) %>%

    # Add formula information
    when(
      length(parse_formula(formula))>0 ~ left_join(.,
                                                   .data %>%
                                                     select(!!.sample, parse_formula(formula) ) %>%
                                                     distinct(),
                                                   by = quo_name(.sample)
                                                  ),
      ~ (.)
    ) %>%

    # Return
    sccomp_glm_data_frame_counts(
      formula = formula,
      .sample = !!.sample,
      .cell_group = !!.cell_group,
      .count = count,
      formula_variability = formula_variability,
      my_glm_model = my_glm_model,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive =  percent_false_positive,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      variance_association = variance_association,
      exclude_priors = exclude_priors,
      cores = cores,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )
}

sccomp_glm_data_frame_counts = function(.data,
                                        formula,
                                        .sample,
                                        .cell_group,
                                        .count,
                                        formula_variability = ~ 1,
                                        my_glm_model,

                                        # Secondary arguments
                                        prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                                        percent_false_positive = 5,
                                        check_outliers = TRUE,
                                        approximate_posterior_inference = "outlier_detection",
                                        test_composition_above_logit_fold_change = 0.2,
                                        verbose = FALSE,
                                        noise_model = "multi_beta_binomial",
                                        variance_association = FALSE,
                                        exclude_priors = FALSE,
                                        cores = 4,
                                        mcmc_seed = sample(1e5, 1),
                                        max_sampling_iterations = 20000,
                                        pass_fit = TRUE) {

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
      formula_variability = formula_variability,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      variance_association = variance_association,
      exclude_priors = exclude_priors,
      cores = cores,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    ) %>%
    add_attr(.sample, ".sample") %>%
    add_attr(.cell_group, ".cell_group")
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
#'     approximate_posterior_inference = "all",
#'     check_outliers = FALSE,
#'     cores = 1
#'   ) |>
#'
#'   replicate_data()
#'
replicate_data <- function(.data,
                           number_of_draws = 1,
                           mcmc_seed = sample(1e5, 1)) {
  UseMethod("replicate_data", .data)
}

#' @export
#'
replicate_data.data.frame = function(.data,
                                     number_of_draws = 1,
                                     mcmc_seed = sample(1e5, 1)){


  # Select model based on noise model
  my_model = attr(.data, "noise_model") %>% when(
    (.) == "multi_beta_binomial" ~ stanmodels$glm_multi_beta_binomial_generate_date,
    (.) == "dirichlet_multinomial" ~ get_model_from_data("model_glm_dirichlet_multinomial_generate_quantities.rds", glm_dirichlet_multinomial_generate_quantities)
  )

  model_input = attr(.data, "model_input")
  .sample = attr(.data, ".sample")
  .cell_group = attr(.data, ".cell_group")

  fit_matrix = as.matrix(attr(.data, "fit") )

  # Generate quantities
  rstan::gqs(
    my_model,
    draws =  fit_matrix[sample(seq_len(nrow(fit_matrix)), size=number_of_draws),, drop=FALSE],
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
#' @param .estimate_object The result of sccomp_glm execution. This is used for sampling from real-data properties.
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
#'     approximate_posterior_inference = "all",
#'     check_outliers = FALSE,
#'     cores = 1
#'   )
#'
#' # Set coefficients for cell_types. In this case all coefficients are 0 for simplicity.
#' counts_obj = counts_obj |> mutate(b_0 = 0, b_1 = 0)
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
                       mcmc_seed = sample(1e5, 1)) {
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
                                    mcmc_seed = sample(1e5, 1)){


  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .coefficients = enquo(.coefficients)

  #Check column class
  check_if_columns_right_class(.data, !!.sample, !!.cell_group)

  model_data = attr(.estimate_object, "model_input")

  # Select model based on noise model
  my_model = attr(.estimate_object, "noise_model") %>% when(
    (.) == "multi_beta_binomial" ~ stanmodels$glm_multi_beta_binomial_simulate_data,
    (.) == "dirichlet_multinomial" ~ get_model_from_data("model_glm_dirichlet_multinomial_generate_quantities.rds", glm_dirichlet_multinomial_generate_quantities),
    (.) == "logit_normal_multinomial" ~ get_model_from_data("glm_multinomial_logit_linear_simulate_data.stan", readr::read_file("dev/stan_models/glm_multinomial_logit_linear_simulate_data.stan"))

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
    fit %>%
    parse_generated_quantities(number_of_draws = number_of_draws) %>%

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

simulate_multinomial_logit_linear = function(model_input, sd = 0.51){

  mu = model_input$X %*% model_input$beta

  proportions =
    rnorm(length(mu), mu, sd) %>%
    matrix(nrow = nrow(model_input$X)) %>%
    boot::inv.logit()
    apply(1, function(x) x/sum(x)) %>%
    t()

  rownames(proportions) = rownames(model_input$X)
  colnames(proportions) = colnames(model_input$beta )
}


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
#' @param prior_mean_variable_association A list of the form list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)). Where for each parameter, we specify mean and standard deviation. This is used to incorporate prior knowledge about the mean/variability association of cell-type proportions.
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
#'     approximate_posterior_inference = "all",
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
                       formula_variability = ~ 1,
                       # Secondary arguments
                       prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                       percent_false_positive = 5,
                       check_outliers = TRUE,
                       approximate_posterior_inference = "outlier_detection",
                       test_composition_above_logit_fold_change = 0.2,
                       verbose = FALSE,
                       noise_model = "multi_beta_binomial",
                       variance_association = FALSE,
                       exclude_priors = FALSE,
                       cores = detectCores(),
                       mcmc_seed = sample(1e5, 1),
                       max_sampling_iterations = 20000,
                       pass_fit = TRUE) {
  UseMethod("sccomp_glm", .data)
}

#' @export
sccomp_glm.Seurat = function(.data,
                             formula = ~ 1 ,
                             .sample,
                             .cell_group,
                             .count = NULL,
                             formula_variability = ~ 1,
                             # Secondary arguments
                             prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                             percent_false_positive = 5,
                             check_outliers = TRUE,
                             approximate_posterior_inference = "outlier_detection",
                             test_composition_above_logit_fold_change = 0.2,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             variance_association = FALSE,
                             exclude_priors = FALSE,
                             cores = detectCores(),
                             mcmc_seed = sample(1e5, 1),
                             max_sampling_iterations = 20000,
                             pass_fit = TRUE) {

  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")

  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)

  .data[[]] %>%
    sccomp_glm(
      formula = formula,!!.sample,!!.cell_group,
      formula_variability = formula_variability,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      noise_model = noise_model,
      variance_association = variance_association,
      exclude_priors = exclude_priors,
      cores = cores,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )


}

#' @export
sccomp_glm.SingleCellExperiment = function(.data,
                                           formula = ~ 1 ,
                                           .sample,
                                           .cell_group,
                                           .count = NULL,
                                           formula_variability = ~ 1,

                                           # Secondary arguments
                                           prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                                           percent_false_positive = 5,
                                           check_outliers = TRUE,
                                           approximate_posterior_inference = "outlier_detection",
                                           test_composition_above_logit_fold_change = 0.2,
                                           verbose = FALSE,
                                           noise_model = "multi_beta_binomial",
                                           variance_association = FALSE,
                                           exclude_priors = FALSE,
                                           cores = detectCores(),
                                           mcmc_seed = sample(1e5, 1),
                                           max_sampling_iterations = 20000,
                                           pass_fit = TRUE) {

  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")


  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)

  .data %>%
    colData() %>%
    sccomp_glm(
      formula = formula,!!.sample,!!.cell_group,
      formula_variability = formula_variability,
      check_outliers = check_outliers,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      approximate_posterior_inference = approximate_posterior_inference,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      noise_model = noise_model,
      variance_association = variance_association,
      exclude_priors = exclude_priors,
      cores = cores,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )


}

#' @export
sccomp_glm.DFrame = function(.data,
                             formula ,
                             .sample,
                             .cell_group,
                             .count = NULL,
                             formula_variability = ~ 1,

                             # Secondary arguments
                             prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                             percent_false_positive = 5,
                             check_outliers = TRUE,
                             approximate_posterior_inference = "outlier_detection",
                             test_composition_above_logit_fold_change = 0.2,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             variance_association = FALSE,
                             exclude_priors = FALSE,
                             cores = detectCores(),
                             mcmc_seed = sample(1e5, 1),
                             max_sampling_iterations = 20000,
                             pass_fit = TRUE) {

  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")


  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)

  .data %>%
    as.data.frame %>%
    sccomp_glm(
      formula = formula,!!.sample,!!.cell_group,
      formula_variability = formula_variability,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      noise_model = noise_model,
      variance_association = variance_association,
      exclude_priors = exclude_priors,
      cores = cores,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )
}

#' @importFrom purrr when
#' @export
sccomp_glm.data.frame = function(.data,
                                 formula = ~ 1 ,
                                 .sample,
                                 .cell_group,
                                 .count = NULL,
                                 formula_variability = ~ 1,

                                 # Secondary arguments
                                 prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                                 percent_false_positive =  5,
                                 check_outliers = TRUE,
                                 approximate_posterior_inference = "outlier_detection",
                                 test_composition_above_logit_fold_change = 0.2,
                                 verbose = FALSE,
                                 noise_model = "multi_beta_binomial",
                                 variance_association = FALSE,
                                 exclude_priors = FALSE,
                                 cores = detectCores(),
                                 mcmc_seed = sample(1e5, 1),
                                 max_sampling_iterations = 20000,
                                 pass_fit = TRUE) {

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
        formula_variability = formula_variability,
        prior_mean_variable_association = prior_mean_variable_association,
        percent_false_positive = percent_false_positive ,
        check_outliers = check_outliers,
        approximate_posterior_inference = approximate_posterior_inference,
        test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
        verbose = verbose,
        my_glm_model = my_glm_model,
        variance_association = variance_association,
        exclude_priors = exclude_priors,
        cores = cores,
        mcmc_seed = mcmc_seed,
        max_sampling_iterations = max_sampling_iterations,
        pass_fit = pass_fit
      ),

      # If the dataframe does includes counts
      ~ sccomp_glm_data_frame_counts(
        .data,
        formula = formula,
        !!.sample,
        !!.cell_group,
        !!.count,
        formula_variability = formula_variability,
        prior_mean_variable_association = prior_mean_variable_association,
        percent_false_positive = percent_false_positive ,
        check_outliers = check_outliers,
        approximate_posterior_inference = approximate_posterior_inference,
        test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
        verbose = verbose,
        my_glm_model = my_glm_model,
        variance_association = variance_association,
        exclude_priors = exclude_priors,
        cores = cores,
        mcmc_seed = mcmc_seed,
        max_sampling_iterations = max_sampling_iterations,
        pass_fit = pass_fit
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
                                     formula_variability = ~ 1,
                                     my_glm_model,

                                     # Secondary arguments
                                     prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                                     percent_false_positive =  5,
                                     check_outliers = TRUE,
                                     approximate_posterior_inference = "outlier_detection",
                                     test_composition_above_logit_fold_change = 0.2,
                                     verbose = FALSE,
                                     noise_model = "multi_beta_binomial",
                                     variance_association = FALSE,
                                     exclude_priors = FALSE,
                                     cores = 4,
                                     mcmc_seed = sample(1e5, 1),
                                     max_sampling_iterations = 20000,
                                     pass_fit = TRUE ) {

  # See https://community.rstudio.com/t/how-to-make-complete-nesting-work-with-quosures-and-tidyeval/16473
  # See https://github.com/tidyverse/tidyr/issues/506


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
    class_list_to_counts(!!.sample, !!.cell_group) %>%

    # Add formula information
    when(
      length(parse_formula(formula))>0 ~ left_join(.,
                                                   .data %>%
                                                     select(!!.sample, parse_formula(formula) ) %>%
                                                     distinct(),
                                                   by = quo_name(.sample)
                                                  ),
      ~ (.)
    ) %>%

    # Return
    sccomp_glm_data_frame_counts(
      formula = formula,
      .sample = !!.sample,
      .cell_group = !!.cell_group,
      .count = count,
      formula_variability = formula_variability,
      my_glm_model = my_glm_model,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive =  percent_false_positive,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      variance_association = variance_association,
      exclude_priors = exclude_priors,
      cores = cores,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )
}

sccomp_glm_data_frame_counts = function(.data,
                                        formula,
                                        .sample,
                                        .cell_group,
                                        .count,
                                        formula_variability = ~ 1,
                                        my_glm_model,

                                        # Secondary arguments
                                        prior_mean_variable_association = list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(5.06983, 8.549324)),
                                        percent_false_positive = 5,
                                        check_outliers = TRUE,
                                        approximate_posterior_inference = "outlier_detection",
                                        test_composition_above_logit_fold_change = 0.2,
                                        verbose = FALSE,
                                        noise_model = "multi_beta_binomial",
                                        variance_association = FALSE,
                                        exclude_priors = FALSE,
                                        cores = 4,
                                        mcmc_seed = sample(1e5, 1),
                                        max_sampling_iterations = 20000,
                                        pass_fit = TRUE) {

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
      formula_variability = formula_variability,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      variance_association = variance_association,
      exclude_priors = exclude_priors,
      cores = cores,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    ) %>%
    add_attr(.sample, ".sample") %>%
    add_attr(.cell_group, ".cell_group") %>%
    add_attr(parse_formula(formula), "covariates" )
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
#'     approximate_posterior_inference = "all",
#'     check_outliers = FALSE,
#'     cores = 1
#'   ) |>
#'
#'   replicate_data()
#'
replicate_data <- function(.data,
                           number_of_draws = 1,
                           mcmc_seed = sample(1e5, 1)) {
  UseMethod("replicate_data", .data)
}

#' @export
#'
replicate_data.data.frame = function(.data,
                                     number_of_draws = 1,
                                     mcmc_seed = sample(1e5, 1)){


  # Select model based on noise model
  my_model = attr(.data, "noise_model") %>% when(
    (.) == "multi_beta_binomial" ~ stanmodels$glm_multi_beta_binomial_generate_date,
    (.) == "dirichlet_multinomial" ~ get_model_from_data("model_glm_dirichlet_multinomial_generate_quantities.rds", glm_dirichlet_multinomial_generate_quantities)
  )

  model_input = attr(.data, "model_input")
  .sample = attr(.data, ".sample")
  .cell_group = attr(.data, ".cell_group")

  fit_matrix = as.matrix(attr(.data, "fit") )

  # Generate quantities
  rstan::gqs(
    my_model,
    draws =  fit_matrix[sample(seq_len(nrow(fit_matrix)), size=number_of_draws),, drop=FALSE],
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
#' @param .estimate_object The result of sccomp_glm execution. This is used for sampling from real-data properties.
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
#'     approximate_posterior_inference = "all",
#'     check_outliers = FALSE,
#'     cores = 1
#'   )
#'
#' # Set coefficients for cell_types. In this case all coefficients are 0 for simplicity.
#' counts_obj = counts_obj |> mutate(b_0 = 0, b_1 = 0)
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
                       variability_multiplier = 5,
                       number_of_draws = 1,
                       mcmc_seed = sample(1e5, 1)) {
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
                                    variability_multiplier = 5,
                                    number_of_draws = 1,
                                    mcmc_seed = sample(1e5, 1)){


  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .coefficients = enquo(.coefficients)

  #Check column class
  check_if_columns_right_class(.data, !!.sample, !!.cell_group)

  model_data = attr(.estimate_object, "model_input")

  # Select model based on noise model
  my_model = attr(.estimate_object, "noise_model") %>% when(
    (.) == "multi_beta_binomial" ~ stanmodels$glm_multi_beta_binomial_simulate_data,
    (.) == "dirichlet_multinomial" ~ get_model_from_data("model_glm_dirichlet_multinomial_generate_quantities.rds", glm_dirichlet_multinomial_generate_quantities),
    (.) == "logit_normal_multinomial" ~ get_model_from_data("glm_multinomial_logit_linear_simulate_data.stan", readr::read_file("dev/stan_models/glm_multinomial_logit_linear_simulate_data.stan"))

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
    data = model_input %>% c(list(variability_multiplier = variability_multiplier)),
    seed = mcmc_seed
  )

  parsed_fit =
    fit %>%
    parse_generated_quantities(number_of_draws = number_of_draws) %>%

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

simulate_multinomial_logit_linear = function(model_input, sd = 0.51){

  mu = model_input$X %*% model_input$beta

  proportions =
    rnorm(length(mu), mu, sd) %>%
    matrix(nrow = nrow(model_input$X)) %>%
    boot::inv.logit()
    apply(1, function(x) x/sum(x)) %>%
    t()

  rownames(proportions) = rownames(model_input$X)
  colnames(proportions) = colnames(model_input$beta )
}


#' plot_summary
#'
#' @description This function plots a summary of the results of the model.
#'
#'
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | covariate columns | Pvalue column | a significance column
#' @return A `ggplot`
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
#'     approximate_posterior_inference = "all",
#'     check_outliers = FALSE,
#'     cores = 1
#'   )
#'
#' estimate |> plot_summary()
#'
plot_summary <- function(.data, .cell_group) {

    multipanel_theme =
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(size=0.5),
      panel.grid.major = element_line(size = 0.1),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.background = element_blank(),
      axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 7),
      axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 7),
      panel.spacing.x=unit(0.1, "lines"),
      axis.text.x = element_text(size=6),
      axis.text.y = element_text(size=6),
      strip.text.x = element_text(size = 7),
      strip.text.y = element_text(size = 7),

      # legend
      legend.key.size = unit(5, 'mm'),
      legend.title = element_text(size=7),
      legend.text = element_text(size=6),

      # Avoid text clipping for facets. Currently not merged remotes::install_github("tidyverse/ggplot2#4223")
      # strip.clip = "off",

      # Title
      plot.title = element_text(size=7),

      axis.line.x = element_line(size=0.2),
      axis.line.y = element_line(size=0.2),
      axis.ticks.x = element_line(size=0.2),
      axis.ticks.y = element_line(size=0.2)
    )

  dropLeadingZero <- function(l){  stringr::str_replace(l, '0(?=.)', '') }

  S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
  IS_sqrt <- function(x){x^2*sign(x)}
  S_sqrt_trans <- function() scales::trans_new("S_sqrt",S_sqrt,IS_sqrt)


  .cell_group = enquo(.cell_group)

  if("v_effect" %in% colnames(.data)){
  # mean-variance association
plot_associations =
  .data %>%

  # Filter where I did not inferred the variance
  filter(!is.na(v_effect)) %>%

  # Plot
  ggplot(aes(c_effect, v_effect, label=!!.cell_group)) +
  geom_vline(xintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
  geom_hline(yintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
  geom_errorbar(aes(xmin=`c_lower`, xmax=`c_upper`, color=`c_FDR`<0.025, alpha=`c_FDR`<0.025), size=0.2) +
  geom_errorbar(aes(ymin=`v_lower`, ymax=`v_upper`, color=`v_FDR`<0.025, alpha=`v_FDR`<0.025), size=0.2) +

  geom_point(size=0.2)  +
  annotate("text", x = 0, y = -3.5, label = "Variable", size=2) +
  annotate("text", x = 5, y = 0, label = "Abundant", size=2, angle=270) +
  scale_color_manual(values = c("#D3D3D3", "#E41A1C")) +
  scale_alpha_manual(values = c(0.4, 1)) +
  facet_wrap(~parameter) +
  multipanel_theme +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  }

if("fit" %in% names(attributes(.data))){

  calc_boxplot_stat <- function(x) {
    coef <- 1.5
    n <- sum(!is.na(x))
    # calculate quantiles
    stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
    names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
    iqr <- diff(stats[c(2, 4)])
    # set whiskers
    outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
    if (any(outliers)) {
      stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
    }
    return(stats)
  }


  data_proportion =
    .data %>%
    pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) %>%
    unnest(count_data) %>%
    with_groups(sample, ~ mutate(.x, proportion = (count)/sum(count)) )

  factor_of_interest = .data %>% attr("covariates") %>% .[1]
  simulated_proportion =
    .data %>%
    replicate_data( number_of_draws = 100) %>%
    left_join(data_proportion %>% distinct(!!as.symbol(factor_of_interest), sample, !!.cell_group))


  ggplot() +

    stat_summary(
      aes(!!as.symbol(factor_of_interest), (generated_proportions)),
      fun.data = calc_boxplot_stat, geom="boxplot",
      fatten = 0.5, lwd=0.2,
      data =
        simulated_proportion %>%

        # Filter uanitles because of limits
        inner_join( data_proportion %>% distinct(!!as.symbol(factor_of_interest), !!.cell_group)) ,
      color="blue"

    ) +

    geom_boxplot(
      aes(!!as.symbol(factor_of_interest), proportion,  group=!!as.symbol(factor_of_interest)), # fill=Effect),
      outlier.shape = NA,
      data = data_proportion |> filter(!outlier), fatten = 0.5, lwd=0.5,
    ) +
    geom_jitter(
      aes(!!as.symbol(factor_of_interest), proportion, shape=outlier,  group=!!as.symbol(factor_of_interest)),
      data = data_proportion,
      position=position_jitterdodge(jitter.height = 0, jitter.width = 0.2),
      size = 0.5
    ) +

    # geom_boxplot(
    #   aes(Condition, generated_proportions),
    #   outlier.shape = NA, alpha=0.2,
    #   data = simulated_proportion, fatten = 0.5, size=0.5,
    # ) +
    # geom_jitter(aes(Condition, generated_proportions), color="black" ,alpha=0.2, size = 0.2, data = simulated_proportion) +

    facet_wrap(
      vars(!!.cell_group) ,# forcats::fct_reorder(!!.cell_group, abs(Effect), .desc = TRUE, na.rm=TRUE),
      scales = "free_y", nrow = 4
    ) +
    #scale_color_manual(values = c("black", "#e11f28")) +
    #scale_fill_manual(values = c("white", "#E2D379")) +
    scale_fill_distiller(palette = "Spectral", na.value = "white") +
    #scale_color_distiller(palette = "Spectral") +

    scale_y_continuous(trans="S_sqrt", labels = dropLeadingZero) +
    #scale_y_continuous(labels = dropLeadingZero, trans="logit") +
    xlab("Biological condition") +
    ylab("Cell-group proportion") +
    guides(color="none", alpha="none", size="none") +
    labs(fill="Compositional difference") +
    multipanel_theme +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x =  element_text(angle=20, hjust = 1))

}

}
