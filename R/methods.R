

#' sccomp_glm main
#'
#' @description The function for linear modelling takes as input a table of cell counts with three columns containing a cell-group identifier, sample identifier, integer count and the factors (continuous or discrete). The user can define a linear model with an input R formula, where the first factor is the factor of interest. Alternatively, sccomp accepts single-cell data containers (Seurat, SingleCellExperiment44, cell metadata or group-size). In this case, sccomp derives the count data from cell metadata.
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
#' @param .data A tibble including a cell_group name column | sample name column | read counts column (optional depending on the input class) | factor columns.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment.
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment. In most cases, if differentially variability is of interest, the formula should only include the factor of interest as a large anount of data is needed to define variability depending to each factors.
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_group identifier
#' @param .count A column name as symbol. The cell_group abundance (read count). Used only for data frame count output. The variable in this column should be of class integer.
#'
#' @param contrasts A vector of character strings. For example if your formula is `~ 0 + treatment` and the factor treatment has values `yes` and `no`, your contrast could be constrasts = c("treatmentyes - treatmentno").
#' @param prior_mean_variable_association A list of the form list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)). Where for intercept and slope parameters, we specify mean and standard deviation, while for standard deviation, we specify shape and rate. This is used to incorporate prior knowledge about the mean/variability association of cell-type proportions.
#' @param check_outliers A boolean. Whether to check for outliers before the fit.
#' @param bimodal_mean_variability_association A boolean. Whether to model the mean-variability as bimodal, as often needed in the case of single-cell RNA sequencing data, and not usually for CyTOF and microbiome data. The plot summary_plot()$credible_intervals_2D can be used to assess whether the bimodality should be modelled.
#' @param enable_loo A boolean. Enable model comparison by the R package LOO. This is helpful when you want to compare the fit between two models, for example, analogously to ANOVA, between a one factor model versus a interceot-only model.
#'
#' @param percent_false_positive A real between 0 and 100 non included. This used to identify outliers with a specific false positive rate.
#' @param exclude_priors A boolean. Whether to run a prior-free model, for benchmarking purposes.
#' @param use_data A booelan. Whether to sun the model data free. This can be used for prior predictive check.
#' @param max_sampling_iterations An integer. This limit the maximum number of iterations in case a large dataset is used, for limiting the computation time.
#' @param pass_fit A boolean. Whether to pass the Stan fit as attribute in the output. Because the Stan fit can be very large, setting this to FALSE can be used to lower the memory imprint to save the output.
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' @param verbose A boolean. Prints progression.
#' @param noise_model A character string. The two noise models available are multi_beta_binomial (default) and dirichlet_multinomial.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#'
#' @return A nested tibble `tbl`, with the following columns
#' \itemize{
#'   \item cell_group - column including the cell groups being tested
#'   \item parameter - The parameter being estimated, from the design matrix dscribed with the input formula_composition and formula_variability
#'   \item factor - The factor in the formula corresponding to the covariate, if exists (e.g. it does not exist in case og Intercept or contrasts, which usually are combination of parameters)
#'
#'   \item c_lower - lower (2.5%) quantile of the posterior distribution for a composition (c) parameter.
#'   \item c_effect - mean of the posterior distribution for a composition (c) parameter.
#'   \item c_upper - upper (97.5%) quantile of the posterior distribution fo a composition (c)  parameter.
#'   \item c_pH0 - Probability of the null hypothesis (no difference) for  a composition (c). This is not a p-value.
#'   \item c_FDR - False-discovery rate of the null hypothesis (no difference) for  a composition (c).
#'   \item c_n_eff - Effective sample size - the number of independent draws in the sample, the higher the better (mc-stan.org/docs/2_25/cmdstan-guide/stansummary.html).
#'   \item c_R_k_hat - R statistic, a measure of chain equilibrium, should be within 0.05 of 1.0 (mc-stan.org/docs/2_25/cmdstan-guide/stansummary.html).
#'
#'   \item v_lower - Lower (2.5%) quantile of the posterior distribution for a variability (v) parameter
#'   \item v_effect - Mean of the posterior distribution for a variability (v) parameter
#'   \item v_upper - Upper (97.5%) quantile of the posterior distribution for a variability (v) parameter
#'   \item v_pH0 - Probability of the null hypothesis (no difference) for a variability (v). This is not a p-value.
#'   \item v_FDR - False-discovery rate of the null hypothesis (no difference), for a variability (v).
#'   \item v_n_eff - Effective sample size for a variability (v) parameter - the number of independent draws in the sample, the higher the better (mc-stan.org/docs/2_25/cmdstan-guide/stansummary.html).
#'   \item v_R_k_hat - R statistic for a variability (v) parameter, a measure of chain equilibrium, should be within 0.05 of 1.0 (mc-stan.org/docs/2_25/cmdstan-guide/stansummary.html).
#'
#'   \item count_data Nested input count data.
#'
#' }
#'
#' @examples
#'
#' data("counts_obj")
#'
#' estimate =
#'   sccomp_glm(
#'   counts_obj ,
#'    ~ type,
#'    ~1,
#'    sample,
#'    cell_group,
#'    count,
#'     check_outliers = FALSE,
#'     cores = 1
#'   )
#'
#' @export
#'
#'
sccomp_glm <- function(.data,
                       formula_composition = ~ 1 ,
                       formula_variability = ~ 1,
                       .sample,
                       .cell_group,
                       .count = NULL,

                       # Secondary arguments
                       contrasts = NULL,
                       prior_mean_variable_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                       check_outliers = TRUE,
                       bimodal_mean_variability_association = FALSE,
                       enable_loo = FALSE,

                       # Tertiary arguments
                       cores = detectCores(),
                       percent_false_positive = 5,
                       approximate_posterior_inference = "none",
                       test_composition_above_logit_fold_change = 0.2,
                       verbose = FALSE,
                       noise_model = "multi_beta_binomial",
                       exclude_priors = FALSE,
                       use_data = TRUE,
                       mcmc_seed = sample(1e5, 1),
                       max_sampling_iterations = 20000,
                       pass_fit = TRUE) {
  UseMethod("sccomp_glm", .data)
}

#' @export
sccomp_glm.Seurat = function(.data,
                             formula_composition = ~ 1 ,
                             formula_variability = ~ 1,
                             .sample,
                             .cell_group,
                             .count = NULL,

                             # Secondary arguments
                             contrasts = NULL,
                             prior_mean_variable_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                             check_outliers = TRUE,
                             bimodal_mean_variability_association = FALSE,
                             enable_loo = FALSE,

                             # Tertiary arguments
                             cores = detectCores(),
                             percent_false_positive = 5,
                             approximate_posterior_inference = "none",
                             test_composition_above_logit_fold_change = 0.2,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             exclude_priors = FALSE,
                             use_data = TRUE,
                             mcmc_seed = sample(1e5, 1),
                             max_sampling_iterations = 20000,
                             pass_fit = TRUE) {

  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")

  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)

  .data[[]] %>%
    sccomp_glm(
      formula_composition = formula_composition,
      formula_variability = formula_variability,

      !!.sample,!!.cell_group,
      contrasts = contrasts,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      noise_model = noise_model,
      exclude_priors = exclude_priors,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      enable_loo = enable_loo,
      use_data = use_data,
      cores = cores,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )


}

#' @export
sccomp_glm.SingleCellExperiment = function(.data,
                                           formula_composition = ~ 1 ,
                                           formula_variability = ~ 1,
                                           .sample,
                                           .cell_group,
                                           .count = NULL,

                                           # Secondary arguments
                                           contrasts = NULL,
                                           prior_mean_variable_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                                           check_outliers = TRUE,
                                           bimodal_mean_variability_association = FALSE,
                                           enable_loo = FALSE,

                                           # Tertiary arguments
                                           cores = detectCores(),
                                           percent_false_positive = 5,
                                           approximate_posterior_inference = "none",
                                           test_composition_above_logit_fold_change = 0.2,
                                           verbose = FALSE,
                                           noise_model = "multi_beta_binomial",
                                           exclude_priors = FALSE,
                                           use_data = TRUE,
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
      formula_composition = formula_composition,
      formula_variability = formula_variability,

      !!.sample,!!.cell_group,
      check_outliers = check_outliers,
      contrasts = contrasts,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      approximate_posterior_inference = approximate_posterior_inference,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      noise_model = noise_model,
      exclude_priors = exclude_priors,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      enable_loo = enable_loo,
      use_data = use_data,
      cores = cores,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )


}

#' @export
sccomp_glm.DFrame = function(.data,
                             formula_composition = ~ 1 ,
                             formula_variability = ~ 1,
                             .sample,
                             .cell_group,
                             .count = NULL,

                             # Secondary arguments
                             contrasts = NULL,
                             prior_mean_variable_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                             check_outliers = TRUE,
                             bimodal_mean_variability_association = FALSE,
                             enable_loo = FALSE,

                             # Tertiary arguments
                             cores = detectCores(),
                             percent_false_positive = 5,
                             approximate_posterior_inference = "none",
                             test_composition_above_logit_fold_change = 0.2,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             exclude_priors = FALSE,
                             use_data = TRUE,
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
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      !!.sample,!!.cell_group,
      contrasts = contrasts,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      noise_model = noise_model,
      exclude_priors = exclude_priors,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      enable_loo = enable_loo,
      use_data = use_data,
      cores = cores,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )
}

#' @importFrom purrr when
#' @export
sccomp_glm.data.frame = function(.data,
                                 formula_composition = ~ 1 ,
                                 formula_variability = ~ 1,
                                 .sample,
                                 .cell_group,
                                 .count = NULL,

                                 # Secondary arguments
                                 contrasts = NULL,
                                 prior_mean_variable_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                                 check_outliers = TRUE,
                                 bimodal_mean_variability_association = FALSE,
                                 enable_loo = FALSE,

                                 # Tertiary arguments
                                 cores = detectCores(),
                                 percent_false_positive = 5,
                                 approximate_posterior_inference = "none",
                                 test_composition_above_logit_fold_change = 0.2,
                                 verbose = FALSE,
                                 noise_model = "multi_beta_binomial",
                                 exclude_priors = FALSE,
                                 use_data = TRUE,
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
        formula_composition = formula_composition,
        formula_variability = formula_variability,

        !!.sample,
        !!.cell_group,
        contrasts = contrasts,
        prior_mean_variable_association = prior_mean_variable_association,
        percent_false_positive = percent_false_positive ,
        check_outliers = check_outliers,
        approximate_posterior_inference = approximate_posterior_inference,
        test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
        verbose = verbose,
        my_glm_model = my_glm_model,
        exclude_priors = exclude_priors,
        bimodal_mean_variability_association = bimodal_mean_variability_association,
        enable_loo = enable_loo,
        use_data = use_data,
        cores = cores,
        mcmc_seed = mcmc_seed,
        max_sampling_iterations = max_sampling_iterations,
        pass_fit = pass_fit
      ),

      # If the dataframe does includes counts
      ~ sccomp_glm_data_frame_counts(
        .data,
        formula_composition = formula_composition,
        formula_variability = formula_variability,

        !!.sample,
        !!.cell_group,
        !!.count,
        contrasts = contrasts,
        prior_mean_variable_association = prior_mean_variable_association,
        percent_false_positive = percent_false_positive ,
        check_outliers = check_outliers,
        approximate_posterior_inference = approximate_posterior_inference,
        test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
        verbose = verbose,
        my_glm_model = my_glm_model,
        exclude_priors = exclude_priors,
        bimodal_mean_variability_association = bimodal_mean_variability_association,
        enable_loo = enable_loo,
        use_data = use_data,
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
                                     formula_composition = ~ 1 ,
                                     formula_variability = ~ 1,
                                     .sample,
                                     .cell_group,
                                     .count = NULL,

                                     # Secondary arguments
                                     contrasts = NULL,
                                     prior_mean_variable_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                                     percent_false_positive =  5,
                                     check_outliers = TRUE,
                                     approximate_posterior_inference = "none",
                                     test_composition_above_logit_fold_change = 0.2,
                                     verbose = FALSE,
                                     my_glm_model,
                                     exclude_priors = FALSE,
                                     bimodal_mean_variability_association = FALSE,
                                     enable_loo = FALSE,
                                     use_data = TRUE,
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
    parse_formula(formula_composition)
  ))

  # Check if any column is NA or null
  check_if_any_NA(.data, c(
    quo_name(.sample),
    quo_name(.cell_group),
    parse_formula(formula_composition)
  ))

  .grouping_for_random_intercept = parse_formula_random_intercept(formula_composition) |> pull(grouping) |> unique()

  # Make counts
  .data %>%
    class_list_to_counts(!!.sample, !!.cell_group) %>%

    # Add formula_composition information
    when(
      length(parse_formula(formula_composition))>0 ~
        left_join(.,
                  .data %>%
                    select(!!.sample, parse_formula(formula_composition), any_of(.grouping_for_random_intercept)) %>%
                    distinct(),
                  by = quo_name(.sample)
        ),
      ~ (.)
    ) %>%

    # Return
    sccomp_glm_data_frame_counts(
      formula_composition = formula_composition,
      formula_variability = formula_variability,

      .sample = !!.sample,
      .cell_group = !!.cell_group,
      .count = count,
      my_glm_model = my_glm_model,
      contrasts = contrasts,
      #.grouping_for_random_intercept = !! .grouping_for_random_intercept,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive =  percent_false_positive,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      exclude_priors = exclude_priors,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      enable_loo = enable_loo,
      use_data = use_data,
      cores = cores,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )
}

sccomp_glm_data_frame_counts = function(.data,
                                        formula_composition = ~ 1 ,
                                        formula_variability = ~ 1,
                                        .sample,
                                        .cell_group,
                                        .count = NULL,

                                        # Secondary arguments
                                        contrasts = NULL,
                                        #.grouping_for_random_intercept = NULL,
                                        prior_mean_variable_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                                        percent_false_positive = 5,
                                        check_outliers = TRUE,
                                        approximate_posterior_inference = "none",
                                        test_composition_above_logit_fold_change = 0.2,
                                        verbose = FALSE,
                                        my_glm_model ,
                                        exclude_priors = FALSE,
                                        bimodal_mean_variability_association = FALSE,
                                        enable_loo = FALSE,
                                        use_data = TRUE,
                                        cores = 4,
                                        mcmc_seed = sample(1e5, 1),
                                        max_sampling_iterations = 20000,
                                        pass_fit = TRUE) {

  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)
  #.grouping_for_random_intercept = enquo(.grouping_for_random_intercept)


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
    parse_formula(formula_composition)
  ))

  # Check if any column is NA or null
  check_if_any_NA(.data, c(
    quo_name(.sample),
    quo_name(.cell_group),
    quo_name(.count),
    parse_formula(formula_composition)
  ))

  # Return
  .data %>%
    my_glm_model(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      .sample = !!.sample,
      .cell_group = !!.cell_group,
      .count = !!.count,
      contrasts = contrasts,
      #.grouping_for_random_intercept = !! .grouping_for_random_intercept,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      exclude_priors = exclude_priors,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      enable_loo = enable_loo,
      use_data = use_data,
      cores = cores,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change,
      verbose = verbose,
      seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    ) %>%
    add_attr(.sample, ".sample") %>%
    add_attr(.cell_group, ".cell_group") %>%
    add_attr(.count, ".count") %>%
    add_attr(parse_formula(formula_composition), "factors" )
}


#' test_contrasts
#'
#' @description This function test ocntrasts from a sccomp result.
#'
#'
#' @param .data A tibble. The result of sccomp_glm.
#' @param contrasts A vector of character strings. For example if your formula is `~ 0 + treatment` and the factor treatment has values `yes` and `no`, your contrast could be "constrasts = c(treatmentyes - treatmentno)".
#' @param percent_false_positive A real between 0 and 100 non included. This used to identify outliers with a specific false positive rate.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#'
#' @return A nested tibble `tbl` with cell_group-wise statistics
#'
#' @export
#'
#' @examples
#'
#' data("counts_obj")
#'
#'   estimates =
#'   sccomp_glm(
#'   counts_obj ,
#'    ~ 0 + type, ~1,  sample, cell_group, count,
#'     check_outliers = FALSE,
#'     cores = 1
#'   ) |>
#'
#'   test_contrasts("typecancer - typebenign")
#'
test_contrasts <- function(.data,
                           contrasts = NULL,
                           percent_false_positive = 5,
                           test_composition_above_logit_fold_change = 0.2) {
  UseMethod("test_contrasts", .data)
}


#' @export
#'
test_contrasts.data.frame = function(.data,
                                     contrasts = NULL,
                                     percent_false_positive = 5,
                                     test_composition_above_logit_fold_change = 0.2){


  .sample = .data |>  attr(".sample")
  .cell_group = .data |>  attr(".cell_group")
  .count = .data |>  attr(".count")
  check_outliers = .data |>  attr("check_outliers")
  model_input = .data |> attr("model_input")
  truncation_df2 =  .data |>  attr("truncation_df2")

  # Abundance
  abundance_CI =
    get_abundance_contrast_draws(.data, contrasts) %>%

    # If my constrasts do not match my model. I have to do something more elegant.
    when(
      "parameter" %in% colnames(.) ~  draws_to_statistics(
        .,
        percent_false_positive/100,
        test_composition_above_logit_fold_change,
        !!.cell_group,
        "c_"
      ),
    ~ (.)
  )


  # Variability
  variability_CI =
    get_variability_contrast_draws(.data, contrasts) |>
    draws_to_statistics(
      percent_false_positive/100,
      test_composition_above_logit_fold_change,
      !!.cell_group,
      "v_"
    )

  # Merge and parse
  abundance_CI |>

    # Add ALPHA
    left_join(variability_CI) |>
    suppressMessages() |>

    # Add easy to understand factor labels
    left_join(
      model_input$factor_parameter_dictionary |>

        # If I don't have factors (~1)
        when(
          !"factor" %in% colnames(.) ~ tibble(`factor`=character(), design_matrix_col = character()),
          ~ (.)
        ) |>

        select(`factor`, design_matrix_col),
      by = c("parameter" = "design_matrix_col" )
    ) %>%
    select(parameter, `factor`, everything()) %>%

    select(!!.cell_group, everything(), -M) %>%

    # Add outlier
    when(
      check_outliers ~ (.) %>%
        left_join(
          truncation_df2 |>
            select(-c(M, N, .variable, mean, se_mean, sd, n_eff, one_of("R_hat", "k_hat"))) |>
            suppressWarnings() |>
            nest(count_data = -!!.cell_group),
          by = quo_name(.cell_group)
        ),
      ~ (.) %>% left_join(
        truncation_df2  |>
          nest(count_data = -!!.cell_group),
        by = quo_name(.cell_group)
      )
    )
}

#' sccomp_replicate
#'
#' @description This function replicates counts from a real-world dataset.
#'
#'
#' @param fit The result of sccomp_glm.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment. In most cases, if differentially variability is of interest, the formula should only include the factor of interest as a large anount of data is needed to define variability depending to each factors. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
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
#' if(.Platform$OS.type == "unix")
#'   sccomp_glm(
#'   counts_obj ,
#'    ~ type, ~1,  sample, cell_group, count,
#'     approximate_posterior_inference = "all",
#'     check_outliers = FALSE,
#'     cores = 1
#'   ) |>
#'
#'   sccomp_replicate()
#'
sccomp_replicate <- function(fit,
                             formula_composition = NULL,
                             formula_variability = NULL,
                             number_of_draws = 1,
                             mcmc_seed = sample(1e5, 1)) {
  UseMethod("sccomp_replicate", fit)
}

#' @export
#'
sccomp_replicate.data.frame = function(fit,
                                       formula_composition = NULL,
                                       formula_variability = NULL,
                                       number_of_draws = 1,
                                       mcmc_seed = sample(1e5, 1)){

  .sample = attr(fit, ".sample")
  .cell_group = attr(fit, ".cell_group")

  rng =
    replicate_data(
      fit,
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      number_of_draws = number_of_draws,
      mcmc_seed = mcmc_seed
    )

  model_input = attr(fit, "model_input")

  # mean generated
  rng |>

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

    select(-N, -M) |>
    select(!!.cell_group, !!.sample, everything())

}

#' sccomp_predict
#'
#' @description This function replicates counts from a real-world dataset.
#'
#'
#' @param fit The result of sccomp_glm.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#' @param new_data A sample-wise data frame including the column that represent the factors in your formula. If you want to predict proportions for 10 samples, there should be 10 rows. T
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
#' if(.Platform$OS.type == "unix")
#'   sccomp_glm(
#'   counts_obj ,
#'    ~ type, ~1,  sample, cell_group, count,
#'     approximate_posterior_inference = "all",
#'     check_outliers = FALSE,
#'     cores = 1
#'   ) |>
#'
#'   sccomp_predict()
#'
sccomp_predict <- function(fit,
                           formula_composition = NULL,
                           new_data = NULL,
                           number_of_draws = 500,
                           mcmc_seed = sample(1e5, 1)) {
  UseMethod("sccomp_predict", fit)
}

#' @export
#'
sccomp_predict.data.frame = function(fit,
                                     formula_composition = NULL,
                                     new_data = NULL,
                                     number_of_draws = 500,
                                     mcmc_seed = sample(1e5, 1)){


  model_input = attr(fit, "model_input")
  .sample = attr(fit, ".sample")
  .cell_group = attr(fit, ".cell_group")

  if(new_data |> is.null())
    new_data =
      fit |>
      select(count_data) |>
      unnest(count_data) |>
      distinct() |>
      .subset(!!.sample)

  sample_names = new_data |> distinct(!!.sample) |> pull(!!.sample)

  rng =
    replicate_data(
      fit,
      formula_composition = formula_composition,
      new_data = new_data,
      number_of_draws = number_of_draws,
      mcmc_seed = mcmc_seed
    )

  # mean generated
  rng |>
    summary_to_tibble("mu", "M", "N") |>
    select(M, N, proportion_mean = mean, proportion_lower = `2.5%`, proportion_upper = `97.5%`) |>

    # Get sample name
    nest(data = -N) %>%
    arrange(N) %>%
    mutate(!!.sample := sample_names) %>%
    unnest(data) %>%

    # get cell type name
    nest(data = -M) %>%
    mutate(!!.cell_group := colnames(model_input$y)) %>%
    unnest(data) %>%

    select(-N, -M) |>
    select(!!.cell_group, !!.sample, everything())

}


#' remove_unwanted_variation
#'
#' @description This function uses the model to remove unwanted variation from a dataset using the estimated of the model. For example if you fit your data with this formula `~ factor_1 + factor_2` and use this formula to remove unwanted variation `~ factor_1`, the `factor_2` will be factored out.
#'
#'
#' @param .data A tibble. The result of sccomp_glm.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment. In most cases, if differentially variability is of interest, the formula should only include the factor of interest as a large anount of data is needed to define variability depending to each factors. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#'
#' @return A nested tibble `tbl` with cell_group-wise statistics
#'
#' @export
#'
#' @examples
#'
#' data("counts_obj")
#'
#'   estimates = sccomp_glm(
#'   counts_obj ,
#'    ~ type, ~1,  sample, cell_group, count,
#'     approximate_posterior_inference = "all",
#'     check_outliers = FALSE,
#'     cores = 1
#'   )
#'
#'   remove_unwanted_variation(estimates)
#'
remove_unwanted_variation <- function(.data,
                                      formula_composition = ~1,
                                      formula_variability = NULL) {
  UseMethod("remove_unwanted_variation", .data)
}

#' @export
#'
remove_unwanted_variation.data.frame = function(.data,
                                                formula_composition = ~1,
                                                formula_variability = NULL){



  model_input = attr(.data, "model_input")
  .sample = attr(.data, ".sample")
  .cell_group = attr(.data, ".cell_group")
  .grouping_for_random_intercept = attr(.data, ".grouping_for_random_intercept")
  .count = attr(.data, ".count")

  fit_matrix = as.matrix(attr(.data, "fit") )

  message("sccomp says: calculating residuals")

  # Residuals
  residuals =
    .data |>
    sccomp_predict(
      number_of_draws = min(dim(fit_matrix)[1], 500)
    ) |>
    distinct(!!.sample, !!.cell_group, proportion_mean) |>
    mutate( proportion_mean =
              proportion_mean |>
             #compress_zero_one() |>
             boot::logit()
    )|>
  # Join counts
  left_join(
    .data |>
      attr("model_input") %$%
      y |>
      as_tibble(rownames = quo_name(.sample)) |>
      pivot_longer(-!!.sample, names_to = quo_name(.cell_group), values_to = quo_name(.count)) |>
      with_groups(!!.sample,  ~ .x |> mutate(observed_proportion := !!.count / sum(!!.count ))) |>

      with_groups(!!.sample,  ~ .x |>  mutate(exposure := sum(!!.count))  ) |>

      mutate(observed_proportion =
               observed_proportion |>
               compress_zero_one() |>
               boot::logit()
      ),
    by = c(quo_name(.sample), quo_name(.cell_group))
  ) |>
  mutate(logit_residuals = observed_proportion - proportion_mean) |>
  select(!!.sample, !!.cell_group, logit_residuals, exposure)


  message("sccomp says: regressing out unwanted factors")

  # Generate quantities
  .data |>
    sccomp_predict(
      formula_composition = formula_composition,
      number_of_draws = min(dim(fit_matrix)[1], 500)
    ) |>
    distinct(!!.sample, !!.cell_group, proportion_mean) |>
    mutate(proportion_mean =
             proportion_mean |>
             # compress_zero_one() |>
             boot::logit()
    ) |>
    left_join(residuals,  by = c(quo_name(.sample), quo_name(.cell_group))) |>
    mutate(adjusted_proportion = proportion_mean + logit_residuals) |>
    mutate(adjusted_proportion = adjusted_proportion |> boot::inv.logit()) |>
    with_groups(!!.sample,  ~ .x |> mutate(adjusted_proportion := adjusted_proportion / sum(adjusted_proportion ))) |>

    # Recostituite counts
    mutate(adjusted_counts = adjusted_proportion * exposure) |>

    select(!!.sample, !!.cell_group, adjusted_proportion, adjusted_counts, logit_residuals)



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
#' @param .data A tibble including a cell_group name column | sample name column | read counts column | factor columns | Pvalue column | a significance column
#' @param .estimate_object The result of sccomp_glm execution. This is used for sampling from real-data properties.
#' @param formula_composition A formula. The sample formula used to perform the differential cell_group abundance analysis
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment. In most cases, if differentially variability is of interest, the formula should only include the factor of interest as a large anount of data is needed to define variability depending to each factors.
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_group identifier
#' @param .coefficients The column names for coefficients, for example, c(b_0, b_1)
#' @param variability_multiplier A real scalar. This can be used for artificially increasing the variability of the simulation for benchmarking purposes.
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
#'  sccomp_glm(
#'  counts_obj ,
#'   ~ type, ~1,  sample, cell_group, count,
#'    approximate_posterior_inference = "all",
#'    check_outliers = FALSE,
#'    cores = 1
#'  )
#'
#' # Set coefficients for cell_groups. In this case all coefficients are 0 for simplicity.
#' counts_obj = counts_obj |> mutate(b_0 = 0, b_1 = 0)

#' # Simulate data
#' simulate_data(counts_obj, estimate, ~type, ~1, sample, cell_group, c(b_0, b_1))
#'
simulate_data <- function(.data,
                          .estimate_object,
                          formula_composition,
                          formula_variability = NULL,
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
#' @importFrom readr read_file
#'
simulate_data.data.frame = function(.data,
                                    .estimate_object,
                                    formula_composition,
                                    formula_variability = NULL,
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
    (.) == "logit_normal_multinomial" ~ get_model_from_data("glm_multinomial_logit_linear_simulate_data.stan", read_file("~/PostDoc/sccomp/dev/stan_models/glm_multinomial_logit_linear_simulate_data.stan"))

  )

  model_input =
    .data %>%
    nest(data___ = -!!.sample) %>%
    mutate(.exposure = sample(model_data$exposure, size = n(), replace = TRUE )) %>%
    unnest(data___) %>%
    data_simulation_to_model_input(
      formula_composition,
      #formula_variability,
      !!.sample, !!.cell_group, .exposure, !!.coefficients
    )

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


#' plot_summary
#'
#' @description This function plots a summary of the results of the model.
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @importFrom tidyr unite
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr with_groups
#'
#' @param .data A tibble including a cell_group name column | sample name column | read counts column | factor columns | Pvalue column | a significance column
#' @param significance_threshold A real. FDR threshold for labelling significant cell-groups.
#'
#' @return A `ggplot`
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
#'    ~ type, ~1, sample, cell_group, count,
#'     approximate_posterior_inference = "all",
#'     check_outliers = FALSE,
#'     cores = 1
#'   )
#'
#' # estimate |> plot_summary()
#'
plot_summary <- function(.data, significance_threshold = 0.025) {

    multipanel_theme =
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(linewidth=0.5),
      panel.grid.major = element_line(linewidth = 0.1),
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

      axis.line.x = element_line(linewidth=0.2),
      axis.line.y = element_line(linewidth=0.2),
      axis.ticks.x = element_line(linewidth=0.2),
      axis.ticks.y = element_line(linewidth=0.2)
    )

  .cell_group = attr(.data, ".cell_group")
  .count = attr(.data, ".count")
  .sample = attr(.data, ".sample")

plots = list()


data_proportion =
  .data %>%

  # Otherwise does not work
  select(-`factor`) %>%

  pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) %>%
  unnest(count_data) %>%
  with_groups(!!.sample, ~ mutate(.x, proportion = (!!.count)/sum(!!.count)) ) |>

  # If I don't have outliers add them
  when(!"outlier" %in% colnames(.) ~ mutate(., outlier = FALSE), ~ (.))


# Boxplot
plots$boxplot =

  # Select non numerical types
  .data %>%
  filter(!is.na(`factor`)) |>
   distinct(`factor`) %>%
    pull(`factor`) |>

  map(
    ~ plot_boxplot(
      .data,
      data_proportion,
      .x,
      !!.cell_group,
      !!.sample,
      significance_threshold = significance_threshold,
      multipanel_theme
    ) +
      ggtitle(sprintf("Grouped by %s (for multi-factor models, associations could be hardly observable with unidimensional data stratification)", .x))
  )

# 1D intervals
plots$credible_intervals_1D = plot_1d_intervals(.data, !!.cell_group, significance_threshold = significance_threshold, multipanel_theme)

# 2D intervals
if("v_effect" %in% colnames(.data) && (.data |> filter(!is.na(v_effect)) |> nrow()) > 0)  plots$credible_intervals_2D = plot_2d_intervals(.data, !!.cell_group, significance_threshold = significance_threshold, multipanel_theme)

plots

}

