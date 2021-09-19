#' multi_beta_binomial main
#'
#' @description This function runs the data modelling and statistical test for the hypothesis that a cell_type includes outlier biological replicate.
#'
#' @importFrom tibble as_tibble
#' @import dplyr
#' @importFrom tidyr spread
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom purrr map2
#' @importFrom purrr map_int
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom purrr map
#' @importFrom tibble rowid_to_column
#' @importFrom purrr map_lgl
#' @importFrom dplyr case_when
#' @importFrom rlang :=
#'
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | covariate columns | Pvaue column | a significance column
#' @param formula A formula. The sample formula used to perform the differential cell_type abundance analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_type A column name as symbol. The cell_type identifier
#' @param .count A column name as symbol. The cell_type abundance (read count)
#' @param check_outliers A boolean. Whether to check for outliers before the fit.
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param variance_association A boolean. Whether the variance should be associated to the factor of interest.
#' @param verbose A boolean. Prints progression.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param seed An integer. Used for development and testing purposes
#'
#'
#' @noRd
#'
#' @return A List object
#'
#'
estimate_multi_beta_binomial_glm = function(.data,
                                            formula = ~ 1,
                                            .sample,
                                            .cell_type,
                                            .count,
                                            prior_mean_variable_association,
                                            percent_false_positive = 5,
                                            check_outliers = FALSE,
                                            approximate_posterior_inference = TRUE,
                                            variance_association = FALSE,
                                            cores = detectCores(), # For development purpose,
                                            seed = sample(1:99999, size = 1),
                                            verbose = FALSE
) {
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  CI = 1 - (percent_false_positive/100)

  # Produce data list
  covariate_names = parse_formula(formula)

  data_for_model =
    .data %>%
    data_to_spread ( formula, !!.sample, !!.cell_type, !!.count) %>%
    data_spread_to_model_input(
      formula, !!.sample, !!.cell_type, !!.count,
      variance_association = variance_association,
      truncation_ajustment = 1.1,
      approximate_posterior_inference = approximate_posterior_inference
    )

  # Pior
  data_for_model$prior_prec_intercept = prior_mean_variable_association$intercept
  data_for_model$prior_prec_slope  = prior_mean_variable_association$slope
  data_for_model$prior_prec_sd = prior_mean_variable_association$standard_deviation

  # Original - old
  # prec_sd ~ normal(0,2);
  # prec_coeff ~ normal(0,5);


  # Start first fit
  message("sccomp says: outlier identification first pass - step 1/3 [ETA: ~20s]")


  fit =
    data_for_model %>%

    # Run the first discovery phase with permissive false discovery rate
    fit_model(stanmodels$glm_multi_beta_binomial, cores= cores,  quantile = CI,  approximate_posterior_inference = approximate_posterior_inference, verbose = verbose, seed = seed)

  rng =  rstan::gqs(
    stanmodels$glm_multi_beta_binomial_generate_date,
    #rstan::stan_model("inst/stan/glm_multi_beta_binomial_generate_date.stan"),
    draws =  as.matrix(fit),
    data = data_for_model
  )

  if(!check_outliers){

    generated_quantities =
      parse_generated_quantities(rng) %>%

      # Get sample name
      nest(data = -N) %>%
      arrange(N) %>%
      mutate(!!.sample := rownames(data_for_model$y)) %>%
      unnest(data) %>%

      # get cell type name
      nest(data = -M) %>%
      mutate(!!.cell_type := colnames(data_for_model$y)) %>%
      unnest(data) %>%

      select(-N, -M) %>%
      nest(generated_data = -c(!!.sample, !!.cell_type))

    list(
      fit = fit,
      data_for_model = data_for_model,
      truncation_df2 =  .data %>% left_join(generated_quantities, by = c(quo_name(.sample), quo_name(.cell_type)))
    )

  }

  else{


    # Detect outliers
    truncation_df =
      .data %>%
      left_join(
        summary_to_tibble(rng, "counts", "N", "M", probs = c(0.05, 0.95)) %>%
          nest(data = -N) %>%
          mutate(!!.sample := rownames(data_for_model$y)) %>%
          unnest(data) %>%
          nest(data = -M) %>%
          mutate(!!.cell_type := colnames(data_for_model$y)) %>%
          unnest(data) ,

        by = c(quo_name(.sample), quo_name(.cell_type))
      ) %>%

      # Add truncation
      mutate(   truncation_down = `5%`,   truncation_up =  `95%`) %>%

      # Add outlier stats
      mutate( outlier = !(!!.count >= `5%` & !!.count <= `95%`) ) %>%
      nest(data = -M) %>%
      mutate(contains_outliers = map_lgl(data, ~ .x %>% filter(outlier) %>% nrow() %>% `>` (0))) %>%
      unnest(data) %>%

      mutate(
        truncation_down = case_when( outlier ~ -1, TRUE ~ truncation_down),
        truncation_up = case_when(outlier ~ -1, TRUE ~ truncation_up),
      )

    # Add censoring
    data_for_model$is_truncated = 1
    data_for_model$truncation_up = truncation_df %>% select(N, M, truncation_up) %>% spread(M, truncation_up) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)
    data_for_model$truncation_down = truncation_df %>% select(N, M, truncation_down) %>% spread(M, truncation_down) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)

    message("sccomp says: outlier identification second pass - step 2/3 [ETA: ~60s]")

    my_quantile_step_2 = 1 - (0.1 / data_for_model$N)

    # This step gets the credible interval to control for within-category false positive rate
    # We want a category-wise false positive rate of 0.1, and we have to correct for how many samples we have in each category
    CI_step_2 = (1-my_quantile_step_2) / 2 * 2


    fit2 =
      data_for_model %>%
      fit_model(stanmodels$glm_multi_beta_binomial, cores = cores, quantile = my_quantile_step_2,  approximate_posterior_inference = approximate_posterior_inference, verbose = verbose, seed = seed)
    #fit_model(stan_model("inst/stan/glm_multi_beta_binomial.stan"), chains= 4, output_samples = 500, approximate_posterior_inference = FALSE, verbose = TRUE)

    rng2 =  rstan::gqs(
      stanmodels$glm_multi_beta_binomial_generate_date,
      #rstan::stan_model("inst/stan/glm_multi_beta_binomial_generate_date.stan"),
      draws =  as.matrix(fit2),
      data = data_for_model
    )

    # Detect outliers
    truncation_df2 =
      .data %>%
      left_join(
        summary_to_tibble(rng2, "counts", "N", "M", probs = c(CI_step_2, 0.5, 1-CI_step_2)) %>%

          # !!! THIS COMMAND RELIES ON POSITION BECAUSE IT'S NOT TRIVIAL TO MATCH
          # !!! COLUMN NAMES BASED ON LIMITED PRECISION AND/OR PERIODICAL QUANTILES
          rename(
            .lower := !!as.symbol(colnames(.)[7]) ,
            .median = `50%`,
            .upper := !!as.symbol(colnames(.)[9])
          ) %>%
          nest(data = -N) %>%
          mutate(!!.sample := rownames(data_for_model$y)) %>%
          unnest(data) %>%
          nest(data = -M) %>%
          mutate(!!.cell_type := colnames(data_for_model$y)) %>%
          unnest(data) ,

        by = c(quo_name(.sample), quo_name(.cell_type))
      ) %>%

      # Add truncation
      mutate(   truncation_down = .lower,   truncation_up =  .upper) %>%

      # Add outlier stats
      mutate( outlier = !(!!.count >= .lower & !!.count <= .upper) ) %>%
      nest(data = -M) %>%
      mutate(contains_outliers = map_lgl(data, ~ .x %>% filter(outlier) %>% nrow() %>% `>` (0))) %>%
      unnest(data) %>%

      mutate(
        truncation_down = case_when( outlier ~ -1, TRUE ~ truncation_down),
        truncation_up = case_when(outlier ~ -1, TRUE ~ truncation_up),
      )


    data_for_model$truncation_up = truncation_df2 %>% select(N, M, truncation_up) %>% spread(M, truncation_up) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)
    data_for_model$truncation_down = truncation_df2 %>% select(N, M, truncation_down) %>% spread(M, truncation_down) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)

    message("sccomp says: outlier-free model fitting - step 3/3 [ETA: ~20s]")

    fit3 =
      data_for_model %>%
      # Run the first discovery phase with permissive false discovery rate
      fit_model(stanmodels$glm_multi_beta_binomial, cores = cores, quantile = CI,  approximate_posterior_inference = approximate_posterior_inference, verbose = verbose, seed = seed)
    #fit_model(stan_model("inst/stan/glm_multi_beta_binomial.stan"), chains= 4, output_samples = 500)

    rng3 =  rstan::gqs(
      stanmodels$glm_multi_beta_binomial_generate_date,
      #rstan::stan_model("inst/stan/glm_multi_beta_binomial_generate_date.stan"),
      draws =  as.matrix(fit),
      data = data_for_model
    )

    generated_quantities3 =
      parse_generated_quantities(rng3) %>%

      # Get sample name
      nest(data = -N) %>%
      arrange(N) %>%
      mutate(!!.sample := rownames(data_for_model$y)) %>%
      unnest(data) %>%

      # get cell type name
      nest(data = -M) %>%
      mutate(!!.cell_type := colnames(data_for_model$y)) %>%
      unnest(data) %>%

      select(-N, -M) %>%
      nest(generated_data = -c(!!.sample, !!.cell_type))


    list(
      fit = fit3,
      data_for_model = data_for_model,
      truncation_df2 =
        truncation_df2 %>%
        left_join(generated_quantities3, by = c(quo_name(.sample), quo_name(.cell_type)))
    )


  }


}

#' multi_beta_binomial main
#'
#' @description This function runs the data modelling and statistical test for the hypothesis that a cell_type includes outlier biological replicate.
#'
#' @importFrom tibble as_tibble
#' @import dplyr
#' @importFrom tidyr spread
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom purrr map2
#' @importFrom purrr map_int
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom purrr map
#' @importFrom tibble rowid_to_column
#' @importFrom purrr map_lgl
#' @importFrom dplyr case_when
#'
#' @param fit The fit object
#' @param data_for_model Parsed data
#' @param check_outliers A boolean
#' @param truncation_df2 Truncation data frame
#'
#' @noRd
#'
#' @return A nested tibble `tbl` with cell_type-wise information: `sample wise data` | plot | `ppc samples failed` | `exposure deleterious outliers`
#'
#'
hypothesis_test_multi_beta_binomial_glm = function( .sample,
                                                    .cell_type, .count, fit, data_for_model, percent_false_positive, check_outliers,  truncation_df2 = NULL) {

  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  parsed_fit =
    fit %>%
    parse_fit(data_for_model, .)

  minimum_logit_fold_change = 0.2

  parsed_fit %>%
    beta_to_CI(false_positive_rate = percent_false_positive/100 ) %>%

    # Hypothesis test
    when(
      ncol(data_for_model$X) > 1 ~
        mutate(
          .,
          significant =

            !!as.symbol(sprintf(".lower_%s", colnames(data_for_model$X)[2])) > minimum_logit_fold_change |
            !!as.symbol(sprintf(".upper_%s", colnames(data_for_model$X)[2])) < -minimum_logit_fold_change
        ) %>%

        # add probability
       # left_join( get_probability_non_zero(parsed_fit), by="M" ),
      ~ (.)
    )






}



#' multi_beta_binomial main
#'
#' @description This function runs the data modelling and statistical test for the hypothesis that a cell_type includes outlier biological replicate.
#'
#' @importFrom tibble as_tibble
#' @import dplyr
#' @importFrom tidyr spread
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom purrr map2
#' @importFrom purrr map_int
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom purrr map
#' @importFrom tibble rowid_to_column
#' @importFrom purrr map_lgl
#' @importFrom dplyr case_when
#' @importFrom rlang :=
#'
#' @keywords internal
#' @noRd
#'
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | covariate columns | Pvaue column | a significance column
#' @param formula A formula. The sample formula used to perform the differential cell_type abundance analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_type A column name as symbol. The cell_type identifier
#' @param .count A column name as symbol. The cell_type abundance (read count)
#' @param check_outliers A boolean. Whether to check for outliers before the fit.
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param variance_association A boolean. Whether the variance should be associated to the factor of interest.
#' @param verbose A boolean. Prints progression.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param seed An integer. Used for development and testing purposes
#'
#' @return A nested tibble `tbl` with cell_type-wise information: `sample wise data` | plot | `ppc samples failed` | `exposure deleterious outliers`
#'
#'
multi_beta_binomial_glm = function(.data,
                                   formula = ~ 1,
                                   .sample,
                                   .cell_type,
                                   .count,
                                   prior_mean_variable_association,
                                   percent_false_positive = 5,
                                   check_outliers = FALSE,
                                   approximate_posterior_inference = TRUE,
                                   variance_association = FALSE,
                                   cores = detectCores(), # For development purpose,
                                   seed = sample(1:99999, size = 1),
                                   verbose = FALSE
) {
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)


  result_list =
    estimate_multi_beta_binomial_glm(
    .data = .data,
    formula = formula,
    .sample = !!.sample,
    .cell_type = !!.cell_type,
    .count = !!.count,
    prior_mean_variable_association = prior_mean_variable_association,
    percent_false_positive = percent_false_positive,
    check_outliers = check_outliers,
    approximate_posterior_inference = approximate_posterior_inference,
    variance_association = variance_association,
    cores = cores, # For development purpose,
    seed = seed,
    verbose = verbose
  )


  hypothesis_test_multi_beta_binomial_glm(
    .sample = !!.sample,
    .cell_type = !!.cell_type,
    .count = !!.count,
    result_list$fit,
    result_list$data_for_model,
    percent_false_positive,
    check_outliers,
    result_list$truncation_df2
  ) %>%

    # Join the precision
    left_join(get_mean_precision(result_list$fit), by="M") %>%


    # Clean
    select(-M) %>%
    mutate(!!.cell_type := result_list$data_for_model$y %>% colnames()) %>%
    select(!!.cell_type, everything()) %>%

    # Join generated data
    # left_join(result_list$generated_quantities, by="M") %>%

    # Add outlier
    when(
      check_outliers ~ (.) %>%
        left_join(
          result_list$truncation_df2 %>%
            select(-c(M, N, .variable, mean, se_mean, sd, n_eff, Rhat)) %>%
            nest(count_data = -!!.cell_type),
          by = quo_name(.cell_type)
        ),
      ~ (.) %>% left_join(result_list$truncation_df2 %>% nest(count_data = -!!.cell_type),  by = quo_name(.cell_type))
    ) %>%

    # Attach association mean concentration
    add_attr(get_mean_precision_association(result_list$fit), "mean_concentration_association")

}

#' @importFrom stats model.matrix
glm_multi_beta_binomial = function(input_df, formula, .sample){

  covariate_names = parse_formula(formula)
  .sample = enquo(.sample)

  sampling(stanmodels$glm_multi_beta_binomial,
           data = list(
             N = input_df %>% nrow(),
             M = input_df %>% select(-!!.sample, -covariate_names, -exposure) %>% ncol(),
             exposure = input_df$exposure,
             y = input_df %>% select(-covariate_names, -exposure) %>% as_matrix(rownames = !!.sample),
             X = input_df %>% select(!!.sample, covariate_names) %>% model.matrix(formula, data=.),
             C = length(covariate_names)
           ),
           cores = 4
  )

}

get_mean_precision = function(fit){
  fit %>%
    summary_to_tibble("alpha", "C", "M") %>%
    select( M, mean, `2.5%` , `97.5%`) %>%
    nest(concentration = -M)
}

get_mean_precision_association = function(fit){
  c(
    fit %>%
      summary("prec_coeff") %$%
      summary %>%
      .[,1] ,

    fit %>%
      summary("prec_sd") %$%
      summary %>%
      .[,1]
  )
}
