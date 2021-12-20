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
#' @param prior_mean_variable_association A list of the form list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(0.4527292, 0.3318759)). Where for each parameter, we specify mean and standard deviation. This is used to incorporate prior knowledge about the mean/variability association of cell-type proportions.
#' @param percent_false_positive A real between 0 and 100. It is the aimed percent of cell types being a false positive. For example, percent_false_positive_genes = 1 provide 1 percent of the calls for significant changes that are actually not significant.
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
                                            formula_variability = ~ 1,
                                            prior_mean_variable_association,
                                            percent_false_positive = 5,
                                            check_outliers = FALSE,
                                            approximate_posterior_inference = "outlier_detection",
                                            variance_association = FALSE,
                                            cores = detectCores(), # For development purpose,
                                            seed = sample(1e5, 1),
                                            verbose = FALSE,
                                            exclude_priors = FALSE,
                                            max_sampling_iterations = 20000
) {
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  CI = 1 - (percent_false_positive/100)

  # Produce data list
  covariate_names = parse_formula(formula)

  # Original - old
  # prec_sd ~ normal(0,2);
  # prec_coeff ~ normal(0,5);

  # If we are NOT checking outliers
  if(!check_outliers){

    message("sccomp says: estimation [ETA: ~20s]")

    data_for_model =
      .data %>%
      data_to_spread ( formula, !!.sample, !!.cell_type, !!.count) %>%
      data_spread_to_model_input(
        formula, !!.sample, !!.cell_type, !!.count,
        variance_association = variance_association,
        truncation_ajustment = 1.1,
        approximate_posterior_inference = approximate_posterior_inference == "all",
        formula_variability = formula_variability
      )

    # Pior
    data_for_model$prior_prec_intercept = prior_mean_variable_association$intercept
    data_for_model$prior_prec_slope  = prior_mean_variable_association$slope
    data_for_model$prior_prec_sd = prior_mean_variable_association$standard_deviation
    data_for_model$exclude_priors = exclude_priors

    # Check that design matrix is not too big
    if(ncol(data_for_model$X)>20)
      warning("sccomp says: the design matrix has more than 20 columns. Possibly some numerical covariates are erroneously of type character/factor.")

    fit =
      data_for_model %>%

      # Run the first discovery phase with permissive false discovery rate
      fit_model(
        stanmodels$glm_multi_beta_binomial,
        cores= cores,
        quantile = CI,
        approximate_posterior_inference = approximate_posterior_inference == "all",
        verbose = verbose,
        seed = seed,
        max_sampling_iterations = max_sampling_iterations,
        pars = c("beta", "alpha", "prec_coeff","prec_sd",   "alpha_normalised")
      )

    list(
      fit = fit,
      data_for_model = data_for_model,
      truncation_df2 =  .data
    )

  }

  # If we are checking outliers
  else{

    message("sccomp says: outlier identification first pass - step 1/3 [ETA: ~20s]")

    # Force variance NOT associated with mean for stringency of outlier detection
    data_for_model =
      .data %>%
      data_to_spread ( formula, !!.sample, !!.cell_type, !!.count) %>%
      data_spread_to_model_input(
        formula, !!.sample, !!.cell_type, !!.count,
        variance_association = FALSE,
        truncation_ajustment = 1.1,
        approximate_posterior_inference = approximate_posterior_inference %in% c("outlier_detection", "all"),
        formula_variability = ~1
      )

    # Pior
    data_for_model$prior_prec_intercept = prior_mean_variable_association$intercept
    data_for_model$prior_prec_slope  = prior_mean_variable_association$slope
    data_for_model$prior_prec_sd = prior_mean_variable_association$standard_deviation
    data_for_model$exclude_priors = exclude_priors


    fit =
      data_for_model %>%

      # Run the first discovery phase with permissive false discovery rate
      fit_model(
        stanmodels$glm_multi_beta_binomial,
        cores= cores,
        quantile = CI,
        approximate_posterior_inference = approximate_posterior_inference %in% c("outlier_detection", "all"),
        verbose = verbose,
        seed = seed,
        max_sampling_iterations = max_sampling_iterations,
        pars = c("beta", "alpha", "prec_coeff","prec_sd",   "alpha_normalised")
      )

    rng =  rstan::gqs(
      stanmodels$glm_multi_beta_binomial_generate_date,
      #rstan::stan_model("inst/stan/glm_multi_beta_binomial_generate_date.stan"),
      draws =  as.matrix(fit),
      data = data_for_model
    )

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

    # Allow variance association
    data_for_model =
      .data %>%
      data_to_spread ( formula, !!.sample, !!.cell_type, !!.count) %>%
      data_spread_to_model_input(
        formula, !!.sample, !!.cell_type, !!.count,
        variance_association = variance_association,
        truncation_ajustment = 1.1,
        approximate_posterior_inference = approximate_posterior_inference %in% c("outlier_detection", "all"),
        formula_variability = formula_variability
      )

    # Pior
    data_for_model$prior_prec_intercept = prior_mean_variable_association$intercept
    data_for_model$prior_prec_slope  = prior_mean_variable_association$slope
    data_for_model$prior_prec_sd = prior_mean_variable_association$standard_deviation
    data_for_model$exclude_priors = exclude_priors

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
      fit_model(
        stanmodels$glm_multi_beta_binomial,
        cores = cores,
        quantile = my_quantile_step_2,
        approximate_posterior_inference = approximate_posterior_inference %in% c("outlier_detection", "all"),
        verbose = verbose,
        seed = seed,
        max_sampling_iterations = max_sampling_iterations,
        pars = c("beta", "alpha", "prec_coeff", "prec_sd",   "alpha_normalised")
      )

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
      fit_model(
        stanmodels$glm_multi_beta_binomial,
        cores = cores,
        quantile = CI,
        approximate_posterior_inference = approximate_posterior_inference %in% c("all"),
        verbose = verbose, seed = seed,
        max_sampling_iterations = max_sampling_iterations,
        pars = c("beta", "alpha", "prec_coeff","prec_sd",   "alpha_normalised")
      )

    #fit_model(stan_model("inst/stan/glm_multi_beta_binomial.stan"), chains= 4, output_samples = 500)

    list(
      fit = fit3,
      data_for_model = data_for_model,
      truncation_df2 =
        truncation_df2
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
                                                    .cell_type,
                                                    .count,
                                                    fit,
                                                    data_for_model,
                                                    percent_false_positive,
                                                    check_outliers,
                                                    truncation_df2 = NULL,
                                                    variance_association = FALSE,
                                                    test_composition_above_logit_fold_change ) {

  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  do_test = ncol(data_for_model$X) > 1

  # parsed_beta =
  #   fit %>%
  #   summary_to_tibble("beta", "C", "M") %>%
  #   left_join(tibble(C=seq_len(ncol(data_for_model$X)), C_name = colnames(data_for_model$X)), by = "C") %>%
  #   nest(beta_posterior_1 = -M)



  # Parse fit
  false_positive_rate = percent_false_positive/100
  factor_of_interest = data_for_model$X %>% colnames() %>% .[-1]
  median_factor_of_interest = sprintf(".median_%s", factor_of_interest)
  effect_column_name = sprintf("composition_effect_%s", factor_of_interest)

  beta_CI =
    fit %>%
    summary_to_tibble("beta", "C", "M", probs = c(false_positive_rate/2,  0.5,  1-(false_positive_rate/2))) %>%

    # Drop columns I dont need
    select(1, 2, 3, 7, 8, 9) %>%

    # Rename column to match %
    rename(
      .lower := !!as.symbol(sprintf("%s%s", false_positive_rate/2*100, "%")) ,
      .median := !!as.symbol(sprintf("%s%s", "50", "%")) ,
      .upper := !!as.symbol(sprintf("%s%s", (1-(false_positive_rate/2))*100, "%")) ,
    ) %>%
    left_join(tibble(C=seq_len(ncol(data_for_model$X)), C_name = colnames(data_for_model$X)), by = "C") %>%
    select(-C) %>%
    pivot_wider(names_from = C_name, values_from=c(.lower , .median ,  .upper)) %>%

    # Create main effect if exists
    when(
      !is.na(factor_of_interest) %>% any() ~ nest(., composition_CI = -c(M, one_of(median_factor_of_interest))) %>%
        setNames(gsub(".median", "composition_effect", colnames(.))),
      ~ nest(., composition_CI = -c(M))
    )

  if(variance_association) {
    variability_effect_column_name = sprintf("variability_effect_%s", factor_of_interest) %>% as.symbol()


    alpha_CI =
      fit %>%
      summary_to_tibble("alpha_normalised", "C", "M", probs = c(false_positive_rate/2,  0.5,  1-(false_positive_rate/2))) %>%

      # Drop columns I dont need
      select(1, 2, 3, 7, 8, 9) %>%

      # Rename column to match %
      rename(
        .lower := !!as.symbol(sprintf("%s%s", false_positive_rate/2*100, "%")) ,
        .median := !!as.symbol(sprintf("%s%s", "50", "%")) ,
        .upper := !!as.symbol(sprintf("%s%s", (1-(false_positive_rate/2))*100, "%")) ,
      ) %>%
      left_join(tibble(C=seq_len(ncol(data_for_model$X)), C_name = colnames(data_for_model$X)), by = "C") %>%
      select(-C) %>%
      pivot_wider(names_from = C_name, values_from=c(.lower , .median ,  .upper)) %>%

      # Create main effect if exists
      when(
        !is.na(factor_of_interest) %>% any() ~ nest(., variability_CI = -c(M, one_of(median_factor_of_interest))) %>%
          setNames(gsub(".median", "variability_effect", colnames(.))),
        ~ nest(., variability_CI = -c(M))
      )
  }

  beta_CI %>%

    # Add probability if do_test
    when(
      do_test ~ left_join(
        .,
        get_probability_non_zero(fit, "beta", prefix="composition", test_above_logit_fold_change = test_composition_above_logit_fold_change) %>%
          setNames(c("M", sprintf("composition_pH0_%s", factor_of_interest))),
        by="M"
        ),
      ~ (.)
    ) %>%

    # Add ALPHA
    when(do_test & variance_association ~ left_join(.,  alpha_CI ,  by="M"), ~(.)) %>%

    # ADD CI alpha
    when(
      do_test & variance_association ~ left_join(
        .,
        get_probability_non_zero(fit, "alpha_normalised", prefix="variability", test_above_logit_fold_change = 0) %>%
          setNames(c("M", sprintf("variability_pH0_%s", factor_of_interest))),
        by="M"
      ),
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
                                   seed = sample(1e5, 1),
                                   verbose = FALSE,
                                   exclude_priors = FALSE,
                                   test_composition_above_logit_fold_change,
                                   max_sampling_iterations = 20000,
                                   pass_fit = TRUE
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
      verbose = verbose,
      exclude_priors = exclude_priors,
      max_sampling_iterations = max_sampling_iterations
    )


  hypothesis_test_multi_beta_binomial_glm(
    .sample = !!.sample,
    .cell_type = !!.cell_type,
    .count = !!.count,
    result_list$fit,
    result_list$data_for_model,
    percent_false_positive,
    result_list$truncation_df2,
    variance_association = variance_association,
    test_composition_above_logit_fold_change = test_composition_above_logit_fold_change
  ) %>%

    # Join the precision
    left_join(get_mean_precision(result_list$fit, result_list$data_for_model), by="M") %>%


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
    add_attr(get_mean_precision_association(result_list$fit), "mean_concentration_association") %>%
    when(pass_fit ~ add_attr(., result_list$fit, "fit"), ~ (.)) %>%
    add_attr(result_list$data_for_model, "model_input")

}

#' @importFrom stats model.matrix
# glm_multi_beta_binomial = function(input_df, formula, .sample){
#
#   covariate_names = parse_formula(formula)
#   .sample = enquo(.sample)
#
#   sampling(stanmodels$glm_multi_beta_binomial,
#            data = list(
#              N = input_df %>% nrow(),
#              M = input_df %>% select(-!!.sample, -covariate_names, -exposure) %>% ncol(),
#              exposure = input_df$exposure,
#              y = input_df %>% select(-covariate_names, -exposure) %>% as_matrix(rownames = !!.sample),
#              X = input_df %>% select(!!.sample, covariate_names) %>% model.matrix(formula, data=.),
#              C = length(covariate_names)
#            ),
#            cores = 4
#   )
#
# }

get_mean_precision = function(fit, data_for_model){
  fit %>%
    summary_to_tibble("alpha", "C", "M") %>%

    # Just grub intercept of alpha
    filter(C==1) %>%
    select( M, mean, `2.5%` , `97.5%`) %>%
    nest(concentration = -M)

  # WRONG ATTEMPT TO PLOT ALPHA FROM MULTIPLE COVARIATES
  # fit %>%
  #   draws_to_tibble_x_y("alpha", "C", "M") %>%
  #   nest(data = -c(M)) %>%
  #
  #   # Add Design matrix
  #   left_join(
  #     as_tibble(data_for_model$X, rownames="M") %>%
  #       nest(X = -M) %>%
  #       mutate(M = as.integer(M)),
  #     by="M"
  #   ) %>%
  #   mutate(concentration = map2(
  #     data, X,
  #     ~ .x %>%
  #           select(.draw, C, .value) %>%
  #           spread(C, .value) %>%
  #           as_matrix(rownames=".draw") %*%
  #
  #         ( .y %>% as_matrix() %>% t() ) %>%
  #       quantile(c(0.025, 0.5, 0.975)) %>%
  #       enframe() %>%
  #       spread(name, value) %>%
  #       rename(mean = `50%`)
  #   )) %>%
  #   select(-data, -X)
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
