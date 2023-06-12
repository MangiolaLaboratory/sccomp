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
#' @importFrom purrr map_chr
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom purrr map
#' @importFrom tibble rowid_to_column
#' @importFrom purrr map_lgl
#' @importFrom dplyr case_when
#' @importFrom rlang :=
#' @importFrom rlang quo_is_symbolic
#' @importFrom readr write_file
#'
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | factor columns | Pvaue column | a significance column
#' @param formula_composition A formula. The sample formula used to perform the differential cell_group abundance analysis
#' @param formula_variability A formula. The sample formula used to perform the differential cell_group variability analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_type identifier
#' @param .count A column name as symbol. The cell_type abundance (read count)
#'
#' @param prior_mean_variable_association A list of the form list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(0.4527292, 0.3318759)). Where for each parameter, we specify mean and standard deviation. This is used to incorporate prior knowledge about the mean/variability association of cell-type proportions.
#' @param percent_false_positive A real between 0 and 100. It is the aimed percent of cell types being a false positive. For example, percent_false_positive_genes = 1 provide 1 percent of the calls for significant changes that are actually not significant.
#' @param check_outliers A boolean. Whether to check for outliers before the fit.
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param enable_loo A boolean. Enable model comparison by the R package LOO. This is helpful when you want to compare the fit between two models, for example, analogously to ANOVA, between a one factor model versus a interceot-only model.
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
                                            formula_composition = ~ 1,
                                            formula_variability = ~ 1,
                                            .sample,
                                            .cell_group,
                                            .count,

                                            # Secondary parameters
                                            contrasts = NULL,
                                            #.grouping_for_random_intercept = NULL,
                                            prior_mean_variable_association,
                                            percent_false_positive = 5,
                                            check_outliers = FALSE,
                                            approximate_posterior_inference = "all",
                                            enable_loo = FALSE,
                                            cores = detectCores(), # For development purpose,
                                            seed = sample(1e5, 1),
                                            verbose = FALSE,
                                            exclude_priors = FALSE,
                                            bimodal_mean_variability_association = bimodal_mean_variability_association,
                                            use_data = use_data,
                                            max_sampling_iterations = 20000
) {
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)

  # Credible interval
  CI = 1 - (percent_false_positive/100)

  # Produce data list
  factor_names = parse_formula(formula_composition)

  # Random intercept
  random_intercept_elements = parse_formula_random_intercept(formula_composition)

  # If no random intercept fake it
  if(nrow(random_intercept_elements)>0){
    .grouping_for_random_intercept = random_intercept_elements |> pull(grouping) |> unique() |>   map(~ quo(!! sym(.x)))

  } else{
    .grouping_for_random_intercept = list(quo(!! sym("random_intercept")))
    .data = .data |> mutate(!!.grouping_for_random_intercept[[1]] := "1")
  }

  # rstan_options(threads_per_chain = floor(cores/chains))
  # Load model
  if(file.exists("glm_multi_beta_binomial_cmdstanr.rds"))
    model = readRDS("glm_multi_beta_binomial_cmdstanr.rds")
  else {
    write_file(glm_multi_beta_binomial, "glm_multi_beta_binomial_cmdstanr.stan")
    model = cmdstan_model( "glm_multi_beta_binomial_cmdstanr.stan", cpp_options = list(stan_threads = TRUE) )
    model  %>% saveRDS("glm_multi_beta_binomial_cmdstanr.rds")
  }
  
  # Original - old
  # prec_sd ~ normal(0,2);
  # prec_coeff ~ normal(0,5);

  # If we are NOT checking outliers
  if(!check_outliers){

    message("sccomp says: estimation")

    data_for_model =
      .data %>%
      data_to_spread ( formula_composition, !!.sample, !!.cell_group, !!.count, .grouping_for_random_intercept) %>%
      data_spread_to_model_input(
        formula_composition, !!.sample, !!.cell_group, !!.count,
        truncation_ajustment = 1.1,
        approximate_posterior_inference = approximate_posterior_inference == "all",
        formula_variability = formula_variability,
        contrasts = contrasts,
        bimodal_mean_variability_association = bimodal_mean_variability_association,
        use_data = use_data,
        random_intercept_elements
      )

    # Print design matrix
    message(sprintf("sccomp says: the composition design matrix has columns: %s", data_for_model$X %>% colnames %>% paste(collapse=", ")))
    message(sprintf("sccomp says: the variability design matrix has columns: %s", data_for_model$Xa %>% colnames %>% paste(collapse=", ")))

    # Pior
    data_for_model$prior_prec_intercept = prior_mean_variable_association$intercept
    data_for_model$prior_prec_slope  = prior_mean_variable_association$slope
    data_for_model$prior_prec_sd = prior_mean_variable_association$standard_deviation
    data_for_model$exclude_priors = exclude_priors
    data_for_model$enable_loo = TRUE & enable_loo

    # # Check that design matrix is not too big
    # if(ncol(data_for_model$X)>20)
    #   message("sccomp says: the design matrix has more than 20 columns. Possibly some numerical factors are erroneously of type character/factor.")

    fit =
      data_for_model %>%

      # Run the first discovery phase with permissive false discovery rate
      fit_model(
        model,
        cores= cores,
        quantile = CI,
        approximate_posterior_inference = approximate_posterior_inference == "all",
        verbose = verbose,
        seed = seed,
        max_sampling_iterations = max_sampling_iterations,
        pars = c("beta", "alpha", "prec_coeff","prec_sd",   "alpha_normalised", "beta_random_intercept", "log_lik")
      )

    list(
      fit = fit,
      data_for_model = data_for_model,
      truncation_df2 =  .data
    )

  }

  # If we are checking outliers
  else{

    message("sccomp says: outlier identification first pass - step 1/3")

    # Force variance NOT associated with mean for stringency of outlier detection
    data_for_model =
      .data %>%
      data_to_spread ( formula_composition, !!.sample, !!.cell_group, !!.count, .grouping_for_random_intercept) %>%
      data_spread_to_model_input(
        formula_composition, !!.sample, !!.cell_group, !!.count,
        truncation_ajustment = 1.1,
        approximate_posterior_inference = approximate_posterior_inference %in% c("outlier_detection", "all"),
        formula_variability = ~1,
        contrasts = contrasts,
        bimodal_mean_variability_association = bimodal_mean_variability_association,
        use_data = use_data,
        random_intercept_elements
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
        model,
        cores= cores,
        quantile = CI,
        approximate_posterior_inference = approximate_posterior_inference %in% c("outlier_detection", "all"),
        verbose = verbose,
        seed = seed,
        max_sampling_iterations = max_sampling_iterations,
        pars = c("beta", "alpha", "prec_coeff","prec_sd",   "alpha_normalised", "beta_random_intercept")
      )

    # Load model
    if(file.exists("glm_multi_beta_binomial_generate_cmdstanr.rds"))
    	mod_rng = readRDS("glm_multi_beta_binomial_generate_cmdstanr.rds")
    else {
    	write_file(glm_multi_beta_binomial_generate, "glm_multi_beta_binomial_generate_cmdstanr.stan")
    	mod_rng = cmdstan_model( "glm_multi_beta_binomial_generate_cmdstanr.stan" )
    	mod_rng  %>% saveRDS("glm_multi_beta_binomial_generate_cmdstanr.rds")
    }
    
    rng = mod_rng$generate_quantities(
    	fit,

      # This is for the new data generation with selected factors to do adjustment
      data = data_for_model |> c(list(

        # Add subset of coefficients
        length_X_which = ncol(data_for_model$X),
        length_XA_which = ncol(data_for_model$XA),
        X_which = seq_len(ncol(data_for_model$X)) |> as.array(),
        XA_which = seq_len(ncol(data_for_model$Xa)) |> as.array(),

        # Random intercept
        length_X_random_intercept_which = ncol(data_for_model$X_random_intercept),
        X_random_intercept_which = seq_len(ncol(data_for_model$X_random_intercept)) |> as.array(),
        create_intercept = FALSE
      )),
    	parallel_chains = ifelse(data_for_model$is_vb, 1, fit$num_chains())
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
          mutate(!!.cell_group := colnames(data_for_model$y)) %>%
          unnest(data) ,

        by = c(quo_name(.sample), quo_name(.cell_group))
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
      data_to_spread ( formula_composition, !!.sample, !!.cell_group, !!.count, .grouping_for_random_intercept) %>%
      data_spread_to_model_input(
        formula_composition, !!.sample, !!.cell_group, !!.count,
        truncation_ajustment = 1.1,
        approximate_posterior_inference = approximate_posterior_inference %in% c("outlier_detection", "all"),
        formula_variability = formula_variability,
        contrasts = contrasts,
        bimodal_mean_variability_association = bimodal_mean_variability_association,
        use_data = use_data,
        random_intercept_elements
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
    data_for_model$truncation_not_idx = (data_for_model$truncation_down >= 0) %>% t() %>% as.vector()  %>% which()
    data_for_model$TNS = length(data_for_model$truncation_not_idx)
    
    # Introduced with cmdstanr
    data_for_model$truncation_not_matrix = (data_for_model$truncation_down >= 0)
    data_for_model$truncation_df = 
      data_for_model$truncation_down |> 
      st(0) |> 
      as_tibble() |> 
      rowid_to_column("N") |>
      pivot_longer(-N, names_to = "M", values_to = "do_truncate") |> 
      mutate(N = as.integer(N), M = as.integer(M)) |> 
      arrange(N, M) |> 
      filter(do_truncate) |> 
      select(N, M) |>
      as.matrix()
    data_for_model$truncation_df_length = nrow(data_for_model$truncation_df)
    data_for_model$truncation_matrix_idx_length = data_for_model$truncation_not_matrix |> apply(1, function(x) x %>% `!` |> which()  |> length())
    data_for_model$truncation_not_matrix_idx_length = data_for_model$truncation_not_matrix |> apply(1, function(x) x |> which() |> length())
    

    message("sccomp says: outlier identification second pass - step 2/3")

    my_quantile_step_2 = 1 - (0.1 / data_for_model$N)

    # This step gets the credible interval to control for within-category false positive rate
    # We want a category-wise false positive rate of 0.1, and we have to correct for how many samples we have in each category
    CI_step_2 = (1-my_quantile_step_2) / 2 * 2


    fit2 =
      data_for_model %>%
      fit_model(
        model,
        cores = cores,
        quantile = my_quantile_step_2,
        approximate_posterior_inference = approximate_posterior_inference %in% c("outlier_detection", "all"),
        verbose = verbose,
        seed = seed,
        max_sampling_iterations = max_sampling_iterations,
        pars = c("beta", "alpha", "prec_coeff", "prec_sd",   "alpha_normalised", "beta_random_intercept")
      )

    #fit_model(stan_model("inst/stan/glm_multi_beta_binomial.stan"), chains= 4, output_samples = 500, approximate_posterior_inference = FALSE, verbose = TRUE)

    rng2 = mod_rng$generate_quantities(
      fit2,
      data = data_for_model |> c(list(

        # Add subset of coefficients
        length_X_which = ncol(data_for_model$X),
        length_XA_which = ncol(data_for_model$XA),
        X_which = seq_len(ncol(data_for_model$X)) |> as.array(),
        XA_which = seq_len(ncol(data_for_model$Xa)) |> as.array(),

        # Random intercept
        length_X_random_intercept_which = ncol(data_for_model$X_random_intercept),
        X_random_intercept_which = seq_len(ncol(data_for_model$X_random_intercept)) |> as.array(),
        create_intercept = FALSE

      )),
      parallel_chains = ifelse(data_for_model$is_vb, 1, fit$num_chains())
    )

    # Detect outliers
    truncation_df2 =
      .data %>%
      left_join(
        summary_to_tibble(rng2, "counts", "N", "M", probs = c(CI_step_2, 0.5, 1-CI_step_2)) %>%

          # !!! THIS COMMAND RELIES ON POSITION BECAUSE IT'S NOT TRIVIAL TO MATCH
          # !!! COLUMN NAMES BASED ON LIMITED PRECISION AND/OR PERIODICAL QUANTILES
          rename(
            .lower := !!as.symbol(colnames(.)[4]) ,
            .median = `50%`,
            .upper := !!as.symbol(colnames(.)[6])
          ) %>%
          nest(data = -N) %>%
          mutate(!!.sample := rownames(data_for_model$y)) %>%
          unnest(data) %>%
          nest(data = -M) %>%
          mutate(!!.cell_group := colnames(data_for_model$y)) %>%
          unnest(data) ,

        by = c(quo_name(.sample), quo_name(.cell_group))
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
        truncation_up = case_when(outlier ~ -1, TRUE ~ truncation_up)
      )


    data_for_model$truncation_up = truncation_df2 %>% select(N, M, truncation_up) %>% spread(M, truncation_up) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)
    data_for_model$truncation_down = truncation_df2 %>% select(N, M, truncation_down) %>% spread(M, truncation_down) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)
    data_for_model$enable_loo = TRUE & enable_loo

    message("sccomp says: outlier-free model fitting - step 3/3")

    # Print design matrix
    message(sprintf("sccomp says: the composition design matrix has columns: %s", data_for_model$X %>% colnames %>% paste(collapse=", ")))
    message(sprintf("sccomp says: the variability design matrix has columns: %s", data_for_model$Xa %>% colnames %>% paste(collapse=", ")))

    fit3 =
      data_for_model %>%
      # Run the first discovery phase with permissive false discovery rate
      fit_model(
        model,
        cores = cores,
        quantile = CI,
        approximate_posterior_inference = approximate_posterior_inference %in% c("all"),
        verbose = verbose, seed = seed,
        max_sampling_iterations = max_sampling_iterations,
        pars = c("beta", "alpha", "prec_coeff","prec_sd",   "alpha_normalised", "beta_random_intercept", "log_lik")
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
#' @importFrom rlang :=
#'
#' @keywords internal
#' @noRd
#'
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | factor columns | Pvaue column | a significance column
#' @param formula_composition A formula. The sample formula used to perform the differential cell_group abundance analysis
#' @param formula_variability A formula. The sample formula used to perform the differential cell_group variability analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_type identifier
#' @param .count A column name as symbol. The cell_type abundance (read count)
#' @param check_outliers A boolean. Whether to check for outliers before the fit.
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param enable_loo A boolean. Enable model comparison by the R package LOO. This is helpful when you want to compare the fit between two models, for example, analogously to ANOVA, between a one factor model versus a interceot-only model.
#' @param verbose A boolean. Prints progression.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param seed An integer. Used for development and testing purposes
#'
#' @return A nested tibble `tbl` with cell_type-wise information: `sample wise data` | plot | `ppc samples failed` | `exposure deleterious outliers`
#'
#'
multi_beta_binomial_glm = function(.data,
                                   formula_composition = ~ 1,
                                   formula_variability = ~1,
                                   .sample,
                                   .cell_group,
                                   .count,

                                   # Secondary parameters
                                   contrasts = NULL,
                                   #.grouping_for_random_intercept = NULL,
                                   prior_mean_variable_association,
                                   percent_false_positive = 5,
                                   check_outliers = FALSE,
                                   approximate_posterior_inference = TRUE,
                                   enable_loo = FALSE,
                                   cores = detectCores(), # For development purpose,
                                   seed = sample(1e5, 1),
                                   verbose = FALSE,
                                   exclude_priors = FALSE,
                                   bimodal_mean_variability_association = FALSE,
                                   use_data = TRUE,
                                   test_composition_above_logit_fold_change,
                                   max_sampling_iterations = 20000,
                                   pass_fit = TRUE
) {

  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)
  #.grouping_for_random_intercept = enquo(.grouping_for_random_intercept)
  #contrasts = contrasts |> enquo() |> quo_names()

  estimates_list =
    estimate_multi_beta_binomial_glm(
      .data = .data,
      formula_composition = formula_composition,
      .sample = !!.sample,
      .cell_group = !!.cell_group,
      .count = !!.count,
      formula_variability = formula_variability,
      contrasts = contrasts,
      #.grouping_for_random_intercept = !!.grouping_for_random_intercept,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      enable_loo = enable_loo,
      cores = cores, # For development purpose,
      seed = seed,
      verbose = verbose,
      exclude_priors = exclude_priors,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      use_data = use_data,
      max_sampling_iterations = max_sampling_iterations
    )

  # Create a dummy tibble
  tibble() |>
    # Attach association mean concentration
    add_attr(estimates_list$fit, "fit") %>%
    add_attr(estimates_list$data_for_model, "model_input") |>
    add_attr(estimates_list$truncation_df2, "truncation_df2") |>
    add_attr(.sample, ".sample") |>
    add_attr(.cell_group, ".cell_group") |>
    add_attr(.count, ".count") |>
    add_attr(check_outliers, "check_outliers") |>
    add_attr(formula_composition, "formula_composition") |>
    add_attr(formula_variability, "formula_variability") |>

    test_contrasts(
      contrasts = contrasts,
      percent_false_positive = percent_false_positive,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change
    )


}

#' @importFrom stats model.matrix
get_mean_precision = function(fit, data_for_model){
  fit %>%
    summary_to_tibble("alpha", "C", "M") %>%

    # Just grub intercept of alpha
    filter(C==1) %>%
    select( M, mean, `2.5%` , `97.5%`) %>%
    nest(concentration = -M)

  # WRONG ATTEMPT TO PLOT ALPHA FROM MULTIPLE factorS
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
  	fit$summary("prec_coeff")$mean ,
  	fit$summary("prec_sd")$mean
  )
}
