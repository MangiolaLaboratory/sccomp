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
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | covariate columns | Pvaue column | a significance column
#' @param formula A formula. The sample formula used to perform the differential cell_type abundance analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_type A column name as symbol. The cell_type identifier
#' @param .count A column name as symbol. The cell_type abundance (read count)
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param seed An integer. Used for development and testing purposes
#'
#' @return A nested tibble `tbl` with cell_type-wise information: `sample wise data` | plot | `ppc samples failed` | `exposure deleterious outliers`
#'
#' @export
#'
multi_beta_binomial_glm = function(.data,
                                   formula = ~ 1,
                                   .sample,
                                   .cell_type,
                                   .count,
                                   check_outliers = FALSE,
                                   approximate_posterior_inference = T,
                                   variance_association = F,
                                   cores = detectCores(), # For development purpose,
                                   seed = sample(1:99999, size = 1),
                                   verbose = FALSE
) {
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)


  # Produce data list
  covariate_names = parse_formula(formula)

  data_for_model =
    .data %>%
    data_to_spread ( formula, !!.sample, !!.cell_type, !!.count) %>%
    data_spread_to_model_input(
      formula, !!.sample, !!.cell_type, !!.count,
      variance_association = variance_association,
      truncation_ajustment = 1.1
    )

  if(!check_outliers){

    data_for_model %>%
      # Run the first discovery phase with permissive false discovery rate
      fit_model(stanmodels$glm_multi_beta_binomial, cores= cores,  quantile = 0.95,  approximate_posterior_inference = approximate_posterior_inference, verbose = verbose, seed = seed) %>%
      parse_fit(data_for_model, .) %>%
      beta_to_CI( ) %>%

      # Join filtered
      mutate(
        significant =
          !!as.symbol(sprintf(".lower_%s", colnames(data_for_model$X)[2])) *
          !!as.symbol(sprintf(".upper_%s", colnames(data_for_model$X)[2])) > 0
      ) %>%

      # Clesn
      select(-M) %>%
      mutate(!!.cell_type := data_for_model$y %>% colnames()) %>%
      select(!!.cell_type, everything())


  }

  else{

    message("sccomp says: outlier identification first pass - step 1/3 [ETA: ~20s]")

    fit =
      data_for_model %>%
      fit_model(stanmodels$glm_multi_beta_binomial, cores= cores,  quantile = 0.95,  approximate_posterior_inference = approximate_posterior_inference, verbose = verbose, seed = seed)
    #fit_model(stan_model("inst/stan/glm_multi_beta_binomial.stan"), chains= 4, output_samples = 500,  approximate_posterior_inference = approximate_posterior_inference)

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

    # Add censoring
    data_for_model$is_truncated = 1
    data_for_model$truncation_up = truncation_df %>% select(N, M, truncation_up) %>% spread(M, truncation_up) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)
    data_for_model$truncation_down = truncation_df %>% select(N, M, truncation_down) %>% spread(M, truncation_down) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)

    message("sccomp says: outlier identification second pass - step 2/3 [ETA: ~60s]")

    my_quantile_step_2 = 1 - (0.1 / data_for_model$N)
    CI_step_2 = (1-my_quantile_step_2) / 2 * 2


    fit2 =
      data_for_model %>%
      fit_model(stanmodels$glm_multi_beta_binomial, cores = cores, quantile = my_quantile_step_2,  approximate_posterior_inference = approximate_posterior_inference, verbose = verbose, seed = seed)
    #fit_model(stan_model("inst/stan/glm_multi_beta_binomial.stan"), chains= 4, output_samples = 500, approximate_posterior_inference = F, verbose = T)

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
        summary_to_tibble(rng2, "counts", "N", "M", probs = c(CI_step_2, 1-CI_step_2)) %>%
          rename(
            .lower := !!as.symbol(paste0(CI_step_2*100,"%")) ,
            .upper := !!as.symbol(paste0((1-CI_step_2)*100,"%"))
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

    my_quantile_step_3 = 0.95

    fit3 =
      data_for_model %>%
      # Run the first discovery phase with permissive false discovery rate
      fit_model(stanmodels$glm_multi_beta_binomial, cores = cores, quantile = my_quantile_step_3,  approximate_posterior_inference = approximate_posterior_inference, verbose = verbose, seed = seed)
    #fit_model(stan_model("inst/stan/glm_multi_beta_binomial.stan"), chains= 4, output_samples = 500)

    fit3 %>%
      parse_fit(data_for_model, .) %>%
      beta_to_CI( ) %>%

      # Join filtered
      mutate(
        significant =
          !!as.symbol(sprintf(".lower_%s", colnames(data_for_model$X)[2])) *
          !!as.symbol(sprintf(".upper_%s", colnames(data_for_model$X)[2])) > 0
      ) %>%

      # Clesn
      select(-M) %>%
      mutate(!!.cell_type := data_for_model$y %>% colnames()) %>%
      select(!!.cell_type, everything()) %>%

      # join outliers
      left_join(
        truncation_df2 %>%
          select(!!.sample, !!.cell_type, outlier,.lower, .upper) %>%
          nest(outliers = -!!.cell_type),
        by = quo_name(.cell_type)
      )

  }


}


#' @export
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
