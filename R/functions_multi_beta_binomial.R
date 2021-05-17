#' multi_beta_binomial main
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
                                   cores = detect_cores(), # For development purpose,
                                   seed = sample(1:99999, size = 1)
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
    data_spread_to_model_input(formula, !!.sample, !!.cell_type, !!.count)




  if(!check_outliers){

    data_for_model %>%
      # Run the first discovery phase with permissive false discovery rate
      fit_model(stanmodels$glm_multi_beta_binomial, chains= 4) %>%
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

    fit =
      data_for_model %>%
      # Run the first discovery phase with permissive false discovery rate
       fit_model(stanmodels$glm_multi_beta_binomial, chains= 4, output_samples = 500)
      #fit_model(stan_model("inst/stan/glm_multi_beta_binomial.stan"), chains= 4, output_samples = 500)

    rng =  rstan::gqs(
      #stanmodels$generated_quantities,
      rstan::stan_model("inst/stan/glm_multi_beta_binomial_generate_date.stan"),
      draws =  as.matrix(fit),
      data = data_for_model
    )

    # Detect outliers
    truncation_df =
      .data %>%
      left_join(
        summary_to_tibble(rng, "counts", "N", "M") %>%
          nest(data = -N) %>%
          mutate(!!.sample := rownames(data_for_model$y)) %>%
          unnest(data) %>%
          nest(data = -M) %>%
          mutate(!!.cell_type := colnames(data_for_model$y)) %>%
          unnest(data) ,

        by = c("sample", "cell_type")
      ) %>%
      mutate(
        truncation_down = if_else(!!.count > `2.5%` & !!.count < `97.5%`, `2.5%`, -1),
        truncation_up = if_else(!!.count > `2.5%` & !!.count < `97.5%`, `97.5%`, -1)
      )

    # Add censoring
    data_for_model$is_truncated = 1
    data_for_model$truncation_up = truncation_df %>% select(N, M, truncation_up) %>% spread(M, truncation_up) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)
    data_for_model$truncation_down = truncation_df %>% select(N, M, truncation_down) %>% spread(M, truncation_down) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)

    fit_2 =
      data_for_model %>%
      # Run the first discovery phase with permissive false discovery rate
     # fit_model(stanmodels$glm_multi_beta_binomial, chains= 4, output_samples = 500)
    fit_model(stan_model("inst/stan/glm_multi_beta_binomial.stan"), chains= 4, output_samples = 500)

    rng =  rstan::gqs(
      #stanmodels$generated_quantities,
      rstan::stan_model("inst/stan/glm_multi_beta_binomial_generate_date.stan"),
      draws =  as.matrix(fit),
      data = data_for_model
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
             y = input_df %>% select(-covariate_names, -exposure) %>% nanny::as_matrix(rownames = !!.sample),
             X = input_df %>% select(!!.sample, covariate_names) %>% model.matrix(formula, data=.),
             C = length(covariate_names)
           ),
           cores = 4
  )

}
