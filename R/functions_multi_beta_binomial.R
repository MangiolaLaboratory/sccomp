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

  # Produce data list
  covariate_names = parse_formula(formula)
  X = .data %>% select(!!.sample, covariate_names) %>% model.matrix(formula, data=.)
  cell_cluster_names = .data %>% select(-!!.sample, -covariate_names, -exposure) %>% colnames()

  data_for_model =
    list(
      N = .data %>% nrow(),
      M = .data %>% select(-!!.sample, -covariate_names, -exposure) %>% ncol(),
      exposure = .data$exposure,
      y = .data %>% select(-covariate_names, -exposure) %>% nanny::as_matrix(rownames = !!.sample),
      X = X,
      C = ncol(X)
    )


  if(!check_outliers){

    data_for_model %>%
      fit_model_and_parse(
        stanmodels$glm_multi_beta_binomial,
        chains = 4
      ) %>%

      # Join filtered
      mutate(
        significant =
          !!as.symbol(sprintf(".lower_%s", colnames(X)[2])) *
          !!as.symbol(sprintf(".upper_%s", colnames(X)[2])) > 0
      ) %>%

      # Clesn
      select(-M) %>%
      mutate(.cell = cell_cluster_names) %>%
      select(.cell, everything())


  }

  else{

    covariate_names = parse_formula(formula)

    sampling(stanmodels$glm_multi_beta_binomial,
             data = list(
               N = .data %>% nrow(),
               M = .data %>% select(-!!.sample, -covariate_names, -exposure) %>% ncol(),
               exposure = .data$exposure,
               y = .data %>% select(-covariate_names, -exposure) %>% nanny::as_matrix(rownames = !!.sample),
               X = .data %>% select(!!.sample, covariate_names) %>% model.matrix(formula, data=.),
               C = length(covariate_names)
             ),
             cores = 4
    )

    x =  rstan::gqs(
      #stanmodels$generated_quantities,
      rstan::stan_model("inst/stan/glm_multi_beta_binomial_generate_date.stan"),
      draws =  as.matrix(fit),
      data = list(
        N = .data %>% nrow(),
        M = .data %>% select(-!!.sample, -covariate_names, -exposure) %>% ncol(),
        exposure = .data$exposure,
        X = .data %>% select(!!.sample, covariate_names) %>% model.matrix(formula, data=.),
        C = length(covariate_names)
      )
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
