#' dirichlet_multinomial_glm main
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
dirichlet_multinomial_glm = function(.data,
                                     formula = ~ 1,
                                     .sample,
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
        stanmodels$glm_dirichlet_multinomial,
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

    .data_1 =
      .data %>%
      fit_model_and_parse_out_no_missing_data(!!.count, formula, X, exposure, iteration = 1, chains = 4)

    .data_2 =
      .data_1 %>%
      select(-contains("posterior")) %>%
      fit_model_and_parse_out_missing_data(!!.count, formula, X, exposure, iteration = 2)

    .data_2 %>%

      # Join filtered
      mutate(significant = map_lgl(
        beta_quantiles_2,
        ~ .x$`2.5%` * .x$`97.5%` > 0
      )) %>%

      #Join unfiltered
      mutate(significant_pre_filtering = map_lgl(
        beta_quantiles_1,
        ~ .x$`2.5%` * .x$`97.5%` > 0
      )) %>%

      # Define outlier
      rename(outlier = outlier_2 ) %>%

      # Clean
      select(-N, -M, -contains("posterior"))

  }


}
