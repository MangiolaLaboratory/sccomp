#' multi_beta_glm main
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
#'
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | covariate columns | Pvaue column | a significance column
#' @param formula A formula. The sample formula used to perform the differential cell_type abundance analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param check_outliers A boolean. Whether to check for outliers before the fit.
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param seed An integer. Used for development and testing purposes
#'
#' @return A nested tibble `tbl` with cell_type-wise information: `sample wise data` | plot | `ppc samples failed` | `exposure deleterious outliers`
#'
#'
multi_beta_glm = function(.data,
                          formula = ~ 1,
                          .sample,
                          check_outliers = FALSE,
                          approximate_posterior_inference = TRUE,
                          cores = detect_cores(), # For development purpose,
                          seed = sample(1e5, 1)
) {
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  if(!check_outliers){

    .data %>%
      glm_multi_beta(formula, !!.sample)



  }

  else{


  }


}

#' @importFrom stats model.matrix
glm_multi_beta = function(input_df, formula, .sample){

  covariate_names = parse_formula(formula)
  .sample = enquo(.sample)

  # Load model
  if(file.exists("glm_multi_beta.rds"))
    model_glm_multi_beta = readRDS("glm_multi_beta.rds")
  else {
    model_glm_multi_beta = stan_model(model_code = model_glm_multi_beta)
    model_glm_multi_beta %>% saveRDS("glm_multi_beta.rds")

  }


  sampling(model_glm_multi_beta,
           data = list(
             N = input_df %>% nrow(),
             M = input_df %>% select(-!!.sample, -any_of(covariate_names)) %>% ncol(),
             y = input_df %>% select(-any_of(covariate_names)) %>% as_matrix(rownames = !!.sample),
             X = input_df %>% select(!!.sample, any_of(covariate_names)) %>% model.matrix(formula, data=.)
           ),
           cores = 4
  )

}
