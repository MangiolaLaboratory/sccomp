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
multi_beta_glm = function(.data,
                          formula = ~ 1,
                          .sample,
                          check_outliers = FALSE,
                          approximate_posterior_inference = T,
                          cores = detect_cores(), # For development purpose,
                          seed = sample(1:99999, size = 1)
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


#' @export
glm_multi_beta = function(input_df, formula, .sample){

  covariate_names = parse_formula(formula)
  .sample = enquo(.sample)

  sampling(stanmodels$glm_multi_beta,
           data = list(
             N = input_df %>% nrow(),
             M = input_df %>% select(-!!.sample, -covariate_names) %>% ncol(),
             y = input_df %>% select(-covariate_names) %>% as_matrix(rownames = !!.sample),
             X = input_df %>% select(!!.sample, covariate_names) %>% model.matrix(formula, data=.)
           ),
           cores = 4
  )

}
