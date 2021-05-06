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
#' @return A nested tibble `tbl` with cell_type-wise information: `sample wise data` | plot | `ppc samples failed` | `tot deleterious outliers`
#'
#' @export
#'
dirichlet_multinomial_glm = function(.data,
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

  # Check that the dataset is squared
  #if(my_df %>% distinct(!!.sample, !!.cell_type) %>% count(!!.cell_type) %>% count(n) %>% nrow %>% `>` (1))
  #  stop("The input data frame does not represent a rectangular structure. Each cell_type must be present in all samples.")

  .data_parsed =
    .data %>%
    mutate(
      N = as.factor(!!.sample) %>% as.integer,
      M = as.factor(!!.cell_type) %>% as.integer
    )

  exposure =
    .data_parsed %>%
    group_by(N) %>%
    summarise(s = sum(count)) %>%
    arrange(N) %>%
    pull(s)

  # Create design matrix
  X =  model.matrix(object = formula,   data =
                      .data_parsed %>%
                      select(N, parse_formula(formula)) %>%
                      distinct() %>%
                      arrange(N)
  )

  if(!check_outliers){

    .data_parsed_0 = .data_parsed %>% count_in_beta(!!.count, formula, X, exposure, iteration = 1, chains = 4)

    .data_parsed_0 %>%

      # Join filtered
      mutate(significant = !!as.symbol(sprintf(".lower_%s", colnames(X)[2])) * !!as.symbol(sprintf(".upper_%s", colnames(X)[2])) > 0 ) %>%

      # Clean
      select(-N, -M, -contains("posterior"))

  }

  else{

    .data_parsed_1 =
      .data_parsed %>%
      count_in_beta_out_no_missing_data(!!.count, formula, X, exposure, iteration = 1, chains = 4)

    .data_parsed_2 =
      .data_parsed_1 %>%
      select(-contains("posterior")) %>%
      count_in_beta_out_missing_data(!!.count, formula, X, exposure, iteration = 2)

    .data_parsed_2 %>%

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


#' multi_beta_glm main
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
#' @return A nested tibble `tbl` with cell_type-wise information: `sample wise data` | plot | `ppc samples failed` | `tot deleterious outliers`
#'
#' @export
#'
multi_beta_glm = function(.data,
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

  .data_parsed =
    .data %>%
    mutate(
      N = as.factor(!!.sample) %>% as.integer,
      M = as.factor(!!.cell_type) %>% as.integer
    )


  # Create design matrix
  X =  model.matrix(object = formula,   data =
                      .data_parsed %>%
                      select(N, parse_formula(formula)) %>%
                      distinct() %>%
                      arrange(N)
  )

  if(!check_outliers){

    .data %>%
      select(!!.sample, !!.cell_type, !!.count, parse_formula(formula)) %>%
      spread(!!.cell_type, !!.count) %>%
      glm_multi_beta(formula, !!.sample)



  }

  else{


  }


}


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
#' @return A nested tibble `tbl` with cell_type-wise information: `sample wise data` | plot | `ppc samples failed` | `tot deleterious outliers`
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

  .data_parsed =
    .data %>%
    mutate(
      N = as.factor(!!.sample) %>% as.integer,
      M = as.factor(!!.cell_type) %>% as.integer
    )


  # Create design matrix
  X =  model.matrix(object = formula,   data =
                      .data_parsed %>%
                      select(N, parse_formula(formula)) %>%
                      distinct() %>%
                      arrange(N)
  )

  if(!check_outliers){

    .data %>%
      nest(data = -!!.sample) %>%
      mutate(sample_tot = map_int(data, ~ .x %>% pull(!!.count) %>% sum() )) %>%
      unnest(data) %>%
      select(!!.sample, !!.cell_type, sample_tot, !!.count, parse_formula(formula)) %>%
      spread(!!.cell_type, !!.count) %>%
      glm_multi_beta_binomial(formula, !!.sample)



  }

  else{


  }


}
