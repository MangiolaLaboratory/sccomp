
#' sccomp_glm main
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
sccomp_glm = function(.data,
                      formula = ~ 1,
                      .sample,
                      .cell_type,
                      .count,
                      noise_model = "multi_beta_binomial",
                      check_outliers = FALSE,
                      approximate_posterior_inference = T,
                      cores = detect_cores(), # For development purpose,
                      seed = sample(1:99999, size = 1)
) {
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  # Get factor of interest
  #factor_of_interest = ifelse(parse_formula(formula) %>% length %>% `>` (0), parse_formula(formula)[1], "")

  # Check if columns exist
  check_columns_exist(.data, !!.sample, !!.cell_type, !!.count)

  # Check if any column is NA or null
  check_if_any_NA(.data, !!.sample, !!.cell_type, !!.count, parse_formula(formula))

  # Check if count column is integer class
  # check_if_count_integer(.data, !!.count)

  # Check that the dataset is squared
  #if(my_df %>% distinct(!!.sample, !!.cell_type) %>% count(!!.cell_type) %>% count(n) %>% nrow %>% `>` (1))
  #  stop("The input data frame does not represent a rectangular structure. Each cell_type must be present in all samples.")

  noise_model %>%
    when(
      equals(., "multi_beta_binomial") ~ multi_beta_binomial_glm(
        .data = .data,
        formula = formula,
        .sample = !!.sample,
        .cell_type = !!.cell_type,
        .count = !!.count,
        check_outliers = check_outliers,
        approximate_posterior_inference = approximate_posterior_inference,
        cores = cores,
        # For development purpose,
        seed = seed
      ),

      equals(., "multi_beta") ~ multi_beta_glm(
        .data = .data,
        formula = formula,
        .sample = !!.sample,
        .cell_type = !!.cell_type,
        .count = !!.count,
        check_outliers = check_outliers,
        approximate_posterior_inference = approximate_posterior_inference,
        cores = cores,
        # For development purpose,
        seed = seed
      ),

      equals(., "dirichlet_multinomial") ~ dirichlet_multinomial_glm(
        .data = .data,
        formula = formula,
        .sample = !!.sample,
        .cell_type = !!.cell_type,
        .count = !!.count,
        check_outliers = check_outliers,
        approximate_posterior_inference = approximate_posterior_inference,
        cores = cores,
        # For development purpose,
        seed = seed
      )
    )



}


#' plot_credible interval for theoretical data distributions
#'
#' @description Plot the data along the theoretical data distribution.
#'
#' @importFrom tibble as_tibble
#'
#' @param .data The tibble returned by sccomp_glm
#'
#' @return A tibble with an additional `plot` column
#'
#' @export
#'
plot_credible_intervals = function(.data){

	.cell_type = .data %>% attr("cell_type_column")
	.count = .data %>% attr("abundance_column")
	.sample = .data %>% attr("sample_column")
	formula = .data %>% attr("formula") %>% as.formula()

	.data %>%

		# Create plots for every tested cell_type
		mutate(plot =
					 	pmap(
					 		list(
					 			`sample wise data`,
					 			.cell_type,
					 			# nested data for plot
					 			.count,
					 			# name of value column
					 			.sample,
					 			# name of sample column
					 			parse_formula(formula)[1] # main covariate
					 		),
					 		~ produce_plots(..1, ..2, ..3, ..4, ..5)
					 	))
}


