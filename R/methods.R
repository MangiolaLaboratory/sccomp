
#' sccomp_glm main
#'
#' @description This function runs the data modeling and statistical test for the hypothesis that a cell_type includes outlier biological replicate.
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
#'
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | covariate columns | Pvaue column | a significance column
#' @param formula A formula. The sample formula used to perform the differential cell_type abundance analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_type A column name as symbol. The cell_type identifier
#' @param .count A column name as symbol. The cell_type abunace (read count)
#' @param .significance A column name as symbol. A column with the Pvalue, or other significanc measure (preferred Pvalue over false discovery rate)
#' @param .do_check A column name as symbol. A column with a booean indicating whether a cell_type was identified as differentially abundant
#' @param percent_false_positive_genes A real between 0 and 100. It is the aimed percent of cell_type being a false positive. For example, percent_false_positive_genes = 1 provide 1 percent of the calls for outlier containing cell_types that has actually not outliers.
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param approximate_posterior_analysis A boolean. Whether the calculation of the credible intervals should be done semi-analitically, rather than with pure ampling from the posterior. It confers execution time and memory advantage.
#' @param how_many_negative_controls An integer. How many cell_type from the bottom non-significant should be taken for inferring the mean-overdispersion trend.
#' @param draws_after_tail An integer. How many draws should on average be after the tail, in a way to inform CI.
#' @param save_generated_quantities A boolean. Used for development and testing purposes
#' @param additional_parameters_to_save A character vector. Used for development and testing purposes
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param pass_fit A boolean. Used for development and testing purposes
#' @param do_check_only_on_detrimental A boolean. Whether to test only for detrimental outliers (same direction as the fold change). It allows to test for less cell_type/sample pairs and therefore higher the probability threshold.
#' @param tol_rel_obj A real. Used for development and testing purposes
#' @param just_discovery A boolean. Used for development and testing purposes
#' @param seed An integer. Used for development and testing purposes
#' @param adj_prob_theshold_2 A boolean. Used for development and testing purposes
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
														 percent_false_positive_genes = 1,
														 how_many_negative_controls = 500,

														 approximate_posterior_inference = T,
														 approximate_posterior_analysis = T,
														 draws_after_tail = 10,

														 save_generated_quantities = F,
														 additional_parameters_to_save = c(),  # For development purpose
														 #cores = detect_cores(), # For development purpose,
														 pass_fit = F,
														 do_check_only_on_detrimental = length(parse_formula(formula)) > 0,
														 tol_rel_obj = 0.01,
														 just_discovery = F,
														 seed = sample(1:99999, size = 1),
														 adj_prob_theshold_2 = NULL,
														 return_fit = FALSE
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
	check_if_count_integer(.data, !!.count)

	# Check that the dataset is squared
	#if(my_df %>% distinct(!!.sample, !!.cell_type) %>% count(!!.cell_type) %>% count(n) %>% nrow %>% `>` (1))
	#  stop("The input data frame does not represent a rectangular structure. Each cell_type must be present in all samples.")

	.data_parsed =
	  .data %>%
	  mutate(
	    N = as.factor(!!.sample) %>% as.integer,
	    M = as.factor(!!.cell_type) %>% as.integer
	   )

	.data_wide =
	  .data_parsed %>%
	  select(N, M, !!.count, parse_formula(formula)) %>%
	  distinct() %>%
	  spread(M, !!.count)

	# Create design matrix
	X =  model.matrix(object = formula,   data =.data_wide)
	C = X %>% ncol

	.data_wide_no_covariates = .data_wide %>% select(-parse_formula(formula))

	exposure =
	  .data_wide_no_covariates %>%
	  mutate(sumrow = rowSums(.[2:ncol(.)])) %>%
	  pull(sumrow)

	# Run the first discovery phase with permissive false discovery rate
	fit_discovery  =
	  .data_wide_no_covariates %>%
		do_inference(
			formula,
			approximate_posterior_inference,
			approximate_posterior_analysis = F,
			C,
			X,
			cores,
			additional_parameters_to_save,
			pass_fit = T,
			tol_rel_obj = tol_rel_obj,
			seed = seed,
			output_samples = 5000
		)


	beta_posterior =
	  fit_discovery %>%
	  draws_to_tibble_x_y("beta", "C", "M") %>%
	  filter(C==2) %>%
	  nest(data = -M) %>%
	  mutate(quantiles = map(
	    data,
	    ~ quantile(
	      .x$.value,
	      probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
	    ) %>%
	      enframe() %>%
	      spread(name, value)
	  )) %>%
	  unnest(quantiles)

	# # For building some figure I just need the discovery run, return prematurely
	# if(just_discovery) return(res_discovery %>% filter(.variable == "counts_rng"))


  # Generate theoretical data
	generated_discovery =
	  generate_quantities(
	    fit_discovery ,
	    nrow(.data_wide_no_covariates),
	    ncol(.data_wide_no_covariates)-1,
	    exposure = exposure
	   )

	# Integrate generated quantiles
	.data_parsed_outliers =
	  .data_parsed %>%
	  left_join(generated_discovery) %>%

	  # Add covariate from design
	  left_join(X %>% as.data.frame %>% rowid_to_column("N")) %>%

	  # join CI
	  mutate(outlier_above = !!.count > `97.5%`) %>%
	  mutate(outlier_below = !!.count < `2.5%`) %>%
	  #mutate(ppc = !!.count %>% between(`2.5%`, `97.5%`))  %>%

	  # Join slope
	  left_join(
	    fit_discovery %>%
	      summary_to_tibble("beta", "C", "M") %>%
	      #tidybayes::gather_draws(beta[C, M]) %>%
	      #median_qi() %>%
	      filter(C==2) %>%
	      select(M, slope = `50%`)
	  ) %>%

	  # Mark if on the right of the covariate scale
	  mutate(is_group_right = !!as.symbol(colnames(X)[2]) > mean( !!as.symbol(colnames(X)[2]) )) %>%

  	# Check if outlier might be deleterious for the statistics
  	mutate(
  	  deleterious_outliers =
  	    (outlier_above & slope > 0 & is_group_right)  |
  	    (outlier_below & slope > 0 & !is_group_right) |
  	    (outlier_above & slope < 0 & !is_group_right) |
  	    (outlier_below & slope < 0 & is_group_right)
  	)


	# Columns of counts to be ignored from the inference
	to_exclude =
	  .data_parsed_outliers %>%
		filter(deleterious_outliers) %>%
		distinct(N, M, `2.5%`, `97.5%`)

	to_include =
	  .data_parsed_outliers %>%
	  filter(!deleterious_outliers) %>%
	  distinct(N, M)


	# # Calculate how many potential non NB cell_type I should check
	# how_namy_to_exclude = to_exclude %>% nrow

	# # Get the credible intervals for which account in the truncated NB model
	# truncation_values =
	# 	res_discovery %>%
	# 	filter(`.variable` == "counts_rng") %>%
	# 	distinct(S, G, .lower, .upper) %>%
	# 	mutate(`.lower` = `.lower` %>% as.integer,
	# 				 `.upper` = `.upper` %>% as.integer)

	# # Get the inferred values from first model to possibly use them in the second model as priors
	# prior_from_discovery =
	# 	res_discovery %>%
	# 	filter(`.variable` != "counts_rng") %>%
	# 	select(`.variable`, S, G, mean, sd)

	.data_parsed_inliers =
	  .data_parsed %>%
	  anti_join(to_exclude,by = c("N", "M")) %>%
	  select(.value = count, N, M)

	# Dirichlet with missing data
	fit_imputation =
	  .data_wide_no_covariates %>%
	  do_inference_imputation(
	    formula,
	    approximate_posterior_inference,
	    approximate_posterior_analysis,
	    C,
	    X,
	    cores,
	    additional_parameters_to_save,
	    pass_fit = pass_fit,
	    to_include = to_include,
	    tol_rel_obj = tol_rel_obj,
	    #truncation_compensation = 0.7352941, # Taken by approximation study
	    seed = seed,
	    precision = fit_discovery %>% extract("precision") %$% precision,
	    exposure = exposure
	  )

	beta_posterior_corrected =
	  fit_imputation %>%
	  draws_to_tibble_x_y("counts", "N", "M") %>%
	  rename(.draw_imputation = .draw) %>%
	  nest(data = -c(.chain ,.iteration, .draw_imputation ,.variable)) %>%
	  sample_n(100) %>%
	  mutate(fit = future_map(
	    data,
	    ~ .x %>%
	      anti_join(to_include,by = c("N", "M")) %>%
	      bind_rows(.data_parsed_inliers) %>%
	      spread(M, .value) %>%

	      # Run model
	      do_inference(
	        formula,
	        approximate_posterior_inference,
	        approximate_posterior_analysis = F,
	        C,
	        X,
	        cores,
	        additional_parameters_to_save,
	        pass_fit = T,
	        tol_rel_obj = tol_rel_obj,
	        seed = seed,
	        output_samples = 50,
	        chains = 1
	      ) %>%
	      draws_to_tibble_x_y("beta", "C", "M") %>%
	      filter(C==2)

	  )) %>%
	  select(.draw_imputation, fit) %>%
	  unnest(fit) %>%
	  nest(data = -M) %>%
	  mutate(quantiles = map(
	    data,
	    ~ quantile(
	      .x$.value,
	      probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
	    ) %>%
	      enframe() %>%
	      spread(name, value)
	  )) %>%
	  unnest(quantiles)


	.data_parsed %>%

	  # Join filtered
	  left_join(
	    beta_posterior_corrected %>%
	      mutate(significant =  `2.5%` * `97.5%` > 0) %>%
	      select(-data) ,
	    by="M"
	  ) %>%

	  #Join unfiltered
	  left_join(
	    beta_posterior %>%
	      mutate(significant_pre_filtering =  `2.5%` * `97.5%` > 0) %>%
	      select(-data) %>%
	      nest(quantiles_pre_filtering = c(`2.5%` ,  `25%` ,  `50%`,    `75%`, `97.5%`)),
	    by="M"
	  ) %>%

	  # Attach outliers
	  left_join(
	    to_exclude %>%
	      select(N, M) %>%
	      mutate(outlier = TRUE),
	    by = c("N", "M")
	  ) %>%
	  mutate(outlier = if_else(is.na(outlier), FALSE, outlier)) %>%

	  # Clean
	  select(-N, -M)


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


