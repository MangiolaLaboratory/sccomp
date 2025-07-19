#' Main Function for SCCOMP Estimate
#'
#' @description
#' The `sccomp_estimate` function performs linear modeling on a table of cell counts or proportions,
#' which includes a cell-group identifier, sample identifier, abundance (counts or proportions), and factors
#' (continuous or discrete). The user can define a linear model using an R formula,
#' where the first factor is the factor of interest. Alternatively, `sccomp` accepts
#' single-cell data containers (e.g., Seurat, SingleCellExperiment, cell metadata, or
#' group-size) and derives the count data from cell metadata.
#'
#' @import dplyr
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom rlang quo_is_null
#' @importFrom SingleCellExperiment colData
#' @importFrom parallel detectCores
#' @importFrom rlang inform
#' @importFrom lifecycle is_present
#' @importFrom lifecycle deprecate_warn
#' @importFrom rlang is_symbolic
#' @importFrom rlang quo_get_expr
#' @importFrom rlang enquo
#' @importFrom magrittr not
#'
#' @param .data A tibble including cell_group name column, sample name column,
#'              abundance column (counts or proportions), and factor columns.
#' @param formula_composition A formula describing the model for differential abundance.
#' @param formula_variability A formula describing the model for differential variability.
#' @param sample A column name as a character string for the sample identifier. Replaces the deprecated `.sample`.
#' @param cell_group A column name as a character string for the cell-group identifier. Replaces the deprecated `.cell_group`.
#' @param abundance A column name as a character string for the cell-group abundance, which can be counts (> 0) or proportions (between 0 and 1, summing to 1 across `cell_group`). Replaces the deprecated `.abundance` and `.count`.
#' @param cores Number of cores to use for parallel calculations.
#' @param bimodal_mean_variability_association Logical, whether to model mean-variability as bimodal.
#' @param prior_mean A list specifying prior knowledge about the mean distribution, including intercept and coefficients.
#' @param prior_overdispersion_mean_association A list specifying prior knowledge about mean/variability association.
#' @param percent_false_positive A real number between 0 and 100 for outlier identification.
#' @param inference_method Character string specifying the inference method to use ('pathfinder', 'hmc', or 'variational'). Replaces the deprecated `approximate_posterior_inference` and `variational_inference`.
#' @param .sample_cell_group_pairs_to_exclude A column name indicating sample/cell-group pairs to exclude.
#' @param output_directory A character string specifying the output directory for Stan draws.
#' @param verbose Logical, whether to print progression details.
#' @param enable_loo Logical, whether to enable model comparison using the LOO package.
#' @param noise_model A character string specifying the noise model (e.g., 'multi_beta_binomial').
#' @param exclude_priors Logical, whether to run a prior-free model.
#' @param use_data Logical, whether to run the model data-free.
#' @param mcmc_seed An integer seed for MCMC reproducibility.
#' @param max_sampling_iterations Integer to limit the maximum number of iterations for large datasets.
#' @param pass_fit Logical, whether to include the Stan fit as an attribute in the output.
#' @param sig_figs Number of significant figures to use for Stan model output. Default is 9.
#' @param .count **DEPRECATED**. Use `abundance` instead.
#' @param approximate_posterior_inference **DEPRECATED**. Use `inference_method` instead.
#' @param variational_inference **DEPRECATED**. Use `inference_method` instead.
#' @param .sample **DEPRECATED**. Use `sample` instead.
#' @param .cell_group **DEPRECATED**. Use `cell_group` instead.
#' @param .abundance **DEPRECATED**. Use `abundance` instead.
#' @param ... Additional arguments passed to the `cmdstanr::sample` function.
#'
#' @return A tibble (`tbl`), with the following columns:
#' \itemize{
#'   \item cell_group - The cell groups being tested.
#'   \item parameter - The parameter being estimated from the design matrix described by the input formula_composition and formula_variability.
#'   \item factor - The covariate factor in the formula, if applicable (e.g., not present for Intercept or contrasts).
#'   \item c_lower - Lower (2.5%) quantile of the posterior distribution for a composition (c) parameter.
#'   \item c_effect - Mean of the posterior distribution for a composition (c) parameter.
#'   \item c_upper - Upper (97.5%) quantile of the posterior distribution for a composition (c) parameter.
#'   \item c_pH0 - Probability of the c_effect being smaller or bigger than the `test_composition_above_logit_fold_change` argument.
#'   \item c_FDR - False discovery rate of the c_effect being smaller or bigger than the `test_composition_above_logit_fold_change` argument. False discovery rate for Bayesian models is calculated differently from frequentists models, as detailed in Mangiola et al, PNAS 2023. 
#'   \item c_n_eff - Effective sample size, the number of independent draws in the sample. The higher, the better.
#'   \item c_R_k_hat - R statistic, a measure of chain equilibrium, should be within 0.05 of 1.0.
#'   \item v_lower - Lower (2.5%) quantile of the posterior distribution for a variability (v) parameter.
#'   \item v_effect - Mean of the posterior distribution for a variability (v) parameter.
#'   \item v_upper - Upper (97.5%) quantile of the posterior distribution for a variability (v) parameter.
#'   \item v_pH0 - Probability of the v_effect being smaller or bigger than the `test_composition_above_logit_fold_change` argument.
#'   \item v_FDR - False discovery rate of the v_effect being smaller or bigger than the `test_composition_above_logit_fold_change` argument. False discovery rate for Bayesian models is calculated differently from frequentists models, as detailed in Mangiola et al, PNAS 2023. 
#'   \item v_n_eff - Effective sample size for a variability (v) parameter.
#'   \item v_R_k_hat - R statistic for a variability (v) parameter, a measure of chain equilibrium.
#' }
#'
#' The function also attaches several attributes to the result:
#' \itemize{
#'   \item count_data - The original count data used in the analysis, stored as an attribute for efficient access.
#'   \item model_input - The model input data used for fitting.
#'   \item formula_composition - The formula used for composition modeling.
#'   \item formula_variability - The formula used for variability modeling.
#'   \item fit - The Stan fit object (if pass_fit = TRUE).
#' }
#'
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
#'
#' @examples
#'
#' print("cmdstanr is needed to run this example.")
#' # Note: Before running the example, ensure that the 'cmdstanr' package is installed:
#' # install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))
#'
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists()) {
#'     data("counts_obj")
#'
#'     estimate <- sccomp_estimate(
#'       counts_obj,
#'       ~ type,
#'       ~1,
#'       "sample",
#'       "cell_group",
#'       "count",
#'       cores = 1
#'     )
#'     
#'    # Note! 
#'    # If counts are available, do not use proportion.
#'    # Using proportion ignores the high uncertainty of low counts
#'    
#'    estimate_proportion <- sccomp_estimate(
#'       counts_obj,
#'       ~ type,
#'       ~1,
#'       "sample",
#'       "cell_group",
#'       "proportion",
#'       cores = 1
#'     )
#'     
#'   }
#' }
#'
#' @export
sccomp_estimate <- function(.data,
                            formula_composition = ~1,
                            formula_variability = ~1,
                            
                            sample,
                            cell_group,
                            abundance = NULL,
                            
                            # Secondary arguments
                            cores = detectCores(),
                            bimodal_mean_variability_association = FALSE,
                            percent_false_positive = 5,
                            inference_method = "pathfinder",
                            prior_mean = list(intercept = c(0, 1), coefficients = c(0, 1)),
                            prior_overdispersion_mean_association = list(
                              intercept = c(5, 2),
                              slope = c(0, 0.6),
                              standard_deviation = c(10, 20)
                            ),
                            .sample_cell_group_pairs_to_exclude = NULL,
                            output_directory = "sccomp_draws_files",
                            verbose = TRUE,
                            enable_loo = FALSE,
                            noise_model = "multi_beta_binomial",
                            exclude_priors = FALSE,
                            use_data = TRUE,
                            mcmc_seed = sample_seed(),
                            max_sampling_iterations = 20000,
                            pass_fit = TRUE,
                            sig_figs = 9,
                            ...,
                            
                            # DEPRECATED
                            .count = NULL,
                            approximate_posterior_inference = NULL,
                            variational_inference = NULL,
                            .sample = NULL,
                            .cell_group = NULL,
                            .abundance = NULL) {
  
  # Run the function
  check_and_install_cmdstanr()
  
  UseMethod("sccomp_estimate", .data)
}

#' @importFrom purrr map2_lgl
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr spread
#' @importFrom stats C
#' @importFrom rlang :=
#' @importFrom tibble enframe
#' @importFrom purrr map_int
data_to_spread = function(.data, formula, .sample, .cell_group, .count, .grouping_for_random_effect){
  
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)

  is_proportion = .data |> pull(!!.count) |> max() <= 1
  
  .data = 
    .data |>
    nest(data = -!!.sample) 
  
  # If proportions exposure = 1
  if(is_proportion) .data = .data |> mutate(exposure = 1)
  else
    .data = 
    .data |>
    mutate(exposure = map_int(data, ~ .x |> pull(!!.count) |> sum() )) 
  
  .data_to_spread = 
    .data |>
    unnest(data) |>
    select(!!.sample, !!.cell_group, exposure, !!.count, parse_formula(formula), any_of(.grouping_for_random_effect)) 
  
  # Check if duplicated samples
  if(
    .data_to_spread |> distinct(!!.sample, !!.cell_group) |> nrow() <
    .data_to_spread |> nrow()
  ) stop("sccomp says: You have duplicated .sample IDs in your input dataset. A .sample .cell_group combination must be unique")
  
  .data_to_spread |>
    spread(!!.cell_group, !!.count)
  
  
}

#' @importFrom rlang is_symbolic
#' @export
sccomp_estimate.Seurat <- function(.data,
                                   formula_composition = ~1,
                                   formula_variability = ~1,
                                   
                                   sample,
                                   cell_group,
                                   abundance = NULL,
                                   
                                   # Secondary arguments
                                   cores = detectCores(),
                                   bimodal_mean_variability_association = FALSE,
                                   percent_false_positive = 5,
                                   inference_method = "pathfinder",
                                   prior_mean = list(intercept = c(0, 1), coefficients = c(0, 1)),
                                   prior_overdispersion_mean_association = list(
                                     intercept = c(5, 2),
                                     slope = c(0, 0.6),
                                     standard_deviation = c(10, 20)
                                   ),
                                   .sample_cell_group_pairs_to_exclude = NULL,
                                   output_directory = "sccomp_draws_files",
                                   verbose = TRUE,
                                   enable_loo = FALSE,
                                   noise_model = "multi_beta_binomial",
                                   exclude_priors = FALSE,
                                   use_data = TRUE,
                                   mcmc_seed = sample_seed(),
                                   max_sampling_iterations = 20000,
                                   pass_fit = TRUE,
                                   sig_figs = 9,
                                   ...,
                                   
                                   # DEPRECATED
                                   .count = NULL,
                                   approximate_posterior_inference = NULL,
                                   variational_inference = NULL,
                                   .sample = NULL,
                                   .cell_group = NULL,
                                   .abundance = NULL) {
  
  .count <- enquo(.count)
  .sample <- enquo(.sample)
  .cell_group <- enquo(.cell_group)
  .abundance <- enquo(.abundance)
  
  if (!quo_is_null(.abundance))
    stop("sccomp says: .abundance argument can be used only for data frame input")
  
  if (!quo_is_null(.count))
    stop("sccomp says: .count argument can be used only for data frame input")
  
  
  # Prepare column names
  
  .sample_cell_group_pairs_to_exclude <- enquo(.sample_cell_group_pairs_to_exclude)
  
  .data[[]] |>
    sccomp_estimate(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      
      sample = sample,
      cell_group = cell_group,
      
      # Secondary arguments
      cores = cores,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      percent_false_positive = percent_false_positive,
      inference_method = inference_method,
      prior_mean = prior_mean,
      prior_overdispersion_mean_association = prior_overdispersion_mean_association,
      .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      output_directory = output_directory,
      verbose = verbose,
      enable_loo = enable_loo,
      noise_model = noise_model,
      exclude_priors = exclude_priors,
      use_data = use_data,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit,
      sig_figs = sig_figs,
      ...,
      .count = !!.count,
      approximate_posterior_inference = approximate_posterior_inference,
      variational_inference = variational_inference,
      .sample = !!.sample,
      .cell_group = !!.cell_group,
      .abundance = !!.abundance
    )
}

#' @importFrom rlang is_symbolic
#' @export
sccomp_estimate.SingleCellExperiment <- function(.data,
                                                 formula_composition = ~1,
                                                 formula_variability = ~1,
                                                 sample,
                                                 cell_group,
                                                 abundance = NULL,
                                                 
                                                 # Secondary arguments
                                                 cores = detectCores(),
                                                 bimodal_mean_variability_association = FALSE,
                                                 percent_false_positive = 5,
                                                 inference_method = "pathfinder",
                                                 prior_mean = list(intercept = c(0, 1), coefficients = c(0, 1)),
                                                 prior_overdispersion_mean_association = list(
                                                   intercept = c(5, 2),
                                                   slope = c(0, 0.6),
                                                   standard_deviation = c(10, 20)
                                                 ),
                                                 .sample_cell_group_pairs_to_exclude = NULL,
                                                 output_directory = "sccomp_draws_files",
                                                 verbose = TRUE,
                                                 enable_loo = FALSE,
                                                 noise_model = "multi_beta_binomial",
                                                 exclude_priors = FALSE,
                                                 use_data = TRUE,
                                                 mcmc_seed = sample_seed(),
                                                 max_sampling_iterations = 20000,
                                                 pass_fit = TRUE,
                                                 sig_figs = 9,
                                                 ...,
                                                 
                                                 # DEPRECATED
                                                 .count = NULL,
                                                 approximate_posterior_inference = NULL,
                                                 variational_inference = NULL,
                                                 .sample = NULL,
                                                 .cell_group = NULL,
                                                 .abundance = NULL) {
  
  
  .count <- enquo(.count)
  .sample <- enquo(.sample)
  .cell_group <- enquo(.cell_group)
  .abundance <- enquo(.abundance)
  
  .sample_cell_group_pairs_to_exclude <- enquo(.sample_cell_group_pairs_to_exclude)
  
  if (!quo_is_null(.abundance))
    stop("sccomp says: .abundance argument can be used only for data frame input")
  
  if (!quo_is_null(.count))
    stop("sccomp says: .count argument can be used only for data frame input")
  
  
  
  .data |>
    colData() |>
    
    sccomp_estimate(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      sample = sample,
      cell_group = cell_group,
      
      # Secondary arguments
      cores = cores,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      percent_false_positive = percent_false_positive,
      inference_method = inference_method,
      prior_mean = prior_mean,
      prior_overdispersion_mean_association = prior_overdispersion_mean_association,
      .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      output_directory = output_directory,
      verbose = verbose,
      enable_loo = enable_loo,
      noise_model = noise_model,
      exclude_priors = exclude_priors,
      use_data = use_data,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit,
      sig_figs = sig_figs,
      ...,
      .count = !!.count,
      approximate_posterior_inference = approximate_posterior_inference,
      variational_inference = variational_inference,
      .sample = !!.sample,
      .cell_group = !!.cell_group,
      .abundance = !!.abundance
    )
}

#' @importFrom rlang is_symbolic
#' @export
sccomp_estimate.DFrame <- function(.data,
                                   formula_composition = ~1,
                                   formula_variability = ~1,
                                   sample,
                                   cell_group,
                                   abundance = NULL,
                                   
                                   # Secondary arguments
                                   cores = detectCores(),
                                   bimodal_mean_variability_association = FALSE,
                                   percent_false_positive = 5,
                                   inference_method = "pathfinder",
                                   prior_mean = list(intercept = c(0, 1), coefficients = c(0, 1)),
                                   prior_overdispersion_mean_association = list(
                                     intercept = c(5, 2),
                                     slope = c(0, 0.6),
                                     standard_deviation = c(10, 20)
                                   ),
                                   .sample_cell_group_pairs_to_exclude = NULL,
                                   output_directory = "sccomp_draws_files",
                                   verbose = TRUE,
                                   enable_loo = FALSE,
                                   noise_model = "multi_beta_binomial",
                                   exclude_priors = FALSE,
                                   use_data = TRUE,
                                   mcmc_seed = sample_seed(),
                                   max_sampling_iterations = 20000,
                                   pass_fit = TRUE,
                                   sig_figs = 9,
                                   ...,
                                   
                                   # DEPRECATED
                                   .count = NULL,
                                   approximate_posterior_inference = NULL,
                                   variational_inference = NULL,
                                   .sample = NULL,
                                   .cell_group = NULL,
                                   .abundance = NULL) {
  
  .count <- enquo(.count)
  .sample <- enquo(.sample)
  .cell_group <- enquo(.cell_group)
  .abundance <- enquo(.abundance)
  
  .sample_cell_group_pairs_to_exclude <- enquo(.sample_cell_group_pairs_to_exclude)
  
  .data %>%
    as.data.frame() %>%
    sccomp_estimate(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      sample = sample,
      cell_group = cell_group,
      abundance = abundance,
      
      # Secondary arguments
      cores = cores,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      percent_false_positive = percent_false_positive,
      inference_method = inference_method,
      prior_mean = prior_mean,
      prior_overdispersion_mean_association = prior_overdispersion_mean_association,
      .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      output_directory = output_directory,
      verbose = verbose,
      enable_loo = enable_loo,
      noise_model = noise_model,
      exclude_priors = exclude_priors,
      use_data = use_data,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit, 
      ...,
      .count = !!.count,
      approximate_posterior_inference = approximate_posterior_inference,
      variational_inference = variational_inference,
      .sample = !!.sample,
      .cell_group = !!.cell_group,
      .abundance = !!.abundance
    )
}


#' @importFrom purrr when
#' @importFrom rlang is_symbolic
#' @export
sccomp_estimate.data.frame <- function(.data,
                                       formula_composition = ~1,
                                       formula_variability = ~1,
                                       
                                       sample,
                                       cell_group,
                                       abundance = NULL,
                                       
                                       # Secondary arguments
                                       cores = detectCores(),
                                       bimodal_mean_variability_association = FALSE,
                                       percent_false_positive = 5,
                                       inference_method = "pathfinder",
                                       prior_mean = list(intercept = c(0, 1), coefficients = c(0, 1)),
                                       prior_overdispersion_mean_association = list(
                                         intercept = c(5, 2),
                                         slope = c(0, 0.6),
                                         standard_deviation = c(10, 20)
                                       ),
                                       .sample_cell_group_pairs_to_exclude = NULL,
                                       output_directory = "sccomp_draws_files",
                                       verbose = TRUE,
                                       enable_loo = FALSE,
                                       noise_model = "multi_beta_binomial",
                                       exclude_priors = FALSE,
                                       use_data = TRUE,
                                       mcmc_seed = sample_seed(),
                                       max_sampling_iterations = 20000,
                                       pass_fit = TRUE,
                                       sig_figs = 9,
                                       ...,
                                       
                                       # DEPRECATED
                                       .count = NULL,
                                       approximate_posterior_inference = NULL,
                                       variational_inference = NULL,
                                       .sample = NULL,
                                       .cell_group = NULL,
                                       .abundance = NULL) {
  
  
  .count <- enquo(.count)
  .sample <- enquo(.sample)
  .cell_group <- enquo(.cell_group)
  .abundance <- enquo(.abundance)
  
  .sample_cell_group_pairs_to_exclude <- enquo(.sample_cell_group_pairs_to_exclude)
  
  
  # Handle deprecated inference method arguments
  if (lifecycle::is_present(approximate_posterior_inference) && !is.null(approximate_posterior_inference)) {
    lifecycle::deprecate_warn("1.7.7", "sccomp::sccomp_estimate(approximate_posterior_inference = )", 
                              details = "The argument approximate_posterior_inference is now deprecated. Please use inference_method. By default, inference_method value is inferred from approximate_posterior_inference.")
    inference_method <- ifelse(approximate_posterior_inference == "all", "variational", "hmc")
  }
  
  if (lifecycle::is_present(variational_inference) && !is.null(variational_inference)) {
    lifecycle::deprecate_warn("1.7.11", "sccomp::sccomp_estimate(variational_inference = )", 
                              details = "The argument variational_inference is now deprecated. Please use inference_method. By default, inference_method value is inferred from variational_inference")
    inference_method <- ifelse(variational_inference, "variational", "hmc")
  }
  
  # Handle deprecated column arguments
  if (lifecycle::is_present(.sample) && 
      !rlang::quo_is_null(.sample) && 
      rlang::is_symbolic(.sample) && 
      !is.null(rlang::eval_tidy(.sample))) {
    lifecycle::deprecate_warn("2.1.1", "sccomp::sccomp_estimate(.sample)", 
                              details = ".sample argument (which were tidy evaluations, i.e. .sample = my_sample) have been deprecated in favour of sample (without trailing dot, and which is now a character string, i.e. sample = \"my_sample\")")
    sample = quo_name(.sample)
  } else {
    # Capture the quosure for sample using a symbol
    if(
      sample |>
      is.character() |>
      not() 
    ) {
      stop("sccomp says: sample must be of character type")
    }
  }
  
  if (lifecycle::is_present(.cell_group) && 
      !rlang::quo_is_null(.cell_group) && 
      rlang::is_symbolic(.cell_group) ) {
    lifecycle::deprecate_warn("2.1.1", "sccomp::sccomp_estimate(.cell_group)", 
                              details = ".cell_group argument (which were tidy evaluations, i.e. .cell_group = my_cell_group) have been deprecated in favour of cell_group (without trailing dot, and which is now a character string, i.e. cell_group = \"my_cell_group\")")
    cell_group = quo_name(.cell_group)
  } else {
    if(
      cell_group |>
      is.character() |>
      not() 
    ) {
      stop("sccomp says: cell_group must be of character type")
    }
  }
  
  # Handle deprecated .count argument
  if (lifecycle::is_present(.count) && 
      !rlang::quo_is_null(.count) && 
      rlang::is_symbolic(.count) ) {
    rlang::warn("The argument '.count' is deprecated in favour of abundance (without trailing dot, and which is now a character string, i.e. abundance = \"my_abundance\")")
    abundance <- quo_name(.count)
  }
  
  if (lifecycle::is_present(.abundance) && 
      !rlang::quo_is_null(.abundance) && 
      rlang::is_symbolic(.abundance) ) {
    lifecycle::deprecate_warn("2.1.1", "sccomp::sccomp_estimate(.abundance)", 
                              details = ".abundance argument (which were tidy evaluations, i.e. .abundance = my_abundance) have been deprecated in favour of abundance (without trailing dot, and which is now a character string, i.e. abundance = \"my_abundance\")")
    abundance = quo_name(.abundance)
  } else {
    if(
      !is.null(abundance) &&
      abundance  |>
      is.character() |>
      not() 
    ) {
      stop("sccomp says: abundance must be of character type")
    }
  }
  
  
  
  if (abundance |> is.null())
    res <- sccomp_glm_data_frame_raw(
      .data,
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      
      sample = sample,
      cell_group = cell_group,
      
      # Secondary arguments
      cores = cores,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      percent_false_positive = percent_false_positive,
      inference_method = inference_method,
      prior_mean = prior_mean,
      prior_overdispersion_mean_association = prior_overdispersion_mean_association,
      .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      output_directory = output_directory,
      verbose = verbose,
      enable_loo = enable_loo,
      exclude_priors = exclude_priors,
      use_data = use_data,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit,
      sig_figs = sig_figs,
      ...
    )
  
  else 
    res = sccomp_glm_data_frame_counts(
      .data,
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      
      sample = sample,
      cell_group = cell_group,
      abundance = abundance,
      
      # Secondary arguments
      cores = cores,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      percent_false_positive = percent_false_positive,
      inference_method = inference_method,
      prior_mean = prior_mean,
      prior_overdispersion_mean_association = prior_overdispersion_mean_association,
      .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      output_directory = output_directory,
      verbose = verbose,
      enable_loo = enable_loo,
      exclude_priors = exclude_priors,
      use_data = use_data,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit,
      sig_figs = sig_figs,
      ...
    )
  
  message("sccomp says: to do hypothesis testing run `sccomp_test()`,
  the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
  0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
  Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).")
  
  res |>
    
    # Track input parameters
    add_attr(noise_model, "noise_model")
}

#' @importFrom tidyr complete
#' @importFrom tidyr nesting
#' @importFrom tidyr replace_na
sccomp_glm_data_frame_raw = function(.data,
                                     formula_composition = ~ 1 ,
                                     formula_variability = ~ 1,
                                     
                                     sample,
                                     cell_group,
                                     abundance = NULL,
                                     
                                     
                                     # Secondary arguments
                                     contrasts = NULL,
                                     prior_mean = list(intercept = c(0,1), coefficients = c(0,1)),                        
                                     prior_overdispersion_mean_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                                     percent_false_positive =  5,
                                     check_outliers = TRUE,
                                     variational_inference = NULL,
                                     inference_method = "variational",
                                     test_composition_above_logit_fold_change = 0.1, .sample_cell_group_pairs_to_exclude = NULL,
                                     output_directory = "sccomp_draws_files",
                                     verbose = FALSE,
                                     exclude_priors = FALSE,
                                     bimodal_mean_variability_association = FALSE,
                                     enable_loo = FALSE,
                                     use_data = TRUE,
                                     cores = 4,
                                     mcmc_seed = sample_seed(),
                                     max_sampling_iterations = 20000,
                                     pass_fit = TRUE,
                                     sig_figs = 9,
                                     ...) {
  
  # See https://community.rstudio.com/t/how-to-make-complete-nesting-work-with-quosures-and-tidyeval/16473
  # See https://github.com/tidyverse/tidyr/issues/506
  
  
  # Prepare column same enquo
  .sample = as.symbol(sample)
  .cell_group = as.symbol(cell_group)
  .sample_cell_group_pairs_to_exclude = enquo(.sample_cell_group_pairs_to_exclude)
  
  # Check if columns exist
  check_columns_exist(.data, c(
    !!.sample,
    !!.cell_group,
    all_of(parse_formula(formula_composition))
  ))
  
  # Check if any column is NA or null
  check_if_any_NA(.data, c(
    !!.sample,
    !!.cell_group,
    all_of(parse_formula(formula_composition))
  ))
  
  .grouping_for_random_effect = parse_formula_random_effect(formula_composition) |> pull(grouping) |> unique()
  
  # Make counts
  .data %>%
    
    class_list_to_counts(!!.sample, !!.cell_group) %>%
    
    # Add formula_composition information
    add_formula_columns(.data, !!.sample,formula_composition) |>
    
    # Attach possible exclusion of data points
    left_join(.data %>%
                as_tibble() |>
                select(
                  !!.sample, !!.cell_group, any_of(quo_name(.sample_cell_group_pairs_to_exclude))
                ) %>%
                distinct(),
              by = c(quo_name(.sample), quo_name(.cell_group))
    ) |>
    mutate(
      across(!!.sample_cell_group_pairs_to_exclude, ~replace_na(.x, 0))
    ) |>
    
    # Return
    sccomp_glm_data_frame_counts(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      
      sample = sample,
      cell_group = cell_group,
      abundance = "count",
      
      contrasts = contrasts,
      #.grouping_for_random_effect = !! .grouping_for_random_effect,
      prior_mean = prior_mean,
      prior_overdispersion_mean_association = prior_overdispersion_mean_association,
      percent_false_positive =  percent_false_positive,
      check_outliers = check_outliers,
      inference_method = inference_method,
      exclude_priors = exclude_priors,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      enable_loo = enable_loo,
      use_data = use_data,
      cores = cores,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change, .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      verbose = verbose,
      output_directory = output_directory,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit,
      sig_figs = sig_figs,
      ...
    )
}



sccomp_glm_data_frame_counts = function(.data,
                                        formula_composition = ~ 1 ,
                                        formula_variability = ~ 1,
                                        
                                        sample,
                                        cell_group,
                                        abundance = NULL,
                                        
                                        # Secondary arguments
                                        contrasts = NULL,
                                        #.grouping_for_random_effect = NULL,
                                        prior_mean = list(intercept = c(0,1), coefficients = c(0,1)),                        
                                        prior_overdispersion_mean_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                                        percent_false_positive = 5,
                                        check_outliers = TRUE,
                                        variational_inference = NULL,
                                        inference_method = "variational",
                                        test_composition_above_logit_fold_change = 0.1, .sample_cell_group_pairs_to_exclude = NULL,
                                        output_directory = "sccomp_draws_files",
                                        verbose = FALSE,
                                        exclude_priors = FALSE,
                                        bimodal_mean_variability_association = FALSE,
                                        enable_loo = FALSE,
                                        use_data = TRUE,
                                        cores = 4,
                                        mcmc_seed = sample_seed(),
                                        max_sampling_iterations = 20000,
                                        pass_fit = TRUE,
                                        sig_figs = 9,
                                        ...) {
  
  # Prepare column same enquo
  .sample = as.symbol(sample)
  .cell_group = as.symbol(cell_group)
  .count = as.symbol(abundance)
  
  .sample_cell_group_pairs_to_exclude = enquo(.sample_cell_group_pairs_to_exclude)
  #.grouping_for_random_effect = enquo(.grouping_for_random_effect)
  
  # Check Sample Consistency of Factors
  check_sample_consistency_of_factors(.data, formula_composition, sample, cell_group)
  
  #Check column class
  check_if_columns_right_class(.data, !!.sample, !!.cell_group)
  
  # Check that count is integer
  if(.data %>% pull(!!.count) %>% is("integer")) message(sprintf("sccomp says: %s column is an integer. The sum-constrained beta binomial model will be used", quo_name(.count)))
  else if(.data %>% pull(!!.count) %>% is("integer") |> not() & .data %>% pull(!!.count) |> dplyr::between(0, 1) |> all()) message(sprintf("sccomp says: %s column is a proportion. The sum-constrained beta model will be used. When possible using counts is preferred as the binomial noise component is often dominating for rare groups (e.g. rare cell types).", quo_name(.count)))
  else stop(sprintf("sccomp: `%s` column must be an integer or a proportion", quo_name(.count)))
  
  # Check if columns exist
  check_columns_exist(.data, c(
    !!.sample,
    !!.cell_group,
    !!.count,
    all_of(parse_formula(formula_composition))
    
  ))
  
  # Check that NAs are not in counts
  check_if_NAs_in_count(.data, !!.count)
  
  
  # Make rectangular data
  .data = .data |> make_rectangular_data(!!.sample, !!.cell_group, !!.count, formula_composition)
  
  # Check if test_composition_above_logit_fold_change is 0, as the Bayesian FDR does not allow it
  if(test_composition_above_logit_fold_change <= 0)
    stop("sccomp says: test_composition_above_logit_fold_change should be > 0 for the FDR to be calculated in the Bayesian context (doi: 10.1093/biostatistics/kxw041). Also, testing for > 0 differences avoids significant but meaningless (because of the small magnitude) estimates.")
  
  
  # Check if any column is NA or null
  check_if_any_NA(.data, c(
    quo_name(.sample),
    quo_name(.cell_group),
    quo_name(.count),
    parse_formula(formula_composition)
  ))
  
  # Return
  
  # Credible interval
  CI = 1 - (percent_false_positive/100)
  
  # Produce data list
  factor_names = parse_formula(formula_composition)
  
  # Random intercept
  random_effect_elements = parse_formula_random_effect(formula_composition)
  
  # Variational only if no random intercept
  if(inference_method=="variational" & random_effect_elements |> filter(factor != "(Intercept)") |> nrow() > 0)
    stop("sccomp says: for random effect modelling plese use `inference_method` = \"hmc\", for the full Bayes HMC inference.")
  
  # If no random intercept fake it
  if(nrow(random_effect_elements)>0){
    .grouping_for_random_effect = random_effect_elements |> pull(grouping) |> unique() |>   map(~ quo(!! sym(.x)))
    
  } else{
    .grouping_for_random_effect = list(quo(!! sym("random_effect")))
    .data = .data |> mutate(!!.grouping_for_random_effect[[1]] := "1")
  }
  
  # If .sample_cell_group_pairs_to_exclude 
  if(quo_is_symbolic(.sample_cell_group_pairs_to_exclude)){
    
    # Error if not logical
    if(.data |> pull(!!.sample_cell_group_pairs_to_exclude) |> is("logical") |> not())
      stop(glue("sccomp says: {quo_name(.sample_cell_group_pairs_to_exclude)} must be logical"))
    
    # Error if not consistent to sample/cell group
    if(.data |> count(!!.sample, !!.cell_group, name = "n") |> filter(n>1) |> nrow() |> gt(0))
      stop(glue("sccomp says: {quo_name(.sample_cell_group_pairs_to_exclude)} must be unique with .sample/.cell_group pairs. You might have a .sample/.cell_group pair with both TRUE and FALSE {quo_name(.sample_cell_group_pairs_to_exclude)}."))
    
  } else {
    
    # If no .sample_cell_group_pairs_to_exclude fake it
    .sample_cell_group_pairs_to_exclude = quo(!! sym(".sample_cell_group_pairs_to_exclude"))
    .data = .data |> mutate(!!.sample_cell_group_pairs_to_exclude := FALSE)
  }
  
  
  # Original - old
  # prec_sd ~ normal(0,2);
  # prec_coeff ~ normal(0,5);
  
  message("sccomp says: estimation")
  
  data_for_model =
    .data %>%
    data_to_spread(
      formula_composition, 
      !!.sample, !!.cell_group, !!.count, 
      .grouping_for_random_effect |> map(~ .x |> quo_name() ) |> unlist()
    ) |>
    
    # This emerged with
    # https://github.com/MangiolaLaboratory/sccomp/issues/175#issuecomment-2622749180
    check_if_sample_is_a_unique_identifier(!!.sample) |> 
    
    # Create input for Stan
    data_spread_to_model_input(
      formula_composition, !!.sample, !!.cell_group, !!.count,
      truncation_ajustment = 1.1,
      approximate_posterior_inference = inference_method %in% c("variational", "pathfinder"),
      formula_variability = formula_variability,
      contrasts = contrasts,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      use_data = use_data,
      random_effect_elements
    )
  
  # Print design matrix
  message(sprintf("sccomp says: the composition design matrix has columns: %s", data_for_model$X %>% colnames %>% paste(collapse=", ")))
  message(sprintf("sccomp says: the variability design matrix has columns: %s", data_for_model$Xa %>% colnames %>% paste(collapse=", ")))
  
  # Force outliers, Get the truncation index
  data_for_model$user_forced_truncation_not_idx = 
    .data |> 
    select(!!.sample, !!.cell_group, !!.sample_cell_group_pairs_to_exclude) |>
    left_join( data_for_model$X |> rownames() |> enframe(name="N", value=quo_name(.sample)), by = join_by(!!.sample) ) |>  
    left_join( data_for_model$y |> colnames() |> enframe(name="M", value=quo_name(.cell_group)), by = join_by(!!.cell_group) ) |> 
    select(!!.sample_cell_group_pairs_to_exclude, N, M) |>
    arrange(N, M) |>
    pull(!!.sample_cell_group_pairs_to_exclude) |>
    not() |>
    which()
  
  data_for_model$truncation_not_idx = data_for_model$user_forced_truncation_not_idx
  data_for_model$TNS = length(data_for_model$truncation_not_idx)
  
  # Prior
  data_for_model$prior_prec_intercept = prior_overdispersion_mean_association$intercept
  data_for_model$prior_prec_slope  = prior_overdispersion_mean_association$slope
  data_for_model$prior_prec_sd = prior_overdispersion_mean_association$standard_deviation
  data_for_model$prior_mean_intercept = prior_mean$intercept
  data_for_model$prior_mean_coefficients = prior_mean$coefficients
  data_for_model$exclude_priors = exclude_priors
  data_for_model$enable_loo = enable_loo
  
  # # Check that design matrix is not too big
  # if(ncol(data_for_model$X)>20)
  #   message("sccomp says: the design matrix has more than 20 columns. Possibly some numerical factors are erroneously of type character/factor.")
  
  fit =
    data_for_model %>%
    
    # Run the first discovery phase with permissive false discovery rate
    fit_model(
      "glm_multi_beta_binomial",
      cores= cores,
      quantile = CI,
      inference_method = inference_method,
      verbose = verbose,
      output_directory = output_directory,
      seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pars = c(
        "beta", "alpha", "prec_coeff","prec_sd",   "alpha_normalised", 
        "random_effect", "random_effect_2", 
        "random_effect_sigma", "random_effect_sigma_2", 
        "log_lik"
      ),
      sig_figs = sig_figs,
      ...
    )
  
  # # Make the fit standalone
  # temp_rds_file <- tempfile(fileext = ".rds")
  # fit$save_object(file = temp_rds_file) 
  # fit = readRDS(temp_rds_file)
  # file.remove(temp_rds_file)
  
  estimate_tibble = 
    # Create a dummy tibble
    tibble() |>
    # Attach association mean concentration
    add_attr(fit, "fit") %>%
    add_attr(data_for_model, "model_input") |> 

    # Save data only taking the columns essential, including the ones use for the formulas 
    add_attr(.data |> select(
      !!.sample, !!.cell_group, !!.count, 
      all_of(parse_formula(formula_composition)),
      all_of(random_effect_elements |> pull(grouping) |> unique())
      ), "count_data") |>
    
    add_attr(.sample |> drop_environment(), ".sample") |>
    add_attr(.cell_group |> drop_environment(), ".cell_group") |>
    add_attr(.count |> drop_environment(), ".count") |>
    add_attr(check_outliers, "check_outliers") |>
    add_attr(formula_composition |> drop_environment(), "formula_composition") |>
    add_attr(formula_variability |> drop_environment(), "formula_variability") |>
    add_attr(parse_formula(formula_composition), "factors" ) |> 
    add_attr(inference_method, "inference_method" ) |> 
    
    # Add class to the tbl
    add_class("sccomp_tbl") |> 
    
    # Print estimates
    sccomp_test() |>
    
    # drop hypothesis testing as the estimation exists without probabilities.
    # For hypothesis testing use sccomp_test
    select(-contains("_FDR"), -contains("_pH0")) 
  
  
  if(inference_method %in% c("variational") && max(na.omit(estimate_tibble$c_R_k_hat)) > 4)
    warning("sccomp says: using variational inference, c_R_k_hat resulted too high for some parameters, indicating lack of convergence of the model. We reccomend using inference_method = \"hmc\" to use the state-of-the-art (although slower) HMC sampler.")
  
  estimate_tibble
}

#' Check if data frame is rectangular
#'
#' @description
#' Checks if the input data frame has the same number of samples for each cell group
#' and vice versa. If not, it makes the data frame rectangular by adding zeros
#' for missing sample/cell group pairs.
#'
#' @param .data A tibble containing the data
#' @param .sample Column containing sample identifiers
#' @param .cell_group Column containing cell group identifiers
#' @param .count Column containing count data
#' @param formula_composition Formula for composition
#'
#' @return A rectangular data frame with zeros added for missing combinations
#' @keywords internal
#' @noRd
make_rectangular_data = function(.data, .sample, .cell_group, .count, formula_composition) {
  
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)
  
  if(
    .data |> count(!!.sample) |> distinct(n) |> nrow() > 1 || 
    .data |> count(!!.cell_group) |> distinct(n) |> nrow() > 1 
  ){
    warning(sprintf("sccomp says: the input data frame does not have the same number of `%s`, for all `%s`. We have made it so, adding 0s for the missing sample/feature pairs.", quo_name(.cell_group), quo_name(.sample)))
    
    .data |> 
      # I need renaming trick because complete list(...) cannot accept quosures
      select(!!.sample, !!.cell_group, count := !!.count) |> 
      distinct() |> 
      complete( !!.sample, !!.cell_group, fill = list(count = 0) ) %>%
      rename(!!.count := count) |> 
      mutate(!!.count := as.integer(!!.count)) |> 
      
      # Add formula_composition information
      add_formula_columns(.data, !!.sample, formula_composition)
  } else {
    .data
  }
}

#' @importFrom rlang ensym
#' @noRd
class_list_to_counts = function(.data, .sample, .cell_group){
  
  .sample_for_tidyr = ensym(.sample)
  .cell_group_for_tidyr = ensym(.cell_group)
  
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  
  
  .data %>%
    count(!!.sample,
          !!.cell_group,
          name = "count") %>%
    
    complete(
      !!.sample_for_tidyr,!!.cell_group_for_tidyr,
      fill = list(count = 0)
    ) %>%
    mutate(count = as.integer(count))
}

add_formula_columns = function(.data, .original_data, .sample,  formula_composition){
  
  .sample = enquo(.sample)
  
  formula_elements = parse_formula(formula_composition)
  
  # If no formula return the input
  if(length(formula_elements) == 0) return(.data)
  
  # Get random intercept
  .grouping_for_random_effect = parse_formula_random_effect(formula_composition) |> pull(grouping) |> unique()
  
  data_frame_formula =
    .original_data %>%
    as_tibble() |>
    select( !!.sample, all_of(formula_elements), any_of(.grouping_for_random_effect) ) %>%
    distinct()
  
  .data |>
    left_join(data_frame_formula, by = quo_name(.sample) )
  
}
