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
#'
#' @param .data A tibble including cell_group name column, sample name column,
#'              abundance column (counts or proportions), and factor columns.
#' @param formula_composition A formula describing the model for differential abundance.
#' @param formula_variability A formula describing the model for differential variability.
#' @param sample_column A column name as a character string for the sample identifier. Replaces the deprecated `.sample`.
#' @param cell_group_column A column name as a character string for the cell-group identifier. Replaces the deprecated `.cell_group`.
#' @param abundance_column A column name as a character string for the cell-group abundance, which can be counts (> 0) or proportions (between 0 and 1, summing to 1 across `cell_group_column`). Replaces the deprecated `.abundance` and `.count`.
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
#' @param .count **DEPRECATED**. Use `abundance_column` instead.
#' @param approximate_posterior_inference **DEPRECATED**. Use `inference_method` instead.
#' @param variational_inference **DEPRECATED**. Use `inference_method` instead.
#' @param .sample **DEPRECATED**. Use `sample_column` instead.
#' @param .cell_group **DEPRECATED**. Use `cell_group_column` instead.
#' @param .abundance **DEPRECATED**. Use `abundance_column` instead.
#' @param ... Additional arguments passed to the `cmdstanr::sample` function.
#'
#' @return A tibble (`tbl`) with the following columns:
#' \itemize{
#'   \item cell_group - The cell groups being tested.
#'   \item parameter - The parameter being estimated from the design matrix described by the input `formula_composition` and `formula_variability`.
#'   \item factor - The covariate factor in the formula, if applicable (e.g., not present for Intercept or contrasts).
#'   \item c_lower - Lower (2.5%) quantile of the posterior distribution for a composition (c) parameter.
#'   \item c_effect - Mean of the posterior distribution for a composition (c) parameter.
#'   \item c_upper - Upper (97.5%) quantile of the posterior distribution for a composition (c) parameter.
#'   \item c_pH0 - Probability of the null hypothesis (no difference) for a composition (c). This is not a p-value.
#'   \item c_FDR - False-discovery rate of the null hypothesis for a composition (c).
#'   \item c_n_eff - Effective sample size for a composition (c) parameter.
#'   \item c_R_k_hat - R statistic for a composition (c) parameter, should be within 0.05 of 1.0.
#'   \item v_lower - Lower (2.5%) quantile of the posterior distribution for a variability (v) parameter.
#'   \item v_effect - Mean of the posterior distribution for a variability (v) parameter.
#'   \item v_upper - Upper (97.5%) quantile of the posterior distribution for a variability (v) parameter.
#'   \item v_pH0 - Probability of the null hypothesis for a variability (v).
#'   \item v_FDR - False-discovery rate of the null hypothesis for a variability (v).
#'   \item v_n_eff - Effective sample size for a variability (v) parameter.
#'   \item v_R_k_hat - R statistic for a variability (v) parameter.
#'   \item count_data - Nested input count data.
#' }
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
#'       sample,
#'       cell_group,
#'       count,
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
#'       sample,
#'       cell_group,
#'       proportion,
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
                            
                            sample_column,
                            cell_group_column,
                            abundance_column = NULL,
                            
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
                            mcmc_seed = sample(1e5, 1),
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
  
  
  
  # rlang::inform(
  #   message = "sccomp says: From version 1.7.12 the logit fold change threshold for significance has been changed from 0.2 to 0.1.",
  #   .frequency = "once",
  #   .frequency_id = "new_logit_fold_change_threshold"
  # )
  
  
  # Run the function
  check_and_install_cmdstanr()
  
  UseMethod("sccomp_estimate", .data)
}

#' @importFrom rlang is_symbolic
#' @export
sccomp_estimate.Seurat <- function(.data,
                                   formula_composition = ~1,
                                   formula_variability = ~1,
                                   
                                   sample_column,
                                   cell_group_column,
                                   abundance_column = NULL,
                                   
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
                                   mcmc_seed = sample(1e5, 1),
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
  
  .sample <- enquo(.sample)
  .cell_group <- enquo(.cell_group)
  
  if (!is.null(.abundance))
    stop("sccomp says: .abundance argument can be used only for data frame input")
  
  if (!is.null(.count))
    stop("sccomp says: .count argument can be used only for data frame input")
  
  # DEPRECATION OF approximate_posterior_inference
  if (lifecycle::is_present(approximate_posterior_inference) & !is.null(approximate_posterior_inference)) {
    lifecycle::deprecate_warn("1.7.7", "sccomp::sccomp_estimate(approximate_posterior_inference = )", details = "The argument approximate_posterior_inference is now deprecated. Please use inference_method. By default, inference_method value is inferred from approximate_posterior_inference.")
    
    inference_method <- ifelse(approximate_posterior_inference == "all", "variational", "hmc")
  }
  
  # DEPRECATION OF variational_inference
  if (lifecycle::is_present(variational_inference) & !is.null(variational_inference)) {
    lifecycle::deprecate_warn("1.7.11", "sccomp::sccomp_estimate(variational_inference = )", details = "The argument variational_inference is now deprecated. Please use inference_method. By default, inference_method value is inferred from variational_inference")
    
    inference_method <- ifelse(variational_inference, "variational", "hmc")
  }
  
  # DEPRECATION OF .sample
  if (lifecycle::is_present(.sample) & !quo_is_null(.sample)) {
    lifecycle::deprecate_warn("2.1.1", "sccomp::sccomp_estimate(.sample)", details = ".sample argument (which were tidy evaluations, i.e. .sample = my_sample) have been deprecated in favour of sample_column (without trailing dot, and which is now a character string, i.e. sample_column = \"my_sample\")")
    
    sample_column = quo_name(.sample)
    
  }
  
  # DEPRECATION OF .cell_group
  if (lifecycle::is_present(.cell_group) & !quo_is_null(.cell_group)) {
    lifecycle::deprecate_warn("2.1.1", "sccomp::sccomp_estimate(.cell_group)", details = ".cell_group argument (which were tidy evaluations, i.e. .cell_group = my_cell_group) have been deprecated in favour of cell_group_column (without trailing dot, and which is now a character string, i.e. cell_group_column = \"my_cell_group\")")
    
    cell_group_column = quo_name(.cell_group)
    
  }
  
  
  # Prepare column names
  
  .sample_cell_group_pairs_to_exclude <- enquo(.sample_cell_group_pairs_to_exclude)
  
  .data[[]] |>
    sccomp_estimate(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      
      sample_column = sample_column,
      cell_group_column = cell_group_column,
      
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
      ...
    )
}

#' @importFrom rlang is_symbolic
#' @export
sccomp_estimate.SingleCellExperiment <- function(.data,
                                                 formula_composition = ~1,
                                                 formula_variability = ~1,
                                                 sample_column,
                                                 cell_group_column,
                                                 abundance_column = NULL,
                                                 
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
                                                 mcmc_seed = sample(1e5, 1),
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
  
  # Prepare column names
  .sample <- enquo(.sample)
  .cell_group <- enquo(.cell_group)
  .sample_cell_group_pairs_to_exclude <- enquo(.sample_cell_group_pairs_to_exclude)
  
  
  if (!is.null(.abundance))
    stop("sccomp says: .abundance argument can be used only for data frame input")
  
  if (!is.null(.count))
    stop("sccomp says: .count argument can be used only for data frame input")
  
  # DEPRECATION OF approximate_posterior_inference
  if (lifecycle::is_present(approximate_posterior_inference) & !is.null(approximate_posterior_inference)) {
    lifecycle::deprecate_warn("1.7.7", "sccomp::sccomp_estimate(approximate_posterior_inference = )", details = "The argument approximate_posterior_inference is now deprecated. Please use inference_method. By default, inference_method value is inferred from approximate_posterior_inference.")
    
    inference_method <- ifelse(approximate_posterior_inference == "all", "variational", "hmc")
  }
  
  # DEPRECATION OF variational_inference
  if (lifecycle::is_present(variational_inference) & !is.null(variational_inference)) {
    lifecycle::deprecate_warn("1.7.11", "sccomp::sccomp_estimate(variational_inference = )", details = "The argument variational_inference is now deprecated. Please use inference_method. By default, inference_method value is inferred from variational_inference")
    
    inference_method <- ifelse(variational_inference, "variational", "hmc")
  }
  
  # DEPRECATION OF .sample
  if (lifecycle::is_present(.sample) & !quo_is_null(.sample)) {
    lifecycle::deprecate_warn("2.1.1", "sccomp::sccomp_estimate(.sample)", details = ".sample argument (which were tidy evaluations, i.e. .sample = my_sample) have been deprecated in favour of sample_column (without trailing dot, and which is now a character string, i.e. sample_column = \"my_sample\")")
    
    sample_column = quo_name(.sample)
    
  }
  
  # DEPRECATION OF .cell_group
  if (lifecycle::is_present(.cell_group) & !quo_is_null(.cell_group)) {
    lifecycle::deprecate_warn("2.1.1", "sccomp::sccomp_estimate(.cell_group)", details = ".cell_group argument (which were tidy evaluations, i.e. .cell_group = my_cell_group) have been deprecated in favour of cell_group_column (without trailing dot, and which is now a character string, i.e. cell_group_column = \"my_cell_group\")")
    
    cell_group_column = quo_name(.cell_group)
    
  }
  
  .data |>
    colData() |>
    
    sccomp_estimate(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      sample_column = sample_column,
      cell_group_column = cell_group_column,
      
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
      ...
    )
}

#' @importFrom rlang is_symbolic
#' @export
sccomp_estimate.DFrame <- function(.data,
                                   formula_composition = ~1,
                                   formula_variability = ~1,
                                   sample_column,
                                   cell_group_column,
                                   abundance_column = NULL,
                                   
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
                                   mcmc_seed = sample(1e5, 1),
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
  
  if (!is.null(.abundance))
    stop("sccomp says: .abundance argument can be used only for data frame input")
  
  if (!is.null(.count))
    stop("sccomp says: .count argument can be used only for data frame input")
  
  # Prepare column names
  .sample <- enquo(.sample)
  .cell_group <- enquo(.cell_group)
  .abundance <- enquo(.abundance)
  .sample_cell_group_pairs_to_exclude <- enquo(.sample_cell_group_pairs_to_exclude)
  .count = enquo(.count)
  
  
  # Deprecation of .count
  if (rlang::quo_is_symbolic(.count)) {
    rlang::warn("The argument '.count' is deprecated. Please use '.abundance' instead. This because now `sccomp` cam model both counts and proportions.")
    .abundance <- .count
  }
  
  
  # DEPRECATION OF approximate_posterior_inference
  if (lifecycle::is_present(approximate_posterior_inference) & !is.null(approximate_posterior_inference)) {
    lifecycle::deprecate_warn("1.7.7", "sccomp::sccomp_estimate(approximate_posterior_inference = )", details = "The argument approximate_posterior_inference is now deprecated. Please use inference_method. By default, inference_method value is inferred from approximate_posterior_inference.")
    
    inference_method <- ifelse(approximate_posterior_inference == "all", "variational", "hmc")
  }
  
  # DEPRECATION OF variational_inference
  if (lifecycle::is_present(variational_inference) & !is.null(variational_inference)) {
    lifecycle::deprecate_warn("1.7.11", "sccomp::sccomp_estimate(variational_inference = )", details = "The argument variational_inference is now deprecated. Please use inference_method. By default, inference_method value is inferred from variational_inference")
    
    inference_method <- ifelse(variational_inference, "variational", "hmc")
  }
  
  # DEPRECATION OF .sample
  if (lifecycle::is_present(.sample) & !quo_is_null(.sample)) {
    lifecycle::deprecate_warn("2.1.1", "sccomp::sccomp_estimate(.sample)", details = ".sample argument (which were tidy evaluations, i.e. .sample = my_sample) have been deprecated in favour of sample_column (without trailing dot, and which is now a character string, i.e. sample_column = \"my_sample\")")
    
    sample_column = quo_name(.sample)
    
  }
  
  # DEPRECATION OF .cell_group
  if (lifecycle::is_present(.cell_group) & !quo_is_null(.cell_group)) {
    lifecycle::deprecate_warn("2.1.1", "sccomp::sccomp_estimate(.cell_group)", details = ".cell_group argument (which were tidy evaluations, i.e. .cell_group = my_cell_group) have been deprecated in favour of cell_group_column (without trailing dot, and which is now a character string, i.e. cell_group_column = \"my_cell_group\")")
    
    cell_group_column = quo_name(.cell_group)
    
  }
  
  # DEPRECATION OF .abundance
  if (lifecycle::is_present(.abundance) & !quo_is_null(.abundance)) {
    lifecycle::deprecate_warn("2.1.1", "sccomp::sccomp_estimate(.abundance)", details = ".abundance argument (which were tidy evaluations, i.e. .abundance = my_abundance) have been deprecated in favour of abundance_column (without trailing dot, and which is now a character string, i.e. abundance_column = \"my_abundance\")")
    
    abundance_column = quo_name(.abundance)
    
  }
  
  .data %>%
    as.data.frame() %>%
    sccomp_estimate(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      sample_column = sample_column,
      cell_group_column = cell_group_column,
      abundance_column = abundance_column,
      
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
      pass_fit = pass_fit, ...
    )
}


#' @importFrom purrr when
#' @importFrom rlang is_symbolic
#' @export
sccomp_estimate.data.frame <- function(.data,
                                       formula_composition = ~1,
                                       formula_variability = ~1,
                                       
                                       sample_column,
                                       cell_group_column,
                                       abundance_column = NULL,
                                       
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
                                       mcmc_seed = sample(1e5, 1),
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
  
  
  # Prepare column names
  .sample <- enquo(.sample)
  .cell_group <- enquo(.cell_group)
  .abundance <- enquo(.abundance)
  .count <- enquo(.count)
  .sample_cell_group_pairs_to_exclude <- enquo(.sample_cell_group_pairs_to_exclude)
  
  # Deprecation of .count
  if (rlang::quo_is_symbolic(.count)) {
    rlang::warn("The argument '.count' is deprecated. Please use '.abundance' instead. This because now `sccomp` cam model both counts and proportions.")
    .abundance <- .count
  }
  
  # DEPRECATION OF approximate_posterior_inference
  if (lifecycle::is_present(approximate_posterior_inference) & !is.null(approximate_posterior_inference)) {
    lifecycle::deprecate_warn("1.7.7", "sccomp::sccomp_estimate(approximate_posterior_inference = )", details = "The argument approximate_posterior_inference is now deprecated. Please use inference_method. By default, inference_method value is inferred from approximate_posterior_inference.")
    
    inference_method <- ifelse(approximate_posterior_inference == "all", "variational", "hmc")
  }
  
  # DEPRECATION OF variational_inference
  if (lifecycle::is_present(variational_inference) & !is.null(variational_inference)) {
    lifecycle::deprecate_warn("1.7.11", "sccomp::sccomp_estimate(variational_inference = )", details = "The argument variational_inference is now deprecated. Please use inference_method. By default, inference_method value is inferred from variational_inference")
    
    inference_method <- ifelse(variational_inference, "variational", "hmc")
  }
  
  # DEPRECATION OF .sample
  if (lifecycle::is_present(.sample) & !quo_is_null(.sample)) {
    lifecycle::deprecate_warn("2.1.1", "sccomp::sccomp_estimate(.sample)", details = ".sample argument (which were tidy evaluations, i.e. .sample = my_sample) have been deprecated in favour of sample_column (without trailing dot, and which is now a character string, i.e. sample_column = \"my_sample\")")
    
    sample_column = quo_name(.sample)
    
  }
  
  # DEPRECATION OF .cell_group
  if (lifecycle::is_present(.cell_group) & !quo_is_null(.cell_group)) {
    lifecycle::deprecate_warn("2.1.1", "sccomp::sccomp_estimate(.cell_group)", details = ".cell_group argument (which were tidy evaluations, i.e. .cell_group = my_cell_group) have been deprecated in favour of cell_group_column (without trailing dot, and which is now a character string, i.e. cell_group_column = \"my_cell_group\")")
    
    cell_group_column = quo_name(.cell_group)
    
  }
  
  # DEPRECATION OF .abundance
  if (lifecycle::is_present(.abundance) & !quo_is_null(.abundance)) {
    lifecycle::deprecate_warn("2.1.1", "sccomp::sccomp_estimate(.abundance)", details = ".abundance argument (which were tidy evaluations, i.e. .abundance = my_abundance) have been deprecated in favour of abundance_column (without trailing dot, and which is now a character string, i.e. abundance_column = \"my_abundance\")")
    
    abundance_column = quo_name(.abundance)
    
  }
  
  
  if (abundance_column |> is.null())
    res <- sccomp_glm_data_frame_raw(
      .data,
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      
      sample_column = sample_column,
      cell_group_column = cell_group_column,
      
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
      
      sample_column = sample_column,
      cell_group_column = cell_group_column,
      abundance_column = abundance_column,
      
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