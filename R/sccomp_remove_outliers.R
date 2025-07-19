#' sccomp_remove_outliers main
#'
#' @description 
#' The `sccomp_remove_outliers` function takes as input a table of cell counts with columns for cell-group identifier, sample identifier, integer count, and factors (continuous or discrete). The user can define a linear model using an input R formula, where the first factor is the factor of interest. Alternatively, `sccomp` accepts single-cell data containers (e.g., Seurat, SingleCellExperiment, cell metadata, or group-size) and derives the count data from cell metadata.
#'
#' @import dplyr
#' @importFrom magrittr %$% divide_by multiply_by equals
#' @importFrom rlang quo_is_null quo_is_symbolic inform
#' @importFrom SingleCellExperiment colData
#' @importFrom parallel detectCores
#' @importFrom tidyr unnest nest
#' @importFrom rlang :=
#' @importFrom dplyr select
#' @importFrom tibble column_to_rownames
#'
#' @param .estimate A tibble including a cell_group name column, sample name column, read counts column (optional depending on the input class), and factor columns.
#' @param percent_false_positive A real number between 0 and 100 (not inclusive), used to identify outliers with a specific false positive rate.
#' @param cores Integer, the number of cores to be used for parallel calculations.
#' @param inference_method Character string specifying the inference method to use ('pathfinder', 'hmc', or 'variational').
#' @param output_directory A character string specifying the output directory for Stan draws.
#' @param verbose Logical, whether to print progression details.
#' @param mcmc_seed Integer, used for Markov-chain Monte Carlo reproducibility. By default, a random number is sampled from 1 to 999999.
#' @param max_sampling_iterations Integer, limits the maximum number of iterations in case a large dataset is used, to limit computation time.
#' @param enable_loo Logical, whether to enable model comparison using the R package LOO. This is useful for comparing fits between models, similar to ANOVA.
#' @param sig_figs Number of significant figures to use for Stan model output. Default is 9.
#' @param cache_stan_model A character string specifying the cache directory for compiled Stan models. 
#'                        The sccomp version will be automatically appended to ensure version isolation.
#'                        Default is `sccomp_stan_models_cache_dir` which points to `~/.sccomp_models`.
#' @param approximate_posterior_inference DEPRECATED, use the `variational_inference` argument.
#' @param variational_inference DEPRECATED Logical, whether to use variational Bayes for posterior inference. It is faster and convenient. Setting this argument to `FALSE` runs full Bayesian (Hamiltonian Monte Carlo) inference, which is slower but the gold standard.
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
#'     estimate = sccomp_estimate(
#'       counts_obj,
#'       ~ type,
#'       ~1,
#'       "sample",
#'       "cell_group",
#'       "count",
#'       cores = 1
#'     ) |>
#'     sccomp_remove_outliers(cores = 1)
#'   }
#' }
#'
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
#'
#' @export
sccomp_remove_outliers <- function(.estimate,
                                   percent_false_positive = 5,
                                   cores = detectCores(),
                                   inference_method = .estimate |> attr("inference_method"),
                                   output_directory = "sccomp_draws_files",
                                   verbose = TRUE,
                                   mcmc_seed = sample_seed(),
                                   max_sampling_iterations = 20000,
                                   enable_loo = FALSE,
                                   sig_figs = 9,
                                   cache_stan_model = sccomp_stan_models_cache_dir,
                                   
                                   # DEPRECATED
                                   approximate_posterior_inference = NULL,
                                   variational_inference = NULL,
                                   ...
) {
  if(inference_method == "variational") 
    rlang::inform(
      message = "sccomp says: From version 1.7.7 the model by default is fit with the variational inference method (inference_method = \"variational\"; much faster). For a full Bayesian inference (HMC method; the gold standard) use inference_method = \"hmc\".", 
      .frequency = "once", 
      .frequency_id = "variational_message"
    )
  
  # Run the function
  check_and_install_cmdstanr()
  
  UseMethod("sccomp_remove_outliers", .estimate)
}

#' @importFrom readr write_file
#' @export
sccomp_remove_outliers.sccomp_tbl = function(.estimate,
                                             percent_false_positive = 5,
                                             cores = detectCores(),
                                             inference_method = .estimate |> attr("inference_method"),
                                             output_directory = "sccomp_draws_files",
                                             verbose = TRUE,
                                             mcmc_seed = sample_seed(),
                                             max_sampling_iterations = 20000,
                                             enable_loo = FALSE,
                                             sig_figs = 9,
                                             cache_stan_model = sccomp_stan_models_cache_dir,
                                             
                                             # DEPRECATED
                                             approximate_posterior_inference = NULL,
                                             variational_inference = NULL,
                                             ...
) {
  
  
  # DEPRECATION OF approximate_posterior_inference
  if (is_present(approximate_posterior_inference) & !is.null(approximate_posterior_inference)) {
    deprecate_warn("1.7.7", "sccomp::sccomp_estimate(approximate_posterior_inference = )", details = "The argument approximate_posterior_inference is now deprecated please use variational_inference. By default variational_inference value is inferred from approximate_posterior_inference.")
    inference_method = ifelse(approximate_posterior_inference == "all", "variational","hmc")
  }
  
  # DEPRECATION OF variational_inference
  if (is_present(variational_inference) & !is.null(variational_inference)) {
    deprecate_warn("1.7.11", "sccomp::sccomp_estimate(variational_inference = )", details = "The argument variational_inference is now deprecated please use inference_method. By default inference_method value is inferred from variational_inference")
    
    inference_method = ifelse(variational_inference, "variational","hmc")
  }
  
  # Prepare column same enquo
  .sample = .estimate |>  attr(".sample")
  .cell_group = .estimate |>  attr(".cell_group")
  .count = .estimate |>  attr(".count")
  .sample_cell_group_pairs_to_exclude = .estimate |> attr(".sample_cell_group_pairs_to_exclude")
  
  # Formulae
  formula_composition = .estimate |> attr("formula_composition")
  formula_variability = .estimate |> attr("formula_variability")
  
  noise_model = .estimate |> attr("noise_model")
  
  # Get model input
  data_for_model = .estimate |> attr("model_input")
  
  # Credible interval
  CI = 1 - (percent_false_positive/100)
  
  # Count data
  .data = 
    .estimate |> 
    attr("count_data") |>
    distinct() |> 
    
    # Drop previous outlier estimation for the new one
    select(-any_of(c(
      ".lower" ,  ".median" ,  ".upper" , "Rhat" ,  "truncation_down" , 
      "truncation_up"   , "outlier"    ,   "contains_outliers"
    )))
  
  # Random intercept
  random_effect_elements = .estimate |> attr("formula_composition") |> parse_formula_random_effect()
  
  # Load model
  mod_rng = load_model("glm_multi_beta_binomial_generate_data", threads = cores, cache_dir = cache_stan_model)
  
  rng = mod_rng |> sample_safe(
    generate_quantities_fx,
    attr(.estimate , "fit")$draws(format = "matrix"),
    
    # This is for the new data generation with selected factors to do adjustment
    data = 
      .estimate |>
      attr("model_input") |> 
      c(list(
        
        # Add subset of coefficients
        X_original = data_for_model$X,
        N_original = data_for_model$N,
        length_X_which = ncol(data_for_model$X),
        length_XA_which = ncol(data_for_model$XA),
        X_which = seq_len(ncol(data_for_model$X)) |> as.array(),
        XA_which = seq_len(ncol(data_for_model$Xa)) |> as.array(),
        
        # Random intercept common variable between grouping 1 and 2
        ncol_X_random_eff_new = ncol(data_for_model$X_random_effect) |> c(ncol(data_for_model$X_random_effect_2) ), # I could put this in the intial data
        length_X_random_effect_which = ncol(data_for_model$X_random_effect) |> c(ncol(data_for_model$X_random_effect_2)),
        
        # Grouping 1
        X_random_effect_which = seq_len(ncol(data_for_model$X_random_effect)) |> as.array(),
        
        # Grouping 2 - Random intercept DUPLICATED
        X_random_effect_which_2 = seq_len(ncol(data_for_model$X_random_effect_2)) |> as.array(),
        
        # Initialize unseen random effect variables
        ncol_X_random_eff_unseen = c(0, 0),
        X_random_effect_unseen = matrix(0, nrow = nrow(data_for_model$X), ncol = 0),
        X_random_effect_2_unseen = matrix(0, nrow = nrow(data_for_model$X), ncol = 0),
        
        create_intercept = FALSE
      )),
    
    parallel_chains = ifelse(
      inference_method %in% c("variational", "pathfinder") | 
        attr(.estimate , "fit") |> is("CmdStanPathfinder"),
      1, 
      attr(.estimate , "fit")$num_chains()
    ), 
    threads_per_chain = cores,
    sig_figs = sig_figs
    
    
  )
  
  # Free memory
  rm(.estimate)
  
  # Detect outliers
  truncation_df =
    .data |>
    left_join(
      summary_to_tibble(rng, "counts", "N", "M", probs = c(0.05, 0.95)) |>
        nest(data = -N) |>
        mutate(!!.sample := rownames(data_for_model$X)) |>
        unnest(data) |>
        nest(data = -M) |>
        mutate(!!.cell_group := colnames(data_for_model$y)) |>
        unnest(data) ,
      
      by = c(quo_name(.sample), quo_name(.cell_group))
    ) |>
    
    # Add truncation
    mutate(   truncation_down = `5%`,   truncation_up =  `95%`) |>
    
    # Add outlier stats
    mutate( outlier = !(!!.count >= `5%` & !!.count <= `95%`) ) |>
    nest(data = -M) |>
    mutate(contains_outliers = map_lgl(data, ~ .x |> filter(outlier) |> nrow() > 0)) |>
    unnest(data) |>
    
    mutate(
      truncation_down = case_when( outlier ~ -1, TRUE ~ truncation_down),
      truncation_up = case_when(outlier ~ -1, TRUE ~ truncation_up),
    )
  
  # Add censoring
  data_for_model$is_truncated = 1
  data_for_model$truncation_up = truncation_df |> 
    select(N, M, truncation_up) |> 
    pivot_wider(names_from = M, values_from = truncation_up) |> 
    column_to_rownames("N") |> 
    as.matrix() |> 
    apply(2, as.integer)
  data_for_model$truncation_down = truncation_df |> 
    select(N, M, truncation_down) |> 
    pivot_wider(names_from = M, values_from = truncation_down) |> 
    column_to_rownames("N") |> 
    as.matrix() |> 
    apply(2, as.integer)
  data_for_model$truncation_not_idx = 
    (data_for_model$truncation_down >= 0) |> 
    t() |> 
    as.vector()  |> 
    which() |>
    intersect(data_for_model$user_forced_truncation_not_idx) |>
    sort()
  data_for_model$TNS = length(data_for_model$truncation_not_idx)
  
  data_for_model$truncation_not_idx_minimal = 
    truncation_df |> 
    select(N, M, truncation_down) |> 
    spread(M, truncation_down) |>
    mutate(row = row_number()) |>
    pivot_longer(
      cols = -c(N, row),
      names_to = "columns",
      values_to = "value"
    ) |>
    filter(value == -1) |>
    select(row, columns) |>
    mutate(columns = as.integer(columns)) |> 
    as.matrix()
  
  data_for_model$TNIM = nrow(data_for_model$truncation_not_idx_minimal)
  
  message("sccomp says: outlier identification - step 1/2")
  
  my_quantile_step_2 = 1 - (0.1 / data_for_model$N)
  
  # This step gets the credible interval to control for within-category false positive rate
  # We want a category-wise false positive rate of 0.1, and we have to correct for how many samples we have in each category
  CI_step_2 = (1-my_quantile_step_2) / 2 * 2
  
  fit2 =
    data_for_model |>
    fit_model(
      "glm_multi_beta_binomial",
      cores = cores,
      quantile = my_quantile_step_2,
      inference_method = inference_method,
      output_directory = output_directory,
      verbose = verbose,
      seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pars = c("beta", "alpha", "prec_coeff", "prec_sd", "alpha_normalised", "random_effect", "random_effect_2"),
      sig_figs = sig_figs,
      cache_stan_model = cache_stan_model,
      ...
    )
  
  rng2 = mod_rng |> sample_safe(
    
    generate_quantities_fx,
    fit2$draws(format = "matrix"),
    
    # This is for the new data generation with selected factors to do adjustment
    data = data_for_model |> c(list(
      
      # Add subset of coefficients
      X_original = data_for_model$X,
      N_original = data_for_model$N,
      length_X_which = ncol(data_for_model$X),
      length_XA_which = ncol(data_for_model$XA),
      X_which = seq_len(ncol(data_for_model$X)) |> as.array(),
      XA_which = seq_len(ncol(data_for_model$Xa)) |> as.array(),
      
      # Random intercept common variable between grouping 1 and 2
      ncol_X_random_eff_new = ncol(data_for_model$X_random_effect) |> c(ncol(data_for_model$X_random_effect_2) ), # I could put this in the intial data
      length_X_random_effect_which = ncol(data_for_model$X_random_effect) |> c(ncol(data_for_model$X_random_effect_2)),
      
      # Grouping 1
      X_random_effect_which = seq_len(ncol(data_for_model$X_random_effect)) |> as.array(),
      
      # Grouping 2 - Random intercept DUPLICATED
      X_random_effect_which_2 = seq_len(ncol(data_for_model$X_random_effect_2)) |> as.array(),
      
      # Initialize unseen random effect variables
      ncol_X_random_eff_unseen = c(0, 0),
      X_random_effect_unseen = matrix(0, nrow = nrow(data_for_model$X), ncol = 0),
      X_random_effect_2_unseen = matrix(0, nrow = nrow(data_for_model$X), ncol = 0),
      
      create_intercept = FALSE
      
    )),
    
    parallel_chains = ifelse(inference_method %in% c("variational", "pathfinder"), 1, fit2$num_chains()), 
    threads_per_chain = cores,
    sig_figs = sig_figs
    
  )
  
  rng2_summary = 
    summary_to_tibble(rng2, "counts", "N", "M", probs = c(CI_step_2, 0.5, 1-CI_step_2)) 
  
  column_quantile_names = rng2_summary |> colnames() |> str_subset("\\%") |> _[c(1,3)]
  
  rng2_summary = 
    rng2_summary |>
    
    # !!! THIS COMMAND RELIES ON POSITION BECAUSE IT'S NOT TRIVIAL TO MATCH
    # !!! COLUMN NAMES BASED ON LIMITED PRECISION AND/OR PERIODICAL QUANTILES
    rename(
      .lower := !!as.symbol(column_quantile_names[1]) ,
      .median = `50%`,
      .upper := !!as.symbol(column_quantile_names[2])
    ) |>
    nest(data = -N) |>
    mutate(!!.sample := rownames(data_for_model$X)) |>
    unnest(data) |>
    nest(data = -M) |>
    mutate(!!.cell_group := colnames(data_for_model$y)) |>
    unnest(data) 
  
  # Detect outliers
  truncation_df2 =
    .data |>
    left_join(rng2_summary,
              by = c(quo_name(.sample), quo_name(.cell_group))
    ) |>
    
    # Add truncation
    mutate(   truncation_down = .lower,   truncation_up =  .upper) |>
    
    # Add outlier stats
    mutate( outlier = !(!!.count >= .lower & !!.count <= .upper) ) |>
    nest(data = -M) |>
    mutate(contains_outliers = map_lgl(data, ~ .x |> filter(outlier) |> nrow() > 0)) |>
    unnest(data) |>
    
    mutate(
      truncation_down = case_when( outlier ~ -1, TRUE ~ truncation_down),
      truncation_up = case_when(outlier ~ -1, TRUE ~ truncation_up)
    )
  
  data_for_model$truncation_up = truncation_df2 |> 
    select(N, M, truncation_up) |> 
    pivot_wider(names_from = M, values_from = truncation_up) |> 
    column_to_rownames("N") |> 
    as.matrix() |> 
    apply(2, as.integer)
  data_for_model$truncation_down = truncation_df2 |> 
    select(N, M, truncation_down) |> 
    pivot_wider(names_from = M, values_from = truncation_down) |> 
    column_to_rownames("N") |> 
    as.matrix() |> 
    apply(2, as.integer)
  data_for_model$truncation_not_idx = 
    (data_for_model$truncation_down >= 0) |> 
    t() |> 
    as.vector()  |> 
    which() |>
    intersect(data_for_model$user_forced_truncation_not_idx) |>
    sort()
  data_for_model$TNS = length(data_for_model$truncation_not_idx)
  
  # LOO
  data_for_model$enable_loo = TRUE & enable_loo
  
  message("sccomp says: outlier-free model fitting - step 2/2")
  
  # Print design matrix
  message(sprintf("sccomp says: the composition design matrix has columns: %s", data_for_model$X |> colnames() |> paste(collapse=", ")))
  message(sprintf("sccomp says: the variability design matrix has columns: %s", data_for_model$Xa |> colnames() |> paste(collapse=", ")))
  
  fit3 =
    data_for_model |>
    # Run the first discovery phase with permissive false discovery rate
    fit_model(
      "glm_multi_beta_binomial",
      cores = cores,
      quantile = CI,
      inference_method = inference_method,
      output_directory = output_directory,
      verbose = verbose, 
      seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pars = c("beta", "alpha", "prec_coeff","prec_sd",   "alpha_normalised", "random_effect", "random_effect_2", "log_lik"),
      cache_stan_model = cache_stan_model,
      ...
    )
  
  # # Make the fit standalone
  # temp_rds_file <- tempfile(fileext = ".rds")
  # fit3$save_object(file = temp_rds_file) 
  # fit3 = readRDS(temp_rds_file)
  # file.remove(temp_rds_file)
  
  # Create a dummy tibble
  tibble() |>
    # Attach association mean concentration
    add_attr(fit3, "fit") |>
    add_attr(data_for_model, "model_input") |>
    add_attr(.sample, ".sample") |>
    add_attr(.cell_group, ".cell_group") |>
    add_attr(.count, ".count") |>
    
    add_attr(formula_composition, "formula_composition") |>
    add_attr(formula_variability, "formula_variability") |>
    add_attr(parse_formula(formula_composition), "factors" ) |> 
    add_attr(inference_method, "inference_method" ) |> 
    
    # Add count data as attribute
    add_attr(.data, "count_data") |>
    
    # Add class to the tbl
    add_class("sccomp_tbl") |> 
    
    # Print estimates
    sccomp_test() |>
    
    # drop hypothesis testing as the estimation exists without probabilities.
    # For hypothesis testing use sccomp_test
    select(-contains("_FDR"), -contains("_pH0")) |> 
    
    add_attr(noise_model, "noise_model") 
  
  
}