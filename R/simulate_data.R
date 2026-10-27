#' sccomp_simulate
#'
#' @description Simulate counts from a fitted sccomp model, optionally replacing
#'   the abundance coefficients and the mean-variability association with
#'   user-chosen values. Variability is not taken from the posterior: for every
#'   variability term it is derived from the abundance coefficient it is
#'   associated with, `alpha = beta * slope + intercept`, plus the fitted
#'   residual scatter around that trend (see `variability_scatter`), so
#'   changing the abundance coefficients changes the simulated variability as
#'   the model expects. Random effects and the remaining parameters come from
#'   the posterior draws.
#'
#' @import dplyr
#' @importFrom rlang quo_is_null sym
#' @importFrom parallel detectCores
#'
#' @param fit The result of sccomp_estimate execution. This is used for sampling from real-data properties.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment. It can be a sub-formula of the fitted one.
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment. If NULL, the fitted variability formula is used.
#' @param new_data A tibble including sample-specific columns (sample identifier and factor columns from formula). If coefficients are provided separately, cell_group column is optional. Otherwise, should include cell_group column. If NULL, the data used for the fit is used.
#' @param coefficients A data frame/tibble with cell-type specific coefficients. Must contain a column matching the cell_group column name with every fitted cell group, and columns matching the design matrix column names (e.g., "(Intercept)", "typecancer"). Design columns that are absent are set to zero. For each design column, the coefficients must sum to zero across cell groups. If NULL, the posterior coefficients are used.
#' @param mean_dispersion_slope Optional slope of the mean-variability association. Either a single number, applied to every variability term, or a vector named by variability design column (e.g. `c("(Intercept)" = -0.7)`), in which case unnamed terms keep their posterior slope. If NULL, the posterior slopes are used, or 0 for a model fitted with `exclude_mean_variability_association = TRUE`.
#' @param mean_dispersion_intercept Optional intercept of the mean-variability association, in the same forms as `mean_dispersion_slope`. If NULL, the posterior intercepts are used, except for terms whose slope is overridden, which take their posterior-mean intercept so that only the slope changes.
#' @param .sample A column name as symbol. The sample identifier. If NULL, taken from the fit.
#' @param .cell_group A column name as symbol. The cell_group identifier. If NULL, taken from the fit.
#' @param variability_multiplier A real scalar. This can be used for artificially increasing the variability of the simulation for benchmarking purposes.
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#' @param cores Integer, the number of cores to be used for parallel calculations.
#' @param sig_figs Number of significant figures to use for Stan model output. Default is 9.
#' @param cache_stan_model A character string specifying the cache directory for compiled Stan models. 
#'                        The sccomp version will be automatically appended to ensure version isolation.
#'                        Default is `sccomp_stan_models_cache_dir` which points to `~/.sccomp_models`.
#' @param variability_scatter Logical. If `TRUE` (default), each cell group's variability is drawn around the mean-variability trend with the fitted residual standard deviation (`prec_sd`), as fitted cell groups scatter around it. If `FALSE`, every cell group sits exactly on the trend.
#' @param .coefficients (Deprecated) The column names for coefficients in new_data, for example, c(b_0, b_1). Use the 'coefficients' parameter instead. This parameter is placed last for backward compatibility.
#' 
#' @return The rows of `new_data` (expanded to every cell group if `new_data`
#'   has no cell_group column), with the following columns added:
#' \itemize{
#'   \item \strong{generated_proportions} - A numeric column representing the generated proportions from the simulation.
#'   \item \strong{generated_counts} - An integer column representing the generated cell counts from the simulation.
#'   \item \strong{replicate} - An integer column representing the replicate number for each draw from the posterior distribution.
#' }
#'
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
#'
#' @export
#'
#' @examples
#'
#' # print("cmdstanr is needed to run this example.")
#' # Note: Before running the example, ensure that the 'cmdstanr' package is installed:
#' # install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))
#'
#' # \donttest{
#' #   if (instantiate::stan_cmdstan_exists()) {
#' #     data("counts_obj")
#' #     library(dplyr)
#'
#' #     estimate = sccomp_estimate(
#' #       counts_obj,
#' #       ~ type, ~1, "sample", "cell_group", "count",
#' #       cores = 1
#' #     )
#'
#' #     # Coefficients for every cell group. All zero for simplicity.
#' #     coefficients = counts_obj |>
#' #       distinct(cell_group) |>
#' #       mutate(`(Intercept)` = 0, typecancer = 0)
#'
#' #     # Simulate data
#' #     sccomp_simulate(
#' #       estimate, ~type, ~1,
#' #       new_data = counts_obj |> distinct(sample, type),
#' #       coefficients = coefficients
#' #     )
#' #   }
#' # }
sccomp_simulate <- function(fit,
                            formula_composition,
                            formula_variability = NULL,
                            new_data = NULL,
                            coefficients = NULL,
                            mean_dispersion_slope = NULL,
                            mean_dispersion_intercept = NULL,
                            .sample = NULL,
                            .cell_group = NULL,
                            variability_multiplier = 1,
                            number_of_draws = 1,
                            mcmc_seed = sample_seed(),
                            cores = detectCores(),
                            sig_figs = 9,
                            cache_stan_model = sccomp_stan_models_cache_dir,
                            variability_scatter = TRUE,
                            .coefficients = NULL) {
  
  # Run the function
  check_and_install_cmdstanr()
  
  UseMethod("sccomp_simulate", fit)
}

#' @export
#'
sccomp_simulate.sccomp_tbl = function(fit,
                                      formula_composition,
                                      formula_variability = NULL,
                                      new_data = NULL,
                                      coefficients = NULL,
                                      mean_dispersion_slope = NULL,
                                      mean_dispersion_intercept = NULL,
                                      .sample = NULL,
                                      .cell_group = NULL,
                                      variability_multiplier = 1,
                                      number_of_draws = 1,
                                      mcmc_seed = sample_seed(),
                                      cores = detectCores(),
                                      sig_figs = 9,
                                      cache_stan_model = sccomp_stan_models_cache_dir,
                                      variability_scatter = TRUE,
                                      .coefficients = NULL) {
  
  .coefficients = enquo(.coefficients)
  
  if(!quo_is_null(.coefficients))
    lifecycle::deprecate_warn(
      "2.5.3",
      "sccomp_simulate(.coefficients)",
      details = "sccomp says: please use the `coefficients` argument with a data frame/tibble instead. See ?sccomp_simulate for details."
    )
  
  simulate_counts(
    fit,
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = new_data,
    coefficients = coefficients,
    coefficient_columns = if(quo_is_null(.coefficients)) NULL else quo_names(.coefficients),
    mean_dispersion_slope = mean_dispersion_slope,
    mean_dispersion_intercept = mean_dispersion_intercept,
    .sample = !!enquo(.sample),
    .cell_group = !!enquo(.cell_group),
    variability_multiplier = variability_multiplier,
    number_of_draws = number_of_draws,
    mcmc_seed = mcmc_seed,
    sig_figs = sig_figs,
    cache_stan_model = cache_stan_model,
    variability_scatter = variability_scatter
  )
}

#' Simulate counts and join them onto the requested rows
#'
#' @description Shared by `sccomp_simulate()` and the deprecated
#'   `simulate_data()`, so the latter does not also warn about `.coefficients`.
#'
#' @keywords internal
#' @noRd
simulate_counts = function(fit,
                           formula_composition,
                           formula_variability = NULL,
                           new_data = NULL,
                           coefficients = NULL,
                           coefficient_columns = NULL,
                           mean_dispersion_slope = NULL,
                           mean_dispersion_intercept = NULL,
                           .sample = NULL,
                           .cell_group = NULL,
                           variability_multiplier = 1,
                           number_of_draws = 1,
                           mcmc_seed = sample_seed(),
                           sig_figs = 9,
                           cache_stan_model = sccomp_stan_models_cache_dir,
                           variability_scatter = TRUE) {
  
  .sample = column_quo_or_default(enquo(.sample), attr(fit, ".sample"))
  .cell_group = column_quo_or_default(enquo(.cell_group), attr(fit, ".cell_group"))
  
  if(is.null(formula_variability)) formula_variability = attr(fit, "formula_variability")
  if(is.null(new_data)) new_data = attr(fit, "count_data")
  
  rng = simulate_data_rng(
    fit,
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = new_data,
    coefficients = coefficients,
    coefficient_columns = coefficient_columns,
    mean_dispersion_slope = mean_dispersion_slope,
    mean_dispersion_intercept = mean_dispersion_intercept,
    .sample = !!.sample,
    .cell_group = !!.cell_group,
    variability_multiplier = variability_multiplier,
    number_of_draws = number_of_draws,
    mcmc_seed = mcmc_seed,
    sig_figs = sig_figs,
    cache_stan_model = cache_stan_model,
    variability_scatter = variability_scatter
  )
  
  # Without a cell_group column, every sample is simulated for every cell group
  if(!quo_name(.cell_group) %in% colnames(new_data))
    new_data = new_data |> cross_join(tibble(!!.cell_group := rng$cell_group_names))
  
  generated =
    rng$fit |>
    parse_generated_quantities(number_of_draws = rng$number_of_draws) |>
    nest(data = -N) |>
    arrange(N) |>
    mutate(!!.sample := rng$sample_names) |>
    unnest(data) |>
    nest(data = -M) |>
    arrange(M) |>
    mutate(!!.cell_group := rng$cell_group_names) |>
    unnest(data) |>
    select(-N, -M)
  
  new_data |>
    left_join(generated, by = c(quo_name(.sample), quo_name(.cell_group)))
}

#' DEPRECATED: simulate_data
#'
#' @description This function is DEPRECATED. Please use \code{\link{sccomp_simulate}} instead.
#'
#' @param .data A tibble including a cell_group name column | sample name column | read counts column | factor columns | Pvalue column | a significance column
#' @param .estimate_object The result of sccomp_estimate execution. This is used for sampling from real-data properties.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_group identifier
#' @param .coefficients The column names for coefficients, for example, c(b_0, b_1). Use the `coefficients` argument of \code{sccomp_simulate()} instead.
#' @param variability_multiplier A real scalar. This can be used for artificially increasing the variability of the simulation for benchmarking purposes.
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#' @param cores Integer, the number of cores to be used for parallel calculations.
#' @param sig_figs Number of significant figures to use for Stan model output. Default is 9.
#' @param cache_stan_model A character string specifying the cache directory for compiled Stan models. 
#'                        The sccomp version will be automatically appended to ensure version isolation.
#'                        Default is `sccomp_stan_models_cache_dir` which points to `~/.sccomp_models`.
#'
#' @return See \code{\link{sccomp_simulate}}.
#'
#' @export
#'
simulate_data <- function(.data,
                          .estimate_object,
                          formula_composition,
                          formula_variability = NULL,
                          .sample = NULL,
                          .cell_group = NULL,
                          .coefficients = NULL,
                          variability_multiplier = 1,
                          number_of_draws = 1,
                          mcmc_seed = sample_seed(),
                          cores = detectCores(),
                          sig_figs = 9,
                          cache_stan_model = sccomp_stan_models_cache_dir) {
  
  lifecycle::deprecate_warn(
    "2.5.3",
    "sccomp::simulate_data()",
    "sccomp::sccomp_simulate()"
  )
  
  check_and_install_cmdstanr()
  
  .coefficients = enquo(.coefficients)
  
  simulate_counts(
    .estimate_object,
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = .data,
    coefficient_columns = if(quo_is_null(.coefficients)) NULL else quo_names(.coefficients),
    .sample = !!enquo(.sample),
    .cell_group = !!enquo(.cell_group),
    variability_multiplier = variability_multiplier,
    number_of_draws = number_of_draws,
    mcmc_seed = mcmc_seed,
    sig_figs = sig_figs,
    cache_stan_model = cache_stan_model
  )
}

#' Use a column argument, or fall back to the one recorded on the fit
#'
#' @param column_quo Quosure of the argument as supplied, a symbol or a string
#' @param default Quosure recorded on the fit
#'
#' @keywords internal
#' @noRd
column_quo_or_default = function(column_quo, default) {
  if(!quo_is_null(column_quo)) return(quo(!!sym(quo_name(column_quo))))
  if(is.null(default)) stop("sccomp says: the column could not be found in the fit; please provide it.")
  default
}

#' Generate simulated counts from a fitted sccomp model
#'
#' @description
#' Builds the simulation design with `prepare_replicate_data()`, turns the
#' optional coefficient and mean-variability overrides into Stan data, and
#' runs `glm_multi_beta_binomial_simulate_data` over posterior draws.
#'
#' @param coefficients Cell-type by design-column table, or NULL
#' @param coefficient_columns Deprecated form: names of columns of `new_data`
#'   holding one coefficient per design column, or NULL
#'
#' @return A list with the generated-quantities `fit`, the `sample_names` and
#'   `cell_group_names` indexing its `N` and `M` dimensions, and the
#'   `number_of_draws` actually used.
#'
#' @keywords internal
#' @noRd
simulate_data_rng = function(fit,
                             formula_composition,
                             formula_variability,
                             new_data,
                             coefficients = NULL,
                             coefficient_columns = NULL,
                             mean_dispersion_slope = NULL,
                             mean_dispersion_intercept = NULL,
                             .sample,
                             .cell_group,
                             variability_multiplier = 1,
                             number_of_draws = 1,
                             mcmc_seed = sample_seed(),
                             sig_figs = 9,
                             cache_stan_model = sccomp_stan_models_cache_dir,
                             variability_scatter = TRUE) {
  
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  
  model_input = attr(fit, "model_input")
  cell_group_names = colnames(model_input$y)
  
  if(!quo_name(.sample) %in% colnames(new_data))
    stop(sprintf("sccomp says: `new_data` must contain the sample column `%s`.", quo_name(.sample)))
  
  # Simulated cell groups must be fitted ones, since their parameters come from the fit
  if(quo_name(.cell_group) %in% colnames(new_data)) {
    unknown_cell_groups = new_data |> pull(!!.cell_group) |> as.character() |> unique() |> setdiff(cell_group_names)
    if(length(unknown_cell_groups) > 0)
      stop(sprintf(
        "sccomp says: the cell groups %s are not present in the fitted model (M = %d).",
        paste(unknown_cell_groups, collapse = ", "), length(cell_group_names)
      ))
  }
  
  sample_data =
    new_data |>
    select(-any_of(c(quo_name(.cell_group), coefficient_columns))) |>
    .subset(!!.sample)
  
  prepared_data = prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect_slots = lapply(seq_len(N_RE_SLOTS), function(k)
      model_input[[paste0("X_random_effect_", k)]]),
    .sample = !!.sample,
    .cell_group = !!.cell_group,
    .count = !!attr(fit, ".count"),
    formula_composition = formula_composition,
    original_formula_composition = attr(fit, "formula_composition"),
    formula_variability = formula_variability,
    new_data = sample_data,
    original_count_data = attr(fit, "count_data") |> .subset(!!attr(fit, ".sample")),
    smooth_results = get_smooth_results(fit)
  )
  
  if(prepared_data$create_intercept)
    stop("sccomp says: simulating with an intercept from an intercept-free fit is not supported. Please use a composition formula with the same intercept as the fitted model.")
  
  design_columns = colnames(prepared_data$X)
  
  beta_provided =
    if(!is.null(coefficients))
      coefficients_to_beta(coefficients, quo_name(.cell_group), design_columns, cell_group_names)
    else if(!is.null(coefficient_columns))
      coefficient_columns_to_beta(new_data, coefficient_columns, quo_name(.cell_group), design_columns, cell_group_names)
  
  if(!is.null(beta_provided)) check_coefficients_sum_to_zero(beta_provided)
  
  variability_terms = colnames(model_input$Xa)
  slope = expand_variability_term_values(mean_dispersion_slope, variability_terms, "mean_dispersion_slope")
  intercept = expand_variability_term_values(mean_dispersion_intercept, variability_terms, "mean_dispersion_intercept")
  
  # A new slope alone keeps the intercept fixed at its posterior mean
  fix_intercept = slope$provided == 1 & intercept$provided == 0
  if(any(fix_intercept)) {
    posterior_mean_intercept = attr(fit, "fit")$summary("prec_intercept_1", "mean")$mean
    intercept$values[fix_intercept] = posterior_mean_intercept[fix_intercept]
    intercept$provided[fix_intercept] = 1L
  }
  
  # Samples take exposures drawn from the fitted ones
  exposure = sample(model_input$exposure, size = prepared_data$N, replace = TRUE) |> as.integer() |> as.array()
  
  stan_data = c(
    model_input[c(
      "M", "C", "A",
      "bimodal_mean_variability_association", "variability_to_composition_map",
      "is_random_effect", "ncol_X_random_eff", "n_groups", "how_many_factors_in_random_design",
      paste0("group_factor_indexes_for_covariance_", seq_len(N_RE_SLOTS))
    )],
    prepared_data[c(
      "N", "X", "Xa", "X_which", "XA_which",
      "ncol_X_random_eff_new", "ncol_X_random_eff_unseen",
      paste0("X_random_effect_", seq_len(N_RE_SLOTS)),
      paste0("X_random_effect_", seq_len(N_RE_SLOTS), "_unseen"),
      paste0("X_random_effect_which_", seq_len(N_RE_SLOTS))
    )],
    list(
      exclude_mean_variability_association = as.integer(isTRUE(as.logical(model_input$exclude_mean_variability_association))),
      exposure = exposure,
      length_X_which = length(prepared_data$X_which),
      length_XA_which = length(prepared_data$XA_which),
      length_X_random_effect_which =
        map_int(seq_len(N_RE_SLOTS), ~ length(prepared_data[[paste0("X_random_effect_which_", .x)]])),
      user_provided_beta = as.integer(!is.null(beta_provided)),
      beta_provided =
        if(is.null(beta_provided)) matrix(0, length(design_columns), length(cell_group_names))
        else beta_provided,
      user_provided_slope = slope$provided,
      slope_provided = slope$values,
      user_provided_intercept = intercept$provided,
      intercept_provided = intercept$values,
      variability_multiplier = variability_multiplier,
      variability_scatter = as.integer(isTRUE(variability_scatter))
    )
  )
  
  fit_obj = attr(fit, "fit")
  
  # generate_quantities needs exactly the fitted parameter block
  parameter_variables = fit_obj$metadata()$model_params
  draws_matrix =
    if(!is.null(parameter_variables) && length(parameter_variables) > 0)
      fit_obj$draws(variables = parameter_variables, format = "matrix")
    else
      fit_obj$draws(format = "matrix")
  
  number_of_draws = min(number_of_draws, get_output_samples(fit_obj), nrow(draws_matrix))
  draws_matrix = draws_matrix[sample(seq_len(nrow(draws_matrix)), size = number_of_draws), , drop = FALSE]
  
  mod_rng = load_model("glm_multi_beta_binomial_simulate_data", cache_dir = cache_stan_model)
  
  list(
    fit = mod_rng |> sample_safe(
      generate_quantities_fx,
      draws_matrix,
      data = stan_data,
      seed = mcmc_seed,
      threads_per_chain = 1,
      sig_figs = sig_figs
    ),
    sample_names = rownames(prepared_data$X),
    cell_group_names = cell_group_names,
    number_of_draws = number_of_draws
  )
}

#' Coefficient matrix from a cell-type by design-column table
#'
#' @return Matrix with one row per design column and one column per fitted
#'   cell group. Design columns absent from the table are zero.
#'
#' @keywords internal
#' @noRd
coefficients_to_beta = function(coefficients, cell_group_column, design_columns, cell_group_names) {
  
  if(!cell_group_column %in% colnames(coefficients))
    stop(sprintf("sccomp says: `coefficients` must contain the cell group column `%s`.", cell_group_column))
  
  coefficient_cell_groups = coefficients[[cell_group_column]] |> as.character()
  
  if(anyDuplicated(coefficient_cell_groups))
    stop("sccomp says: `coefficients` must have one row per cell group.")
  
  check_coefficient_cell_groups(coefficient_cell_groups, cell_group_names)
  
  row_order = match(cell_group_names, coefficient_cell_groups)
  
  beta = matrix(0, length(design_columns), length(cell_group_names), dimnames = list(design_columns, cell_group_names))
  for(column in intersect(design_columns, colnames(coefficients)))
    beta[column, ] = coefficients[[column]][row_order]
  
  beta
}

#' Coefficient matrix from coefficient columns of `new_data` (deprecated form)
#'
#' @description Columns are matched to design columns by name, and otherwise by
#'   position, as `simulate_data()` always did.
#'
#' @keywords internal
#' @noRd
coefficient_columns_to_beta = function(new_data, coefficient_columns, cell_group_column, design_columns, cell_group_names) {
  
  if(!cell_group_column %in% colnames(new_data))
    stop(sprintf("sccomp says: `.coefficients` requires the cell group column `%s` in the data.", cell_group_column))
  
  coefficients =
    new_data |>
    distinct(across(all_of(c(cell_group_column, coefficient_columns)))) |>
    mutate(across(all_of(cell_group_column), as.character))
  
  if(anyDuplicated(coefficients[[cell_group_column]]))
    stop("sccomp says: the `.coefficients` columns must have one value per cell group.")
  
  check_coefficient_cell_groups(coefficients[[cell_group_column]], cell_group_names)
  
  row_order = match(cell_group_names, coefficients[[cell_group_column]])
  
  beta = matrix(0, length(design_columns), length(cell_group_names), dimnames = list(design_columns, cell_group_names))
  for(i in seq_along(design_columns)) {
    column =
      if(design_columns[i] %in% coefficient_columns) design_columns[i]
      else if(i <= length(coefficient_columns)) coefficient_columns[i]
      else NA
    if(!is.na(column)) beta[i, ] = coefficients[[column]][row_order]
  }
  
  beta
}

#' @keywords internal
#' @noRd
check_coefficient_cell_groups = function(coefficient_cell_groups, cell_group_names) {
  
  unknown = setdiff(coefficient_cell_groups, cell_group_names)
  if(length(unknown) > 0)
    stop(sprintf(
      "sccomp says: the coefficients include cell groups that are not present in the fitted model (M = %d): %s.",
      length(cell_group_names), paste(unknown, collapse = ", ")
    ))
  
  missing = setdiff(cell_group_names, coefficient_cell_groups)
  if(length(missing) > 0)
    stop(sprintf(
      "sccomp says: the coefficients must cover every fitted cell group, since they sum to zero across them. Missing: %s.",
      paste(missing, collapse = ", ")
    ))
}

#' Coefficients of a compositional model must sum to zero across cell groups
#'
#' @param beta Matrix with one row per design column
#'
#' @keywords internal
#' @noRd
check_coefficients_sum_to_zero = function(beta) {
  
  tolerance = .Machine$double.eps * 100 * max(max(abs(beta)), 1) * ncol(beta)
  row_sums = rowSums(beta)
  invalid = which(abs(row_sums) > tolerance)
  
  if(length(invalid) > 0)
    stop(sprintf(
      paste(
        "sccomp says: user-provided coefficients for design column '%s' do not sum to zero (sum = %s).",
        "Coefficients in compositional models must sum to zero across cell groups.",
        "Please normalise your coefficients, e.g. by subtracting their mean."
      ),
      rownames(beta)[invalid[1]], signif(row_sums[invalid[1]], 6)
    ))
}

#' Expand a mean-variability override to one value per variability term
#'
#' @param value NULL, a single number for every term, or a vector named by
#'   variability design column
#' @param variability_terms Variability design column names of the fit
#' @param argument Argument name, for error messages
#'
#' @return A list with `values` and `provided` (0/1), one entry per term
#'
#' @keywords internal
#' @noRd
expand_variability_term_values = function(value, variability_terms, argument) {
  
  A = length(variability_terms)
  values = rep(0, A)
  provided = rep(0L, A)
  
  if(!is.null(value)) {
    
    if(!is.numeric(value) || anyNA(value))
      stop(sprintf("sccomp says: `%s` must be numeric.", argument))
    
    if(is.null(names(value))) {
      if(length(value) != 1)
        stop(sprintf(
          "sccomp says: `%s` must be a single number, or a vector named by variability term (%s).",
          argument, paste(variability_terms, collapse = ", ")
        ))
      values[] = value
      provided[] = 1L
    } else {
      unknown = setdiff(names(value), variability_terms)
      if(length(unknown) > 0)
        stop(sprintf(
          "sccomp says: `%s` names terms that are not in the variability design: %s. Available terms: %s.",
          argument, paste(unknown, collapse = ", "), paste(variability_terms, collapse = ", ")
        ))
      idx = match(names(value), variability_terms)
      values[idx] = unname(value)
      provided[idx] = 1L
    }
  }
  
  list(values = as.array(values), provided = as.array(provided))
}
