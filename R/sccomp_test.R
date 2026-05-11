#' sccomp_test
#'
#' @description This function test contrasts from a sccomp result.
#'
#' @param .data A tibble. The result of sccomp_estimate.
#' @param contrasts A vector of character strings. For example if your formula is `~ 0 + treatment` and the factor treatment has values `yes` and `no`, your contrast could be "constrasts = c(treatmentyes - treatmentno)".
#' @param percent_false_positive A real between 0 and 100 non included. This used to identify outliers with a specific false positive rate.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' @param pass_fit A boolean. Whether to pass the Stan fit as attribute in the output. Because the Stan fit can be very large, setting this to FALSE can be used to lower the memory imprint to save the output.
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
#'   \item c_FDR - False-discovery rate of the null hypothesis (no difference) for a composition (c).
#'   \item c_n_eff - Effective sample size - the number of independent draws in the sample, the higher the better.
#'   \item c_R_k_hat - R statistic, a measure of chain equilibrium, should be within 0.05 of 1.0.
#'   \item v_lower - Lower (2.5%) quantile of the posterior distribution for a variability (v) parameter.
#'   \item v_effect - Mean of the posterior distribution for a variability (v) parameter.
#'   \item v_upper - Upper (97.5%) quantile of the posterior distribution for a variability (v) parameter.
#'   \item v_pH0 - Probability of the null hypothesis (no difference) for a variability (v).
#'   \item v_FDR - False-discovery rate of the null hypothesis (no difference) for a variability (v).
#'   \item v_n_eff - Effective sample size for a variability (v) parameter.
#'   \item v_R_k_hat - R statistic for a variability (v) parameter.
#'   \item count_data - Nested input count data.
#' }
#'
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
#'
#' @export
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
#'     estimates = sccomp_estimate(
#'       counts_obj,
#'       ~ 0 + type, ~1, "sample", "cell_group", "count",
#'       cores = 1
#'     ) |>
#'     sccomp_test("typecancer - typebenign")
#'   }
#' }
#'
sccomp_test <- function(.data,
                        contrasts = NULL,
                        percent_false_positive = 5,
                        test_composition_above_logit_fold_change = 0.1,
                        pass_fit = TRUE) {
  
  # Run the function
  check_and_install_cmdstanr()
  
  UseMethod("sccomp_test", .data)
}


#'
#' @importFrom dplyr any_of
#'
#' @export
#'
#'
sccomp_test.sccomp_tbl = function(.data,
                                  contrasts = NULL,
                                  percent_false_positive = 5,
                                  test_composition_above_logit_fold_change = 0.1,
                                  pass_fit = TRUE){
  

  .sample = .data |>  attr(".sample")
  .cell_group = .data |>  attr(".cell_group")
  .count = .data |>  attr(".count")
  model_input = .data |> attr("model_input")
  truncation_df2 =  .data |>  attr("truncation_df2")
  inference_method = .data |>  attr("inference_method")
  
  # sccomp_test always computes pH0/FDR from posterior draws.
  result  <- 
    get_abundance_contrast_draws(.data, contrasts) |>
      draws_to_statistics(
        percent_false_positive / 100,
        test_composition_above_logit_fold_change,
        !!.cell_group,
        "c_"
      )
  
  variability_draws <- get_variability_contrast_draws(.data, contrasts)

  # If I have variability draws for those contrasts, compute the variability CI
  if ("parameter" %in% colnames(variability_draws)) {
    result = 
    result |>
    left_join(
      variability_draws |>
        draws_to_statistics(
          percent_false_positive / 100,
          test_composition_above_logit_fold_change,
          !!.cell_group,
          "v_"
        )
    )
  }

  
  factor_parameter_dictionary =
    model_input$factor_parameter_dictionary |>
    select(`factor`, design_matrix_col)
  
  # Merge and parse
  result =
   result |>
    
    # Keep factor labels in the result schema; plotting/subsetting helpers rely on this.
    left_join(factor_parameter_dictionary, by = c("parameter" = "design_matrix_col")) |>
    select(!!.cell_group, parameter, `factor`, everything(), -M) 
      
  if(pass_fit)
    result =
    result |>
    add_attr(.data |> attr("fit") , "fit")
  
  result |>
    
    # TEMPORARILY DROPPING KHAT
    # select(-contains("n_eff"), -contains("_hat")) |> 
    
    add_attr(test_composition_above_logit_fold_change, "test_composition_above_logit_fold_change") |>
    
    # Attach association mean concentration
    add_attr(.data |> attr("model_input") , "model_input") |>
    add_attr(.data |> attr("truncation_df2"), "truncation_df2") |>
    add_attr(.data |> attr("noise_model") , "noise_model") |>
    
    add_attr(.sample, ".sample") |>
    add_attr(.cell_group, ".cell_group") |>
    add_attr(.count, ".count") |>
    
    add_attr(.data |> attr("formula_composition"), "formula_composition") |>
    add_attr(.data |> attr("formula_variability"), "formula_variability") |>
    add_attr(inference_method, "inference_method" ) |> 
    
    # Add count data as attribute
    add_attr(.data |> attr("count_data"), "count_data") |>
    add_attr(.data |> attr("outliers"), "outliers") |>
    
  
    # Add class to the tbl
    add_class("sccomp_tbl") 
}

#' Build per-parameter posterior summaries (mean, equal-tailed intervals, rhat, ESS)
#' for default reporting from `sccomp_estimate()` / `sccomp_remove_outliers()` without
#' loading full posterior draws. Does not compute pH0 / FDR (use [sccomp_test()]).
#'
#' `design_colnames` may be a subset of `full_design_colnames`; in that case only those
#' Stan elements are passed to `fit$summary()`.
#'
#' @keywords internal
#' @noRd
summarise_stan_matrix_for_estimate <- function(
    fit,
    model_input,
    stan_parameter,
    parameter_names,
    probs,
    prefix) {
  # Map Stan matrix column index M back to user-facing cell-group labels.
  cell_group_levels <- colnames(model_input$y)

  # Query Stan summaries for this parameter family (beta/random effects).
  summ <- summary_to_tibble(fit, stan_parameter, "C", "M", probs = probs)

  # Identify lower/upper quantile column names from the summary output.
  quantile_columns <- names(summ)[grepl("%$", names(summ))]
  ordered_quantile_columns <- quantile_columns |>
    tibble::enframe(value = "qcol") |>
    dplyr::mutate(pct = as.numeric(gsub("%$", "", qcol, perl = TRUE))) |>
    dplyr::arrange(pct) |>
    dplyr::pull(qcol)
  lower_quantile_column <- ordered_quantile_columns[[1]]
  upper_quantile_column <- ordered_quantile_columns[[2]]

  # Provide a stable mapping from C indices to design-matrix parameter names.
  par_by_C <-
    tibble::tibble(C = seq_along(parameter_names), parameter = parameter_names) |>
    dplyr::distinct()

  # Attach parameter labels to each C/M summary row.
  joined <- summ |>
    dplyr::left_join(par_by_C, by = "C")

  # Keep output column names configurable for c_/v_ style reuse.
  lower_col <- paste0(prefix, "lower")
  effect_col <- paste0(prefix, "effect")
  upper_col <- paste0(prefix, "upper")
  rhat_col <- paste0(prefix, "rhat")
  ess_bulk_col <- paste0(prefix, "ess_bulk")
  ess_tail_col <- paste0(prefix, "ess_tail")

  # Return compact, user-facing summary table with diagnostics.
  joined |>
    dplyr::transmute(
      cell_group = cell_group_levels[M],
      M,
      parameter,
      !!lower_col := .data[[lower_quantile_column]],
      !!effect_col := mean,
      !!upper_col := .data[[upper_quantile_column]],
      !!rhat_col := rhat,
      !!ess_bulk_col := ess_bulk,
      !!ess_tail_col := ess_tail
    )
}

#' Build variability summaries from R-side alpha normalisation draws.
#'
#' This mirrors `summarise_stan_matrix_for_estimate()` but computes quantiles for
#' alpha_normalised from derived draws, while borrowing convergence diagnostics
#' (rhat/ESS) from the corresponding `alpha` elements.
#'
#' @keywords internal
#' @noRd
summarise_alpha_normalised_for_estimate <- function(
    fit,
    model_input,
    design_colnames,
    probs,
    cell_group_colname,
    prefix) {
  # Map Stan matrix column index M back to user-facing cell-group labels.
  cell_group_levels <- colnames(model_input$y)

  # Match variability design columns to Stan C indices.
  C_idx <- match(design_colnames, colnames(model_input$XA))
  if (anyNA(C_idx)) {
    stop(
      "sccomp says: each variability design column must appear in model_input$XA.",
      call. = FALSE
    )
  }

  n_M <- length(cell_group_levels)
  g <- expand.grid(C = C_idx, M = seq_len(n_M), stringsAsFactors = FALSE)
  alpha_subset <- sprintf("alpha[%d,%d]", g$C, g$M)

  # Compute alpha_normalised summaries from derived R-side draws.
  draws_summary <- compute_alpha_normalised_draws(
    fit = fit,
    model_input = model_input,
    alpha_variable_subset = alpha_subset
  ) |>
    dplyr::group_by(C, M) |>
    dplyr::summarise(
      effect = mean(.value),
      lower_quantile = stats::quantile(.value, probs = probs[[1]], na.rm = TRUE),
      upper_quantile = stats::quantile(.value, probs = probs[[2]], na.rm = TRUE),
      .groups = "drop"
    )

  # Pull diagnostics from the corresponding base alpha parameters in Stan.
  alpha_diagnostics <- summary_to_tibble(fit, alpha_subset, "C", "M", probs = probs) |>
    dplyr::select(C, M, rhat, ess_bulk, ess_tail) |>
    dplyr::distinct()

  # Combine derived summaries, diagnostics, and design parameter labels.
  joined <- draws_summary |>
    dplyr::left_join(alpha_diagnostics, by = c("C", "M")) |>
    dplyr::left_join(
      tibble::tibble(C = C_idx, parameter = design_colnames) |>
        dplyr::distinct(),
      by = "C"
    )

  # Keep output column names configurable for v_ style summaries.
  lower_col <- paste0(prefix, "lower")
  effect_col <- paste0(prefix, "effect")
  upper_col <- paste0(prefix, "upper")
  rhat_col <- paste0(prefix, "rhat")
  ess_bulk_col <- paste0(prefix, "ess_bulk")
  ess_tail_col <- paste0(prefix, "ess_tail")

  # Return compact, user-facing summary table with diagnostics.
  joined |>
    dplyr::transmute(
      !!rlang::sym(cell_group_colname) := cell_group_levels[M],
      M,
      parameter,
      !!lower_col := lower_quantile,
      !!effect_col := effect,
      !!upper_col := upper_quantile,
      !!rhat_col := rhat,
      !!ess_bulk_col := ess_bulk,
      !!ess_tail_col := ess_tail
    )
}

#' @keywords internal
#' @noRd
sccomp_summarise_posterior_for_estimate <- function(
    fit,
    model_input,
    .cell_group,
    percent_false_positive = 5) {
  check_and_install_cmdstanr()

  cg <- rlang::quo_name(.cell_group)
  fp <- percent_false_positive / 100
  probs <- c(fp / 2, 1 - fp / 2)

  abundance_parts <- list(
    summarise_stan_matrix_for_estimate(
      fit = fit,
      model_input = model_input,
      stan_parameter = "beta",
      parameter_names = colnames(model_input$X),
      probs = probs,
      prefix = "c_"
    )
  )
  if (model_input$n_random_eff > 0) {
    abundance_parts <- c(
      abundance_parts,
      list(summarise_stan_matrix_for_estimate(
        fit = fit,
        model_input = model_input,
        stan_parameter = "random_effect",
        parameter_names = colnames(model_input$X_random_effect),
        probs = probs,
        prefix = "c_"
      ))
    )
  }
  if (model_input$n_random_eff > 1) {
    abundance_parts <- c(
      abundance_parts,
      list(summarise_stan_matrix_for_estimate(
        fit = fit,
        model_input = model_input,
        stan_parameter = "random_effect_2",
        parameter_names = colnames(model_input$X_random_effect_2),
        probs = probs,
        prefix = "c_"
      ))
    )
  }

  abundance <- dplyr::bind_rows(abundance_parts) |>
    dplyr::rename(!!cg := cell_group)

  variability <- summarise_alpha_normalised_for_estimate(
    fit = fit,
    model_input = model_input,
    design_colnames = colnames(model_input$XA),
    probs = probs,
    cell_group_colname = cg,
    prefix = "v_"
  ) |>
    dplyr::mutate(
      v_lower = -v_lower,
      v_effect = -v_effect,
      v_upper = -v_upper
    ) |>
    dplyr::rename(v_lower = v_upper, v_upper = v_lower)

  factor_parameter_dictionary <-
    model_input$factor_parameter_dictionary |>
    dplyr::select(`factor`, design_matrix_col)

  result <-
    abundance |>
    dplyr::left_join(variability, by = c(cg, "M", "parameter")) |>
    dplyr::left_join(
      factor_parameter_dictionary,
      by = c("parameter" = "design_matrix_col")
    ) |>
    dplyr::select(parameter, `factor`, dplyr::everything()) |>
    dplyr::select(!!.cell_group, dplyr::everything(), -M)

  result
}

#' Identify direct covariate contrasts for the summary-only shortcut.
#'
#' For now, the shortcut is enabled only when each contrast maps to exactly one
#' covariate and all covariates are present in the composition design matrix `X`.
#' Otherwise the draw-based contrast arithmetic is used.
#'
#' @keywords internal
#' @noRd
sccomp_identify_covariate_contrasts <- function(contrasts, model_input) {
  if (is.null(contrasts)) {
    return(NULL)
  }
  if (is.null(names(contrasts))) {
    names(contrasts) <- contrasts
  }

  contrast_mapping <- tibble::tibble(
    contrast = names(contrasts),
    design_param = vapply(
      contrasts,
      function(s) {
        extracted <- contrasts_to_parameter_list(s)
        if (length(extracted) != 1) return(NA_character_)
        extracted[[1]]
      },
      character(1)
    )
  )

  if (anyNA(contrast_mapping$design_param)) return(NULL)

  if (!all(contrast_mapping$design_param %in% colnames(model_input$X))) return(NULL)

  variab_ok <- all(contrast_mapping$design_param %in% colnames(model_input$XA))

  list(
    par_name = "beta",
    full_design = colnames(model_input$X),
    contrast_mapping = contrast_mapping,
    variab_ok = variab_ok
  )
}

#' Build Stan indexed-variable subset from contrast terms.
#'
#' @param contrasts Character vector of contrast expressions. Can be `NULL`.
#' @param design_columns Character vector of design-matrix column names for the target parameter block.
#' @param stan_parameter Character scalar Stan parameter name (for example, `"beta"` or `"alpha"`).
#' @param model_input Model input list containing `y`, used to infer the `M` index size.
#'
#' @keywords internal
#' @noRd
build_stan_parameter_subset <- function(contrasts, design_columns, stan_parameter, model_input) {
  design_columns <- as.character(design_columns)

  if (is.null(contrasts)) {
    return(
      tibble::tibble(
        parameter = design_columns,
        variable = rep(stan_parameter, length(design_columns))
      )
    )
  }

  candidate_terms <- contrasts_to_parameter_list(contrasts)
  if (is.null(candidate_terms) || length(candidate_terms) == 0) {
    return(tibble::tibble(parameter = character(), variable = character()))
  }

  candidate_terms <- unique(candidate_terms)
  candidate_terms <- candidate_terms[!is.na(candidate_terms) & candidate_terms != ""]
  if (length(candidate_terms) == 0) {
    return(tibble::tibble(parameter = character(), variable = character()))
  }

  matched <- intersect(candidate_terms, design_columns)
  if (length(matched) == 0) {
    return(tibble::tibble(parameter = character(), variable = character()))
  }

  n_M <- ncol(model_input$y)
  C_idx <- match(matched, design_columns)
  g <- expand.grid(C = C_idx, M = seq_len(n_M), stringsAsFactors = FALSE)

  tibble::tibble(
    parameter = matched[g$C],
    variable = sprintf("%s[%d,%d]", stan_parameter, g$C, g$M)
  )
}

# this can be helpful if we want to draw PCA with uncertainty
get_abundance_contrast_draws = function(.data, contrasts = NULL){
  
  # Define the variables as NULL to avoid CRAN NOTES
  X <- NULL
  .value <- NULL
  X_random_effect <- NULL
  .variable <- NULL
  y <- NULL
  M <- NULL
  khat <- NULL
  parameter <- NULL
  n_eff <- NULL
  R_k_hat <- NULL
  
  
  .cell_group = .data |>  attr(".cell_group")
  model_input <- .data |> attr("model_input")
  cell_index_map <-
    model_input %$%
    y %>%
    colnames() |>
    enframe(name = "M", value  = quo_name(.cell_group))
  
  # Beta
  beta_covariates = model_input %$% X |> colnames()

  beta_subset <- build_stan_parameter_subset(
    contrasts = contrasts,
    design_columns = beta_covariates,
    stan_parameter = "beta",
    model_input = model_input
  )
  beta_parameters <- beta_subset |> dplyr::pull("parameter") |> unique()
  beta_variable_subset <- beta_subset |> dplyr::pull("variable") |> unique()
  

    draws =
      .data |>
      attr("fit") %>%
      draws_to_tibble_x_y( beta_variable_subset, "C", "M") |>
    left_join(
      beta_covariates |> enframe(name = "C", value = "parameters_name"),
      by = "C"
    )  |>
    select(-C) |>
    pivot_wider(names_from = parameters_name, values_from = .value) |>
    select(-.variable)

  
  
  # Random effect
  
  random_effect_covariates = model_input %$% X_random_effect |> colnames()
  beta_random_effect_subset <- build_stan_parameter_subset(
    contrasts = contrasts,
    design_columns = random_effect_covariates,
    stan_parameter = "random_effect",
    model_input = model_input
  )
  beta_random_effect_parameters <- beta_random_effect_subset |> dplyr::pull("parameter") |> unique()
  beta_random_effect_variables <- beta_random_effect_subset |> dplyr::pull("variable") |> unique()
  
  if(
    .data |> attr("model_input") %$% n_random_eff > 0 &&
    (
      contrasts |> is.null() || 
      length(beta_random_effect_parameters) > 0
    )  
  ){
    
    
    beta_random_effect =
      .data |>
      attr("fit") %>%
      draws_to_tibble_x_y(beta_random_effect_variables,  "C",  "M"  ) 
    
    # Add last component
    other_group_random_effect = 
      beta_random_effect |> 
      with_groups(c(C, .chain, .iteration, .draw, .variable ), ~ .x |> summarise(.value = sum(.value))) |> 
      mutate(.value = -.value, M = beta_random_effect |> pull(M) |> max() + 1)
    
    
    beta_random_effect = 
      beta_random_effect |> 
      bind_rows( other_group_random_effect )
    
    
    # Reshape
    # Speed up if I have contrasts
    if(!contrasts |> is.null())
      beta_random_effect = 
      beta_random_effect |> 
      left_join(
        random_effect_covariates |> enframe(name = "C", value = "parameters_name"),
        by = "C"
      )  |> 
      filter(parameters_name %in% beta_random_effect_parameters) |> 
      select(-C) |> 
      pivot_wider(names_from = parameters_name, values_from = .value)
    
    else
      beta_random_effect = 
      beta_random_effect |>
      pivot_wider(names_from = C, values_from = .value) %>%
      setNames(colnames(.)[1:5] |> c(random_effect_covariates))

    # If I don't have fix nor 1st level random effect
    if(draws |> nrow() == 0)
      draws = select(beta_random_effect, -.variable)
    else 
      draws = draws |> 
      left_join(select(beta_random_effect, -.variable),
                by = c("M", ".chain", ".iteration", ".draw")
      )
    
  }  else {
    random_effect_covariates = ""
  }
  
  # Second random effect. IN THE FUTURE THIS WILL BE VECTORISED TO ARBUTRARY GRI+OUING
  random_effect_covariates_2 = model_input %$% X_random_effect_2 |> colnames()
  beta_random_effect_subset_2 <- build_stan_parameter_subset(
    contrasts = contrasts,
    design_columns = random_effect_covariates_2,
    stan_parameter = "random_effect_2",
    model_input = model_input
  )
  beta_random_effect_parameters_2 <- beta_random_effect_subset_2 |> dplyr::pull("parameter") |> unique()
  beta_random_effect_variables_2 <- beta_random_effect_subset_2 |> dplyr::pull("variable") |> unique()
  
  if(
    .data |> attr("model_input") %$% n_random_eff > 1 &&
    (
      contrasts |> is.null() || 
      length(beta_random_effect_parameters_2) > 0
    )
  ){
    
    beta_random_effect_2 =
      .data |>
      attr("fit") %>%
      draws_to_tibble_x_y(  beta_random_effect_variables_2,   "C",    "M"   ) 
    
    # Add last component
    other_group_random_effect = 
      beta_random_effect_2 |> 
      with_groups(c(C, .chain, .iteration, .draw, .variable ), ~ .x |> summarise(.value = sum(.value))) |> 
      mutate(.value = -.value, M = beta_random_effect_2 |> pull(M) |> max() + 1)
    
    beta_random_effect_2 = 
      beta_random_effect_2 |> 
      bind_rows( other_group_random_effect )
    
    # Reshape
    # Speed up if I have contrasts
    if(!contrasts |> is.null())
      beta_random_effect_2 = 
      beta_random_effect_2 |> 
      left_join(
        random_effect_covariates_2 |> enframe(name = "C", value = "parameters_name"),
        by = "C"
      )  |> 
      filter(parameters_name %in% beta_random_effect_parameters_2) |> 
      select(-C) |> 
      pivot_wider(names_from = parameters_name, values_from = .value)
    
    else
      beta_random_effect_2 = 
      beta_random_effect_2 |>
      pivot_wider(names_from = C, values_from = .value) %>%
      setNames(colnames(.)[1:5] |> c(random_effect_covariates_2))
    
    # If I don't have fix nor 1st level random effect
    if(draws |> nrow() == 0)
      draws = select(beta_random_effect_2, -.variable)
    else 
      draws = draws |> 
      left_join(select(beta_random_effect_2, -.variable),
                by = c("M", ".chain", ".iteration", ".draw")
      )
  } else {
    random_effect_covariates_2 = ""
  }
  
  
  # If I have constrasts calculate
  if(!is.null(contrasts))
    draws = 
    draws |> 
    mutate_from_expr_list(contrasts, ignore_errors = FALSE) |>
    select(- any_of(c(beta_covariates, random_effect_covariates) |> setdiff(contrasts)) )
  

  draws = draws |>
    pivot_longer(-c(1:4), names_to = "parameter", values_to = ".value") |>
    
    # Reorder because pivot long is bad
    mutate(parameter = parameter |> fct_relevel(colnames(draws)[-c(1:4)])) |>
    arrange(parameter)

     # Add cell name
  draws = draws |> 
    left_join(cell_index_map, by = "M" ) %>%
    select(!!.cell_group, everything())

  draws
  
}

#' @importFrom forcats fct_relevel
#' @noRd
get_variability_contrast_draws = function(.data, contrasts){
  
  # Define the variables as NULL to avoid CRAN NOTES
  XA <- NULL
  .value <- NULL
  y <- NULL
  M <- NULL
  khat <- NULL
  parameter <- NULL
  n_eff <- NULL
  R_k_hat <- NULL
  
  .cell_group = .data |>  attr(".cell_group")
  model_input <- .data |> attr("model_input")
  cell_index_map <-
    model_input %$%
    y %>%
    colnames() |>
    enframe(name = "M", value  = quo_name(.cell_group))
  
  variability_covariates = model_input %$% XA |> colnames()
  alpha_subset <- build_stan_parameter_subset(
    contrasts = contrasts,
    design_columns = variability_covariates,
    stan_parameter = "alpha",
    model_input = model_input
  )
  alpha_parameters <- alpha_subset |> dplyr::pull("parameter") |> unique()
  alpha_variable_subset <- alpha_subset |> dplyr::pull("variable") |> unique()
  if (length(alpha_variable_subset) == 0) alpha_variable_subset <- NULL
  if (!is.null(contrasts) && length(alpha_parameters) == 0) {
    return(cell_index_map |> dplyr::select(!!.cell_group, M))
  }
  
  draws =
    compute_alpha_normalised_draws(
      fit = .data |> attr("fit"),
      model_input = model_input,
      alpha_variable_subset = alpha_variable_subset
    ) |>
    
    # We want variability, not concentration
    mutate(.value = -.value) 
  
  # Reshape
  # Speed up if I have contrasts
  if(!contrasts |> is.null())
    draws = 
    draws |> 
    left_join(
      variability_covariates |> enframe(name = "C", value = "parameters_name"),
      by = "C"
    )  |> 
    select(-C) |> 
    pivot_wider(names_from = parameters_name, values_from = .value) |> 
    select( -.variable) 
  
  else
    draws =
    draws |>  
    pivot_wider(names_from = C, values_from = .value) %>%
    setNames(colnames(.)[1:5] |> c(variability_covariates)) |>
    select( -.variable) 
  
  # If I have constrasts calculate
  if (!is.null(contrasts)) 
    draws <- mutate_from_expr_list(draws, contrasts, ignore_errors = TRUE)
  
  # If no contrasts of interest just return an empty data frame
  if(ncol(draws) <= 4) return(cell_index_map |> dplyr::select(!!.cell_group, M) |> distinct())
  
  draws = draws |>
    pivot_longer(-c(1:4), names_to = "parameter", values_to = ".value") |>

    # Reorder because pivot long is bad
    mutate(parameter = parameter |> fct_relevel(colnames(draws)[-c(1:4)])) |>
    arrange(parameter)

  draws = draws |> 
    left_join(cell_index_map, by = "M") %>%
    select(!!.cell_group, everything())
  
  draws
  
}

#' Mutate Data Frame Based on Expression List
#'
#' @noRd
add_missing_contrast_names = function(formula_expr){
  contrast_names = names(formula_expr)
  if (is.null(contrast_names)) {
    contrast_names = formula_expr
  } else {
    missing_names = is.na(contrast_names) | contrast_names == ""
    contrast_names[missing_names] = formula_expr[missing_names]
  }

  make.unique(contrast_names)
}

#' @description
#' `mutate_from_expr_list` takes a data frame and a list of formula expressions, 
#' and mutates the data frame based on these expressions. It allows for ignoring 
#' errors during the mutation process.
#'
#' @param x A data frame to be mutated.
#' @param formula_expr A named list of formula expressions used for mutation.
#' @param ignore_errors Logical flag indicating whether to ignore errors during mutation.
#'
#' @return A mutated data frame with added or modified columns based on `formula_expr`.
#'
#' @details
#' The function performs various checks and transformations on the formula expressions,
#' ensuring that the specified transformations are valid and can be applied to the data frame.
#' It supports advanced features like handling special characters in column names and intelligent
#' parsing of formulas.
#'
#' @importFrom purrr map2_dfc
#' @importFrom tibble add_column
#' @importFrom tidyselect last_col
#' @importFrom dplyr mutate
#' @importFrom stringr str_subset
#' 
#' @noRd
#' 
mutate_from_expr_list = function(x, formula_expr, ignore_errors = TRUE){
  
  names(formula_expr) = add_missing_contrast_names(formula_expr)
  
  # Creating a named vector where the names are the strings to be replaced
  # and the values are empty strings
  contrasts_elements = contrasts_to_parameter_list(formula_expr, drop_back_quotes = FALSE)
  
  # Check if all elements of contrasts are in the parameter
  parameter_names = x |> colnames()
  
  # Check is backquoted are not used
  require_back_quotes = !contrasts_elements |>  str_remove_all("`") |> contains_only_valid_chars_for_column() 
  has_left_back_quotes = contrasts_elements |>  str_detect("^`") 
  has_right_back_quotes = contrasts_elements |>  str_detect("`$") 
  if_true_not_good = require_back_quotes & !(has_left_back_quotes & has_right_back_quotes)
  
  if(any(if_true_not_good))
    warning(sprintf("sccomp says: for columns which have special characters e.g. %s, you need to use surrounding backquotes ``.", paste(contrasts_elements[!if_true_not_good], sep=", ")))
  
  # Check if columns exist
  contrasts_not_in_the_model = 
    contrasts_elements |> 
    str_remove_all("`") |> 
    setdiff(parameter_names)
  
  contrasts_not_in_the_model = contrasts_not_in_the_model[contrasts_not_in_the_model!=""]
  
  if(length(contrasts_not_in_the_model) > 0 & !ignore_errors)
    warning(sprintf("sccomp says: These components of your contrasts are not present in the model as parameters: %s. Factors including special characters, e.g. \"(Intercept)\" require backquotes e.g. \"`(Intercept)`\" ", paste(contrasts_not_in_the_model, sep = ", ")))
  
  # Calculate
  if(ignore_errors) my_mutate = mutate_ignore_error
  else my_mutate = mutate
  
  map2_dfc(
    formula_expr,
    names(formula_expr),
    ~  x |>
      my_mutate(!!.y := eval(rlang::parse_expr(.x))) |>
      # mutate(!!column_name := eval(rlang::parse_expr(.x))) |>
      select(any_of(.y))
  )  |>
    
    # I could drop this to just result contrasts
    add_column(x |> select(M, .chain, .iteration, .draw), .before = 1)
  
}

#' @importFrom dplyr pull
#' @importFrom posterior as_draws_df summarise_draws rhat ess_bulk ess_tail
#' @noRd
draws_to_statistics = function(draws, false_positive_rate, test_composition_above_logit_fold_change, .cell_group, prefix = ""){
  
  # Define the variables as NULL to avoid CRAN NOTES
  M <- NULL
  .chain <- NULL
  .iteration <- NULL
  .draw <- NULL
  .value <- NULL
  parameter <- NULL
  bigger_zero <- NULL
  smaller_zero <- NULL
  lower <- NULL
  effect <- NULL
  upper <- NULL
  pH0 <- NULL
  FDR <- NULL
  n_eff <- NULL
  R_k_hat <- NULL
  
  .cell_group = enquo(.cell_group)
  lower_prob = false_positive_rate / 2
  upper_prob = 1 - lower_prob
  
  draw_summaries =
    draws |>
    group_by(!!.cell_group, M, parameter) |>
    group_modify(
      ~ {
        draws_df = 
          .x |>
          select(.chain, .iteration, .draw, value = .value) |> 
          as_draws_df() |>
          summarise_draws(
            lower = function(.x) stats::quantile(.x, probs = lower_prob, na.rm = TRUE),
            effect = function(.x) mean(.x, na.rm = TRUE),
            upper = function(.x) stats::quantile(.x, probs = upper_prob, na.rm = TRUE),
            rhat,
            ess_bulk,
            ess_tail
          ) |>
          select(-variable)

        # Standardise CI column names when posterior returns percentage names (e.g. `5%`, `97.5%`).
        quantile_cols = names(draws_df)[grepl("%$", names(draws_df))]
        names(draws_df)[names(draws_df) %in% quantile_cols] = c("lower", "upper")

        draws_df
      }
    ) |>
    ungroup()

  probability_stats =
    draws %>%
    group_by(!!.cell_group, M, parameter) %>%
    summarise(
      bigger_zero = sum(.value > test_composition_above_logit_fold_change, na.rm = TRUE),
      smaller_zero = sum(.value < -test_composition_above_logit_fold_change, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  draws =
    draw_summaries |>
    left_join(probability_stats, by = c(quo_name(.cell_group), "M", "parameter")) |>
    # Calculate probability non 0
    mutate(pH0 =  (1 - (pmax(bigger_zero, smaller_zero) / n))) |>
    with_groups(parameter, ~ mutate(.x, FDR = get_FDR(pH0))) |>
    
    select(!!.cell_group, M, parameter, lower, effect, upper, pH0, FDR, any_of(c("rhat", "ess_bulk", "ess_tail"))) |>
    suppressWarnings()
  
  # Setting up names separately because |> is not flexible enough
  draws |>
    setNames(c(colnames(draws)[1:3], sprintf("%s%s", prefix, colnames(draws)[4:ncol(draws)])))
}

mutate_ignore_error = function(x, ...){
  tryCatch(
    {  x |> mutate(...) },
    error=function(cond) {  x  }
  )
}

#' @importFrom dplyr cummean
#' @noRd
get_FDR = function(x){
  
  # Define the variables as NULL to avoid CRAN NOTES
  value <- NULL
  name <- NULL
  FDR <- NULL
  
  
  enframe(x) %>%
    arrange(value) %>%
    mutate(FDR = cummean(value)) %>%
    arrange(name) %>%
    pull(FDR)
}