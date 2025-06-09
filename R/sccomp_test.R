

#' sccomp_test
#'
#' @description This function test contrasts from a sccomp result.
#'
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
#'   \item count_data - Nested input count data.
#' }#'
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
  
  # Abundance
  abundance_CI =
    get_abundance_contrast_draws(.data, contrasts)
  
  # If my contrasts do not match my model. I have to do something more elegant.
  if ("parameter" %in% colnames(abundance_CI))
    abundance_CI =
    abundance_CI |>
    draws_to_statistics(
      percent_false_positive/100,
      test_composition_above_logit_fold_change,
      !!.cell_group,
      "c_"
    )
  
  variability_CI =
    get_variability_contrast_draws(.data, contrasts)
  
  # Variability
  if ("parameter" %in% colnames(variability_CI))
    variability_CI = 
    variability_CI |>
    draws_to_statistics(
      percent_false_positive / 100,
      test_composition_above_logit_fold_change,
      !!.cell_group,
      "v_"
    )
  
  # If I don't have factors (~1)
  if (!"factor" %in% colnames(model_input$factor_parameter_dictionary))
    factor_parameter_dictionary = tibble(`factor` = character(), design_matrix_col = character())
  else
    factor_parameter_dictionary =
    model_input$factor_parameter_dictionary |>
    select(`factor`, design_matrix_col)
  
  # Merge and parse
  result =
    abundance_CI |>
    
    # Add ALPHA
    left_join(variability_CI) |>
    suppressMessages() |>
    
    # Add easy to understand factor labels
    left_join(factor_parameter_dictionary,
              by = c("parameter" = "design_matrix_col")) |>
    select(parameter, `factor`, everything()) |>
    
    select(!!.cell_group, everything(),-M)
  
  
  result =
    result |>
    left_join(
      truncation_df2 |>
        select(-any_of(c(
          "M",
          "N",
          ".variable",
          "mean",
          "se_mean",
          "sd",
          "n_eff",
          "R_hat",
          "k_hat",
          "Rhat",
          ".lower", ".median", ".upper"
        ))) |>
        nest(count_data = -!!.cell_group),
      by = quo_name(.cell_group)
    )
  
  # result =
  #   result |>
  # 
  #   # Add back attributes
  #   add_attr(
  #     .data |> attr("fit") |> get_mean_precision_association(),
  #     "mean_concentration_association"
  #   )
  
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
    
    # Add class to the tbl
    add_class("sccomp_tbl") 
}

# this can be helpful if we want to draw PCA with uncertainty
get_abundance_contrast_draws = function(.data, contrasts){
  
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
  
  # Beta
  beta_factor_of_interest = .data |> attr("model_input") %$% X |> colnames()
  # beta =
  #   .data |>
  #   attr("fit") %>%
  #   draws_to_tibble_x_y("beta", "C", "M") |>
  #   pivot_wider(names_from = C, values_from = .value) %>%
  #   setNames(colnames(.)[1:5] |> c(beta_factor_of_interest))
  
  if(contrasts |> is.null())
    draws = 
    .data |>
    attr("fit") %>%
    draws_to_tibble_x_y("beta", "C", "M") |> 
    pivot_wider(names_from = C, values_from = .value) %>%
    setNames(colnames(.)[1:5] |> c(beta_factor_of_interest))
  
  else if(
    (beta_factor_of_interest %in% contrasts_to_parameter_list(contrasts)) |> which() |> length() > 0
  )
    
    draws =
    .data |>
    attr("fit") %>%
    draws_to_tibble_x_y("beta", "C", "M") |> 
    left_join(
      beta_factor_of_interest |> enframe(name = "C", value = "parameters_name"),
      by = "C"
    )  |> 
    filter(parameters_name %in% contrasts_to_parameter_list(contrasts)) |> 
    select(-C) |> 
    pivot_wider(names_from = parameters_name, values_from = .value)
  
  else 
    draws = tibble()
  
  
  # Abundance
  draws = draws |> select(-.variable)
  
  
  # Random effect
  
  beta_random_effect_factor_of_interest = .data |> attr("model_input") %$% X_random_effect |> colnames()
  
  if(
    .data |> attr("model_input") %$% n_random_eff > 0 &&
    (
      contrasts |> is.null() || 
      (beta_random_effect_factor_of_interest %in% contrasts_to_parameter_list(contrasts)) |> which() |> length() > 0
    )  
  ){
    
    
    beta_random_effect =
      .data |>
      attr("fit") %>%
      draws_to_tibble_x_y("random_effect", "C", "M") 
    
    # Add last component
    other_group_random_effect = 
      beta_random_effect |> 
      with_groups(c(C, .chain, .iteration, .draw, .variable ), ~ .x |> summarise(.value = sum(.value))) |> 
      mutate(.value = -.value, M = beta_random_effect |> pull(M) |> max() + 1)
    
    # I HAVE TO REGULARISE THE LAST COMPONENT
    mean_of_the_sd_of_the_point_estimates = 
      beta_random_effect |> 
      group_by(M, C) |> 
      summarise(point_estimate = mean(.value)) |> 
      group_by(M) |> 
      summarise(sd_of_point_estimates = sd(point_estimate)) |> 
      pull(sd_of_point_estimates) |> 
      mean()
    
    other_sd_of_the_point_estimates = 
      other_group_random_effect |> 
      group_by(M, C) |> 
      summarise(point_estimate = mean(.value)) |> 
      group_by(M) |> 
      summarise(sd_of_point_estimates = sd(point_estimate)) |> 
      pull(sd_of_point_estimates)
    
    other_group_random_effect = 
      other_group_random_effect |> 
      mutate(.value = .value / (other_sd_of_the_point_estimates / mean_of_the_sd_of_the_point_estimates))
    
    
    beta_random_effect = 
      beta_random_effect |> 
      bind_rows( other_group_random_effect )
    
    # mutate(is_treg = cell_type =="treg") |>
    #   nest(data = -is_treg) |>
    #   mutate(data = map2(
    #     data, is_treg,
    #     ~ {
    #       if(.y) .x |> mutate(c_effect = c_effect/5 )
    #       else(.x)
    #     }
    #   )) |>
    #   unnest(data) |>
    
    
    # Reshape
    # Speed up if I have contrasts
    if(!contrasts |> is.null())
      beta_random_effect = 
      beta_random_effect |> 
      left_join(
        beta_random_effect_factor_of_interest |> enframe(name = "C", value = "parameters_name"),
        by = "C"
      )  |> 
      filter(parameters_name %in% contrasts_to_parameter_list(contrasts)) |> 
      select(-C) |> 
      pivot_wider(names_from = parameters_name, values_from = .value)
    
    else
      beta_random_effect = 
      beta_random_effect |>
      pivot_wider(names_from = C, values_from = .value) %>%
      setNames(colnames(.)[1:5] |> c(beta_random_effect_factor_of_interest))
    
    # If I don't have fix nor 1st level random effect
    if(draws |> nrow() == 0)
      draws = select(beta_random_effect, -.variable)
    else 
      draws = draws |> 
      left_join(select(beta_random_effect, -.variable),
                by = c("M", ".chain", ".iteration", ".draw")
      )
    
  }  else {
    beta_random_effect_factor_of_interest = ""
  }
  
  # Second random effect. IN THE FUTURE THIS WILL BE VECTORISED TO ARBUTRARY GRI+OUING
  beta_random_effect_factor_of_interest_2 = .data |> attr("model_input") %$% X_random_effect_2 |> colnames()
  
  if(
    .data |> attr("model_input") %$% n_random_eff > 1 &&
    (
      contrasts |> is.null() || 
      (beta_random_effect_factor_of_interest_2 %in% contrasts_to_parameter_list(contrasts)) |> which() |> length() > 0
    )
  ){
    
    beta_random_effect_2 =
      .data |>
      attr("fit") %>%
      draws_to_tibble_x_y("random_effect_2", "C", "M") 
    
    # Add last component
    other_group_random_effect = 
      beta_random_effect_2 |> 
      with_groups(c(C, .chain, .iteration, .draw, .variable ), ~ .x |> summarise(.value = sum(.value))) |> 
      mutate(.value = -.value, M = beta_random_effect_2 |> pull(M) |> max() + 1)
    
    # I HAVE TO REGULARISE THE LAST COMPONENT
    mean_of_the_sd_of_the_point_estimates = 
      beta_random_effect_2 |> 
      group_by(M, C) |> 
      summarise(point_estimate = mean(.value)) |> 
      group_by(M) |> 
      summarise(sd_of_point_estimates = sd(point_estimate)) |> 
      pull(sd_of_point_estimates) |> 
      mean()
    
    other_sd_of_the_point_estimates = 
      other_group_random_effect |> 
      group_by(M, C) |> 
      summarise(point_estimate = mean(.value)) |> 
      group_by(M) |> 
      summarise(sd_of_point_estimates = sd(point_estimate)) |> 
      pull(sd_of_point_estimates)
    
    other_group_random_effect = 
      other_group_random_effect |> 
      mutate(.value = .value / (other_sd_of_the_point_estimates / mean_of_the_sd_of_the_point_estimates))
    
    
    beta_random_effect_2 = 
      beta_random_effect_2 |> 
      bind_rows( other_group_random_effect )
    
    # Reshape
    # Speed up if I have contrasts
    if(!contrasts |> is.null())
      beta_random_effect_2 = 
      beta_random_effect_2 |> 
      left_join(
        beta_random_effect_factor_of_interest_2 |> enframe(name = "C", value = "parameters_name"),
        by = "C"
      )  |> 
      filter(parameters_name %in% contrasts_to_parameter_list(contrasts)) |> 
      select(-C) |> 
      pivot_wider(names_from = parameters_name, values_from = .value)
    
    else
      beta_random_effect_2 = 
      beta_random_effect_2 |>
      pivot_wider(names_from = C, values_from = .value) %>%
      setNames(colnames(.)[1:5] |> c(beta_random_effect_factor_of_interest_2))
    
    # If I don't have fix nor 1st level random effect
    if(draws |> nrow() == 0)
      draws = select(beta_random_effect_2, -.variable)
    else 
      draws = draws |> 
      left_join(select(beta_random_effect_2, -.variable),
                by = c("M", ".chain", ".iteration", ".draw")
      )
  } else {
    beta_random_effect_factor_of_interest_2 = ""
  }
  
  
  # If I have constrasts calculate
  if(!is.null(contrasts))
    draws = 
    draws |> 
    mutate_from_expr_list(contrasts, ignore_errors = FALSE) |>
    select(- any_of(c(beta_factor_of_interest, beta_random_effect_factor_of_interest) |> setdiff(contrasts)) ) 
  
  # Add cell name
  draws = draws |> 
    left_join(
      .data |>
        attr("model_input") %$%
        y %>%
        colnames() |>
        enframe(name = "M", value  = quo_name(.cell_group)),
      by = "M"
    ) %>%
    select(!!.cell_group, everything())
  
  
  # If no contrasts of interest just return an empty data frame
  if(ncol(draws)==5) return(draws |> distinct(M, !!.cell_group))
  
  # Get convergence
  convergence_df =
    .data |>
    attr("fit") |>
    summary_to_tibble("beta", "C", "M") |>
    
    # Add cell name
    left_join(
      .data |>
        attr("model_input") %$%
        y %>%
        colnames() |>
        enframe(name = "M", value  = quo_name(.cell_group)),
      by = "M"
    ) |>
    
    # factor names
    left_join(
      beta_factor_of_interest |>
        enframe(name = "C", value = "parameter"),
      by = "C"
    )
  
  # if ("Rhat" %in% colnames(convergence_df)) {
  #   convergence_df <- rename(convergence_df, R_k_hat = Rhat)
  # } else if ("khat" %in% colnames(convergence_df)) {
  #   convergence_df <- rename(convergence_df, R_k_hat = khat)
  # }
  
  
  convergence_df =
    convergence_df |> 
    select(!!.cell_group, parameter, any_of(c("n_eff", "R_k_hat", "rhat", "ess_bulk", "ess_tail"))) |>
    suppressWarnings()
  
  draws |>
    pivot_longer(-c(1:5), names_to = "parameter", values_to = ".value") |>
    
    # Attach convergence if I have no contrasts
    left_join(convergence_df, by = c(quo_name(.cell_group), "parameter")) |>
    
    # Reorder because pivot long is bad
    mutate(parameter = parameter |> fct_relevel(colnames(draws)[-c(1:5)])) |>
    arrange(parameter)
  
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
  
  variability_factor_of_interest = .data |> attr("model_input") %$% XA |> colnames()
  
  draws =
    
    .data |>
    attr("fit") %>%
    draws_to_tibble_x_y("alpha_normalised", "C", "M") |>
    
    # We want variability, not concentration
    mutate(.value = -.value) 
  
  # Reshape
  # Speed up if I have contrasts
  if(!contrasts |> is.null())
    draws = 
    draws |> 
    left_join(
      variability_factor_of_interest |> enframe(name = "C", value = "parameters_name"),
      by = "C"
    )  |> 
    filter(parameters_name %in% contrasts_to_parameter_list(contrasts)) |> 
    select(-C) |> 
    pivot_wider(names_from = parameters_name, values_from = .value) |> 
    select( -.variable) 
  
  else
    draws =
    draws |>  
    pivot_wider(names_from = C, values_from = .value) %>%
    setNames(colnames(.)[1:5] |> c(variability_factor_of_interest)) |>
    select( -.variable) 
  
  # If I have constrasts calculate
  if (!is.null(contrasts)) 
    draws <- mutate_from_expr_list(draws, contrasts, ignore_errors = TRUE)
  
  draws =  draws |>
    
    # Add cell name
    left_join(
      .data |> attr("model_input") %$%
        y %>%
        colnames() |>
        enframe(name = "M", value  = quo_name(.cell_group)),
      by = "M"
    ) %>%
    select(!!.cell_group, everything())
  
  # If no contrasts of interest just return an empty data frame
  if(ncol(draws)==5) return(draws |> distinct(M, !!.cell_group))
  
  # Get convergence
  convergence_df =
    .data |>
    attr("fit") |>
    summary_to_tibble("alpha_normalised", "C", "M") |>
    
    # Add cell name
    left_join(
      .data |>
        attr("model_input") %$%
        y %>%
        colnames() |>
        enframe(name = "M", value  = quo_name(.cell_group)),
      by = "M"
    ) |>
    
    # factor names
    left_join(
      variability_factor_of_interest |>
        enframe(name = "C", value = "parameter"),
      by = "C"
    )
  
  convergence_df =
    convergence_df |> 
    select(!!.cell_group, parameter, any_of(c("n_eff", "R_k_hat", "rhat", "ess_bulk", "ess_tail"))) |>
    suppressWarnings()
  
  
  draws |>
    pivot_longer(-c(1:5), names_to = "parameter", values_to = ".value") |>
    
    # Attach convergence if I have no contrasts
    left_join(convergence_df, by = c(quo_name(.cell_group), "parameter")) |>
    
    # Reorder because pivot long is bad
    mutate(parameter = parameter |> fct_relevel(colnames(draws)[-c(1:5)])) |>
    arrange(parameter)
  
}

#' Mutate Data Frame Based on Expression List
#'
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
  
  if(formula_expr |> names() |> is.null())
    names(formula_expr) = formula_expr
  
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

draws_to_statistics = function(draws, false_positive_rate, test_composition_above_logit_fold_change, .cell_group, prefix = ""){
  
  # Define the variables as NULL to avoid CRAN NOTES
  M <- NULL
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
  
  draws =
    draws %>%
    group_by(!!.cell_group, M, parameter, rhat, ess_bulk, ess_tail) %>%
    summarise(
      lower = quantile(.value, false_positive_rate / 2),
      effect = quantile(.value, 0.5),
      upper = quantile(.value, 1 - (false_positive_rate / 2)),
      bigger_zero = sum(.value > test_composition_above_logit_fold_change),
      smaller_zero = sum(.value < -test_composition_above_logit_fold_change),
      # R_k_hat = unique(R_k_hat),
      # n_eff = unique(n_eff),
      n = n(),
      .groups = "drop"  # To ungroup the output if needed
    ) |> 
    
    # Calculate probability non 0
    mutate(pH0 =  (1 - (pmax(bigger_zero, smaller_zero) / n))) |>
    with_groups(parameter, ~ mutate(.x, FDR = get_FDR(pH0))) |>
    
    select(!!.cell_group, M, parameter, lower, effect, upper, pH0, FDR, any_of(c("n_eff", "R_k_hat", "rhat", "ess_bulk", "ess_tail"))) |>
    suppressWarnings()
  
  # Setting up names separately because |> is not flexible enough
  draws |>
    setNames(c(colnames(draws)[1:3], sprintf("%s%s", prefix, colnames(draws)[4:ncol(draws)])))
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