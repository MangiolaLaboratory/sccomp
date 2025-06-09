

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
#'       ~ 0 + type, ~1, sample, cell_group, count,
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