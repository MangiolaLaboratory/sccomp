#' sccomp_predict
#'
#' @description This function replicates counts from a real-world dataset.
#'
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
#'
#' @param fit The result of sccomp_estimate.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#' @param new_data A sample-wise data frame including the column that represent the factors in your formula. If you want to predict proportions for 10 samples, there should be 10 rows. T
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#' @param summary_instead_of_draws Return the summary values (i.e. mean and quantiles) of the predicted proportions, or return single draws. Single draws can be helful to better analyse the uncertainty of the prediction.
#' @param robust A logical. If TRUE, use robust statistics (median and median absolute deviation) instead of classical statistics (mean and standard deviation) for the summary calculations.
#'
#' @return A tibble (`tbl`) with the following columns:
#' \itemize{
#'   \item \strong{cell_group} - A character column representing the cell group being tested.
#'   \item \strong{sample} - A factor column representing the sample name for which the predictions are made.
#'   \item \strong{proportion_mean} - A numeric column representing the predicted mean (or median when robust=TRUE) proportions from the model.
#'   \item \strong{proportion_lower} - A numeric column representing the lower bound (2.5%) of the 95% credible interval for the predicted proportions.
#'   \item \strong{proportion_upper} - A numeric column representing the upper bound (97.5%) of the 95% credible interval for the predicted proportions.
#' }
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
#'   if (instantiate::stan_cmdstan_exists() && .Platform$OS.type == "unix") {
#'     data("counts_obj")
#'
#'     sccomp_estimate(
#'       counts_obj,
#'       ~ type, ~1, "sample", "cell_group", "count",
#'       cores = 1
#'     ) |>
#'     sccomp_predict()
#'   }
#' }
#' 
#'
sccomp_predict <- function(fit,
                           formula_composition = NULL,
                           new_data = NULL,
                           number_of_draws = 500,
                           mcmc_seed = sample_seed(),
                           summary_instead_of_draws = TRUE,
                           robust = FALSE) {
  
  
  # Run the function
  check_and_install_cmdstanr()
  
  
  UseMethod("sccomp_predict", fit)
}

#' @export
#'
sccomp_predict.sccomp_tbl = function(fit,
                                     formula_composition = NULL,
                                     new_data = NULL,
                                     number_of_draws = 500,
                                     mcmc_seed = sample_seed(),
                                     summary_instead_of_draws = TRUE,
                                     robust = FALSE){
  
  
  
  model_input = attr(fit, "model_input")
  .sample = attr(fit, ".sample")
  .cell_group = attr(fit, ".cell_group")
  
  rng =
    replicate_data(
      fit,
      formula_composition = formula_composition,
      formula_variability = ~ 1,
      new_data = new_data,
      number_of_draws = number_of_draws,
      mcmc_seed = mcmc_seed
    )
  
  # New data
  if(new_data |> is.null())
    new_data =
    fit |>
    attr("count_data") |>
    distinct() |> 
    .subset(!!.sample)
  
  # If seurat
  else if(new_data |> is("Seurat")) 
    new_data = new_data[[]] 
  
  sample_names =
    new_data |>
    distinct(!!.sample) |> 
    pull(!!.sample)
  
  # mean generated
  if(summary_instead_of_draws){
    
    prediction_df = 
      rng |>
      summary_to_tibble("mu", "M", "N", robust = robust) |>
      select(M, N, 
             
             # This is not great, consider making proportion_mean more general such as proportion_point_estimate
             proportion_mean = any_of(c("mean", "median")), 
             proportion_lower = `2.5%`, 
             proportion_upper = `97.5%`) 
  }

  else
    prediction_df = 
    rng |>
    draws_to_tibble_x_y("mu", "M", "N") |> 
    select(M, N, proportion = .value, .draw) 
  
  prediction_df = 
    prediction_df |>
    
    # Get sample name
    nest(data = -N) |>
    arrange(N) |>
    mutate(!!.sample := sample_names) |>
    unnest(data) |>
    
    # get cell type name
    nest(data = -M) |>
    mutate(!!.cell_group := colnames(model_input$y)) |>
    unnest(data) |>
    
    select(-N, -M) |>
    select(!!.cell_group, !!.sample, everything())
  
  new_data |> 
    left_join(prediction_df, by = join_by(!!.sample))
}