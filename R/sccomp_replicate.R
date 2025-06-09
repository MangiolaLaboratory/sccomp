


#' sccomp_replicate
#'
#' @description This function replicates counts from a real-world dataset.
#'
#'
#' @param fit The result of sccomp_estimate.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment. In most cases, if differentially variability is of interest, the formula should only include the factor of interest as a large anount of data is needed to define variability depending to each factors. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#'
#' @return A tibble `tbl` with cell_group-wise statistics
#'
#' @return A tibble (`tbl`), with the following columns:
#' \itemize{
#'   \item \strong{cell_group} - A character column representing the cell group being tested.
#'   \item \strong{sample} - A factor column representing the sample name from which data was generated.
#'   \item \strong{generated_proportions} - A numeric column representing the proportions generated from the model.
#'   \item \strong{generated_counts} - An integer column representing the counts generated from the model.
#'   \item \strong{replicate} - An integer column representing the replicate number, where each row corresponds to a different replicate of the data.
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
#'       ~ type, ~1, sample, cell_group, count,
#'       cores = 1
#'     ) |>
#'     sccomp_replicate()
#'   }
#' }
sccomp_replicate <- function(fit,
                             formula_composition = NULL,
                             formula_variability = NULL,
                             number_of_draws = 1,
                             mcmc_seed = sample(1e5, 1)) {
  
  # Run the function
  check_and_install_cmdstanr()
  
  UseMethod("sccomp_replicate", fit)
}

#' @export
#'
sccomp_replicate.sccomp_tbl = function(fit,
                                       formula_composition = NULL,
                                       formula_variability = NULL,
                                       number_of_draws = 1,
                                       mcmc_seed = sample(1e5, 1)){
  
  .sample = attr(fit, ".sample")
  .cell_group = attr(fit, ".cell_group")
  
  sample_names =
    fit |>
    select(count_data) |>
    unnest(count_data) |>
    distinct(!!.sample) |> 
    pull(!!.sample)
  
  rng =
    replicate_data(
      fit,
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      number_of_draws = number_of_draws,
      mcmc_seed = mcmc_seed
    )
  
  model_input = attr(fit, "model_input")
  
  # mean generated
  rng |>
    
    # Parse
    parse_generated_quantities(number_of_draws = number_of_draws) |>
    
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
  
}