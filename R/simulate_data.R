#' simulate_data
#'
#' @description This function simulates data from a fitted model.
#'
#' @import dplyr
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom rlang quo_is_null
#' @importFrom SingleCellExperiment colData
#' @importFrom parallel detectCores
#'
#' @param .data A tibble including a cell_group name column | sample name column | read counts column | factor columns | Pvalue column | a significance column
#' @param .estimate_object The result of sccomp_estimate execution. This is used for sampling from real-data properties.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_group identifier
#' @param .coefficients The column names for coefficients, for example, c(b_0, b_1)
#' @param variability_multiplier A real scalar. This can be used for artificially increasing the variability of the simulation for benchmarking purposes.
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()#' @param cores Integer, the number of cores to be used for parallel calculations.
#' @param cores Integer, the number of cores to be used for parallel calculations.
#' @param sig_figs Number of significant figures to use for Stan model output. Default is 9.
#' 
#' @return A tibble (`tbl`) with the following columns:
#' \itemize{
#'   \item \strong{sample} - A character column representing the sample name.
#'   \item \strong{type} - A factor column representing the type of the sample.
#'   \item \strong{phenotype} - A factor column representing the phenotype in the data.
#'   \item \strong{count} - An integer column representing the original cell counts.
#'   \item \strong{cell_group} - A character column representing the cell group identifier.
#'   \item \strong{b_0} - A numeric column representing the first coefficient used for simulation.
#'   \item \strong{b_1} - A numeric column representing the second coefficient used for simulation.
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
#' print("cmdstanr is needed to run this example.")
#' # Note: Before running the example, ensure that the 'cmdstanr' package is installed:
#' # install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))
#'
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists()) {
#'     data("counts_obj")
#'     library(dplyr)
#'
#'     estimate = sccomp_estimate(
#'       counts_obj,
#'       ~ type, ~1, "sample", "cell_group", "count",
#'       cores = 1
#'     )
#'
#'     # Set coefficients for cell_groups. In this case all coefficients are 0 for simplicity.
#'     counts_obj = counts_obj |> mutate(b_0 = 0, b_1 = 0)
#'
#'     # Simulate data
#'     simulate_data(counts_obj, estimate, ~type, ~1, sample, cell_group, c(b_0, b_1))
#'   }
#' }
#'
simulate_data <- function(.data,
                          .estimate_object,
                          formula_composition,
                          formula_variability = NULL,
                          .sample = NULL,
                          .cell_group = NULL,
                          .coefficients = NULL,
                          variability_multiplier = 5,
                          number_of_draws = 1,
                          mcmc_seed = sample_seed(),
                          cores = detectCores(),
                          sig_figs = 9) {
  
  # Run the function
  check_and_install_cmdstanr()
  
  UseMethod("simulate_data", .data)
}

#' @export
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom SingleCellExperiment counts
#' @importFrom purrr map_dbl
#' @importFrom purrr pmap
#' @importFrom readr read_file
#'
simulate_data.tbl = function(.data,
                             .estimate_object,
                             formula_composition,
                             formula_variability = NULL,
                             .sample = NULL,
                             .cell_group = NULL,
                             .coefficients = NULL,
                             variability_multiplier = 5,
                             number_of_draws = 1,
                             mcmc_seed = sample_seed(),
                             cores = detectCores(),
                             sig_figs = 9) {
  
  
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .coefficients = enquo(.coefficients)
  
  #Check column class
  check_if_columns_right_class(.data, !!.sample, !!.cell_group)
  
  model_data = attr(.estimate_object, "model_input")
  
  # # Select model based on noise model
  # if(attr(.estimate_object, "noise_model") == "multi_beta_binomial") my_model = stanmodels$glm_multi_beta_binomial_simulate_data
  # else if(attr(.estimate_object, "noise_model") == "dirichlet_multinomial") my_model = get_model_from_data("model_glm_dirichlet_multinomial_generate_quantities.rds", glm_dirichlet_multinomial_generate_quantities)
  # else if(attr(.estimate_object, "noise_model") == "logit_normal_multinomial") my_model = get_model_from_data("glm_multinomial_logit_linear_simulate_data.stan", read_file("~/PostDoc/sccomp/dev/stan_models/glm_multinomial_logit_linear_simulate_data.stan"))
  
  
  data_for_model =
    .data |>
    nest(data___ = -!!.sample) |>
    mutate(.exposure = sample(model_data$exposure, size = n(), replace = TRUE )) |>
    unnest(data___) |>
    data_simulation_to_model_input(
      formula_composition,
      #formula_variability,
      !!.sample, !!.cell_group, .exposure, !!.coefficients
    )
  names(data_for_model)  = names(data_for_model) |> stringr::str_c("_simulated")
  
  # Drop data from old input
  original_data = .estimate_object |> attr("model_input")
  original_data = original_data[(names(original_data) %in% c("C", "M", "A", "ncol_X_random_eff", "is_random_effect", "how_many_factors_in_random_design"  ))]
  
  # [1]  5.6260004 -0.6940178
  # prec_sd  = 0.816423129
  
  mod_rng = load_model("glm_multi_beta_binomial_simulate_data", threads = cores)
  
  fit = mod_rng |> sample_safe(
    generate_quantities_fx,
    attr(.estimate_object , "fit")$draws(format = "matrix"),
    
    # This is for the new data generation with selected factors to do adjustment
    data = data_for_model |> c(original_data) |> c(list(variability_multiplier = variability_multiplier)),
    seed = mcmc_seed,
    parallel_chains = attr(.estimate_object , "fit")$metadata()$threads_per_chain, 
    threads_per_chain = cores,
    sig_figs = sig_figs
  )
  
  parsed_fit =
    fit |>
    parse_generated_quantities(number_of_draws = number_of_draws) |>
    
    # Get sample name
    nest(data = -N) |>
    arrange(N) |>
    mutate(!!.sample := rownames(data_for_model$X_simulated)) |>
    unnest(data) |>
    
    # get cell type name
    nest(data = -M) |>
    mutate(!!.cell_group := colnames(data_for_model$beta_simulated)) |>
    unnest(data) |>
    
    select(-N, -M)
  
  .data |>
    left_join(
      parsed_fit,
      by = c(quo_name(.sample), quo_name(.cell_group))
    )
  
  
}

#' @importFrom purrr when
#' @importFrom stats model.matrix
#'
#' @keywords internal
#' @noRd
#'
data_simulation_to_model_input =
  function(.data, formula, .sample, .cell_type, .exposure, .coefficients, truncation_ajustment = 1, approximate_posterior_inference ){
    
    # Define the variables as NULL to avoid CRAN NOTES
    sd <- NULL
    . <- NULL
    
    
    # Prepare column same enquo
    .sample = enquo(.sample)
    .cell_type = enquo(.cell_type)
    .exposure = enquo(.exposure)
    .coefficients = enquo(.coefficients)
    
    factor_names = parse_formula(formula)
    
    sample_data =
      .data %>%
      select(!!.sample, any_of(factor_names)) %>%
      distinct() %>%
      arrange(!!.sample)
    X =
      sample_data %>%
      model.matrix(formula, data=.) %>%
      apply(2, function(x) {
        
        if(sd(x)==0 ) x
        else x |> scale(scale=FALSE)
        
      } ) %>%
      {
        .x = (.)
        rownames(.x) = sample_data %>% pull(!!.sample)
        .x
      }
    
    if(factor_names == "1") XA = X[,1, drop=FALSE]
    else XA = X[,c(1,2), drop=FALSE]
    
    XA = XA |> 
      as_tibble()  |> 
      distinct()
    
    cell_cluster_names =
      .data %>%
      distinct(!!.cell_type) %>%
      arrange(!!.cell_type) %>%
      pull(!!.cell_type)
    
    coefficients =
      .data %>%
      select(!!.cell_type, !!.coefficients) %>%
      unnest(!!.coefficients) %>%
      distinct() %>%
      arrange(!!.cell_type) %>%
      as_matrix(rownames = quo_name(.cell_type)) %>%
      t()
    
    list(
      N = .data %>% distinct(!!.sample) %>% nrow(),
      M = .data %>% distinct(!!.cell_type) %>% nrow(),
      exposure = .data %>% distinct(!!.sample, !!.exposure) %>% arrange(!!.sample) %>% pull(!!.exposure),
      X = X,
      XA = XA,
      C = ncol(X),
      A =  ncol(XA),
      beta = coefficients
    )
    
  }
