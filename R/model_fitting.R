#' @importFrom readr write_file
fit_model = function(
    data_for_model, model_name, censoring_iteration = 1, cores = detectCores(), quantile = 0.95,
    warmup_samples = 300, approximate_posterior_inference = NULL, inference_method, verbose = TRUE,
    seed , pars = c("beta", "alpha", "prec_coeff","prec_sd"), output_samples = NULL, chains=NULL, max_sampling_iterations = 20000, 
    output_directory = "sccomp_draws_files",
    sig_figs = 9,
    ...
)
{
  
  # # if analysis approximated
  # # If posterior analysis is approximated I just need enough
  # how_many_posterior_draws_practical = ifelse(approximate_posterior_analysis, 1000, how_many_posterior_draws)
  # additional_parameters_to_save = additional_parameters_to_save %>% c("lambda_log_param", "sigma_raw") %>% unique
  
  
  # Find number of draws
  draws_supporting_quantile = 50
  if(is.null(output_samples)){
    
    output_samples =
      (draws_supporting_quantile/((1-quantile)/2)) %>% # /2 because I have two tails
      max(4000) 
    
    if(output_samples > max_sampling_iterations) {
      # message("sccomp says: the number of draws used to defined quantiles of the posterior distribution is capped to 20K.") # This means that for very low probability threshold the quantile could become unreliable. We suggest to limit the probability threshold between 0.1 and 0.01")
      output_samples = max_sampling_iterations
      
    }}
  
  # Find optimal number of chains
  if(is.null(chains))
    chains =
      find_optimal_number_of_chains(
        how_many_posterior_draws = output_samples,
        warmup = warmup_samples,
        parallelisation_start_penalty = 100
      ) %>%
      min(cores)
  
  # chains = 3
  
  init_list=list(
    prec_coeff = c(5,0),
    prec_sd = 1,
    alpha = matrix(c(rep(5, data_for_model$M), rep(0, (data_for_model$A-1) *data_for_model$M)), nrow = data_for_model$A, byrow = TRUE),
    beta_raw = matrix(0, data_for_model$C , data_for_model$M) ,
    mix_p = 0.1 
  )
  
  if(data_for_model$n_random_eff>0){
    init_list$random_effect_raw = matrix(0, data_for_model$ncol_X_random_eff[1]  , data_for_model$M)
    init_list$random_effect_sigma_raw = matrix(0, data_for_model$M , data_for_model$how_many_factors_in_random_design[1])
    init_list$sigma_correlation_factor = array(0, dim = c(
      data_for_model$M, 
      data_for_model$how_many_factors_in_random_design[1], 
      data_for_model$how_many_factors_in_random_design[1]
    ))
    
    # init_list$random_effect_sigma_mu = 0.5 |> as.array()
    # init_list$random_effect_sigma_sigma = 0.2 |> as.array()
    init_list$zero_random_effect = rep(0, size = 1) |> as.array()
    
  } 
  
  if(data_for_model$n_random_eff>1){
    init_list$random_effect_raw_2 = matrix(0, data_for_model$ncol_X_random_eff[2]  , data_for_model$M)
    init_list$random_effect_sigma_raw_2 = matrix(0, data_for_model$M , data_for_model$how_many_factors_in_random_design[2])
    init_list$sigma_correlation_factor_2 = array(0, dim = c(
      data_for_model$M, 
      data_for_model$how_many_factors_in_random_design[2], 
      data_for_model$how_many_factors_in_random_design[2]
    ))
    
  } 
  
  init = map(1:chains, ~ init_list) %>%
    setNames(as.character(1:chains))
  
  #output_directory = "sccomp_draws_files"
  dir.create(output_directory, showWarnings = FALSE)
  
  # Fit
  mod = load_model(model_name, threads = cores)
  
  # Avoid 0 proportions
  if(data_for_model$is_proportion && min(data_for_model$y_proportion)==0){
    warning("sccomp says: your proportion values include 0. Assuming that 0s derive from a precision threshold (e.g. deconvolution), 0s are converted to the smaller non 0 proportion value.")
    data_for_model$y_proportion[data_for_model$y_proportion==0] =
      min(data_for_model$y_proportion[data_for_model$y_proportion>0])
  }
  
  if(inference_method == "hmc"){
    
 # tryCatch({
    mod$sample(
        data = data_for_model ,
        chains = chains,
        parallel_chains = chains,
        threads_per_chain = ceiling(cores/chains),
        iter_warmup = warmup_samples,
        iter_sampling = as.integer(output_samples /chains),
        #refresh = ifelse(verbose, 1000, 0),
        seed = seed,
        save_warmup = FALSE,
        init = init,
        output_dir = output_directory,
        show_messages = verbose,
        sig_figs = sig_figs,
        show_exceptions = verbose,
        ...
      ) 
      
    # },
    # error = function(e) {
    # 
    #   # I don't know why thi is needed nd why the model sometimes is not compliled correctly
    #   if(e |> as.character() |>  str_detect("Model not compiled"))
    #     model = load_model(model_name, force=TRUE, threads = cores)
    #   else
    #     stop(e)
    # 
    # })

    
  } else{
    
    if(inference_method=="pathfinder") init = pf
    else if(inference_method=="variational") init = list(init_list)
    
    vb_iterative(
      mod,
      output_samples = output_samples ,
      iter = 10000,
      tol_rel_obj = 0.01,
      data = data_for_model, refresh = ifelse(verbose, 1000, 0),
      seed = seed,
      output_dir = output_directory,
      init = init,
      inference_method = inference_method, 
      cores = cores,
      psis_resample = FALSE, 
      verbose = verbose,
      sig_figs = sig_figs,
      show_exceptions = FALSE,
      ...
    ) 
    
  }
  
  
  
}

get_model_from_data = function(file_compiled_model, model_code){
  if(file.exists(file_compiled_model))
    readRDS(file_compiled_model)
  else {
    model_generate = stan_model(model_code = model_code)
    model_generate  %>% saveRDS(file_compiled_model)
    model_generate
    
  }
}

#' Load, Compile, and Cache a Stan Model
#'
#' This function attempts to load a precompiled Stan model using the `instantiate` package.
#' If the model is not found in the cache or force recompilation is requested, it will locate
#' the Stan model file within the `sccomp` package, compile it using `cmdstanr`, and save the 
#' compiled model to the cache directory for future use.
#'
#' @param name A character string representing the name of the Stan model (without the `.stan` extension).
#' @param cache_dir A character string representing the path to the cache directory where compiled models are saved. 
#' Defaults to `sccomp_stan_models_cache_dir`.
#' @param force A logical value. If `TRUE`, the model will be recompiled even if it exists in the cache. 
#' Defaults to `FALSE`.
#' @param threads An integer specifying the number of threads to use for compilation. 
#' Defaults to `1`.
# 
#' @return A compiled Stan model object from `cmdstanr`.
#' 
#' @importFrom instantiate stan_package_model
#' @importFrom instantiate stan_package_compile
#' 
#' @noRd
#' 
#' @examples
#' \donttest{
#'   model <- load_model("glm_multi_beta_binomial_", "~/cache", force = FALSE, threads = 1)
#' }
load_model <- function(name, cache_dir = sccomp_stan_models_cache_dir, force=FALSE, threads = 1) {
  
  
  # tryCatch({
  #   # Attempt to load a precompiled Stan model using the instantiate package
  #   instantiate::stan_package_model(
  #     name = name,
  #     package = "sccomp"
  #   )
  # }, error = function(e) {
  # Try to load the model from cache
  
  # RDS compiled model
  cache_dir |> dir.create(showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(cache_dir, paste0(name, ".rds"))
  
  # .STAN raw model
  stan_model_path <- system.file("stan", paste0(name, ".stan"), package = "sccomp")
  
  if (file.exists(cache_file) & !force) {
    message("Loading model from cache...")
    return(readRDS(cache_file))
  }
  
  # If loading the precompiled model fails, find the Stan model file within the package
  message("Precompiled model not found. Compiling the model...")
  
  # Compile the Stan model using cmdstanr with threading support enabled
  instantiate::stan_package_compile(
    stan_model_path, 
    cpp_options = list(stan_threads = TRUE),
    force_recompile = TRUE, 
    threads = threads, 
    dir = system.file("stan", package = "sccomp")
  )
  mod = instantiate::stan_package_model(
    name = name, 
    package = "sccomp", 
    compile = TRUE,
    cpp_options = list(stan_threads = TRUE)
  ) |> suppressWarnings()
  
  # Save the compiled model object to cache
  saveRDS(mod, file = cache_file)
  message("Model compiled and saved to cache successfully.")
  
  return(mod)
  # })
  
}

#' Check and Install cmdstanr and CmdStan
#'
#' This function checks if the `cmdstanr` package (version 0.9.0 or higher) and CmdStan are installed. 
#' If they are not installed, it installs them automatically in non-interactive sessions
#' or asks for permission to install them in interactive sessions.
#' 
#' The function requires cmdstanr version 0.9.0 or higher for support of the new `sum_to_zero_vector` type.
#'
#' @importFrom instantiate stan_cmdstan_exists
#' @importFrom rlang check_installed
#' @importFrom rlang abort
#' @importFrom rlang check_installed
#' @return NULL
#' 
#' @noRd
check_and_install_cmdstanr <- function() {
  
  # Check if cmdstanr is installed
  # from https://github.com/wlandau/instantiate/blob/33989d74c26f349e292e5efc11c267b3a1b71d3f/R/utils_assert.R#L114
  
  # tryCatch(
    rlang::check_installed(
      pkg = "cmdstanr",
      reason = paste(
        "The {cmdstanr} package is required in order to install",
        "CmdStan and run Stan models. Please install it manually using",
        "install.packages(pkgs = \"cmdstanr\",",
        "repos = c(\"https://mc-stan.org/r-packages/\", getOption(\"repos\"))"
      ),

      # I have to see if Bioconductor is compatible with this
      action = function(...) install.packages(..., repos = c('https://stan-dev.r-universe.dev', 'https://cloud.r-project.org'))
    )
  #   ,
  #   error = function(e) {
  #     clear_stan_model_cache()
  #     stan_error(conditionMessage(e))
  #   }
  # )
  
  # Check cmdstanr version - requires 0.9.0 or higher for sum_to_zero_vector support
  if (packageVersion("cmdstanr") < "0.9.0") {
    stop(
      "cmdstanr version 0.9.0 or higher is required for sum_to_zero_vector support.\n\n",
      "Current version: ", packageVersion("cmdstanr"), "\n\n",
      "Please update cmdstanr using:\n",
      "install.packages(pkgs = \"cmdstanr\", repos = c(\"https://mc-stan.org/r-packages/\", getOption(\"repos\")))\n\n",
      "Or install the latest development version:\n",
      "remotes::install_github(\"stan-dev/cmdstanr\")"
    )
  }
  
  # Check cmdstanr version - requires 0.9.0 or higher for sum_to_zero_vector support
  if (packageVersion("cmdstanr") < "0.9.0") {
    stop(
      "cmdstanr version 0.9.0 or higher is required for sum_to_zero_vector support.\n\n",
      "Current version: ", packageVersion("cmdstanr"), "\n\n",
      "Please update cmdstanr using:\n",
      "install.packages(pkgs = \"cmdstanr\", repos = c(\"https://mc-stan.org/r-packages/\", getOption(\"repos\")))\n\n",
      "Or install the latest development version:\n",
      "remotes::install_github(\"stan-dev/cmdstanr\")"
    )
  }
  
  # Check if CmdStan is installed
  if (!stan_cmdstan_exists()) {
    
    clear_stan_model_cache()
    
    stop(
      "cmdstan is required to proceed.\n\n",
      "You can install CmdStan by running the following command:\n",
      "cmdstanr::check_cmdstan_toolchain(fix = TRUE)\n",
      "cmdstanr::install_cmdstan()\n",
      "This will install the latest version of CmdStan. For more information, visit:\n",
      "https://mc-stan.org/users/interfaces/cmdstan"
    )
  }
}

#' vb_iterative
#'
#' @description Runs iteratively variational bayes until it suceeds
#'
#'
#' @keywords internal
#' @noRd
#'
#' @param model A Stan model
#' @param output_samples An integer of how many samples from posteriors
#' @param iter An integer of how many max iterations
#' @param tol_rel_obj A real
#' @param additional_parameters_to_save A character vector
#' @param data A data frame
#' @param seed An integer
#' @param sig_figs Number of significant figures to use for Stan model output. Default is 9.
#' @param ... List of paramaters for vb function of Stan
#'
#' @return A Stan fit object
#'
vb_iterative = function(model,
                        output_samples,
                        iter,
                        tol_rel_obj,
                        additional_parameters_to_save = c(),
                        data,
                        output_dir = output_dir,
                        seed, 
                        init = "random",
                        inference_method,
                        cores = 1, 
                        verbose = TRUE,
                        psis_resample = FALSE,
                        sig_figs = 9,
                        ...) {
  res = NULL
  i = 0
  while  (is.null(res) & i < 5) {
    res = tryCatch({
      
      if(inference_method=="pathfinder")
        my_res = model |> 
          sample_safe(
            pathfinder_fx,
            data = data,
            tol_rel_obj = tol_rel_obj,
            output_dir = output_dir,
            seed = seed+i,
            # init = init,
            num_paths=50, 
            num_threads = cores,
            single_path_draws = output_samples / 50 ,
            max_lbfgs_iters=100, 
            history_size = 100, 
            show_messages = verbose,
            psis_resample = psis_resample,
            sig_figs = sig_figs,
            ...
          )
      
      else if(inference_method=="variational")
        my_res = model |> 
          sample_safe(
            variational_fx,
            data = data,
            output_samples = output_samples,
            iter = iter,
            tol_rel_obj = tol_rel_obj,
            output_dir = output_dir,
            seed = seed+i,
            init = init,
            show_messages = verbose,
            threads = cores,
            sig_figs = sig_figs,
            ...
          )
      
      boolFalse <- TRUE
      return(my_res)
    },
    error = function(e) {
      if(e$message |> str_detect("The Stan file used to create the `CmdStanModel` object does not exist\\.")) {
        clear_stan_model_cache()
        model <<-  load_model(model_name, force=TRUE, threads = cores)
      }
      else writeLines(sprintf("Further attempt with Variational Bayes: %s", e))     
      
      return(NULL)
    },
    finally = {
    })
    i = i + 1
  }
  
  if(is.null(res)) stop(sprintf("sccomp says: variational Bayes did not converge after %s attempts. Please use variational_inference = FALSE for a HMC fitting.", i))
  
  return(res)
}

#' Choose the number of chains baed on how many draws we need from the posterior distribution
#' Because there is a fix cost (warmup) to starting a new chain,
#' we need to use the minimum amount that we can parallelise
#' @param how_many_posterior_draws A real number of posterior draws needed
#' @param max_number_to_check A sane upper plateau
#'
#' @keywords internal
#' @noRd
#'
#' @return A Stan fit object
find_optimal_number_of_chains = function(how_many_posterior_draws = 100,
                                         max_number_to_check = 100, warmup = 200, parallelisation_start_penalty = 100) {
  
  
  
  # Define the variables as NULL to avoid CRAN NOTES
  chains <- NULL
  
  
  chains_df =
    tibble(chains = seq_len(max_number_to_check)) %>%
    mutate(tot = (how_many_posterior_draws / chains) + warmup + (parallelisation_start_penalty * chains))
  
  d1 <- diff(chains_df$tot) / diff(seq_len(nrow(chains_df))) # first derivative
  abs(d1) %>% order() %>% .[1] # Find derivative == 0
  
  
}



sample_safe <- function(model, fx, ...) {
  tryCatch(
    {
      # Attempt to run the arbitrary fx with additional arguments
      result <- fx(model, ...)
      return(result)
    },
    error = function(e) {
      # Check if the error message is "model not compiled"
      if (grepl("Model not compiled", e$message)) {
        message("Model not compiled. Compiling.")
        # Load the model with force=TRUE
        model <- model$compile(cpp_options = list(stan_threads = TRUE)) # load_model("glm_multi_beta_binomial_generate_data", force = TRUE)
        # Retry the arbitrary fx with additional arguments
        result <- fx(model, ...)
        return(result)
      } else {
        # If a different error, rethrow the error
        stop(e)
      }
    }
  )
}


sample_fx <- function(model, ...) {
  model$sample(...)
}

pathfinder_fx <- function(model, ...) {
  model$pathfinder(...)
}

variational_fx <- function(model, ...) {
  model$variational(...)
}

generate_quantities_fx <- function(model, ...) {
  model$generate_quantities(...)
}