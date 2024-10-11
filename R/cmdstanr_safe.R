

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