




#' @importFrom stats model.matrix
get_mean_precision = function(fit, data_for_model){
  
  # Define the variables as NULL to avoid CRAN NOTES
  M <- NULL
  `2.5%` <- NULL
  `97.5%` <- NULL
  
  fit %>%
    summary_to_tibble("alpha", "C", "M") %>%

    # Just grub intercept of alpha
    filter(C==1) %>%
    select( M, mean, `2.5%` , `97.5%`) %>%
    nest(concentration = -M)

  # WRONG ATTEMPT TO PLOT ALPHA FROM MULTIPLE factorS
  # fit %>%
  #   draws_to_tibble_x_y("alpha", "C", "M") %>%
  #   nest(data = -c(M)) %>%
  #
  #   # Add Design matrix
  #   left_join(
  #     as_tibble(data_for_model$X, rownames="M") %>%
  #       nest(X = -M) %>%
  #       mutate(M = as.integer(M)),
  #     by="M"
  #   ) %>%
  #   mutate(concentration = map2(
  #     data, X,
  #     ~ .x %>%
  #           select(.draw, C, .value) %>%
  #           spread(C, .value) %>%
  #           as_matrix(rownames=".draw") %*%
  #
  #         ( .y %>% as_matrix() %>% t() ) %>%
  #       quantile(c(0.025, 0.5, 0.975)) %>%
  #       enframe() %>%
  #       spread(name, value) %>%
  #       rename(mean = `50%`)
  #   )) %>%
  #   select(-data, -X)
}

get_mean_precision_association = function(fit){
  
  # fit$summary(variables = "prec_coeff", "mean", ~quantile(.x, probs = c(0.05, 0.95),  na.rm=TRUE))
  
  c(
    fit$summary("prec_coeff")[,2] |> set_names("prec_coeff") ,
    fit$summary("prec_sd")[,2] |> set_names("prec_sd")
  )
}
