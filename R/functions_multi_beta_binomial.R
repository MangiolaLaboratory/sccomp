




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

}

get_mean_precision_association = function(fit){
  
  # fit$summary(variables = "prec_coeff", "mean", ~quantile(.x, probs = c(0.05, 0.95),  na.rm=TRUE))
  
  c(
    fit$summary("prec_coeff")[,2] |> set_names("prec_coeff") ,
    fit$summary("prec_sd")[,2] |> set_names("prec_sd")
  )
}
