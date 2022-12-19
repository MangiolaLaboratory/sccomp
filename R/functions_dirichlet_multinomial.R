#' dirichlet_multinomial_glm main
#'
#' @description This function runs the data modelling and statistical test for the hypothesis that a cell_type includes outlier biological replicate.
#'
#' @importFrom tibble as_tibble
#' @import dplyr
#' @importFrom tidyr spread
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom purrr map2
#' @importFrom purrr map_int
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom purrr map
#' @importFrom tibble rowid_to_column
#' @importFrom purrr map_lgl
#' @importFrom rlang :=
#' @importFrom rstan stan_model
#'
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | factor columns | Pvaue column | a significance column
#' @param formula A formula. The sample formula used to perform the differential cell_type abundance analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_type A column name as symbol. The cell_type identifier
#' @param .count A column name as symbol. The cell_type abundance (read count)
#' @param prior_mean_variable_association A list of the form list(intercept = c(4.436925, 1.304049), slope = c(-0.73074903,  0.06532897), standard_deviation = c(0.4527292, 0.3318759)). Where for each parameter, we specify mean and standard deviation. This is used to incorporate prior knowledge about the mean/variability association of cell-type proportions.
#' @param percent_false_positive A real between 0 and 100. It is the aimed percent of cell types being a false positive. For example, percent_false_positive_genes = 1 provide 1 percent of the calls for significant changes that are actually not significant.
#' @param check_outliers A boolean. Whether to check for outliers before the fit.
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param verbose A boolean. Prints progression.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param seed An integer. Used for development and testing purposes
#'
#' @noRd
#'
#'
#' @return A nested tibble `tbl` with cell_type-wise information: `sample wise data` | plot | `ppc samples failed` | `exposure deleterious outliers`
#'
#'
dirichlet_multinomial_glm = function(.data,
                                     formula = ~ 1,
                                     .sample,
                                     .cell_type,
                                     .count,
                                     prior_mean_variable_association = NULL,
                                     percent_false_positive = 5,
                                     check_outliers = FALSE,
                                     approximate_posterior_inference = "none",
                                     variance_association = FALSE,
                                     test_composition_above_logit_fold_change = NULL,
                                     verbose = TRUE,
                                     exclude_priors = FALSE,
                                     cores = detect_cores(), # For development purpose,
                                     seed = sample(1e5, 1),
                                     max_sampling_iterations = NULL,
                                     pass_fit = TRUE,
                                     formula_variability  = NULL
) {
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  data_for_model =
    data_to_spread (.data, formula, !!.sample, !!.cell_type, !!.count) %>%
    data_spread_to_model_input(
      formula, !!.sample, !!.cell_type, !!.count,
      approximate_posterior_inference= approximate_posterior_inference == "all"
    )

  false_positive_rate = percent_false_positive/100


  # Load model
  if(file.exists("model_glm_dirichlet_multinomial.rds"))
    model_glm_dirichlet_multinomial = readRDS("model_glm_dirichlet_multinomial.rds")
  else {
    model_glm_dirichlet_multinomial = stan_model(model_code = glm_dirichlet_multinomial)
    model_glm_dirichlet_multinomial  %>% saveRDS("model_glm_dirichlet_multinomial.rds")

  }

  # Produce data list
  if(!check_outliers){

   fit =
      data_for_model %>%
      # Run the first discovery phase with permissive false discovery rate
      fit_model(model_glm_dirichlet_multinomial, censoring_iteration, cores = cores, chains= 4, pars = c("beta", "precision"), seed = seed, approximate_posterior_inference= approximate_posterior_inference == "all", verbose=verbose)

    parsed_fit =
      fit %>%
      parse_fit(data_for_model, ., censoring_iteration = 1, chains)

    return_df =
      parsed_fit %>%
      {
        # Add precision as attribute
        add_attr(.,
            fit %>% extract(., "precision") %$% precision,
                 "precision"
        )
      } %>%
      beta_to_CI(censoring_iteration = 1, false_positive_rate = false_positive_rate, factor_of_interest = data_for_model$X %>% colnames() %>% .[2] ) %>%

      # add probability
      left_join( get_probability_non_zero_OLD(parsed_fit, prefix = "composition"), by="M" ) %>%

      # Clean
      select(-M) %>%
      mutate(!!.cell_type := data_for_model$y %>% colnames()) %>%
      select(!!.cell_type, everything())

    count_data =
      .data %>%
      nest(count_data = -!!.cell_type)

    return_df =
      return_df %>%
      left_join(  count_data, by = quo_name(.cell_type)  )

  }

  else{

    .data_1 =
      .data %>%
      fit_model_and_parse_out_no_missing_data(data_for_model, model_glm_dirichlet_multinomial, formula, !!.sample, !!.cell_type, !!.count, iteration = 1, chains = 4, seed = seed, approximate_posterior_inference = approximate_posterior_inference != "none")

    return_df =
      .data_1 %>%
      select(-contains("posterior")) %>%
      fit_model_and_parse_out_missing_data(
        model_glm_dirichlet_multinomial,
        formula,
        !!.sample,
        !!.cell_type,
        !!.count,
        iteration = 2,
        seed = seed,
        approximate_posterior_inference = approximate_posterior_inference == "all",
        false_positive_rate = false_positive_rate,
        variance_association = variance_association
      )

  fit = attr(return_df, "fit")
  }

  return_df %>%

    # Attach association mean concentration
    add_attr(fit, "fit") %>%
    add_attr(data_for_model, "model_input")


}

#' @importFrom rlang :=
#' @importFrom purrr map2_lgl
#' @importFrom stats setNames
fit_model_and_parse_out_no_missing_data = function(.data, data_for_model, model_glm_dirichlet_multinomial, formula, .sample, .cell_type, .count, iteration = 1, chains, seed, approximate_posterior_inference){


  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  model_generate = get_model_from_data("model_glm_dirichlet_multinomial_generate_quantities.rds", glm_dirichlet_multinomial_generate_quantities)


  # Run the first discovery phase with permissive false discovery rate
  fit_and_generated  = fit_and_generate_quantities(data_for_model, model_glm_dirichlet_multinomial, model_generate, iteration, chains= 4, output_samples = 5000, seed= seed, approximate_posterior_inference = approximate_posterior_inference)

  # Integrate
  .data %>%
    left_join(tibble(.sample = rownames(data_for_model$y), N = seq_len(nrow(data_for_model$y))), by = ".sample" %>% setNames(quo_name(.sample))) %>%
    left_join(tibble(.cell_type = colnames(data_for_model$y), M = seq_len(ncol(data_for_model$y))), by = ".cell_type" %>% setNames(quo_name(.cell_type))) %>%

    # Add factor from design
    left_join(fit_and_generated, by = c("M", "N")) %>%

    # Add theoretical data quantiles
    mutate(!!as.symbol(sprintf("generated_data_quantiles_%s", iteration)) := map(
      !!as.symbol(sprintf("generated_data_posterior_%s", iteration)),
      ~ quantile(
        .x$generated_quantity,
        probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
      ) %>%
        enframe() %>%
        spread(name, value)
    )) %>%

    # # Attach beta
    # mutate(!!as.symbol(sprintf("beta_quantiles_%s", iteration)) := map(
    #   !!as.symbol(sprintf("beta_posterior_%s", iteration)),
    #   ~ quantile(
    #     .x$.value,
    #     probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    #   ) %>%
    #     enframe() %>%
    #     spread(name, value)
    # ))   %>%

  #     # Join slope
  # left_join(  beta_posterior %>% select(M, !!as.symbol(sprintf("slope_%s", iteration)) := `50%`)  ) %>%

  mutate(!!as.symbol(sprintf("outlier_%s", iteration)) := map2_lgl(
    !!.count, !!as.symbol(sprintf("generated_data_quantiles_%s", iteration)),
    ~ .x < .y$`5%` | .x > .y$`95%`)
  ) %>%

    # Add precision as attribute
    add_attr( attr(fit_and_generated, "precision"), "precision" )



}

#' @importFrom rlang :=
#' @importFrom stats C
#' @importFrom rstan sflist2stanfit
#' @importFrom rstan Rhat
fit_model_and_parse_out_missing_data = function(.data, model_glm_dirichlet_multinomial, formula, .sample, .cell_type, .count, iteration, seed, approximate_posterior_inference, false_positive_rate,
                                                variance_association){


  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  # Produce data list
  data_for_model =
    .data %>%
    data_to_spread (formula, !!.sample, !!.cell_type, !!.count) %>%
    data_spread_to_model_input(formula, !!.sample, !!.cell_type, !!.count, approximate_posterior_inference = approximate_posterior_inference)

  # .count = enquo(.count)
  #
  # .data_wide =
  #   .my_data %>%
  #   select(N, M, !!.count, parse_formula(formula)) %>%
  #   distinct() %>%
  #   spread(M, !!.count)

  # .data_wide_no_factors = .data_wide %>% select(-parse_formula(formula))

  to_exclude =
    .data %>%
    filter(!!as.symbol(sprintf("outlier_%s", iteration - 1))  ) %>%
    distinct(N, M)

  to_include =
    .data %>%
    filter(!!as.symbol(sprintf("outlier_%s", iteration - 1)) %>% `!`  ) %>%
    distinct(N, M)

  # To mix with imputed data
  .data_parsed_inliers =
    .data %>%
    anti_join(to_exclude,by = c("N", "M")) %>%
    select(.value = count, N, M)

  # Dirichlet with missing data

  fit_imputation =
    data_for_model %>%
    do_inference_imputation(
      approximate_posterior_inference = approximate_posterior_inference,
      approximate_posterior_analysis = FALSE,
      cores = cores,
      additional_parameters_to_save = additional_parameters_to_save,
      pass_fit = pass_fit,
      to_include = to_include,
      tol_rel_obj = tol_rel_obj,
      #truncation_compensation = 0.7352941, # Taken by approximation study
      seed = seed,
      precision = .data %>% attr("precision")
    )

  # Check if package is installed, otherwise install
  # if (find.package("furrr", quiet = TRUE) %>% length %>% equals(0)) {
  #   message("Installing furrr")
  #   install.packages("furrr", repos = "https://cloud.r-project.org")
  # }

  if (!requireNamespace("furrr")) stop("sccomp says: please install furrr to use this function.")

  model_generate = get_model_from_data("model_glm_dirichlet_multinomial_generate_quantities.rds", glm_dirichlet_multinomial_generate_quantities)


  beta_posterior_corrected =
    fit_imputation %>%
    draws_to_tibble_x_y("counts", "N", "M") %>%
    rename(.draw_imputation = .draw) %>%
    nest(data = -c(.chain ,.iteration, .draw_imputation ,.variable)) %>%
    sample_n(100) %>%
    mutate(fit_list = furrr::future_map(
      data,
      ~ {
        data_for_model$y =
          .x %>%
          anti_join(to_include,by = c("N", "M")) %>%
          bind_rows(.data_parsed_inliers) %>%
          spread(M, .value) %>%
          as_matrix(rownames = "N")

        # Run model
        fit_and_generate_quantities(data_for_model, model_glm_dirichlet_multinomial, model_generate, iteration, chains=1, output_samples = 200, seed = seed, approximate_posterior_inference = approximate_posterior_inference)
      }
    )) %>%

    # Add precision
    mutate(precision = map(
      fit_list,
      ~ attr(.x, "precision")
    )) %>%
    mutate(fit = map(
      fit_list,
      ~ attr(.x, "fit")
    ))

  fit = beta_posterior_corrected %>% pull(fit) %>% sflist2stanfit()

  beta_summary =
    summary_to_tibble(
    fit, "beta", "C","M",
    probs=c(false_positive_rate/2,  0.5,  1-(false_positive_rate/2))
  ) %>%
    select(-.variable, -n_eff, -`Rhat`, -mean, -se_mean,  -  sd) %>%
    setNames(c(colnames(.)[c(1,2)], ".lower", ".median", ".upper")) %>%
    left_join(
      tibble(
        C=seq_len(ncol(data_for_model$X)),
        C_name = colnames(data_for_model$X)
      ),
    by = "C") %>%
    select(-C) %>%
    pivot_wider(names_from = C_name, values_from=c(.lower , .median ,  .upper))

# Detect outliers
  outlier_df =
    beta_posterior_corrected %>%

    select(fit_list) %>%
    unnest(fit_list) %>%
    .nest_subset(data = -c(N, M)) %>%

    # Merge posterior data
    mutate(!!as.symbol(sprintf(
      "generated_data_posterior_%s", iteration
    )) := map(data,
              ~ .x %>%
                select(!!as.symbol(
                  sprintf("generated_data_posterior_%s", iteration)
                )) %>%
                unnest(!!as.symbol(
                  sprintf("generated_data_posterior_%s", iteration)
                )))) %>%

    # Add theoretical data quantiles
    mutate(!!as.symbol(sprintf(
      "generated_data_quantiles_%s", iteration
    )) := map(
      !!as.symbol(sprintf(
        "generated_data_posterior_%s", iteration
      )),
      ~ quantile(
        .x$generated_quantity,
        probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
      ) %>%
        enframe() %>%
        spread(name, value)
    )) %>%

    select(-data) %>%

    right_join(.data %>% select(!!.count, N, M)) %>%
    unnest(!!as.symbol(
      sprintf("generated_data_quantiles_%s", iteration))) %>%


    mutate(!!as.symbol(sprintf("outlier_%s", iteration)) := !!.count < `5%` | !!.count > `95%`) %>%

    select(N, M, outlier = !!as.symbol(sprintf("outlier_%s", iteration)), .lower = `5%`, .median = `50%`, .upper = `95%`)

count_data =
  .data %>%
  left_join(outlier_df, by = c("N", "M") ) %>%
  select(-M, -N) %>%
  nest(count_data = -!!.cell_type)

beta_summary %>%

  left_join(.data %>% select(!!.cell_type, M), by = "M") %>%
  select(-M) %>%
  left_join(  count_data, by = quo_name(.cell_type)  ) %>%
  select(!!.cell_type, everything()) %>%

    # Add precision as attribute
    add_attr(
      beta_posterior_corrected %>%
        select(precision) %>%
        unnest(precision) %>%
        as.matrix(),
      "precision" ) %>%
    add_attr(  fit,  "fit" )

}

#' @importFrom rlang :=
fit_and_generate_quantities = function(data_for_model, model, model_generate, censoring_iteration, chains, output_samples = 2000, seed, approximate_posterior_inference){


  # fit_discovery  = fit_model(data_for_model, model,  iteration, chains= 4, output_samples = output_samples)

  # Run the first discovery phase with permissive false discovery rate
  fit  = fit_model(
    data_for_model, model, censoring_iteration, chains= chains,
    output_samples = output_samples, verbose = TRUE,
    seed = seed, pars = c("beta", "precision", "alpha"),
    approximate_posterior_inference= approximate_posterior_inference
  )


  # fitted = parse_fit(data_for_model, fit, censoring_iteration = censoring_iteration, chains) %>%
  #   {
  #     # Add precision as attribute
  #     add_attr(.,
  #              fit %>% extract(., "precision") %$% precision,
  #              "precision"
  #     )
  #   }

  # # For building some figure I just need the discovery run, return prematurely
  # if(just_discovery) return(res_discovery %>% filter(.variable == "counts_rng"))

  # Generate theoretical data
  generated_discovery = generate_quantities(fit,  data_for_model, model_generate)

  # Integrate
  data_for_model$X %>%
    as.data.frame %>%
    as_tibble() %>%
    rowid_to_column("N") %>%

    # Drop values for X
    select(N) %>%

    # Add theoretical data posteiror
    left_join(
      generated_discovery %>%
        nest(!!as.symbol(sprintf("generated_data_posterior_%s", censoring_iteration)) := -c(M, N)),
      by="N"
    ) %>%

    # # Attach beta posterior
    # left_join(fitted,  by="M") %>%

    # label_deleterious_outliers()

    # Add precision as attribute
    add_attr(
      fit %>% extract("precision") %$% precision,
      "precision"
    ) %>%
    add_attr(fit, "fit")

}

do_inference_imputation = function(.data,
                                   approximate_posterior_inference = FALSE,
                                   approximate_posterior_analysis = FALSE,
                                   cores,
                                   additional_parameters_to_save,
                                   to_include = tibble(N = integer(), M = integer()),
                                   truncation_compensation = 1,
                                   save_generated_quantities = FALSE,
                                   inits_fx = "random",
                                   prior_from_discovery = tibble(`.variable` = character(),
                                                                 mean = numeric(),
                                                                 sd = numeric()),
                                   pass_fit = FALSE,
                                   tol_rel_obj = 0.01,
                                   write_on_disk = FALSE,
                                   seed,
                                   precision) {




  # # if analysis approximated
  # # If posterior analysis is approximated I just need enough
  # how_many_posterior_draws_practical = ifelse(approximate_posterior_analysis, 1000, how_many_posterior_draws)
  # additional_parameters_to_save = additional_parameters_to_save %>% c("lambda_log_param", "sigma_raw") %>% unique

  # Correct for 0 prop ##############################
  ###################################################

  #https://www.rdocumentation.org/packages/DirichletReg/versions/0.3-0/topics/DR_data
  fix_zeros = function(proportions){
    ( proportions*(nrow(proportions)-1) + (1/ncol(proportions)) ) / nrow(proportions)
  }
  # proportions[proportions==0] = min(proportions[proportions>0])
  # proportions = proportions/apply(proportions,1,sum)

  ###################################################
  ###################################################

  # Convert to log ratios
  .data$y  =
    .data$y %>%
    divide_by(rowSums(.data$y )) %>%
    fix_zeros() %>%
    apply(1, function(x)  x %>% boot::logit() %>% scale(scale = FALSE) %>% as.numeric) %>%
    t()

  # fit =
  #   vb_iterative(
  #     #stanmodels$glm_imputation,
  #     stan_model("inst/stan/glm_imputation.stan"),
  #     output_samples = 2000,
  #     iter = 5000,
  #     tol_rel_obj = 0.01,
  #     data = list(
  #       N = nrow(.data_clr),
  #       M = ncol(.data_clr)-1,
  #       y = .data_clr %>% dplyr::select(-N),
  #       X = X
  #     )
  #   )

  # Load model
  if(file.exists("model_glm_dirichlet_multinomial_imputation.rds"))
    model_imputation = readRDS("model_glm_dirichlet_multinomial_imputation.rds")
  else {
    model_imputation = stan_model(model_code = glm_dirichlet_multinomial_imputation)
    model_imputation  %>% saveRDS("model_glm_dirichlet_multinomial_imputation.rds")

  }

  sampling(
    model_imputation,
    #stan_model("glm_dirichlet_multinomial.stan"),
    data = .data %>%
      c(list(
        precision = precision,
        to_include= to_include,
        how_namy_to_include = to_include %>% nrow,
        I = precision %>% nrow
      )),
    cores = 4
    #, iter = 5000, warmup = 300
  )

  # # Plot results
  # fit %>%
  #   draws_to_tibble_x_y("beta", "C", "M")
  #
  #   tidybayes::gather_draws(beta[C, M]) %>%
  #   median_qi() %>%
  #   filter(C==2) %>%
  #   bind_cols(cell_type = colnames(.data)[-1]) %>%
  #   ggplot(aes(forcats::fct_reorder(cell_type, .value), .value)) +
  #   geom_point() +
  #   geom_errorbar(aes(ymin = .lower, ymax =.upper)) +
  #   geom_hline(yintercept = 0) +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))






}



#' @importFrom tibble enframe
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom boot logit
#' @importFrom SingleCellExperiment counts
#'
#' @keywords internal
#' @noRd
generate_quantities = function(fit, data_for_model, model_generate){

  rstan::gqs(
    model_generate,
    #rstan::stan_model("inst/stan/generated_quantities.stan"),
    draws =  as.matrix(fit),
    data = data_for_model
  ) %>%

    extract("counts") %$% counts %>%
    as.data.frame() %>%
    as_tibble(rownames = "draw") %>%
    gather(N_M, generated_quantity, -draw) %>%
    nest(data = -N_M) %>%
    separate(N_M, c("N", "M")) %>%
    mutate(N = as.integer(N), M = as.integer(M)) %>%
    unnest(data)


}


