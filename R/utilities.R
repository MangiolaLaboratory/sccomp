

# Greater than
gt = function(a, b){	a > b }

# Smaller than
st = function(a, b){	a < b }

# Negation
not = function(is){	!is }

#' Add attribute to abject
#'
#' @keywords internal
#' @noRd
#'
#'
#' @param var A tibble
#' @param attribute An object
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_attr = function(var, attribute, name) {
  attr(var, name) <- attribute
  var
}


#' Formula parser
#'
#' @param fm A formula
#'
#' @return A character vector
#'
#' @keywords internal
#' @noRd
parse_formula <- function(fm) {
  if (attr(terms(fm), "response") == 1)
    stop("The formula must be of the kind \"~ covariates\" ")
  else
    as.character(attr(terms(fm), "variables"))[-1]
}

#' Get matrix from tibble
#'
#' @import dplyr
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @importFrom purrr as_mapper
#'
#' @param .x A tibble
#' @param .p A boolean
#' @param .f1 A function
#' @param .f2 A function
#'
#' @return A tibble
ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
  switch(.p %>% `!` %>% sum(1),
         as_mapper(.f1)(.x),
         if (.f2 %>% is.null %>% `!`)
           as_mapper(.f2)(.x)
         else
           .x)

}

#' @importFrom tidyr gather
#' @importFrom magrittr set_rownames
#'
#' @keywords internal
#' @noRd
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#'
#' @return A matrix
as_matrix <- function(tbl, rownames = NULL) {
  tbl %>%

    ifelse_pipe(
      tbl %>%
        ifelse_pipe(!is.null(rownames),		~ .x %>% dplyr::select(-contains(rownames))) %>%
        summarise_all(class) %>%
        gather(variable, class) %>%
        pull(class) %>%
        unique() %>%
        `%in%`(c("numeric", "integer")) %>% not() %>% any(),
      ~ {
        warning("to_matrix says: there are NON-numerical columns, the matrix will NOT be numerical")
        .x
      }
    ) %>%
    as.data.frame() %>%

    # Deal with rownames column if present
    ifelse_pipe(!is.null(rownames),
                ~ .x  %>%
                  set_rownames(tbl %>% pull(!!rownames)) %>%
                  select(-!!rownames)) %>%

    # Convert to matrix
    as.matrix()
}

#' vb_iterative
#'
#' @description Runs iteratively variational bayes until it suceeds
#'
#' @importFrom rstan vb
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
                        seed,
                        ...) {
  res = NULL
  i = 0
  while (res %>% is.null | i > 5) {
    res = tryCatch({
      my_res = vb(
        model,
        data = data,
        output_samples = output_samples,
        iter = iter,
        tol_rel_obj = tol_rel_obj,
        seed = seed,
        #pars=c("counts_rng", "exposure_rate", "alpha_sub_1", additional_parameters_to_save),
        ...
      )
      boolFalse <- TRUE
      return(my_res)
    },
    error = function(e) {
      i = i + 1
      writeLines(sprintf("Further attempt with Variational Bayes: %s", e))
      return(NULL)
    },
    finally = {
    })
  }

  return(res)
}

#' function to pass initialisation values
#'
#' @importFrom stats setNames
#' @importFrom stats rnorm
#' @importFrom stats sd
#'
#' @keywords internal
#' @noRd
#'
#' @return A list
inits_fx =
  function () {
    pars =
      res_discovery %>%
      filter(`.variable` != "counts_rng") %>%
      distinct(`.variable`) %>%
      pull(1)

    foreach(
      par = pars,
      .final = function(x)
        setNames(x, pars)
    ) %do% {
      res_discovery %>%
        filter(`.variable` == par) %>%
        mutate(init = rnorm(n(), mean, sd)) %>%
        mutate(init = 0) %>%
        select(`.variable`, S, G, init) %>%
        pull(init)
    }
  }

#' fit_to_counts_rng
#'
#' @importFrom tidyr separate
#' @importFrom tidyr nest
#' @importFrom rstan summary
#'
#' @param fit A fit object
#' @param adj_prob_theshold fit real
#'
#' @keywords internal
#' @noRd
fit_to_counts_rng = function(fit, adj_prob_theshold){

  writeLines(sprintf("executing %s", "fit_to_counts_rng"))

  fit %>%
    rstan::summary("counts_rng",
                   prob = c(adj_prob_theshold, 1 - adj_prob_theshold)) %$%
    summary %>%
    as_tibble(rownames = ".variable") %>%
    separate(.variable,
             c(".variable", "S", "G"),
             sep = "[\\[,\\]]",
             extra = "drop") %>%
    mutate(S = S %>% as.integer, G = G %>% as.integer) %>%
    select(-one_of(c("n_eff", "Rhat", "khat"))) %>%
    rename(`.lower` = (.) %>% ncol - 1,
           `.upper` = (.) %>% ncol)
}

#' draws_to_tibble_x_y
#'
#' @importFrom tidyr pivot_longer
#' @importFrom rstan extract
#' @importFrom rlang :=
#'
#' @param fit A fit object
#' @param par A character vector. The parameters to extract.
#' @param x A character. The first index.
#' @param y A character. The first index.
#'
#' @keywords internal
#' @noRd
draws_to_tibble_x_y = function(fit, par, x, y, number_of_draws = NULL) {

  par_names =
    names(fit) %>% grep(sprintf("%s", par), ., value = TRUE)

  fit %>%
    extract(par, permuted=FALSE) %>%
    as.data.frame %>%
    as_tibble() %>%
    mutate(.iteration = seq_len(n())) %>%

    #when(!is.null(number_of_draws) ~ sample_n(., number_of_draws), ~ (.)) %>%

    pivot_longer(
      names_to = c("dummy", ".chain", ".variable", x, y),
      cols = contains(par),
      names_sep = "\\.|\\[|,|\\]|:",
      names_ptypes = list(
        ".variable" = character()),
      values_to = ".value"
    ) %>%

    # Warning message:
    # Expected 5 pieces. Additional pieces discarded
    suppressWarnings() %>%

    mutate(
      !!as.symbol(x) := as.integer(!!as.symbol(x)),
      !!as.symbol(y) := as.integer(!!as.symbol(y))
    ) %>%
    select(-dummy) %>%
    arrange(.variable, !!as.symbol(x), !!as.symbol(y), .chain) %>%
    group_by(.variable, !!as.symbol(x), !!as.symbol(y)) %>%
    mutate(.draw = seq_len(n())) %>%
    ungroup() %>%
    select(!!as.symbol(x), !!as.symbol(y), .chain, .iteration, .draw ,.variable ,     .value) %>%
    filter(.variable == par)

}

draws_to_tibble_x = function(fit, par, x, number_of_draws = NULL) {

  par_names = names(fit) %>% grep(sprintf("%s", par), ., value = TRUE)

  fit %>%
    rstan::extract(par_names, permuted=FALSE) %>%
    as.data.frame %>%
    as_tibble() %>%
    mutate(.iteration = seq_len(n())) %>%
    pivot_longer(names_to = c("dummy", ".chain", ".variable", x),  cols = contains(par), names_sep = "\\.|\\[|,|\\]|:", values_to = ".value") %>%

    mutate(
      !!as.symbol(x) := as.integer(!!as.symbol(x)),
    ) %>%

    select(-dummy) %>%
    arrange(.variable, !!as.symbol(x), .chain) %>%
    group_by(.variable, !!as.symbol(x)) %>%
    mutate(.draw = seq_len(n())) %>%
    ungroup() %>%
    select(!!as.symbol(x), .chain, .iteration, .draw ,.variable ,     .value)

}

#' @importFrom tidyr separate
#' @importFrom purrr when
#' @importFrom rstan summary
#'
#' @param fit A fit object
#' @param par A character vector. The parameters to extract.
#' @param x A character. The first index.
#' @param y A character. The first index.
#' @param probs A numrical vector. The quantiles to extract.
#'
#' @keywords internal
#' @noRd
summary_to_tibble = function(fit, par, x, y = NULL, probs = c(0.025, 0.25, 0.50, 0.75, 0.975)) {

  par_names = names(fit) %>% grep(sprintf("%s", par), ., value = TRUE)

  # Avoid bug
  if(fit@stan_args[[1]]$method %>% is.null) fit@stan_args[[1]]$method = "hmc"

  fit %>%
    rstan::summary(par_names, probs = probs) %$%
    summary %>%
    as_tibble(rownames = ".variable") %>%
    when(
      is.null(y) ~ (.) %>% tidyr::separate(col = .variable,  into = c(".variable", x, y), sep="\\[|,|\\]", convert = TRUE, extra="drop"),
      ~ (.) %>% tidyr::separate(col = .variable,  into = c(".variable", x, y), sep="\\[|,|\\]", convert = TRUE, extra="drop")
    ) %>%
    filter(.variable == par)

}





#' @importFrom rlang :=
label_deleterious_outliers = function(.my_data){

  .my_data %>%

    # join CI
    mutate(outlier_above = !!.count > `95%`) %>%
    mutate(outlier_below = !!.count < `5%`) %>%

    # Mark if on the right of the covariate scale
    mutate(is_group_right = !!as.symbol(colnames(X)[2]) > mean( !!as.symbol(colnames(X)[2]) )) %>%

    # Check if outlier might be deleterious for the statistics
    mutate(
      !!as.symbol(sprintf("deleterious_outlier_%s", iteration)) :=
        (outlier_above & slope > 0 & is_group_right)  |
        (outlier_below & slope > 0 & !is_group_right) |
        (outlier_above & slope < 0 & !is_group_right) |
        (outlier_below & slope < 0 & is_group_right)
    ) %>%

    select(-outlier_above, -outlier_below, -is_group_right)

}

fit_model = function(
  data_for_model, model, censoring_iteration = 1, cores = detectCores(), quantile = 0.95,
  warmup_samples = 300, approximate_posterior_inference = TRUE, verbose = FALSE,
  seed , pars = c("beta", "alpha", "prec_coeff","prec_sd"), output_samples = NULL, chains=NULL, max_sampling_iterations = 20000
)
{


  # # if analysis approximated
  # # If posterior analysis is approximated I just need enough
  # how_many_posterior_draws_practical = ifelse(approximate_posterior_analysis, 1000, how_many_posterior_draws)
  # additional_parameters_to_save = additional_parameters_to_save %>% c("lambda_log_param", "sigma_raw") %>% unique


  # Find number of draws
  draws_supporting_quantile = 50
  if(is.null(output_samples))
    output_samples =
      (draws_supporting_quantile/((1-quantile)/2)) %>% # /2 because I have two tails
      max(4000) %>%

      # If it's bigger than 20K CAP because it would get too extreme
      when(
        (.) > max_sampling_iterations ~ {
          # warning("sccomp says: the number of draws used to defined quantiles of the posterior distribution is capped to 20K. This means that for very low probability threshold the quantile could become unreliable. We suggest to limit the probability threshold between 0.1 and 0.01")
          max_sampling_iterations
        },
        (.)
      )

  # Find optimal number of chains
  if(is.null(chains))
    chains =
      find_optimal_number_of_chains(
        how_many_posterior_draws = output_samples,
        warmup = warmup_samples
      ) %>%
      min(cores)

  init_list=list(
    prec_coeff = c(5,0),
    prec_sd = 1,
    alpha = matrix(c(rep(5, data_for_model$M), rep(0, (data_for_model$A-1) *data_for_model$M)), nrow = data_for_model$A, byrow = TRUE)
  )

  init = map(1:chains, ~ init_list) %>%
    setNames(as.character(1:chains))

  # Fit
  if(!approximate_posterior_inference)
    sampling(
      model,
      data = data_for_model,
      chains = chains,
      cores = chains,
      iter = as.integer(output_samples /chains) + warmup_samples,
      warmup = warmup_samples,
      refresh = ifelse(verbose, 1000, 0),
      seed = seed,
      pars = pars,
      save_warmup = FALSE,
      init = init
    ) %>%
      suppressWarnings()

  else
    vb_iterative(
      model,
      output_samples = output_samples ,
      iter = 10000,
      tol_rel_obj = 0.01,
      data = data_for_model, refresh = ifelse(verbose, 1000, 0),
      seed = seed,
      pars = pars
    ) %>%
      suppressWarnings()


}


#' @importFrom purrr map2_lgl
#' @importFrom tidyr pivot_wider
#' @importFrom rstan extract
#' @importFrom rlang :=
#'
#' @keywords internal
#' @noRd
parse_fit = function(data_for_model, fit, censoring_iteration = 1, chains){

  fit %>%
    draws_to_tibble_x_y("beta", "C", "M") %>%
    left_join(tibble(C=seq_len(ncol(data_for_model$X)), C_name = colnames(data_for_model$X)), by = "C") %>%
    nest(!!as.symbol(sprintf("beta_posterior_%s", censoring_iteration)) := -M)

}

#' @importFrom purrr map2_lgl
#' @importFrom tidyr pivot_wider
#' @importFrom stats C
#' @importFrom rlang :=
#'
#' @keywords internal
#' @noRd
beta_to_CI = function(fitted, censoring_iteration = 1, false_positive_rate, factor_of_interest){

  effect_column_name = sprintf("composition_effect_%s", factor_of_interest) %>% as.symbol()

  fitted %>%
    unnest(!!as.symbol(sprintf("beta_posterior_%s", censoring_iteration))) %>%
    nest(data = -c(M, C, C_name)) %>%
    # Attach beta
    mutate(!!as.symbol(sprintf("beta_quantiles_%s", censoring_iteration)) := map(
      data,
      ~ quantile(
        .x$.value,
        probs = c(false_positive_rate/2,  0.5,  1-(false_positive_rate/2))
      ) %>%
        enframe() %>%
        mutate(name = c(".lower", ".median", ".upper")) %>%
        spread(name, value)
    )) %>%
    unnest(!!as.symbol(sprintf("beta_quantiles_%s", censoring_iteration))) %>%
    select(-data, -C) %>%
    pivot_wider(names_from = C_name, values_from=c(.lower , .median ,  .upper)) %>%

    # Create main effect if exists
    when(
      !is.na(factor_of_interest) ~ mutate(., !!effect_column_name := !!as.symbol(sprintf(".median_%s", factor_of_interest))) %>%
        nest(composition_CI = -c(M, !!effect_column_name)),
      ~ nest(., composition_CI = -c(M))
    )




}

#' @importFrom purrr map2_lgl
#' @importFrom tidyr pivot_wider
#' @importFrom stats C
#' @importFrom rlang :=
#'
#' @keywords internal
#' @noRd
alpha_to_CI = function(fitted, censoring_iteration = 1, false_positive_rate, factor_of_interest){

  effect_column_name = sprintf("variability_effect_%s", factor_of_interest) %>% as.symbol()

  fitted %>%
    unnest(!!as.symbol(sprintf("alpha_%s", censoring_iteration))) %>%
    nest(data = -c(M, C, C_name)) %>%
    # Attach beta
    mutate(!!as.symbol(sprintf("alpha_quantiles_%s", censoring_iteration)) := map(
      data,
      ~ quantile(
        .x$.value,
        probs = c(false_positive_rate/2,  0.5,  1-(false_positive_rate/2))
      ) %>%
        enframe() %>%
        mutate(name = c(".lower", ".median", ".upper")) %>%
        spread(name, value)
    )) %>%
    unnest(!!as.symbol(sprintf("alpha_quantiles_%s", censoring_iteration))) %>%
    select(-data, -C) %>%
    pivot_wider(names_from = C_name, values_from=c(.lower , .median ,  .upper)) %>%
    mutate(!!effect_column_name := !!as.symbol(sprintf(".median_%s", factor_of_interest))) %>%
    nest(variability_CI = -c(M, !!effect_column_name))



}

#' .formula parser
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom stats terms
#'
#' @param fm a formula
#' @return A character vector
#'
#'
parse_formula <- function(fm) {
  if (attr(terms(fm), "response") == 1)
    stop("tidybulk says: The .formula must be of the kind \"~ covariates\" ")
  else
    as.character(attr(terms(fm), "variables"))[-1]
}

#' @importFrom purrr when
#' @importFrom stats model.matrix
#' @importFrom tidyr expand_grid
#' @importFrom stringr str_detect
#'
#' @keywords internal
#' @noRd
#'
data_spread_to_model_input =
  function(
    .data_spread, formula, .sample, .cell_type, .count,
    variance_association = FALSE, truncation_ajustment = 1, approximate_posterior_inference ,
    formula_variability = ~ 1,
    contrasts = NULL,
    bimodal_mean_variability_association = FALSE,
    use_data = TRUE){

    # Prepare column same enquo
    .sample = enquo(.sample)
    .cell_type = enquo(.cell_type)
    .count = enquo(.count)


    get_design_matrix = function(formula, .data_spread){

      .data_spread %>%
        select(!!.sample, parse_formula(formula)) %>%
        model.matrix(formula, data=.) %>%
        apply(2, function(x)
          x %>% when(
            sd(.)==0 ~ (.),

            # If I only have 0 and 1 for a binomial factor
            unique(.) %>% sort() %>% equals(c(0,1)) %>% all() ~ (.),

            # If continuous
            ~ scale(., scale=FALSE)
          )
        )
    }

    X  = get_design_matrix(formula, .data_spread)

    Xa  = get_design_matrix(formula_variability, .data_spread)

    XA = Xa %>%
      as_tibble() %>%
      distinct()

    A = ncol(XA);
    Ar = nrow(XA);

    covariate_names = parse_formula(formula)
    cell_cluster_names = .data_spread %>% select(-!!.sample, -covariate_names, -exposure) %>% colnames()

    data_for_model =
      list(
        N = .data_spread %>% nrow(),
        M = .data_spread %>% select(-!!.sample, -covariate_names, -exposure) %>% ncol(),
        exposure = .data_spread$exposure,
        y = .data_spread %>% select(-covariate_names, -exposure) %>% as_matrix(rownames = quo_name(.sample)),
        X = X,
        XA = XA,
        Xa = Xa,
        C = ncol(X),
        A = A,
        Ar = Ar,
        truncation_ajustment = truncation_ajustment,
        is_vb = as.integer(approximate_posterior_inference),
        bimodal_mean_variability_association = bimodal_mean_variability_association,
        use_data = use_data
      )

    # Add censoring
    data_for_model$is_truncated = 0
    data_for_model$truncation_up = matrix(rep(-1, data_for_model$M * data_for_model$N), ncol = data_for_model$M)
    data_for_model$truncation_down = matrix(rep(-1, data_for_model$M * data_for_model$N), ncol = data_for_model$M)
    data_for_model$truncation_not_idx = seq_len(data_for_model$M*data_for_model$N)
    data_for_model$TNS = length(data_for_model$truncation_not_idx)

    # Add parameter covariate dictionary
    data_for_model$covariate_parameter_dictionary =
      .data_spread  |>
      select(parse_formula(formula))  |>
      distinct()  |>

      # Drop numerical
      select_if(function(x) !is.numeric(x)) |>
      pivot_longer(everything(), names_to =  "covariate", values_to = "parameter") %>%
      unite("design_matrix_col", c(covariate, parameter), sep="", remove = FALSE)  |>
      select(-parameter) |>
      filter(design_matrix_col %in% colnames(data_for_model$X)) %>%
      distinct()

    # If constrasts is set it is a bit more complicated
    if(! contrasts |> is.null())
      data_for_model$covariate_parameter_dictionary =
        data_for_model$covariate_parameter_dictionary |>
        distinct() |>
        expand_grid(parameter=contrasts) |>
        filter(str_detect(contrasts, design_matrix_col )) |>
        select(-design_matrix_col) |>
        rename(design_matrix_col = parameter) |>
        distinct()

    data_for_model$intercept_in_design = X[,1] |> unique() |> identical(1)

    # How many intercept columns
    data_for_model$A_intercept_columns = when(data_for_model$intercept_in_design, (.) ~ 1, ~ .data_spread |> select(covariate_names[1]) |> distinct() |> nrow() )

    # Return
    data_for_model
  }

data_to_spread = function(.data, formula, .sample, .cell_type, .count){

  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  .data %>%
    nest(data = -!!.sample) %>%
    mutate(exposure = map_int(data, ~ .x %>% pull(!!.count) %>% sum() )) %>%
    unnest(data) %>%
    select(!!.sample, !!.cell_type, exposure, !!.count, parse_formula(formula)) %>%
    spread(!!.cell_type, !!.count)

}

#' @importFrom purrr when
#' @importFrom stats model.matrix
#'
#' @keywords internal
#' @noRd
#'
data_simulation_to_model_input =
  function(.data, formula, .sample, .cell_type, .exposure, .coefficients, variance_association = FALSE, truncation_ajustment = 1, approximate_posterior_inference ){

    # Prepare column same enquo
    .sample = enquo(.sample)
    .cell_type = enquo(.cell_type)
    .exposure = enquo(.exposure)
    .coefficients = enquo(.coefficients)

    covariate_names = parse_formula(formula)

    sample_data =
      .data %>%
      select(!!.sample, covariate_names) %>%
      distinct() %>%
      arrange(!!.sample)
    X =
      sample_data %>%
      model.matrix(formula, data=.) %>%
      apply(2, function(x) x %>% when(sd(.)==0 ~ (.), ~ scale(., scale=FALSE))) %>%
      {
        .x = (.)
        rownames(.x) = sample_data %>% pull(!!.sample)
        .x
      }

    XA = variance_association %>%
      when((.) == FALSE ~ X[,1, drop=FALSE], ~ X[,c(1,2), drop=FALSE]) %>%
      as_tibble() %>%
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



  chains_df =
    tibble(chains = seq_len(max_number_to_check)) %>%
    mutate(tot = (how_many_posterior_draws / chains) + warmup + (parallelisation_start_penalty * chains))

  d1 <- diff(chains_df$tot) / diff(seq_len(nrow(chains_df))) # first derivative
  abs(d1) %>% order() %>% .[1] # Find derivative == 0


}


get.elbow.points.indices <- function(x, y, threshold) {
  # From https://stackoverflow.com/questions/41518870/finding-the-elbow-knee-in-a-curve
  d1 <- diff(y) / diff(x) # first derivative
  d2 <- diff(d1) / diff(x[-1]) # second derivative
  indices <- which(abs(d2) > threshold)
  return(indices)
}

#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom stats C
#'
#' @keywords internal
#' @noRd
#'
get_probability_non_zero_OLD = function(.data, prefix = "", test_above_logit_fold_change = 0){

  probability_column_name = sprintf("%s_prob_H0", prefix) %>% as.symbol()

  total_draws = .data %>% pull(2) %>% .[[1]] %>% distinct(.draw) %>% nrow()

  .data %>%
    unnest(2 ) %>%
    filter(C ==2) %>%
    nest(data = -c(M, C_name)) %>%
    mutate(
      bigger_zero = map_int(data, ~ .x %>% filter(.value>test_above_logit_fold_change) %>% nrow),
      smaller_zero = map_int(data, ~ .x %>% filter(.value< -test_above_logit_fold_change) %>% nrow)
    ) %>%
    rowwise() %>%
    mutate(
      !!probability_column_name :=
        1 - (
          max(bigger_zero, smaller_zero) %>%
            #max(1) %>%
            divide_by(total_draws)
        )
    )  %>%
    ungroup() %>%
    select(M, !!probability_column_name)
  # %>%
  # mutate(false_discovery_rate = cummean(prob_non_zero))

}

#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom stats C
#' @importFrom stats setNames
#'
#' @keywords internal
#' @noRd
#'
get_probability_non_zero_ = function(fit, parameter, prefix = "", test_above_logit_fold_change = 0){


  draws = rstan::extract(fit, parameter)[[1]]

  total_draws = dim(draws)[1]


  bigger_zero =
    draws %>%
    apply(2, function(y){
      y %>%
        apply(2, function(x) (x>test_above_logit_fold_change) %>% which %>% length)
    })


  smaller_zero =
    draws %>%
    apply(2, function(y){
      y %>%
        apply(2, function(x) (x< -test_above_logit_fold_change) %>% which %>% length)
    })


  (1 - (pmax(bigger_zero, smaller_zero) / total_draws)) %>%
    as.data.frame() %>%
    rowid_to_column(var = "M")

}

get_probability_non_zero = function(draws, test_above_logit_fold_change = 0, probability_column_name){

  draws %>%
    with_groups(c(M, C_name), ~ .x |> summarise(
      bigger_zero = which(.value>test_above_logit_fold_change) |> length(),
      smaller_zero = which(.value< -test_above_logit_fold_change) |> length(),
      n=n()
    )) |>
    mutate(!!as.symbol(probability_column_name) :=  (1 - (pmax(bigger_zero, smaller_zero) / n)))

}

#' @keywords internal
#' @noRd
#'
parse_generated_quantities = function(rng, number_of_draws = 1){

  draws_to_tibble_x_y(rng, "counts", "N", "M", number_of_draws) %>%
    with_groups(c(.draw, N), ~ .x %>% mutate(generated_proportions = .value/sum(.value))) %>%
    filter(.draw<= number_of_draws) %>%
    rename(generated_counts = .value, replicate = .draw) %>%

    mutate(generated_counts = as.integer(generated_counts)) %>%
    select(M, N, generated_proportions, generated_counts, replicate)

}

#' design_matrix_and_coefficients_to_simulation
#'
#' @description Create simulation from design matrix and coefficient matrix
#'
#' @importFrom dplyr left_join
#' @importFrom tidyr expand_grid
#'
#' @keywords internal
#' @noRd
#'
#' @param design_matrix A matrix
#' @param coefficient_matrix A matrix
#'
#' @return A data frame
#'
#'
#'
design_matrix_and_coefficients_to_simulation = function(

  design_matrix, coefficient_matrix, .estimate_object

){

  design_df = as.data.frame(design_matrix)
  coefficient_df = as.data.frame(coefficient_matrix)

  rownames(design_df) = sprintf("sample_%s", seq_len(nrow(design_df)))
  colnames(design_df) = sprintf("covariate_%s", seq_len(ncol(design_df)))

  rownames(coefficient_df) = sprintf("cell_type_%s", seq_len(nrow(coefficient_df)))
  colnames(coefficient_df) = sprintf("beta_%s", seq_len(ncol(coefficient_df)))

  input_data =
    expand_grid(
      sample = rownames(design_df),
      cell_type = rownames(coefficient_df)
    ) |>
    left_join(design_df |> as_tibble(rownames = "sample") , by = "sample") |>
    left_join(coefficient_df |>as_tibble(rownames = "cell_type"), by = "cell_type")

  simulate_data(.data = input_data,

                .estimate_object = .estimate_object,

                formula = ~ covariate_1 ,
                .sample = sample,
                .cell_group = cell_type,
                .coefficients = c(beta_1, beta_2),
                mcmc_seed = sample(1e5, 1)
  )


}


design_matrix_and_coefficients_to_dir_mult_simulation =function(design_matrix, coefficient_matrix, precision = 100, seed = sample(1:100000, size = 1)){


  # design_df = as.data.frame(design_matrix)
  # coefficient_df = as.data.frame(coefficient_matrix)
  #
  # rownames(design_df) = sprintf("sample_%s", 1:nrow(design_df))
  # colnames(design_df) = sprintf("covariate_%s", 1:ncol(design_df))
  #
  # rownames(coefficient_df) = sprintf("cell_type_%s", 1:nrow(coefficient_df))
  # colnames(coefficient_df) = sprintf("beta_%s", 1:ncol(coefficient_df))

  exposure = 500

  prop.means =
    design_matrix %*%
    t(coefficient_matrix) %>%
    boot::inv.logit()

  extraDistr::rdirmnom(length(design_matrix), exposure, prop.means * precision) %>%
    as_tibble(.name_repair = "unique", rownames = "sample") %>%
    mutate(covariate_1= design_matrix) %>%
    gather(cell_type, generated_counts, -sample, -covariate_1) %>%
    mutate(generated_counts = as.integer(generated_counts))


}

#' @importFrom rlang ensym
class_list_to_counts = function(.data, .sample, .cell_group){

  .sample_for_tidyr = ensym(.sample)
  .cell_group_for_tidyr = ensym(.cell_group)

  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)



  .data %>%
    count(!!.sample,
          !!.cell_group,
          name = "count") %>%

    complete(
      !!.sample_for_tidyr,!!.cell_group_for_tidyr,
      fill = list(count = 0)
    ) %>%
    mutate(count = as.integer(count))
}

#' @importFrom dplyr cummean
get_FDR = function(x){
  enframe(x) %>%
    arrange(value) %>%
    mutate(FDR = cummean(value)) %>%
    arrange(name) %>%
    pull(FDR)
}

#' @importFrom patchwork wrap_plots
#' @importFrom forcats fct_reorder
#' @importFrom tidyr drop_na
plot_1d_intervals = function(.data, .cell_group, significance_threshold= 0.025, my_theme){

  .cell_group = enquo(.cell_group)

  .data |>
    filter(parameter != "(Intercept)") |>

    # Reshape
    pivot_longer(c(contains("c_"), contains("v_")),names_sep = "_" , names_to=c("which", "estimate") ) |>
    drop_na() |>
    pivot_wider(names_from = estimate, values_from = value) |>

    nest(data = -c(parameter, which)) |>
    mutate(plot = pmap(
      list(data, which, parameter),
      ~  ggplot(..1, aes(x=effect, y=fct_reorder(!!.cell_group, effect))) +
        geom_vline(xintercept = 0.2, colour="grey") +
        geom_vline(xintercept = -0.2, colour="grey") +
        geom_errorbar(aes(xmin=lower, xmax=upper, color=FDR<significance_threshold)) +
        geom_point() +
        scale_color_brewer(palette = "Set1") +
        xlab("Credible interval of the slope") +
        ylab("Cell group") +
        ggtitle(sprintf("%s %s", ..2, ..3)) +
        theme(legend.position = "bottom") +
        my_theme
    )) %>%
    pull(plot) |>
    wrap_plots(ncol=2)


}

plot_2d_intervals = function(.data, .cell_group, significance_threshold = 0.025, my_theme){

  .cell_group = enquo(.cell_group)

  # mean-variance association
  .data %>%

    # Filter where I did not inferred the variance
    filter(!is.na(v_effect)) %>%

    # Add labels
    with_groups(
      parameter,
      ~ .x %>%
        arrange(c_FDR) %>%
        mutate(cell_type_label = if_else(row_number()<=3 & c_FDR < significance_threshold & parameter!="(Intercept)", !!.cell_group, ""))
    ) %>%
    with_groups(
      parameter,
      ~ .x %>%
        arrange(v_FDR) %>%
        mutate(cell_type_label = if_else((row_number()<=3 & v_FDR < significance_threshold & parameter!="(Intercept)"), !!.cell_group, cell_type_label))
    ) %>%

    {
      .x = (.)
      # Plot
      ggplot(.x, aes(c_effect, v_effect)) +
        geom_vline(xintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
        geom_hline(yintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
        geom_errorbar(aes(xmin=`c_lower`, xmax=`c_upper`, color=`c_FDR`<significance_threshold, alpha=`c_FDR`<significance_threshold), size=0.2) +
        geom_errorbar(aes(ymin=v_lower, ymax=v_upper, color=`v_FDR`<significance_threshold, alpha=`v_FDR`<significance_threshold), size=0.2) +

        geom_point(size=0.2)  +
        annotate("text", x = 0, y = 3.5, label = "Variable", size=2) +
        annotate("text", x = 5, y = 0, label = "Abundant", size=2, angle=270) +

        geom_text_repel(aes(c_effect, -v_effect, label = cell_type_label), size = 2.5, data = .x %>% filter(cell_type_label!="") ) +

        scale_color_manual(values = c("#D3D3D3", "#E41A1C")) +
        scale_alpha_manual(values = c(0.4, 1)) +
        facet_wrap(~parameter, scales="free") +
        my_theme +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
    }

}

#' @importFrom scales trans_new
#' @importFrom stringr str_replace
#' @importFrom stats quantile
plot_boxplot = function(
    .data, data_proportion, factor_of_interest, .cell_group,
    .sample, significance_threshold = 0.025, my_theme
  ){

  calc_boxplot_stat <- function(x) {
    coef <- 1.5
    n <- sum(!is.na(x))
    # calculate quantiles
    stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
    names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
    iqr <- diff(stats[c(2, 4)])
    # set whiskers
    outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
    if (any(outliers)) {
      stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
    }
    return(stats)
  }

  dropLeadingZero <- function(l){  stringr::str_replace(l, '0(?=.)', '') }

  S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
  IS_sqrt <- function(x){x^2*sign(x)}
  S_sqrt_trans <- function() scales::trans_new("S_sqrt",S_sqrt,IS_sqrt)


  .cell_group = enquo(.cell_group)
  .sample = enquo(.sample)


  if(.data |> attr("contrasts") |> is.null())
    significance_colors =
      .data %>%
      pivot_longer(
        c(contains("c_"), contains("v_")),
        names_pattern = "([cv])_([a-zA-Z0-9]+)",
        names_to = c("which", "stats_name"),
        values_to = "stats_value"
      ) %>%
      filter(stats_name == "FDR") %>%
      filter(parameter != "(Intercept)") %>%
      filter(stats_value < significance_threshold) %>%
      filter(covariate == factor_of_interest) %>%
      unite("name", c(which, parameter), remove = FALSE) %>%
      distinct() %>%
      # Get clean parameter
      mutate(!!as.symbol(factor_of_interest) := str_replace(parameter, sprintf("^%s", covariate), "")) %>%

      with_groups(c(!!.cell_group, !!as.symbol(factor_of_interest)), ~ .x %>% summarise(name = paste(name, collapse = ", ")))

  else
    significance_colors =
      .data %>%
        pivot_longer(
          c(contains("c_"), contains("v_")),
          names_pattern = "([cv])_([a-zA-Z0-9]+)",
          names_to = c("which", "stats_name"),
          values_to = "stats_value"
        ) %>%
        filter(stats_name == "FDR") %>%
        filter(parameter != "(Intercept)") %>%
        filter(stats_value < significance_threshold) %>%
        filter(covariate == factor_of_interest) |>
        mutate(count_data = map(count_data, ~ .x |> select(factor_of_interest) |> distinct())) |>
        unnest(count_data) |>

        # Filter relevant parameters
        mutate( !!as.symbol(factor_of_interest) := as.character(!!as.symbol(factor_of_interest) ) ) |>
        filter(str_detect(parameter, !!as.symbol(factor_of_interest) )) |>

        # Rename
        select(!!.cell_group, !!as.symbol(factor_of_interest), name = parameter) |>

    # Merge contrasts
    with_groups(c(!!.cell_group, !!as.symbol(factor_of_interest)), ~ .x %>% summarise(name = paste(name, collapse = ", ")))



  my_boxplot =  ggplot()

  if("fit" %in% names(attributes(.data))){

    simulated_proportion =
      replicate_data(.data, number_of_draws = 100) %>%
      left_join(data_proportion %>% distinct(!!as.symbol(factor_of_interest), !!.sample, !!.cell_group))

    my_boxplot = my_boxplot +

      stat_summary(
        aes(!!as.symbol(factor_of_interest), (generated_proportions)),
        fun.data = calc_boxplot_stat, geom="boxplot",
        outlier.shape = NA, outlier.color = NA,outlier.size = 0,
        fatten = 0.5, lwd=0.2,
        data =
          simulated_proportion %>%

          # Filter uanitles because of limits
          inner_join( data_proportion %>% distinct(!!as.symbol(factor_of_interest), !!.cell_group)) ,
        color="blue"

      )

    # hideOutliers <- function(x) {
    #   if (x$hoverinfo == 'y') {
    #     x$marker = list(opacity = 0)
    #     x$hoverinfo = NA
    #   }
    #   return(x)
    # }
    #
    # my_boxplot[["x"]][["data"]] <- map(my_boxplot[["x"]][["data"]], ~ hideOutliers(.))

  }

  # Get the exception if no significant cell types. This is not elegant
  if(nrow(significance_colors)==0){
    my_boxplot=
      my_boxplot +

      geom_boxplot(
        aes(!!as.symbol(factor_of_interest), proportion,  group=!!as.symbol(factor_of_interest), fill = NULL), # fill=Effect),
        outlier.shape = NA, outlier.color = NA,outlier.size = 0,
        data =
          data_proportion |>
          mutate(!!as.symbol(factor_of_interest) := as.character(!!as.symbol(factor_of_interest))) %>%
          left_join(significance_colors, by = c(quo_name(.cell_group), factor_of_interest)),
        fatten = 0.5,
        lwd=0.5,
      )
  }

  # If I have significance
  else {
    my_boxplot=
      my_boxplot +

      geom_boxplot(
        aes(!!as.symbol(factor_of_interest), proportion,  group=!!as.symbol(factor_of_interest), fill = name), # fill=Effect),
        outlier.shape = NA, outlier.color = NA,outlier.size = 0,
        data =
          data_proportion |>
          mutate(!!as.symbol(factor_of_interest) := as.character(!!as.symbol(factor_of_interest))) %>%
          left_join(significance_colors, by = c(quo_name(.cell_group), factor_of_interest)),
        fatten = 0.5,
        lwd=0.5,
      )
  }



  my_boxplot +
    geom_jitter(
      aes(!!as.symbol(factor_of_interest), proportion, shape=outlier, color=outlier,  group=!!as.symbol(factor_of_interest)),
      data = data_proportion,
      position=position_jitterdodge(jitter.height = 0, jitter.width = 0.2),
      size = 0.5
    ) +

    # geom_boxplot(
    #   aes(Condition, generated_proportions),
    #   outlier.shape = NA, alpha=0.2,
    #   data = simulated_proportion, fatten = 0.5, size=0.5,
    # ) +
    # geom_jitter(aes(Condition, generated_proportions), color="black" ,alpha=0.2, size = 0.2, data = simulated_proportion) +

    facet_wrap(
      vars(!!.cell_group) ,# forcats::fct_reorder(!!.cell_group, abs(Effect), .desc = TRUE, na.rm=TRUE),
      scales = "free_y",
      nrow = 4
    ) +
    scale_color_manual(values = c("black", "#e11f28")) +
    #scale_fill_manual(values = c("white", "#E2D379")) +
    #scale_fill_distiller(palette = "Spectral", na.value = "white") +
    #scale_color_distiller(palette = "Spectral") +

    scale_y_continuous(trans=S_sqrt_trans(), labels = dropLeadingZero) +
    scale_fill_discrete(na.value = "white") +
    #scale_y_continuous(labels = dropLeadingZero, trans="logit") +
    xlab("Biological condition") +
    ylab("Cell-group proportion") +
    guides(color="none", alpha="none", size="none") +
    labs(fill="Significant difference") +
    my_theme +
    theme(axis.text.x =  element_text(angle=20, hjust = 1))



}

draws_to_statistics = function(draws, contrasts, X, false_positive_rate, test_composition_above_logit_fold_change, prefix = ""){

  factor_of_interest = X %>% colnames()

  if(contrasts |> is.null()){

    draws =
      draws |>
      left_join(tibble(C=seq_len(ncol(X)), parameter = colnames(X)), by = "C") %>%
      select(-C, -.variable)
  }
  else {
    draws =
      draws |>
      pivot_wider(names_from = C, values_from = .value) %>%
      setNames(colnames(.)[1:5] |> c(factor_of_interest)) |>
      mutate_from_expr_list(contrasts) |>
      select(-!!factor_of_interest)

    # If no contrasts of interest just return an empty data frame
    if(ncol(draws)==5) return(draws |> distinct(M))

    draws =
      draws |>
      pivot_longer(-c(1:5), names_to = "parameter", values_to = ".value")

  }

  draws =
    draws |>
    with_groups(c(M, parameter), ~ .x |> summarise(
      lower = quantile(.value, false_positive_rate/2),
      effect = quantile(.value, 0.5),
      upper = quantile(.value, 1-(false_positive_rate/2)),
      bigger_zero = which(.value>test_composition_above_logit_fold_change) |> length(),
      smaller_zero = which(.value< -test_composition_above_logit_fold_change) |> length(),
      n=n()
    )) |>

    # Calculate probability non 0
    mutate(pH0 =  (1 - (pmax(bigger_zero, smaller_zero) / n))) |>
    with_groups(parameter, ~ mutate(.x, FDR = get_FDR(pH0))) |>

    select(M, parameter, lower, effect, upper, pH0, FDR)

  # Setting up names separately because |> is not flexible enough
  draws |>
    setNames(c(colnames(draws)[1:2], sprintf("%s%s", prefix, colnames(draws)[3:ncol(draws)])))
}

enquos_from_list_of_symbols <- function(...) {
  enquos(...)
}

contrasts_to_enquos = function(contrasts){
  contrasts |> enquo() |> quo_names() |> syms() %>% do.call(enquos_from_list_of_symbols, .)
}

#' @importFrom purrr map_dfc
#' @importFrom tibble add_column
#' @importFrom dplyr last_col
mutate_from_expr_list = function(x, formula_expr){
  map_dfc(
    formula_expr,
    ~ x |>
      mutate_ignore_error(!!.x := eval(rlang::parse_expr(.x))) |>
      select(-colnames(x))
  ) |>
    add_column(x, .before = 1)

}

mutate_ignore_error = function(x, ...){
  tryCatch(
    {  x |> mutate(...) },
    error=function(cond) {  x  }
  )
}

simulate_multinomial_logit_linear = function(model_input, sd = 0.51){

  mu = model_input$X %*% model_input$beta

  proportions =
    rnorm(length(mu), mu, sd) %>%
    matrix(nrow = nrow(model_input$X)) %>%
    boot::inv.logit()
  apply(1, function(x) x/sum(x)) %>%
    t()

  rownames(proportions) = rownames(model_input$X)
  colnames(proportions) = colnames(model_input$beta )
}
