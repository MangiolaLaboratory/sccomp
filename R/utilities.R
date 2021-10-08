

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

  par_names = names(fit) %>% grep(sprintf("%s", par), ., value = TRUE)

  fit %>%
    extract(par_names, permuted=FALSE) %>%
    as.data.frame %>%
    as_tibble() %>%
    mutate(.iteration = 1:n()) %>%

    when(!is.null(number_of_draws) ~ sample_n(., number_of_draws), ~ (.)) %>%

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
    mutate(.draw = 1:n()) %>%
    ungroup() %>%
    select(!!as.symbol(x), !!as.symbol(y), .chain, .iteration, .draw ,.variable ,     .value) %>%
    filter(.variable == par)

}

draws_to_tibble_x = function(fit, par, x, number_of_draws = NULL) {

  par_names = names(fit) %>% grep(sprintf("%s", par), ., value = TRUE)

  fit %>%
    extract(par_names, permuted=FALSE) %>%
    as.data.frame %>%
    as_tibble() %>%
    mutate(.iteration = 1:n()) %>%
    pivot_longer(names_to = c("dummy", ".chain", ".variable", x),  cols = contains(par), names_sep = "\\.|\\[|,|\\]|:", names_ptypes = list(".chain" = integer(), ".variable" = character(), "A" = integer(), "C" = integer()), values_to = ".value") %>%
    select(-dummy) %>%
    arrange(.variable, !!as.symbol(x), .chain) %>%
    group_by(.variable, !!as.symbol(x)) %>%
    mutate(.draw = 1:n()) %>%
    ungroup() %>%
    select(!!as.symbol(x), .chain, .iteration, .draw ,.variable ,     .value)

}

#' @importFrom tidyr separate
#' @importFrom purrr when
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
  data_for_model, model, censoring_iteration = 1, cores, quantile = 0.95,
  warmup_samples = 200, approximate_posterior_inference = TRUE, verbose = FALSE,
  seed , pars = c("beta", "alpha", "prec_coeff","prec_sd"), output_samples = NULL, chains=NULL
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
        (.) > 20000 ~ {
          warning("sccomp says: the number of draws used to defined quantiles of the posterior distribution is capped to 20K. This means that for very low probability threshold the quantile could become unreliable. We suggest to limit the probability threshold between 0.1 and 0.01")
          20000
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
      save_warmup = FALSE
    ) %>%
    suppressWarnings()

  else
    vb_iterative(
      model,
      output_samples = output_samples ,
      iter = output_samples,
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
    left_join(tibble(C=1:ncol(data_for_model$X), C_name = colnames(data_for_model$X)), by = "C") %>%
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
    mutate(!!effect_column_name := !!as.symbol(sprintf(".median_%s", factor_of_interest))) %>%
    nest(composition_CI = -c(M, !!effect_column_name))



}

#' @importFrom purrr map2_lgl
#' @importFrom tidyr pivot_wider
#' @importFrom stats C
#' @importFrom rlang :=
#'
#' @keywords internal
#' @noRd
alpha_to_CI = function(fitted, censoring_iteration = 1, false_positive_rate, factor_of_interest){

  effect_column_name = sprintf("heterogeneity_effect_%s", factor_of_interest) %>% as.symbol()

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
    nest(heterogeneity_CI = -c(M, !!effect_column_name))



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
#'
#' @keywords internal
#' @noRd
#'
data_spread_to_model_input =
  function(.data_spread, formula, .sample, .cell_type, .count, variance_association = FALSE, truncation_ajustment = 1, approximate_posterior_inference ){

  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  covariate_names = parse_formula(formula)
  X =
    .data_spread %>%
    select(!!.sample, covariate_names) %>%
    model.matrix(formula, data=.) %>%
    apply(2, function(x) x %>% when(sd(.)==0 ~ (.), ~ scale(., scale=FALSE)))

  XA = variance_association %>%
    when((.) == FALSE ~ X[,1, drop=FALSE], ~ X[,1:2, drop=FALSE]) %>%
    as_tibble() %>%
    distinct()

  A = ncol(XA);
  cell_cluster_names = .data_spread %>% select(-!!.sample, -covariate_names, -exposure) %>% colnames()



  data_for_model =
    list(
      N = .data_spread %>% nrow(),
      M = .data_spread %>% select(-!!.sample, -covariate_names, -exposure) %>% ncol(),
      exposure = .data_spread$exposure,
      y = .data_spread %>% select(-covariate_names, -exposure) %>% as_matrix(rownames = quo_name(.sample)),
      X = X,
      XA = XA,
      C = ncol(X),
      A = A,
      truncation_ajustment = truncation_ajustment,
      is_vb = as.integer(approximate_posterior_inference)
    )

  # Add censoring
  data_for_model$is_truncated = 0
  data_for_model$truncation_up = matrix(rep(-1, data_for_model$M * data_for_model$N), ncol = data_for_model$M)
  data_for_model$truncation_down = matrix(rep(-1, data_for_model$M * data_for_model$N), ncol = data_for_model$M)

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
      when((.) == FALSE ~ X[,1, drop=FALSE], ~ X[,1:2, drop=FALSE]) %>%
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
                                         max_number_to_check = 100, warmup = 200) {

  parallelisation_start_penalty = 60

  chains_df =
    tibble(chains = 1:max_number_to_check) %>%
    mutate(tot = (how_many_posterior_draws / chains) + warmup + (parallelisation_start_penalty * chains))

  d1 <- diff(chains_df$tot) / diff(1:nrow(chains_df)) # first derivative
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
get_probability_non_zero = function(.data, prefix = "", test_above_logit_fold_change = 0){

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

#' @importFrom dplyr left_join
#' @importFrom tidyr expand_grid
#'
#' @export
#'
#'
#'
design_matrix_and_coefficients_to_simulation = function(design_matrix, coefficient_matrix, seed = sample(1:100000, size = 1)){

  design_df = as.data.frame(design_matrix)
  coefficient_df = as.data.frame(coefficient_matrix)

  rownames(design_df) = sprintf("sample_%s", 1:nrow(design_df))
  colnames(design_df) = sprintf("covariate_%s", 1:ncol(design_df))

  rownames(coefficient_df) = sprintf("cell_type_%s", 1:nrow(coefficient_df))
  colnames(coefficient_df) = sprintf("beta_%s", 1:ncol(coefficient_df))

  input_data =
    expand_grid(
    sample = rownames(design_df),
    cell_type = rownames(coefficient_df)
  ) |>
    left_join(design_df |> as_tibble(rownames = "sample") , by = "sample") |>
    left_join(coefficient_df |>as_tibble(rownames = "cell_type"), by = "cell_type")

  simulate_data(.data = input_data,
                formula = ~ covariate_1 ,
                .sample = sample,
                .cell_group = cell_type,
                .coefficients = c(beta_1, beta_2),
                seed = seed
  )


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
