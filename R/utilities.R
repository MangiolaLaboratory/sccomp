# Define global variable
sccomp_stan_models_cache_dir = file.path(path.expand("~"), ".sccomp_models", packageVersion("sccomp"))

# Greater than
gt = function(a, b){	a > b }

# Smaller than
st = function(a, b){	a < b }

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
#' @importFrom stringr str_subset
#' @importFrom magrittr extract2
#' @importFrom stats terms
#'
#' @return A character vector
#'
#' @keywords internal
#' @noRd
parse_formula <- function(fm) {
  stopifnot("The formula must be of the kind \"~ factors\" " = attr(terms(fm), "response") == 0)
  
  as.character(attr(terms(fm), "variables")) |>
    str_subset("\\|", negate = TRUE) %>%
    
    # Does not work the following
    # |>
    # extract2(-1)
    .[-1]
}


#' Formula parser
#'
#' @param fm A formula
#'
#' @importFrom stringr str_subset
#' @importFrom stringr str_split
#' @importFrom stringr str_remove_all
#' @importFrom rlang set_names
#' @importFrom purrr map_dfr
#' @importFrom stringr str_trim
#'
#' @importFrom magrittr extract2
#'
#' @return A character vector
#'
#' @keywords internal
#' @noRd
formula_to_random_effect_formulae <- function(fm) {
  
  # Define the variables as NULL to avoid CRAN NOTES
  formula <- NULL
  
  stopifnot("The formula must be of the kind \"~ factors\" " = attr(terms(fm), "response") == 0)
  
  random_effect_elements =
    as.character(attr(terms(fm), "variables")) |>
    
    # Select random intercept part
    str_subset("\\|")
  
  if(length(random_effect_elements) > 0){
    
    random_effect_elements |>
      
      # Divide grouping from factors
      str_split("\\|") |>
      
      # Set name
      map_dfr(~ .x |> set_names(c("formula", "grouping"))) |>
      
      # Create formula
      mutate(formula = map(formula, ~ formula(glue("~ {.x}")))) |>
      mutate(grouping = grouping |> str_trim())
    
  }
  
  else
    tibble(`formula` = list(), grouping = character())
  
}

#' Formula parser
#'
#' @param fm A formula
#'
#' @importFrom stringr str_subset
#' @importFrom stringr str_split
#' @importFrom stringr str_remove_all
#' @importFrom rlang set_names
#' @importFrom purrr map_dfr
#'
#' @importFrom magrittr extract2
#'
#' @return A character vector
#'
#' @keywords internal
#' @noRd
parse_formula_random_effect <- function(fm) {
  
  # Define the variables as NULL to avoid CRAN NOTES
  formula <- NULL
  
  stopifnot("The formula must be of the kind \"~ factors\" " = attr(terms(fm), "response") == 0)
  
  random_effect_elements =
    as.character(attr(terms(fm), "variables")) |>
    
    # Select random intercept part
    str_subset("\\|")
  
  if(length(random_effect_elements) > 0){
    
    formula_to_random_effect_formulae(fm) |>
      
      # Divide factors
      mutate(factor = map(
        formula,
        ~
          # Attach intercept
          .x |>
          terms() |>
          attr("intercept") |>
          str_replace("^1$", "(Intercept)") |>
          str_subset("0", negate = TRUE) |>
          
          # Attach variables
          c(
            .x |>
              terms() |>
              attr("variables") |>
              as.character() |>
              str_split("\\+") |>
              as.character() %>%
              .[-1]
          )
      )) |>
      unnest(factor)
    
  }
  
  else
    tibble(factor = character(), grouping = character())
  
}

#' @importFrom tidyr gather
#' @importFrom magrittr set_rownames
#' @importFrom tibble deframe
#' @importFrom magrittr not
#'
#' @keywords internal
#' @noRd
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#'
#' @return A matrix
as_matrix <- function(tbl, rownames = NULL) {
  
  # Define the variables as NULL to avoid CRAN NOTES
  variable <- NULL
  
  # Check for non-numerical columns
  has_non_numerical <- tbl %>%
    {
      if (!is.null(rownames)) {
        dplyr::select(., -contains(rownames))
      } else {
        .
      }
    } %>%
    summarise_all(class) %>%
    gather(variable, class) %>%
    pull(class) %>%
    unique() %>%
    `%in%`(c("numeric", "integer")) %>% 
    not() %>% 
    any()
  
  if (has_non_numerical) {
    warning("to_matrix says: there are NON-numerical columns, the matrix will NOT be numerical")
  }
  
  # Convert to data frame
  tbl <- as.data.frame(tbl)
  
  # Handle rownames if present
  if (!is.null(rownames)) {
    tbl <- tbl %>%
      set_rownames(tbl %>% pull(!!rownames)) %>%
      select(-!!rownames)
  }
  
  # Convert to matrix
  as.matrix(tbl)
}

#' draws_to_tibble_x_y
#'
#' @importFrom tidyr pivot_longer
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
  
  # Define the variables as NULL to avoid CRAN NOTES
  dummy <- NULL
  .variable <- NULL
  .chain <- NULL
  .iteration <- NULL
  .draw <- NULL
  .value <- NULL

  
  fit$draws(variables = par, format = "draws_df") %>%
    mutate(.iteration = seq_len(n())) %>%
    
    pivot_longer(
      names_to = "parameter", # c( ".chain", ".variable", x, y),
      cols = contains(par),
      #names_sep = "\\.?|\\[|,|\\]|:",
      # names_ptypes = list(
      #   ".variable" = character()),
      values_to = ".value"
    ) %>%
    tidyr::extract(parameter, c(".chain", ".variable", x, y), "([1-9]+)?\\.?([a-zA-Z0-9_\\.]+)\\[([0-9]+),([0-9]+)") |> 
    
    # Warning message:
    # Expected 5 pieces. Additional pieces discarded
    suppressWarnings() %>%
    
    mutate(
      !!as.symbol(x) := as.integer(!!as.symbol(x)),
      !!as.symbol(y) := as.integer(!!as.symbol(y))
    ) %>%
    arrange(.variable, !!as.symbol(x), !!as.symbol(y), .chain) %>%
    group_by(.variable, !!as.symbol(x), !!as.symbol(y)) %>%
    mutate(.draw = seq_len(n())) %>%
    ungroup() %>%
    select(!!as.symbol(x), !!as.symbol(y), .chain, .iteration, .draw ,.variable ,     .value) %>%
    filter(.variable == par)
  
}


#' @importFrom tidyr separate
#' @importFrom purrr when
#' 
#' @param fit A fit object from a statistical model, from the 'rstan' package.
#' @param par A character vector specifying the parameters to extract from the fit object.
#' @param x A character string specifying the first index in the parameter names.
#' @param y A character string specifying the second index in the parameter names (optional).
#' @param probs A numerical vector specifying the quantiles to extract.
#' 
#' @keywords internal
#' @noRd
summary_to_tibble = function(fit, par, x, y = NULL, probs = c(0.025, 0.25, 0.50, 0.75, 0.975)) {
  
  # Avoid bug
  #if(fit@stan_args[[1]]$method %>% is.null) fit@stan_args[[1]]$method = "hmc"
  
  summary = 
    fit$summary(variables = par, "mean", ~quantile(.x, probs = probs,  na.rm=TRUE), "rhat", "ess_bulk", "ess_tail") %>%
    rename(.variable = variable ) %>%
    
    when(
      is.null(y) ~ (.) %>% tidyr::separate(col = .variable,  into = c(".variable", x, y), sep="\\[|,|\\]", convert = TRUE, extra="drop"),
      ~ (.) %>% tidyr::separate(col = .variable,  into = c(".variable", x, y), sep="\\[|,|\\]", convert = TRUE, extra="drop")
    )
  
  # # summaries are returned only for HMC
  # if(!"n_eff" %in% colnames(summary)) summary = summary |> mutate(n_eff = NA)
  # if(!"R_k_hat" %in% colnames(summary)) summary = summary |> mutate(R_k_hat = NA)
  
  summary
  
}

#' @importFrom purrr map2_lgl
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr spread
#' @importFrom stats C
#' @importFrom rlang :=
#' @importFrom tibble enframe
#'
#' @keywords internal
#' @noRd
beta_to_CI = function(fitted, censoring_iteration = 1, false_positive_rate, factor_of_interest){
  
  # Define the variables as NULL to avoid CRAN NOTES
  M <- NULL
  C_name <- NULL
  .lower <- NULL
  .median <- NULL
  .upper <- NULL
  
  effect_column_name = sprintf("composition_effect_%s", factor_of_interest) %>% as.symbol()
  
  CI = fitted %>%
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
    pivot_wider(names_from = C_name, values_from=c(.lower , .median ,  .upper)) 
  
  # Create main effect if exists
  if(!is.na(factor_of_interest) )
    CI |>
    mutate(!!effect_column_name := !!as.symbol(sprintf(".median_%s", factor_of_interest))) %>%
    nest(composition_CI = -c(M, !!effect_column_name))
  
  else 
    CI |> nest(composition_CI = -c(M))
  
}

#' @importFrom purrr map2_lgl
#' @importFrom tidyr pivot_wider
#' @importFrom stats C
#' @importFrom rlang :=
#'
#' @keywords internal
#' @noRd
alpha_to_CI = function(fitted, censoring_iteration = 1, false_positive_rate, factor_of_interest){
  
  # Define the variables as NULL to avoid CRAN NOTES
  M <- NULL
  C_name <- NULL
  .lower <- NULL
  .median <- NULL
  .upper <- NULL
  
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

#' Get Random Intercept Design 3
#'
#' This function processes the formula composition elements in the data and creates design matrices
#' for random intercept models.
#'
#' @param .data_ A data frame containing the data.
#' @param .sample A quosure representing the sample variable.
#' @param formula_composition A data frame containing the formula composition elements.
#' 
#' @return A data frame with the processed design matrices for random intercept models.
#' 
#' @importFrom glue glue
#' @importFrom magrittr subtract
#' @importFrom purrr map2
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate_all
#' @importFrom dplyr mutate_if
#' @importFrom dplyr as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom rlang enquo
#' @importFrom rlang quo_name
#' @importFrom tidyselect all_of
#' @importFrom readr type_convert
#' @importFrom tibble rownames_to_column
#' @noRd
get_random_effect_design3 = function(
  .data_, formula, grouping, .sample, 
  accept_NA_as_average_effect = FALSE 
){
  
  # Define the variables as NULL to avoid CRAN NOTES
  .sample = enquo(.sample)
  
  mydesign = .data_ |> get_design_matrix(formula, !!.sample, accept_NA_as_average_effect = accept_NA_as_average_effect)
  
  # Create a matrix of group assignments
  group_matrix = .data_ |> 
    select(all_of(grouping)) |> 
    pull(1) |> 
    rep(ncol(mydesign)) |> 
    matrix(ncol = ncol(mydesign))
  
  # Create a mask where design matrix has non-zero values
  mask = mydesign != 0
  
  # Apply mask to group matrix (setting non-masked values to "")
  group_matrix[!mask] = ""
  
  # Handle NAs in group matrix
  if(accept_NA_as_average_effect) {
    # For rows with NA in grouping, set group to "NA"
    group_matrix[is.na(group_matrix)] = "NA"
  }
  
  colnames(group_matrix) = colnames(mydesign)
  rownames(group_matrix) = rownames(mydesign)
  
  # Create the result matrix
  result = matrix(0, nrow = nrow(mydesign), ncol = ncol(mydesign))
  rownames(result) = rownames(mydesign)
  colnames(result) = colnames(mydesign)
  
  # For each column in the design matrix
  for(col in seq_len(ncol(mydesign))) {
    # Get the non-zero values and their corresponding groups
    non_zero = mydesign[, col] != 0
    if(any(non_zero)) {
      # Set the values in the result matrix
      result[non_zero, col] = mydesign[non_zero, col]
    }
  }
  
  # Convert to long format with the correct column names
  result_long = result |>
    as.data.frame() |>
    rownames_to_column("sample") |>
    pivot_longer(-sample, names_to = "factor", values_to = "value") |>
    filter(value != 0) |>
    mutate(
      group___label = paste0(factor, "___", group_matrix[cbind(sample, factor)]),
      group___numeric = as.integer(factor(group___label)),
      factor___numeric = as.integer(factor(factor))
    ) |>
    select(sample, group___label, group___numeric, value)
  
  result_long
}

#' @importFrom glue glue
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr pull
#' @importFrom tidyr where
#' @importFrom rlang enquo
#' @noRd
get_design_matrix = function(.data_spread, formula, .sample, accept_NA_as_average_effect = FALSE){
  
  .sample = enquo(.sample)
  
  .data_spread = .data_spread %>%
    
    select(!!.sample, parse_formula(formula)) |>
    mutate(across(where(is.numeric),  scale)) 

  # Check for NAs in the data
  has_na = any(is.na(.data_spread |> select(parse_formula(formula))))
  
  # If NAs are present and we don't accept them as average effects, throw an error
  if(has_na && !accept_NA_as_average_effect){
    stop("sccomp says: NA values are present in the design matrix factors. Set accept_NA_as_average_effect = TRUE to handle NAs as average effects across factor levels.")
  }

  # Check if we should handle NAs as average effects
  if(accept_NA_as_average_effect && has_na){
    return(get_design_matrix_with_na_handling(.data_spread, formula, !!.sample))
  }
  
  design_matrix =
    .data_spread |>
    model.matrix(formula, data=_)
  
  rownames(design_matrix) = .data_spread |> pull(!!.sample)
  
  design_matrix
}

#' @description
#' This function handles NA values in a special way for design matrices.
#' When NA values are present in categorical variables (factors or characters), they are treated as a separate
#' level "NA" initially, but then this level is removed from the design matrix. For rows with NA values,
#' the function sets equal weights (0.5) across all remaining factor levels. This approach allows for modeling
#' the average effect across all levels of a factor when the specific level is unknown or missing, rather than
#' excluding these observations or treating NA as its own distinct category. This is particularly useful when
#' you want to include observations with missing factor values but don't want to bias the model toward any
#' particular level of that factor.
#' 
#' @param .data_spread A data frame containing the data.
#' @param formula The formula to use for creating the design matrix.
#' @param .sample A quosure representing the sample variable.
#' 
#' @return A design matrix with rows named by the sample variable.
#' 
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr pull
#' @importFrom tidyr where
#' @importFrom rlang enquo
#' @noRd
get_design_matrix_with_na_handling = function(.data_spread, formula, .sample){
  
  .sample = enquo(.sample)
  
  # Convert NA to a factor level "NA" for categorical variables
  data_with_na <- .data_spread %>%
    mutate(
      # 1) Handle factor/character: turn NA into a real "NA" level, last in the ordering
      across(
        where(~ is.factor(.x) || is.character(.x)),
        ~{
          if (any(is.na(.x))) {
            x_char <- as.character(.x)
            x_char[is.na(x_char)] <- "NA"
            f <- factor(x_char)
            non_na_levels <- setdiff(levels(f), "NA")
            new_levels <- c(non_na_levels, "NA")
            factor(f, levels = new_levels)
          } else {
            .x
          }
        }
      ),
      # 2) Handle numeric: replace NA by the (non-NA) mean of that column
      across(
        where(is.numeric),
        ~{
          if (any(is.na(.x))) {
            mean_val <- mean(.x, na.rm = TRUE)
            ifelse(is.na(.x), mean_val, .x)  # Only replace NAs, keep original values
          } else {
            .x
          }
        }
      )
    )
  
  # Create design matrix
  design_matrix = model.matrix(formula, data = data_with_na, na.action = NULL)
  
  # Process rows that had NA values and remove the "NA" columns
  for(col_name in names(data_with_na)) {
    if(is.factor(data_with_na[[col_name]]) || is.character(data_with_na[[col_name]])) {
      # Get unique levels of this factor (excluding NA)
      factor_levels = levels(factor(data_with_na[[col_name]]))
      factor_levels = factor_levels[factor_levels != "NA"]
      
      # Find the "NA" column for this factor
      na_col_pattern = paste0(col_name, "NA$")
      na_col = grep(na_col_pattern, colnames(design_matrix))
      
      if(length(na_col) > 0) {
        # Find rows with NA values using the design matrix
        na_rows = which(design_matrix[, na_col] == 1)
        
        if(length(na_rows) > 0) {
          # Find all columns related to this factor in the design matrix
          factor_cols = c()
          for(level in factor_levels) {
            # Skip the reference level which won't appear in the design matrix
            level_col = grep(paste0("^", col_name, level, "$"), colnames(design_matrix))
            if(length(level_col) > 0) {
              factor_cols = c(factor_cols, level_col)
            }
          }
          
          # If we can't find factor columns by level names, try the generic pattern
          if(length(factor_cols) == 0) {
            factor_pattern = paste0("^", col_name)
            factor_cols = grep(factor_pattern, colnames(design_matrix))
            factor_cols = setdiff(factor_cols, na_col)
          }
          
          # Remove NA column from design matrix
          design_matrix = design_matrix[, -na_col, drop = FALSE]
          
          # Number of visible levels in design matrix (accounting for reference level)
          n_levels = length(factor_cols)
          
          # For multi-level factors, we need to consider the reference level too
          # The reference level is implied when all other level columns are 0
          # So we distribute 1/(n_levels+1) to each visible level
          # Base level gets effect of 1 - sum(other_effects)
          
          if(n_levels > 0) {
            for(row in na_rows) {
              # Equal weight for each level
              design_matrix[row, factor_cols] = 1/(n_levels + 1)
            }
          }
        }
      }
    }
  }

  rownames(design_matrix) = .data_spread |> pull(!!.sample)
  
  design_matrix
}


#' Check Random Intercept Design
#'
#' This function checks the validity of the random intercept design in the data.
#'
#' @param .data A data frame containing the data.
#' @param factor_names A character vector of factor names.
#' @param random_effect_elements A data frame containing the random intercept elements.
#' @param formula The formula used for the model.
#' @param X The design matrix.
#' 
#' @return A data frame with the checked random intercept elements.
#' 
#' @importFrom tidyr nest
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom rlang set_names
#' @importFrom tidyr unite
#' @importFrom purrr map2
#' @importFrom stringr str_subset
#' @importFrom readr type_convert
#' @noRd
check_random_effect_design = function(.data, factor_names, random_effect_elements, formula, X){
  
  # Define the variables as NULL to avoid CRAN NOTES
  factors <- NULL
  groupings <- NULL
  
  
  .data_ = .data
  
  # Loop across groupings
  random_effect_elements |>
    nest(factors = `factor` ) |>
    mutate(checked = map2(
      grouping, factors,
      ~ {
        
        .y = unlist(.y)
        
        # Check that the group column is categorical
        stopifnot("sccomp says: the grouping column should be categorical (not numeric)" =
                    .data_ |>
                    select(all_of(.x)) |>
                    pull(1) |>
                    class() %in%
                    c("factor", "logical", "character")
        )
        
        
        # # Check sanity of the grouping if only random intercept
        # stopifnot(
        #   "sccomp says: the random intercept completely confounded with one or more discrete factors" =
        #     !(
        #       !.y |> equals("(Intercept)") &&
        #         .data_ |> select(any_of(.y)) |> suppressWarnings() |>  pull(1) |> class() %in% c("factor", "character") |> any() &&
        #         .data_ |>
        #         select(.x, any_of(.y)) |>
        #         select_if(\(x) is.character(x) | is.factor(x) | is.logical(x)) |>
        #         distinct() %>%
        #
        #         # TEMPORARY FIX
        #         set_names(c(colnames(.)[1], 'factor___temp')) |>
        #
        #         count(factor___temp) |>
        #         pull(n) |>
        #         equals(1) |>
        #         any()
        #     )
        # )
        
        # # Check if random intercept with random continuous slope. At the moment is not possible
        # # Because it would require I believe a multivariate prior
        # stopifnot(
        #   "sccomp says: continuous random slope is not supported yet" =
        #     !(
        #       .y |> str_subset("1", negate = TRUE) |> length() |> gt(0) &&
        #         .data_ |>
        #         select(
        #           .y |> str_subset("1", negate = TRUE)
        #         ) |>
        #         map_chr(class) %in%
        #         c("integer", "numeric")
        #     )
        # )
        
        # Check if random intercept with random continuous slope. At the moment is not possible
        # Because it would require I believe a multivariate prior
        stopifnot(
          "sccomp says: currently, discrete random slope is only supported in a intercept-free model. For example ~ 0 + treatment + (treatment | group)" =
            !(
              # If I have both random intercept and random discrete slope
              
              .y |> equals("(Intercept)") |> any() &&
                length(.y) > 1 &&
                # If I have random slope and non-intercept-free model
                .data_ |> select(any_of(.y)) |> suppressWarnings() |>  pull(1) |> class() %in% c("factor", "character") |> any()
              
            )
        )
        
        
      }
    ))
  
  random_effect_elements |>
    nest(groupings = grouping ) |>
    mutate(checked = map2(`factor`, groupings, ~{
      # Check the same group spans multiple factors
      stopifnot(
        "sccomp says: the groups in the formula (factor | group) should be present in only one factor, including the intercept" =
          !(
            # If I duplicated groups
            .y |> unlist() |> length() |> gt(1)
            
          )
      )
      
      
    }))
  
  
  
  
}

#' @importFrom purrr when
#' @importFrom stats model.matrix
#' @importFrom tidyr expand_grid
#' @importFrom stringr str_detect
#' @importFrom stringr str_remove_all
#' @importFrom purrr reduce
#' @importFrom purrr map_int
#' @importFrom stats as.formula
#'
#' @keywords internal
#' @noRd
#'
data_spread_to_model_input =
  function(
    .data_spread, formula, .sample, .cell_type, .count,
    truncation_ajustment = 1, approximate_posterior_inference,
    formula_variability = ~ 1,
    contrasts = NULL,
    bimodal_mean_variability_association = FALSE,
    use_data = TRUE,
    random_effect_elements,
    accept_NA_as_average_effect = FALSE){
    
    # Define the variables as NULL to avoid CRAN NOTES
    exposure <- NULL
    design <- NULL
    mat <- NULL
    factor___numeric <- NULL
    mean_idx <- NULL
    design_matrix <- NULL
    minus_sum <- NULL
    group___numeric <- NULL
    idx <- NULL
    group___label <- NULL
    parameter <- NULL
    group <- NULL
    design_matrix_col <- NULL
    
    # Prepare column same enquo
    .sample = enquo(.sample)
    .cell_type = enquo(.cell_type)
    .count = enquo(.count)
    .grouping_for_random_effect =
      random_effect_elements |>
      pull(grouping) |>
      unique() 
    
    if (length(.grouping_for_random_effect)==0 ) .grouping_for_random_effect = "random_effect"
    
    
    X  =
      
      .data_spread |>
      get_design_matrix(
        # Drop random intercept
        formula |>
          as.character() |>
          str_remove_all("\\+ ?\\(.+\\|.+\\)") |>
          paste(collapse="") |>
          as.formula(),
        !!.sample,
        accept_NA_as_average_effect = accept_NA_as_average_effect
      )
    
    Xa  =
      .data_spread |>
      get_design_matrix(
        # Drop random intercept
        formula_variability |>
          as.character() |>
          str_remove_all("\\+ ?\\(.+\\|.+\\)") |>
          paste(collapse="") |>
          as.formula() ,
        !!.sample,
        accept_NA_as_average_effect = accept_NA_as_average_effect
      )
    
    XA = Xa %>%
      as_tibble() %>%
      distinct()
    
    A = ncol(XA);
    Ar = nrow(XA);
    
    factor_names = parse_formula(formula)
    factor_names_variability = parse_formula(formula_variability)
    cell_cluster_names = .data_spread %>% select(-!!.sample, -any_of(factor_names), -exposure, -!!.grouping_for_random_effect) %>% colnames()
    
    # Random intercept
    if(nrow(random_effect_elements)>0 ) {
      
  
      #check_random_effect_design(.data_spread, any_of(factor_names), random_effect_elements, formula, X)
      random_effect_grouping = formula |>
        formula_to_random_effect_formulae() |>
        mutate(design = map2(
          formula, grouping,
          ~ {
            get_random_effect_design3(.data_spread, .x, .y, !!.sample )
          }))
      
      # Actual parameters, excluding for the sum to one parameters
      is_random_effect = 1
      
      random_effect_grouping =
        random_effect_grouping |>
        mutate(design_matrix = map(
          design,
          ~ ..1 |>
            select(!!.sample, group___label, value) |>
            pivot_wider(names_from = group___label, values_from = value) |>
            mutate(across(everything(), ~ .x |> replace_na(0)))
        )) 
      
      
      X_random_effect = 
        random_effect_grouping |> 
        pull(design_matrix) |> 
        _[[1]]  |>  
        as_matrix(rownames = quo_name(.sample))
      
      # Separate NA group column into X_random_effect_unseen
      X_random_effect_unseen = X_random_effect[, colnames(X_random_effect) |> str_detect("___NA$"), drop = FALSE]
      X_random_effect = X_random_effect[, !colnames(X_random_effect) |> str_detect("___NA$"), drop = FALSE]
      
      # For now that stan does not have tuples, I just allow max two levels
      if(random_effect_grouping |> nrow() > 2) stop("sccomp says: at the moment sccomp allow max two groupings")
      # This will be modularised with the new stan
      if(random_effect_grouping |> nrow() > 1){
        X_random_effect_2 =   
        random_effect_grouping |> 
        pull(design_matrix) |> 
        _[[2]] |>  
        as_matrix(rownames = quo_name(.sample))

        # Separate NA group column into X_random_effect_2_unseen  
        X_random_effect_2_unseen = X_random_effect_2[, colnames(X_random_effect_2) |> str_detect("___NA$"), drop = FALSE]
        X_random_effect_2 = X_random_effect_2[, !colnames(X_random_effect_2) |> str_detect("___NA$"), drop = FALSE]
      }

      else X_random_effect_2 =  X_random_effect[,0,drop=FALSE]
      
      n_random_eff = random_effect_grouping |> nrow()
      
      ncol_X_random_eff =
        random_effect_grouping |>
        mutate(n = map_int(design, ~.x |> distinct(group___numeric) |> nrow())) |>
        pull(n) 
      
      if(ncol_X_random_eff |> length() < 2) ncol_X_random_eff[2] = 0
      
      # TEMPORARY
      group_factor_indexes_for_covariance = 
        X_random_effect |> 
        colnames() |> 
        enframe(value = "parameter", name = "order")  |> 
        separate(parameter, c("factor", "group"), "___", remove = FALSE) |> 
        complete(factor, group, fill = list(order=0)) |> 
        select(-parameter) |> 
        pivot_wider(names_from = group, values_from = order)  |> 
        as_matrix(rownames = "factor")
      
      n_groups = group_factor_indexes_for_covariance |> ncol()
      
      # This will be modularised with the new stan
      if(random_effect_grouping |> nrow() > 1)
        group_factor_indexes_for_covariance_2 = 
        X_random_effect_2 |> 
        colnames() |> 
        enframe(value = "parameter", name = "order")  |> 
        separate(parameter, c("factor", "group"), "___", remove = FALSE) |> 
        complete(factor, group, fill = list(order=0)) |> 
        select(-parameter) |> 
        pivot_wider(names_from = group, values_from = order)  |> 
        as_matrix(rownames = "factor")
      else group_factor_indexes_for_covariance_2 = matrix()[0,0, drop=FALSE]
      
      n_groups = n_groups |> c(group_factor_indexes_for_covariance_2 |> ncol())
      
      how_many_factors_in_random_design = list(group_factor_indexes_for_covariance, group_factor_indexes_for_covariance_2) |> map_int(nrow)
      
      
    } else {
      X_random_effect = matrix(rep(1, nrow(.data_spread)))[,0, drop=FALSE]
      X_random_effect_2 = matrix(rep(1, nrow(.data_spread)))[,0, drop=FALSE] # This will be modularised with the new stan
      is_random_effect = 0
      ncol_X_random_eff = c(0,0)
      n_random_eff = 0
      n_groups = c(0,0)
      how_many_factors_in_random_design = c(0,0)
      group_factor_indexes_for_covariance = matrix()[0,0, drop=FALSE]
      group_factor_indexes_for_covariance_2 = matrix()[0,0, drop=FALSE] # This will be modularised with the new stan
    }
    
    
    y = .data_spread %>% select(-any_of(factor_names), -exposure, -!!.grouping_for_random_effect) %>% as_matrix(rownames = quo_name(.sample))
    
    # If proportion ix 0 issue
    is_proportion = y |> as.numeric() |> max()  |> between(0,1) |> all()
    if(is_proportion){
      y_proportion = y
      y = y[0,,drop = FALSE]
    }
    else{
      y = y
      y_proportion = y[0,,drop = FALSE]
    }
    
    data_for_model =
      list(
        N = .data_spread %>% nrow(),
        M = .data_spread %>% select(-!!.sample, -any_of(factor_names), -exposure, -!!.grouping_for_random_effect) %>% ncol(),
        exposure = .data_spread$exposure,
        is_proportion = is_proportion,
        y = y,
        y_proportion = y_proportion,
        X = X,
        XA = XA,
        Xa = Xa,
        C = ncol(X),
        A = A,
        Ar = Ar,
        truncation_ajustment = truncation_ajustment,
        is_vb = as.integer(approximate_posterior_inference),
        bimodal_mean_variability_association = bimodal_mean_variability_association,
        use_data = use_data,
        
        # Random intercept
        is_random_effect = is_random_effect,
        ncol_X_random_eff = ncol_X_random_eff,
        n_random_eff = n_random_eff,
        n_groups  = n_groups,
        X_random_effect = X_random_effect,
        X_random_effect_2 = X_random_effect_2,
        group_factor_indexes_for_covariance = group_factor_indexes_for_covariance,
        group_factor_indexes_for_covariance_2 = group_factor_indexes_for_covariance_2,
        how_many_factors_in_random_design = how_many_factors_in_random_design,
        
        # For parallel chains
        grainsize = 1,
        
        ## LOO
        enable_loo = FALSE
      )
    
    # Add censoring
    data_for_model$is_truncated = 0
    data_for_model$truncation_up = matrix(rep(-1, data_for_model$M * data_for_model$N), ncol = data_for_model$M)
    data_for_model$truncation_down = matrix(rep(-1, data_for_model$M * data_for_model$N), ncol = data_for_model$M)
    data_for_model$truncation_not_idx = seq_len(data_for_model$M*data_for_model$N)
    data_for_model$TNS = length(data_for_model$truncation_not_idx)
    data_for_model$truncation_not_idx_minimal = matrix(c(1,1), nrow = 1)[0,,drop=FALSE]
    data_for_model$TNIM = 0
    
    # Add parameter factor dictionary
    data_for_model$factor_parameter_dictionary = tibble()
    
    if(.data_spread  |> select(any_of(parse_formula(formula))) |> lapply(class) %in% c("factor", "character") |> any())
      data_for_model$factor_parameter_dictionary =
      data_for_model$factor_parameter_dictionary |> bind_rows(
        # For discrete
        .data_spread  |>
          select(any_of(parse_formula(formula)))  |>
          distinct()  |>
          
          # Drop numerical
          select_if(function(x) !is.numeric(x)) |>
          pivot_longer(everything(), names_to =  "factor", values_to = "parameter") %>%
          unite("design_matrix_col", c(`factor`, parameter), sep="", remove = FALSE)  |>
          select(-parameter) |>
          filter(design_matrix_col %in% colnames(data_for_model$X)) %>%
          distinct()
        
      )
    
    # For continuous
    if(.data_spread  |> select(all_of(parse_formula(formula))) |> lapply(class) |> equals("numeric") |> any())
      data_for_model$factor_parameter_dictionary =
      data_for_model$factor_parameter_dictionary |>
      bind_rows(
        tibble(
          design_matrix_col =  .data_spread  |>
            select(all_of(parse_formula(formula)))  |>
            distinct()  |>
            
            # Drop numerical
            select_if(function(x) is.numeric(x)) |>
            names()
        ) |>
          mutate(`factor` = design_matrix_col)
      )
    
    # If constrasts is set it is a bit more complicated
    if(! is.null(contrasts))
      data_for_model$factor_parameter_dictionary =
      data_for_model$factor_parameter_dictionary |>
      distinct() |>
      expand_grid(parameter=contrasts) |>
      filter(str_detect(parameter, design_matrix_col )) |>
      select(-design_matrix_col) |>
      rename(design_matrix_col = parameter) |>
      distinct()
    
    data_for_model$intercept_in_design = X[,1] |> unique() |> identical(1)
    
    
    if (data_for_model$intercept_in_design | length(factor_names_variability) == 0) {
      data_for_model$A_intercept_columns = 1
    } else {
      data_for_model$A_intercept_columns = 
        .data_spread |> 
        select(any_of(factor_names[1])) |> 
        distinct() |> 
        nrow()
    }
    
    
    if (data_for_model$intercept_in_design ) {
      data_for_model$B_intercept_columns = 1
    } else {
      data_for_model$B_intercept_columns = 
        .data_spread |> 
        select(any_of(factor_names[1])) |> 
        distinct() |> 
        nrow()
    }
    
    # Default all grouping known. This is used for data generation to estimate unknown groupings.
    data_for_model$unknown_grouping = c(FALSE, FALSE)
    
    
    # Return
    data_for_model
  }

#' @importFrom purrr map_int
data_to_spread = function(.data, formula, .sample, .cell_type, .count, .grouping_for_random_effect){
  
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  is_proportion = .data |> pull(!!.count) |> max() <= 1
  
  .data = 
    .data |>
    nest(data = -!!.sample) 
  
  # If proportions exposure = 1
  if(is_proportion) .data = .data |> mutate(exposure = 1)
  else
    .data = 
    .data |>
    mutate(exposure = map_int(data, ~ .x |> pull(!!.count) |> sum() )) 
  
  .data_to_spread = 
    .data |>
    unnest(data) |>
    select(!!.sample, !!.cell_type, exposure, !!.count, parse_formula(formula), any_of(.grouping_for_random_effect)) 
  
  # Check if duplicated samples
  if(
    .data_to_spread |> distinct(!!.sample, !!.cell_type) |> nrow() <
    .data_to_spread |> nrow()
  ) stop("sccomp says: You have duplicated .sample IDs in your input dataset. A .sample .cell_group combination must be unique")
  
  .data_to_spread |>
    spread(!!.cell_type, !!.count)
  
  
}

contrasts_to_parameter_list = function(contrasts, drop_back_quotes = TRUE){
  
  if(is.null(names(contrasts)))
    names(contrasts) <- contrasts
  
  contrast_list <-
    contrasts |>
    
    # Remove fractions used in multiplication before or after '*'
    str_remove_all_ignoring_if_inside_backquotes("([0-9]+/[0-9]+ ?\\* ?)|(\\* ?[0-9]+/[0-9]+)") |>
    
    # Remove decimals used in multiplication before or after '*'
    str_remove_all_ignoring_if_inside_backquotes("([-+]?[0-9]+\\.[0-9]+ ?\\* ?)|(\\* ?[-+]?[0-9]+\\.[0-9]+)") |>
    
    # Remove fractions used in divisions
    str_remove_all_ignoring_if_inside_backquotes("/ ?[0-9]+") |>
    
    # Remove standalone numerical constants not inside backquotes
    str_remove_all_ignoring_if_inside_backquotes("\\b[0-9]+\\b") |>
    
    # Split by "+", "-", "*", "/"
    str_split_ignoring_if_inside_backquotes("\\+|-|\\*|/") |>
    unlist() |>
    
    # Remove parentheses and spaces
    str_remove_all_ignoring_if_inside_backquotes("[\\(\\) ]") |>
    
    # Remove empty strings
    {\(x) x[x != ""]}()
  
  if(drop_back_quotes)
    contrast_list <-
    contrast_list |>
    str_remove_all("`")
  
  contrast_list |> unique()
}




mutate_ignore_error = function(x, ...){
  tryCatch(
    {  x |> mutate(...) },
    error=function(cond) {  x  }
  )
}


#' chatGPT - Remove Specified Regex Pattern from Each String in a Vector
#'
#' This function takes a vector of strings and a regular expression pattern.
#' It removes occurrences of the pattern from each string, except where the pattern
#' is found inside backticks. The function returns a vector of cleaned strings.
#'
#' @param text_vector A character vector with the strings to be processed.
#' @param regex A character string containing a regular expression pattern to be removed
#' from the text.
#'
#' @return A character vector with the regex pattern removed from each string.
#' Occurrences of the pattern inside backticks are not removed.
#'
#' @examples
#' texts <- c("A string with (some) parentheses and `a (parenthesis) inside` backticks",
#'            "Another string with (extra) parentheses")
#' cleaned_texts <- str_remove_all_ignoring_if_inside_backquotes(texts, "\\(")
#' print(cleaned_texts)
#' 
#' @noRd
str_remove_all_ignoring_if_inside_backquotes <- function(text_vector, regex) {
  # Nested function to handle regex removal for a single string
  remove_regex_chars <- function(text, regex) {
    inside_backticks <- FALSE
    result <- ""
    skip <- 0
    
    chars <- strsplit(text, "")[[1]]
    for (i in seq_along(chars)) {
      if (skip > 0) {
        skip <- skip - 1
        next
      }
      
      char <- chars[i]
      if (char == "`") {
        inside_backticks <- !inside_backticks
        result <- paste0(result, char)
      } else if (!inside_backticks) {
        # Check the remaining text against the regex
        remaining_text <- paste(chars[i:length(chars)], collapse = "")
        match <- regexpr(regex, remaining_text)
        
        if (attr(match, "match.length") > 0 && match[1] == 1) {
          # Skip the length of the matched text
          skip <- attr(match, "match.length") - 1
          next
        } else {
          result <- paste0(result, char)
        }
      } else {
        result <- paste0(result, char)
      }
    }
    
    return(result)
  }
  
  # Apply the function to each element in the vector
  sapply(text_vector, remove_regex_chars, regex)
}


#' chatGPT - Split Each String in a Vector by a Specified Regex Pattern
#'
#' This function takes a vector of strings and a regular expression pattern. It splits
#' each string based on the pattern, except where the pattern is found inside backticks.
#' The function returns a list, with each element being a vector of the split segments
#' of the corresponding input string.
#'
#' @param text_vector A character vector with the strings to be processed.
#' @param regex A character string containing a regular expression pattern used for splitting
#' the text.
#'
#' @return A list of character vectors. Each list element corresponds to an input string
#' from `text_vector`, split according to `regex`, excluding occurrences inside backticks.
#'
#' @examples
#' texts <- c("A string with, some, commas, and `a, comma, inside` backticks",
#'            "Another string, with, commas")
#' split_texts <- split_regex_chars_from_vector(texts, ",")
#' print(split_texts)
#' 
#' @noRd
str_split_ignoring_if_inside_backquotes <- function(text_vector, regex) {
  # Nested function to handle regex split for a single string
  split_regex_chars <- function(text, regex) {
    inside_backticks <- FALSE
    result <- c()
    current_segment <- ""
    
    chars <- strsplit(text, "")[[1]]
    for (i in seq_along(chars)) {
      char <- chars[i]
      if (char == "`") {
        inside_backticks <- !inside_backticks
        current_segment <- paste0(current_segment, char)
      } else if (!inside_backticks) {
        # Check the remaining text against the regex
        remaining_text <- paste(chars[i:length(chars)], collapse = "")
        match <- regexpr(regex, remaining_text)
        
        if (attr(match, "match.length") > 0 && match[1] == 1) {
          # Add current segment to result and start a new segment
          result <- c(result, current_segment)
          current_segment <- ""
          # Skip the length of the matched text
          skip <- attr(match, "match.length") - 1
          i <- i + skip
        } else {
          current_segment <- paste0(current_segment, char)
        }
      } else {
        current_segment <- paste0(current_segment, char)
      }
    }
    
    # Add the last segment to the result
    result <- c(result, current_segment)
    return(result)
  }
  
  # Apply the function to each element in the vector
  lapply(text_vector, split_regex_chars, regex)
}

#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#'
#' @param v A array of quosures (e.g. c(col_a, col_b))
#'
#' @return A character vector
quo_names <- function(v) {
  
  v = quo_name(quo_squash(v))
  gsub('^c\\(|`|\\)$', '', v) |>
    strsplit(', ') |>
    unlist()
}

#' Add class to abject
#'
#' @keywords internal
#' @noRd
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_class = function(var, name) {
  
  if(!name %in% class(var)) class(var) <- append(class(var),name, after = 0)
  
  var
}


#' Get Output Samples from a Stan Fit Object
#'
#' This function retrieves the number of output samples from a Stan fit object, 
#' supporting different methods (MHC and Variational) based on the available data within the object.
#'
#' @param fit A `stanfit` object, which is the result of fitting a model via Stan.
#' @return The number of output samples used in the Stan model. 
#'         Returns from MHC if available, otherwise from Variational inference.
#' @examples
#' # Assuming 'fit' is a stanfit object obtained from running a Stan model
#' print("samples_count = get_output_samples(fit)")
#'
#' @noRd
#' 
get_output_samples = function(fit){
  
  # Check if the output_samples field is present in the metadata of the fit object
  # This is generally available when the model is fit using MHC (Markov chain Monte Carlo)
  if(!is.null(fit$metadata()$output_samples)) {
    # Return the output_samples from the metadata
    fit$metadata()$output_samples
  }
  
  # If the output_samples field is not present, check for iter_sampling
  # This occurs typically when the model is fit using Variational inference methods
  else if(!is.null(fit$metadata()$iter_sampling)) {
    # Return the iter_sampling from the metadata
    fit$metadata()$iter_sampling
  }
  else
    fit$metadata()$num_psis_draws
}



#' Check and Install cmdstanr and CmdStan
#'
#' This function checks if the `cmdstanr` package and CmdStan are installed. 
#' If they are not installed, it installs them automatically in non-interactive sessions
#' or asks for permission to install them in interactive sessions.
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
  
  stan_error <- function(message = NULL) {
    stan_stop(
      message = message,
      class = c("stan_error", "stan")
    )
  }
  
  stan_stop <- function(message, class) {
    old <- getOption("rlang_backtrace_on_error")
    on.exit(options(rlang_backtrace_on_error = old))
    options(rlang_backtrace_on_error = "none")
    abort(message = message, class = class, call = emptyenv())
  }
  
  tryCatch(
    rlang::check_installed(
      pkg = "cmdstanr",
      reason = paste(
        "The {cmdstanr} package is required in order to install",
        "CmdStan and run Stan models. Please install it manually using",
        "install.packages(pkgs = \"cmdstanr\",",
        "repos = c(\"https://mc-stan.org/r-packages/\", getOption(\"repos\"))"
      )
    ),
    error = function(e) {
      clear_stan_model_cache()
      stan_error(conditionMessage(e))
    }
  )
  
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


drop_environment <- function(obj) {
  # Check if the object has an environment
  if (!is.null(environment(obj))) {
    environment(obj) <- new.env()
  }
  return(obj)
}


#' harmonise_factor_levels
#'
#' @description
#' A helper function to make sure that factor levels in new data match the old data
#'
#' @param new_data A data frame containing potential factor variables
#' @param old_data A data frame containing the reference factor levels
#'
#' @return A data frame with the same dimensions as `new_data`
#' @noRd
#'
harmonise_factor_levels <- function(new_data, old_data) {
  if (!is.data.frame(new_data) || !is.data.frame(old_data)) {
    stop("Both new_data and old_data must be data frames")
  }
  
  # Ensure arguments are in the correct order
  if (!identical(names(formals(harmonise_factor_levels))[1:2], c("new_data", "old_data"))) {
    warning("Arguments appear to be in the wrong order. 'new_data' should be the data to be harmonized, 'old_data' is the reference.")
  }
  
  # Original implementation follows
  f_old <- sapply(old_data, is.factor)
  f_new <- sapply(new_data, is.factor)
  
  # Get the common factor column names
  common_f <- intersect(names(old_data)[f_old], names(new_data)[f_new])
  
  # If there are no common factor columns, return the new_data as is
  if (length(common_f) == 0) {
    return(new_data)
  }
  
  # For each common factor column
  for (col in common_f) {
    # Get all unique levels from both datasets
    all_levels <- unique(c(levels(old_data[[col]]), levels(new_data[[col]])))
    
    # Check if there are new levels in new_data that don't exist in old_data
    new_levels <- setdiff(levels(new_data[[col]]), levels(old_data[[col]]))
    if (length(new_levels) > 0) {
      message(paste("sccomp says: New levels found in column", col, ":", paste(new_levels, collapse = ", ")))
    }
    
    # Set levels of new_data to match old_data, adding any missing levels
    new_data[[col]] <- factor(new_data[[col]], levels = levels(old_data[[col]]))
  }
  
  # Also check for character columns in new_data that are factors in old_data
  char_new <- sapply(new_data, is.character)
  factor_cols_in_old <- names(old_data)[f_old]
  
  char_to_factor <- intersect(names(new_data)[char_new], factor_cols_in_old)
  
  for (col in char_to_factor) {
    # Convert character to factor with the same levels as in old_data
    new_data[[col]] <- factor(new_data[[col]], levels = levels(old_data[[col]]))
    
    # Check if there are values in new_data that don't exist in old_data's levels
    invalid_values <- setdiff(unique(as.character(new_data[[col]])), levels(old_data[[col]]))
    if (length(invalid_values) > 0 && !all(is.na(invalid_values))) {
      warning(paste("sccomp says: Values in column", col, "not found in reference levels:", paste(invalid_values[!is.na(invalid_values)], collapse = ", ")))
    }
  }
  
  return(new_data)
}




