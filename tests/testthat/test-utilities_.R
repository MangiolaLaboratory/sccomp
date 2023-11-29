library(dplyr)
library(stringr)

seurat_obj[[]] = seurat_obj[[]] |>  mutate(cell_group = if_else(cell_group |> str_detect("Mono"), "mono", "other"))

estimate = seurat_obj |> 
  sccomp_estimate(
  formula_composition = ~ type ,
  formula_variability = ~ 1,
  sample, cell_group,
  approximate_posterior_inference = FALSE,
  cores = 1,
  mcmc_seed = 42,
  max_sampling_iterations = 100
) 
 
 fit = estimate |> 
  attr("fit") 

data_for_model = 
  estimate |> attr("model_input")

test_that("gt function works correctly", {
  expect_true(gt(5, 3))
  expect_false(gt(3, 5))
  expect_false(gt(3, 3))
})


test_that("st function works correctly", {
  expect_false(st(5, 3))
  expect_true(st(3, 5))
  expect_false(st(3, 3))
})

test_that("not function works correctly", {
  expect_false(not(TRUE))
  expect_true(not(FALSE))
})

library(tibble)

test_that("add_attr function works correctly", {
  # Create a tibble
  my_tibble <- tibble(x = 1:5, y = letters[1:5])
  # Add an attribute
  my_tibble <- add_attr(my_tibble, "example attribute", "my_attr")
  # Check if the attribute is added
  expect_equal(attr(my_tibble, "my_attr"), "example attribute")
})


library(testthat)
library(magrittr)

test_that("parse_formula function works correctly", {
  # Test with a valid formula
  expect_equal(parse_formula(~ a + b), c("a", "b"))
  
  # Test with a formula containing a response variable
  expect_error(parse_formula(a ~ b), "The formula must be of the kind \"~ factors\"")
  
  # Test with a more complex formula
  expect_equal(parse_formula(~ a + b + c), c("a", "b", "c"))
  
  # Test with a formula that includes interaction terms
  expect_equal(parse_formula(~ a + b + a:b), c("a", "b"))

})


library(testthat)
library(magrittr)
library(dplyr)
library(stringr)
library(purrr)
library(glue)
library(tibble)

test_that("formula_to_random_effect_formulae works correctly", {
  # Test with a formula without random effects
  expect_equal(
    formula_to_random_effect_formulae(~ a + b),
    tibble(`formula` = list(), grouping = character())
  )
  
  # Test with a formula with a single random effect
  result <- formula_to_random_effect_formulae(~ a + (b | c))
  expect_equal(result$grouping, "c")
  expect_equal(deparse(result$formula[[1]]), "~b")
  
  # Test with a formula with multiple random effects
  result <- formula_to_random_effect_formulae(~ a + (b | c) + (d | e))
  expect_equal(result$grouping, c("c", "e"))
  expect_equal(deparse(result$formula[[1]]), "~b")
  expect_equal(deparse(result$formula[[2]]), "~d")
  
  # Test with a formula containing nested random effects
  result <- formula_to_random_effect_formulae(~ a + (b | c/d))
  expect_equal(result$grouping, "c/d")
  expect_equal(deparse(result$formula[[1]]), "~b")
  
  # Add more test cases as needed to cover different scenarios and edge cases
})


library(testthat)
library(magrittr)
library(dplyr)
library(stringr)
library(purrr)
library(glue)
library(tibble)

test_that("parse_formula_random_intercept works correctly", {
  # Test with a formula without random intercepts
  expect_equal(
    parse_formula_random_intercept(~ a + b),
    tibble(factor = character(), grouping = character())
  )
  
  # Test with a formula with a single random intercept
  result <- parse_formula_random_intercept(~ a + (b | c))
  expect_equal(result$factor, c("(Intercept)", "b"))
  expect_equal(result$grouping, c("c", "c"))
  
  # Test with a formula with multiple random intercepts
  result <- parse_formula_random_intercept(~ a + (b | c) + (d | e))
  expect_equal(result$factor, c("(Intercept)", "b", "(Intercept)", "d"))
  expect_equal(result$grouping, c("c", "c", "e", "e"))
  
  # Test with a formula containing nested random intercepts
  result <- parse_formula_random_intercept(~ a + (b | c/d))
  expect_equal(result$factor, c("(Intercept)", "b"))
  expect_equal(result$grouping, c("c/d", "c/d"))
  
  # Add more test cases as needed to cover different scenarios and edge cases
})


library(testthat)
library(rstan)
library(magrittr)

# Assuming 'simple_stan_model' is a simplified Stan model suitable for testing
# You would need to create this model before running the tests

test_that("vb_iterative runs variational inference correctly", {

  # Set parameters for the vb_iterative function
  output_samples <- 1000
  iter <- 100
  tol_rel_obj <- 0.01
  seed <- 12345
  
  # Run vb_iterative function
  result <- vb_iterative(
    model = stanmodels$glm_multi_beta_binomial,
    output_samples = output_samples,
    iter = iter,
    tol_rel_obj = tol_rel_obj,
    data = data_for_model,
    seed = seed
  )
  
  # Check if the result is a Stan fit object
  expect_true("stanfit" %in% class(result))
  
  # Optionally, you can add more checks related to the content of the Stan fit object
  # For example, checking the existence of certain parameters, convergence metrics, etc.
  # expect_true("parameter_name" %in% names(rstan::extract(result)))
  # expect_true(all(rstan::summary(result)$summary[, "Rhat"] < 1.1))
})


library(testthat)
library(rstan)
library(tidyr)
library(dplyr)
library(rlang)

# Assuming 'fit' is a pre-fitted Stan model object relevant to your context
# You would need to create this fit object before running tests

test_that("draws_to_tibble_x_y works correctly", {
  
  # Run the function with a pre-fitted Stan model
  result <- sccomp:::draws_to_tibble_x_y(fit, "beta", "index_x", "index_y")
  
  # Perform some basic checks
  expect_true("tbl" %in% class(result))
  expect_true(all(c("index_x", "index_y", ".chain", ".iteration", ".draw", ".variable", ".value") %in% colnames(result)))
  
  # More specific tests can be added depending on the expected output structure
})

test_that("draws_to_tibble_x works correctly", {

  
  # Run the function with a pre-fitted Stan model
  result <- sccomp:::draws_to_tibble_x(fit, "prec_coeff", "index_x")
  
  # Perform some basic checks
  expect_true("tbl" %in% class(result))
  expect_true(all(c("index_x", ".chain", ".iteration", ".draw", ".variable", ".value") %in% colnames(result)))
  
  # More specific tests can be added depending on the expected output structure
})

# Note: These tests require a working Stan model and appropriate data. 
# Adjust parameters and index names according to the specifics of your Stan model.


library(testthat)
library(rstan)
library(tidyr)
library(dplyr)
library(purrr)

# Assuming 'fit' is a pre-fitted Stan model object relevant to your context
# You would need to create this fit object before running tests

test_that("summary_to_tibble works correctly", {
  # Define parameters to extract
  par <- "beta" # Replace with an actual parameter name from your model
  x <- "index_x" # Replace with actual index names
  y <- "index_y" # Optional: replace or set to NULL if not used
  probs <- c(0.025, 0.25, 0.50, 0.75, 0.975)
  percentage = (probs*100) |> as.character() |>  paste0("%")
  
  # Run the function with a pre-fitted Stan model
  result <- sccomp:::summary_to_tibble(fit, par, x, y, probs)
  
  # Perform some basic checks
  expect_true("tbl" %in% class(result))
  expect_true(all(c(".variable", x, percentage) %in% colnames(result)))  # Adjust based on expected columns
  
  # More specific tests can be added depending on the expected output structure
  # For example:
  # expect_true(all(result$mean >= result$`2.5%`))
  # expect_true(all(result$mean <= result$`97.5%`))
  
  # Note: More detailed tests would require knowledge of the specific output structure
  # and values expected from your Stan model and 'summary_to_tibble' function
})

# Note: This test requires a working Stan model and appropriate data. 
# Adjust parameters and index names according to the specifics of your Stan model.

library(testthat)
library(rstan)
library(purrr)

# Assuming 'simple_stan_model' is a simplified Stan model suitable for testing
# You would need to create this model before running the tests


test_that("fit_model runs correctly", {
  # Set up parameters for the test
  quantile <- 0.95
  seed <- 1234
  pars <- c("beta", "alpha", "prec_coeff", "prec_sd")
  
  # Run the fit_model function
  result <- fit_model(
    data_for_model = data_for_model,
    model = stanmodels$glm_multi_beta_binomial,
    quantile = quantile,
    seed = seed,
    pars = pars
  )
  
  # Check if the result is a Stan fit object
  expect_true(inherits(result, "stanfit"))
  
  # Add more checks as needed based on the specifics of your model and fit_model function
  # For example, checking the existence of certain parameters, convergence metrics, etc.
})

# Note: This test requires a working Stan model and appropriate mock data. Adjust the mock data
# and parameters according to the specifics of your Stan model and data structure.

library(testthat)
library(rstan)
library(tidyr)
library(dplyr)
library(purrr)

# Assuming 'fit' is a pre-fitted Stan model object relevant to your context
# You would need to create this fit object before running tests


test_that("parse_fit works correctly", {
  # Run the function with a pre-fitted Stan model and mock data model
  result <- sccomp:::parse_fit(
    data_for_model = data_for_model,
    fit = fit,
    censoring_iteration = 1,
    chains = mock_chains
  )
  
  # Perform basic checks
  expect_true("tbl" %in% class(result))
  
  # Check if the result contains the expected nested tibble structure
  # You might need to adjust these checks based on your specific requirements
  expect_true(any(str_detect(names(result), "beta_posterior_1")))
  
  # More specific tests can be added depending on the expected output structure
})

# Note: This test requires a working Stan model and appropriate mock data. 
# Adjust the mock data and parameters according to the specifics of your Stan model.



library(testthat)
library(dplyr)
library(tidyr)
library(glue)
library(magrittr)
library(purrr)

# Mock .data_ and .sample based on the expected structure and contents
# Replace with actual structure and contents as per your application
mock_data_ <- tibble(
  sample_id = paste0("sample", 1:10),
  grouping_factor = sample(c("A", "B"), 10, replace = TRUE),
  value1 = rnorm(10),
  value2 = rnorm(10)
)

# Example formula composition for the test
formula_composition <- ~ value1 + value2 + (value1 | grouping_factor)

test_that("get_random_intercept_design2 works correctly", {
  result <- get_random_intercept_design2(mock_data_, sample_id, formula_composition)
  
  # Perform basic structural checks
  expect_true("tbl" %in% class(result))
  expect_true("tbl" %in% class(result$design[[1]]))
  
  # Check for the existence of key columns in the result
  expect_true(all(c("factor", "grouping", "mean_idx", "minus_sum", "group___label", "group___numeric", "factor___numeric") %in% colnames(result$design[[1]])))
  
  # More specific tests can be added based on the expected behavior of the function
  # For example, checking the correct mapping of grouping to factors, etc.
})

# Note: Adjust the mock data and parameters according to the specifics of your application and the expected structure of the input data.



library(testthat)
library(dplyr)
library(tidyr)
library(glue)
library(purrr)
library(magrittr)

# Mock .data_ and .sample based on the expected structure and contents
# Replace with actual structure and contents as per your application
mock_data_ <- tibble(
  sample_id = paste0("sample", 1:10),
  grouping_factor = sample(c("A", "B"), 10, replace = TRUE),
  continuous_factor = rnorm(10),
  discrete_factor = sample(c("X", "Y"), 10, replace = TRUE)
)

# Example random intercept elements for the test
mock_random_intercept_elements <- tibble(
  `factor` = c("(Intercept)", "continuous_factor", "discrete_factor"),
  grouping = c("grouping_factor", "grouping_factor", "grouping_factor")
)

test_that("get_random_intercept_design works correctly", {
  result <- get_random_intercept_design(mock_data_, sample_id, mock_random_intercept_elements)
  
  # Perform basic structural checks
  expect_true("tbl" %in% class(result))
  expect_true("tbl" %in% class(result$design[[1]]))
  
  # Check for the existence of key columns in the result
  expect_true(all(c("mean_idx", "minus_sum", "group___numeric", "factor___numeric") %in% colnames(result$design[[1]])))
  
  # Additional tests based on the expected behavior of the function
  # For example, checking that continuous and discrete factors are handled correctly
  # and that the index calculations are as expected
})

# Note: Adjust the mock data and parameters according to the specifics of your application and the expected structure of the input data.


library(testthat)
library(dplyr)
library(glue)
library(rlang)

# Mock data for the test
mock_data_spread <- tibble(
  sample_id = paste0("sample", 1:10),
  factor1 = rnorm(10),
  factor2 = sample(c("A", "B"), 10, replace = TRUE),
  response = rnorm(10)
)

# Example formula for the test
example_formula <- 

test_that("get_design_matrix works correctly", {
  # Generate the design matrix
  
  get_design_matrix(mock_data_spread, response ~ factor1 + factor2, sample_id) |> 
    expect_error("The formula must be of the kind")
  
  
  design_matrix <- get_design_matrix(mock_data_spread,  ~ factor1 + factor2, sample_id)
  
  # Check if the result is a matrix
  expect_true(is.matrix(design_matrix))
  
  # Check if the row names of the matrix match the sample IDs
  expect_equal(rownames(design_matrix), mock_data_spread$sample_id)
  
  # Check for the correct number of columns (including intercept if present)
  # Adjust the expected number of columns based on your formula
  expect_equal(ncol(design_matrix), 3) # Intercept + 2 factors + response
  
  # Additional checks can be added to verify the content and structure of the design matrix
  # For example, checking if factor levels are correctly set, if numeric columns are scaled, etc.
})



# Note: This test assumes the mock data is correctly structured for your specific use case.
# You may need to adjust the mock data or add additional checks to test for specific error conditions.

library(testthat)



test_that("data_spread_to_model_input processes data correctly", {
  
  # Sample data for testing
  spread_data <- tibble(
    sample = c(
      "10x_6K", "10x_8K", "GSE115189", "SCP345_580", "SCP345_860",
      "SCP424_pbmc1", "SCP424_pbmc2", "SCP591", "SI-GA-E5", "SI-GA-E7",
      "SI-GA-E8", "SI-GA-G6", "SI-GA-G7", "SI-GA-G8", "SI-GA-G9",
      "SI-GA-H1", "SI-GA-H3", "SI-GA-H4", "SRR11038995"
    ),
    exposure = c(
      4691, 7524, 2346, 5761, 6422,
      2670, 2984, 569, 4179, 7339,
      8361, 3374, 2560, 5182, 22124,
      2186, 2403, 2793, 9804
    ),
    type = c(
      "healthy", "healthy", "healthy", "healthy", "healthy",
      "healthy", "healthy", "healthy", "cancer", "cancer",
      "cancer", "cancer", "cancer", "cancer", "cancer",
      "cancer", "cancer", "cancer", "healthy"
    ),
    continuous_covariate = c(
      0.782, 0.0746, 0.919, -0.0561, -0.156,
      -1.47, -0.478, 0.418, -0.0162, 0.821,
      0.594, -2.21, 1.12, -0.0449, 0.944,
      1.51, 0.39, -0.621, -1.99
    ),
    random_intercept = rep(1, 19),
    `B immature` = c(
      257, 1105, 264, 530, 989,
      282, 605, 14, 98, 737,
      923, 182, 38, 142, 1111,
      113, 29, 56, 534
    ),
    `B mem` = c(
      365, 69, 34, 93, 167,
      63, 190, 19, 12, 253,
      184, 12, 1, 75, 54,
      10, 3, 1, 200
    ),
    `CD4 cm high cytokine` = c(
      0, 6, 0, 0, 0,
      0, 3, 0, 17, 15,
      110, 99, 10, 38, 66,
      74, 26, 31, 0
    ),
    `CD4 cm ribosome` = c(
      0, 27, 8, 11, 5,
      2, 9, 19, 106, 87,
      241, 177, 53, 75, 1341,
      51, 55, 60, 35
    ),
    `CD4 cm S100A4` = c(
      1006, 1143, 606, 1144, 1054,
      238, 300, 63, 372, 589,
      1083, 658, 226, 711, 240,
      363, 385, 351, 782
    )
  )
  
  
  result <- data_spread_to_model_input(
    spread_data, ~type + continuous_covariate, sample, cell_group, count,
    random_intercept_elements = tibble(
      factor = character(0),  # Empty character vector
      grouping = character(0) # Empty character vector
    ), 
    approximate_posterior_inference = FALSE 
  )
  
  # Test the structure and contents of the result
  expect_type(result, "list")
  expect_true("N" %in% names(result))
  expect_true("M" %in% names(result))
  # Additional checks for other elements in the result
})

# Note: Adjust the mock data and add more checks as per your specific function requirements.


test_that("data_to_spread transforms data correctly", {
  original_data <- tibble(
    sample = c("Sample1", "Sample2"),
    cell_type = c("Type1", "Type2"),
    count = c(100, 200),
    factor1 = c("A", "B")
  )
  formula <- ~ factor1
  sample_col <- "sample"
  cell_type_col <- "cell_type"
  count_col <- "count"
  grouping_for_random_intercept <- c("Type1")
  
  result <- data_to_spread(
    original_data, formula, sample_col, cell_type_col, count_col,
    grouping_for_random_intercept
  )
  
  # Test the structure and contents of the result
  expect_s3_class(result, "tbl")
  expect_true("sample" %in% colnames(result))
  expect_true("Type1" %in% colnames(result))
  expect_true("Type2" %in% colnames(result))
  # Additional checks for other aspects of the transformed data
})

# Note: Adjust the mock data and add more checks as per your specific function requirements.


library(testthat)

test_that("data_simulation_to_model_input processes data correctly", {

  mock_data = 
    tibble(
      sample = c(
        "10x_6K", "10x_6K", "10x_6K", "10x_6K", "10x_6K",
        "10x_6K", "10x_6K", "10x_6K", "10x_6K", "10x_6K"
      ),
      type = c(
        "benign", "benign", "benign", "benign", "benign",
        "cancer", "cancer", "cancer", "cancer", "cancer"
      ),
      phenotype = c(
        "b_cell_macrophage_precursor_or_follicular_LTB_IGKG_CD79a_MS4A1",
        "B_cell:immature", "B_cell:immature_IGLC3_IGLC2",
        "B_cell:Memory_ITM2C_IGHA1_MZB1_JCHAIN", "Dendritic_CD11_CD1_high_mito",
        "HSC_-G-CSF_LMNA", "HSC_-G-CSF_S100A9_LYZ_S100A8_FCN1_CYBB_LST1",
        "HSC_Bcell_mac_precuror?_IRF8", "Mixed_CD3_close_to_myeloid_IL7R_CD7_CD247_CD14low",
        "Mixed_CD3_close_to_myeloid_TUBA4A"
      ),
      count = c(42, 361, 57, 40, 75, 187, 269, 190, 123, 19),
      cell_group = c("BM", "B1", "B2", "B3", "Dm", "H1", "H2", "H3", "TM1", "TM2"),
      b_0 = rep(0, 10),
      b_1 = rep(0, 10),
      .exposure = rep(9804, 10)
    )
  
  # Wrong formula class
  mock_data |> 
    data_simulation_to_model_input(
     " ~ type", sample, cell_group, .exposure, c(b_0, b_1)
    ) |> 
    expect_error("sccomp says: the formula argument must be of formula class")
  
  result <- 
    mock_data |> 
    data_simulation_to_model_input(
     ~ type, sample, cell_group, .exposure, c(b_0, b_1)
  )
  
  # Test the structure and contents of the result
  expect_type(result, "list")
  expect_true("N" %in% names(result))
  expect_true("M" %in% names(result))
  expect_true("exposure" %in% names(result))
  expect_true("X" %in% names(result))
  expect_true("XA" %in% names(result))
  expect_true("C" %in% names(result))
  expect_true("A" %in% names(result))
  expect_true("beta" %in% names(result))
  
  # Additional checks for the content and structure of each element
  expect_equal(result$N, 1)
  expect_equal(result$M, length(unique(mock_data[["cell_group"]])))
  # More checks...
})

# Note: Adjust the mock data and add more checks as per your specific function requirements.


