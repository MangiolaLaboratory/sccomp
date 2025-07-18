library(dplyr)
library(tidyr)
library(sccomp)
data("seurat_obj")
data("sce_obj")
data("counts_obj")

counts_obj = 
  counts_obj |>
  mutate(count = count+1) |> 
  with_groups("sample", ~ .x |> mutate(proportion = count/sum(count))) 

set.seed(42)

n_iterations = 1000

if (instantiate::stan_cmdstan_exists()){
  
  my_estimate = 
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ continuous_covariate * type ,
      formula_variability = ~ 1,
      "sample", "cell_group",
      
      cores = 1, 
      inference_method = "pathfinder",
      # mcmc_seed = 42,
      max_sampling_iterations = n_iterations, verbose=FALSE
    )
  
  my_estimate_inverse_factor = 
    seurat_obj[[]] |>
    mutate(type = type |> forcats::fct_relevel("healthy")) |> 
    sccomp_estimate(
      formula_composition = ~ continuous_covariate * type ,
      formula_variability = ~ 1,
      "sample", "cell_group",
      
      cores = 1, 
      inference_method = "pathfinder",
      # mcmc_seed = 42,
      max_sampling_iterations = n_iterations, verbose=FALSE
    )
  
  my_estimate_full_bayes = 
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ continuous_covariate * type ,
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1, 
      inference_method = "hmc",
      mcmc_seed = 42,
      max_sampling_iterations = n_iterations, verbose=FALSE
    )
  
  my_estimate_no_intercept = 
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ 0 + type + continuous_covariate,
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      mcmc_seed = 42,
      max_sampling_iterations = n_iterations, verbose=FALSE
    )
  
  my_estimate_random = 
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ type + (type | group__),
      
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      mcmc_seed = 42,     
      max_sampling_iterations = n_iterations, verbose=FALSE
    )
  
  
  
}

test_that("correct columns",{
  skip_cmdstan()
  
  my_estimate = 
    seurat_obj |>
    sccomp_estimate(
      sample = "sampleX", 
      cell_group = "cell_group"
    ) |> 
    expect_error("Can't select columns that don't exist") |> 
    expect_warning("please check typos in your")
  
  my_estimate = 
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ typeX,
      sample = "sample", 
      cell_group = "cell_group"
    ) |> 
    expect_error("Can't subset elements that don't exist") |> 
    expect_warning("please check typos in your")
  
  
})

test_that("Generate data",{
  skip_cmdstan()

  my_estimate |>

    sccomp_replicate() |>
    nrow() |>
    expect_equal(600)

  # With grouping
  my_estimate_random |>

    sccomp_replicate(~ 1 + type) |>
    nrow() |>
    expect_equal(600)


})

test_that("Predict data",{
  skip_cmdstan()
  library(stringr)
  
  new_data_seurat = seurat_obj[, seurat_obj[[]]$sample %in% c("10x_8K", "SI-GA-E5")] 
  new_data_seurat[[]]$sample = new_data_seurat[[]]$sample |> str_replace("SI", "AB") |>  str_replace("10x", "9x") 
  new_data_tibble = new_data_seurat[[]] |> distinct(sample, type, continuous_covariate, group__)
  
  # With new tibble data
  my_estimate |>
    
    sccomp_predict(
      formula_composition = ~ type,
      new_data = new_data_tibble
    ) |>
    nrow() |>
    expect_equal(60)

  
  # With new covariates
  new_data_tibble$continuous_covariate  =  c(1, 2)
  
  my_estimate |>
    
    sccomp_predict(
      formula_composition = ~ continuous_covariate,
      new_data = new_data_tibble
    ) |>
    nrow() |>
    expect_equal(60)
  
  # With random effects
  my_estimate_random |>
    sccomp_predict(~ 1 + type) |>
    nrow() |>
    expect_equal(600)
  
  # With new seurat data if you ever need this
  my_estimate |>
    
    sccomp_predict(
      formula_composition = ~ type,
      new_data = new_data_tibble, 
      number_of_draws = 1
    ) |>
    nrow() |>
    expect_equal(60)
  
  
  # Predict random
  
  my_estimate_random |> 
    sccomp_predict(
      formula_composition = ~ type + (1 | group__),
      new_data = new_data_tibble, 
      number_of_draws = 1
    )
  
  # Test robust parameter
  prediction_robust = 
    my_estimate |>
    sccomp_predict(
      formula_composition = ~ type,
      new_data = new_data_tibble,
      robust = TRUE,
      number_of_draws = 1
    )
  
  # Test that robust prediction works and returns expected columns
  expect_true("proportion_mean" %in% colnames(prediction_robust))
  expect_true("proportion_lower" %in% colnames(prediction_robust))
  expect_true("proportion_upper" %in% colnames(prediction_robust))
  expect_equal(nrow(prediction_robust), 60)
  
})

test_that("outliers",{
  
  skip_cmdstan()
  
  my_estimate |>
    sccomp_remove_outliers(
      cores = 1, 
      max_sampling_iterations = n_iterations,
      inference_method = "hmc",
       verbose=FALSE
    )
  
})

test_that("multilevel multi beta binomial from Seurat",{
  
  skip_cmdstan()
  
 res =
    seurat_obj |>
    ## filter(cell_group %in% c("NK cycling", "B immature")) |>
    sccomp_estimate(
      formula_composition = ~ type + (1 | group__),
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      mcmc_seed = 42,     
      max_sampling_iterations = n_iterations, verbose=FALSE
    )

  # # Check order
  # res |>
  #   filter(parameter == "typehealthy") |>
  #   arrange(desc(abs(c_effect))) |>
  #   slice(1:3) |>
  #   pull(cell_group) |>
  #   sort() |> 
  #   expect_equal(c("CD4 cm high cytokine", "CD4 ribosome", "Mono NKG7 2"  ) |> sort())

  # Check convergence
  # res |>
  #   filter(c_R_k_hat > 4) |>
  #   nrow() |>
  #   expect_equal(0)

  # res |>
  #   filter(parameter == "typecancer - typehealthy") |>
  #   arrange(desc(abs(c_effect))) |>
  #   slice(1:3) |>
  #   pull(cell_group) |>
  #   sort() |>
  #   expect_equal(c("B mem"  ,  "CD4 cm high cytokine" ,"CD4 ribosome"         ))

  # Check convergence
  # res |>
  #   filter(c_R_k_hat > 4) |>
  #   nrow() |>
  #   expect_equal(0)

})

test_that("multilevel nested",{
  
  skip_cmdstan()
  
  library(tidyseurat)
  library(sccomp)

  res =
    seurat_obj |>
    dplyr::left_join(
      tibble(
        sample = c("SI-GA-H1", "SI-GA-H3", "SI-GA-H4", "SI-GA-G6", "SI-GA-G7",
                   "SI-GA-G8", "SI-GA-E5", "SI-GA-G9", "SI-GA-E7", "SI-GA-E8",
                   "GSE115189", "10x_6K", "10x_8K", "SRR11038995", "SRR7244582",
                   "SCP345_580", "SCP345_860", "SCP424_pbmc1", "SCP424_pbmc2", "SCP591"),
        group__ = c("GROUP1", "GROUP1", "GROUP1", "GROUP1", "GROUP1",
                    "GROUP2", "GROUP2", "GROUP2", "GROUP2", "GROUP2",
                    "GROUP3", "GROUP3", "GROUP3", "GROUP3", "GROUP3",
                    "GROUP4", "GROUP4", "GROUP4", "GROUP4", "GROUP4"),
        nested_group = c("GROUP1_Group_1", "GROUP1_Group_2", "GROUP1_Group_1", 
                         "GROUP1_Group_2", "GROUP1_Group_1", "GROUP2_Group_1", 
                         "GROUP2_Group_2", "GROUP2_Group_1", "GROUP2_Group_2", 
                         "GROUP2_Group_1", "GROUP3_Group_1", "GROUP3_Group_2", 
                         "GROUP3_Group_1", "GROUP3_Group_2", "GROUP3_Group_1", 
                         "GROUP4_Group_1", "GROUP4_Group_2", "GROUP4_Group_1", 
                         "GROUP4_Group_2", "GROUP4_Group_1")
      )
    ) |> 
    sccomp_estimate(
      formula_composition = ~ type + (1 | group__) + (1 | nested_group),
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      mcmc_seed = 42,     
      max_sampling_iterations = n_iterations, verbose=FALSE
    ) |> 
    expect_no_error()
  
  res |> 
    sccomp_predict() |> 
    expect_no_error()
  
  res |> 
    sccomp_predict(formula_composition = ~ type + (1 | nested_group)) |> 
    expect_no_error()
  
})


test_that("multilevel multi beta binomial from Seurat with intercept and continuous covariate",{

  skip_cmdstan()
  
  library(sccomp)

  res =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ continuous_covariate + (1 + continuous_covariate | group__),
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      mcmc_seed = 42,   
      max_sampling_iterations = n_iterations,
      inference_method = "hmc", verbose=FALSE
    )

    expect(
      "CD8 em 1" %in%
      (res |>
        filter(parameter == "continuous_covariate") |>
        arrange(desc(abs(c_effect))) |>
        slice(1:3) |>
        pull(cell_group)),
      TRUE
    )


  # Check convergence
  # res |>
  #   filter(c_R_k_hat > 4) |>
  #   nrow() |>
  #   expect_equal(0)


})

# test_that("wrongly-set groups",{
#
#   # library(tidyseurat)
#   # seurat_obj =
#   #   seurat_obj |>
#   #   nest(data = -c("sample", type)) |>
#   #   mutate(group__wrong = c(1,1,1,1,1, 2,2,2,2,2) |> as.character()) |>
#   #   unnest(data)
#
#     expect_error(
#       object =
#         seurat_obj |>
#         ## filter(cell_group %in% c("NK cycling", "B immature")) |>
#         sccomp_estimate(
#           formula_composition = ~ 0 + type + (type | group__wrong),
#           formula_variability = ~ 1,
#           "sample", "cell_group",
#           contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
#           cores = 1,
#           mcmc_seed = 42,       max_sampling_iterations = n_iterations
#         ) ,
#       regexp = "should not be shared"
#     )
#
# })

test_that("multi beta binomial from Seurat",{

  skip_cmdstan()
  
  my_estimate |>
    filter(parameter == "typehealthy") |>
    arrange(desc(abs(c_effect))) |>
    slice(1) |>
    pull(cell_group) |>
    expect_in(c("B mem", "CD4 cm high cytokine"))
  
  my_estimate_full_bayes |>
    filter(parameter == "typehealthy") |>
    arrange(desc(abs(c_effect))) |>
    slice(1) |>
    pull(cell_group) |>
    expect_in(c("B mem", "CD4 cm high cytokine"))
  
  # Check convergence
  # my_estimate |>
  #   filter(c_R_k_hat > 4) |>
  #   nrow() |>
  #   expect_equal(0)

})

test_that("calculate residuals",{

  skip_cmdstan()
  library(dplyr)
  
  my_estimate_random |> 
    sccomp:::sccomp_calculate_residuals() |> 
    pull(residuals) |> 
    max() |> 
    expect_lt(1)
  
})

test_that("remove unwanted effects",{

  skip_cmdstan()
  
  library(tidyseurat)
  
  data =
    seurat_obj |>
    
    # Add batch
    nest(data = -c("sample", type)) |>
    mutate(batch = rep(c(0,1), 10)) |>
    unnest(data)
  
  # Estimate
  estimate =
    data |>
    sccomp_estimate(
      formula_composition = ~ type + batch,
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      mcmc_seed = 43,    
      max_sampling_iterations = n_iterations, verbose = FALSE
    )

})

test_that("multi beta binomial from SCE",{

  skip_cmdstan()
  
    res =
      sce_obj |>
    sccomp_estimate(
      formula_composition = ~ type,
      formula_variability = ~ 1,
      sample = "sample",
      cell_group = "cell_group",
      cores = 1,
      mcmc_seed = 42,      
      max_sampling_iterations = n_iterations, verbose = FALSE
    )

  # res |>
  #   filter(parameter == "typehealthy") |>
  #   arrange(desc(abs(c_effect))) |>
  #   slice(1) |>
  #   pull(cell_group) |>
  #   expect_equal(c("B mem"  ))

  # Check convergence
  # res |>
  #   filter(c_R_k_hat > 4) |>
  #   nrow() |>
  #   expect_equal(0)
})

if (instantiate::stan_cmdstan_exists()){
res_composition =
  seurat_obj[[]] |>
  sccomp_estimate(
    formula_composition = ~ type,
    formula_variability = ~ 1,
    "sample",
    "cell_group",
    cores = 1,
    mcmc_seed = 42,   
    max_sampling_iterations = n_iterations, verbose = FALSE
  )

res_composition_variability =
  seurat_obj[[]] |>
  sccomp_estimate(
    formula_composition = ~ type,
    formula_variability = ~ type,
    "sample",
    "cell_group",
    cores = 1,
    mcmc_seed = 42,    
    max_sampling_iterations = n_iterations, verbose = FALSE
  )
}

test_that("multi beta binomial from metadata",{

  skip_cmdstan()
  
  # res_composition  |>
  #   filter(parameter == "typehealthy") |>
  #   arrange(desc(abs(c_effect))) |>
  #   slice(1) |>
  #   pull(cell_group) |>
  #   expect_equal(c("B mem"  ))

  # Check convergence
  # res_composition |>
  #   filter(c_R_k_hat > 4) |>
  #   nrow() |>
  #   expect_equal(0)

})

test_that("plot test composition",{

  skip_cmdstan()
  
  my_estimate |> 
    sccomp_test() |> 
    plot()

  my_estimate |> 
    plot() |> 
    expect_error("sccomp says: to produce plots, you need to run the function sccomp_test")
  
})

test_that("plot test variability",{

  skip_cmdstan()
  
  res_composition_variability |> 
    sccomp_test() |> 
    plot()


})

# Test for plot_boxplot function
test_that("plot_boxplot function works correctly", {
  
  skip_cmdstan()
  
  my_estimate |> 
    sccomp_test() |> 
    sccomp_boxplot("type", significance_threshold = 0.025) |> 
  expect_s3_class( "ggplot")
  
  # With removal unwanted variation
  my_estimate |> 
    sccomp_test() |> 
    sccomp_boxplot("type", significance_threshold = 0.025, remove_unwanted_effects = TRUE) |> 
    expect_s3_class( "ggplot")
})

test_that("test constrasts",{

  skip_cmdstan()
  
  new_test =
    my_estimate |>
    sccomp_test() 
  
  # Right and wrong contrasts
  my_estimate |> 
    sccomp_test(contrasts = c("a" = "`(Intercept)`", "b" = "typehealthy")) |> 
    expect_s3_class("tbl")
  
  my_estimate |> 
    sccomp_test(contrasts = "(Intercept)") |> 
    expect_error() |> 
    expect_warning("sccomp says: These components of your contrasts")
  
  my_estimate |> 
    sccomp_test(contrasts = "typehealthy_") |> 
    expect_error() |> 
    expect_warning("These components of your contrasts are not present")
  
  
  res = my_estimate_no_intercept |> 
    sccomp_test(contrasts = c("1/2*typecancer - 1/2*typehealthy", "1/2*typehealthy - 1/2*typecancer") )
  
  
  expect_equal(
    res[1,"c_effect"] |> as.numeric(),
    -res[2,"c_effect"] |> as.numeric()
  )
  
  # Interaction
  my_estimate |> 
    sccomp_test(contrasts = c("(1/2*`continuous_covariate:typehealthy` + 1/2*`continuous_covariate:typehealthy`) -  `continuous_covariate:typehealthy`") ) |> 
    expect_s3_class("tbl")
  
  # Wrong interaction
  my_estimate |> 
    sccomp_test(contrasts = c("(1/2*continuous_covariate:typehealthy + 1/2*`continuous_covariate:typehealthy`) -  `continuous_covariate:typehealthy`") ) |> 
    expect_error() |> 
    expect_warning("sccomp says: for columns which have special characters") 
  

})

# test_that("proportions",{
#   
#   skip_cmdstan()
#   
#   counts_obj |>
#       sccomp_estimate(
#       formula_composition = ~ type, 
#       sample = "sample",  
#       cell_group = "cell_group", 
#       abundance = "proportion",
#       cores = 1,
#       mcmc_seed = 42,
#       max_sampling_iterations = n_iterations, verbose = FALSE
#     ) |> 
#       expect_warning("sccomp says: your proportion values include 0")
#  
#   # counts_obj |>
#   #   sccomp_estimate(
#   #     formula_composition = ~ type ,
#   #     sample = "sample",
#   #     cell_group = "cell_group",
#   #     abundance = "proportion",
#   #     cores = 1,
#   #     mcmc_seed = 42,
#   #     max_sampling_iterations = n_iterations
#   #   ) |>
#   #   expect_warning("sccomp says: your proportion values include 0.*")
#   
# })

test_that("sccomp_proportional_fold_change",{
  
  skip_cmdstan()
  
  
  my_estimate |> 
    sccomp_proportional_fold_change(formula_composition = ~  type, from =  "healthy", to = "cancer") |> 
    expect_no_error()
  
  
  my_estimate_inverse_factor |> 
    sccomp_proportional_fold_change(formula_composition = ~  type, from =  "healthy", to = "cancer") |> 
    pull(proportion_fold_change) |> 
    unique() |> 
    length() |> 
    expect_equal(1)
  
})

test_that("plotting for no significance",{
  
  skip_cmdstan()
  
  
  no_significance_df |>
    mutate(count = count |> as.integer()) |> 
    sccomp_estimate(formula_composition = ~ condition,
                    sample = "sample",
                    cell_group = "cell_group",abundance = "count",
                    bimodal_mean_variability_association = TRUE,      
                    cores = 1, max_sampling_iterations = n_iterations, verbose = FALSE
) |>
    sccomp_test() |> 
    sccomp_boxplot("condition") |> 
    expect_warning() |> 
    expect_no_error()
  
  
  
  
})

# Load necessary libraries
library(testthat)
library(stringr)

# Define helper functions (assuming they are already defined as before)

# Begin unit tests
test_that("contrasts_to_parameter_list handles various contrasts correctly", {
  skip_cmdstan()
  # Define a variety of contrasts
  contrasts <- c(
    # Simple contrast without special characters
    "GroupA - GroupB",
    
    # Contrast with variable names containing spaces (requires backquotes)
    "`Variable with spaces` - `Another variable`",
    
    # Contrast with special characters in variable names
    "`Variable#1` + `Variable-2` - `Variable.3`",
    
    # Contrast with multiplication and division
    "0.5 * Treatment / Control",
    
    # Contrast with fractions
    "(1/2) * (`Var1` + `Var2`) - (`Var3` / 2)",
    
    # Contrast with fractions used in multiplication after '*'
    "(`Variable1` + `Variable2`) * 1/2",
    
    # Contrast with nested backquotes and special characters
    "`age_bin_sex_specificSenior___adipose tissue` + `age_bin_sex_specificSenior:sexmale___adipose tissue`",
    
    # Complex contrast provided earlier
    "((`age_bin_sex_specificSenior___adipose tissue` + `age_bin_sex_specificSenior:sexmale___adipose tissue` + `age_bin_sex_specificMiddle Age___adipose tissue` + `age_bin_sex_specificMiddle Age:sexmale___adipose tissue`) / 2) - ((`age_bin_sex_specificYoung Adulthood___adipose tissue` + `age_bin_sex_specificYoung Adulthood:sexmale___adipose tissue`) / 1)"
  )
  
  # Apply the function to extract parameters
  extracted_params <- contrasts_to_parameter_list(contrasts, drop_back_quotes = FALSE)
  
  # Expected parameters after extraction
  expected_params <- c(
    "GroupA", "GroupB",
    "`Variable with spaces`", "`Another variable`",
    "`Variable#1`", "`Variable-2`", "`Variable.3`",
    "Treatment", "Control",
    "`Var1`", "`Var2`", "`Var3`",
    "`Variable1`", "`Variable2`",
    "`age_bin_sex_specificSenior___adipose tissue`", "`age_bin_sex_specificSenior:sexmale___adipose tissue`",
    "`age_bin_sex_specificMiddle Age___adipose tissue`", "`age_bin_sex_specificMiddle Age:sexmale___adipose tissue`",
    "`age_bin_sex_specificYoung Adulthood___adipose tissue`", "`age_bin_sex_specificYoung Adulthood:sexmale___adipose tissue`"
  )

  # Remove duplicates from expected_params
  expected_params <- unique(expected_params)

  # Check if the extracted parameters match the expected parameters
  expect_equal(sort(extracted_params), sort(expected_params))

  # Now, check if backquotes are correctly applied where necessary
  contrasts_elements <- extracted_params

  # Function to check for valid column characters
  contains_only_valid_chars_for_column <- function(string){
    return(all(str_detect(string, "^[A-Za-z\\.][A-Za-z0-9\\._]*$")))
  }

  # # Check if backquotes are required
  # require_back_quotes <- !str_remove_all(contrasts_elements, "`") |> contains_only_valid_chars_for_column()
  #
  # # Check if backquotes are correctly placed
  # has_left_back_quotes <- str_detect(contrasts_elements, "^`")
  # has_right_back_quotes <- str_detect(contrasts_elements, "`$")
  #
  # # Identify elements that require backquotes but don't have them correctly
  # if_true_not_good <- require_back_quotes & !(has_left_back_quotes & has_right_back_quotes)
  #
  # # Expect that all elements either don't require backquotes or have them correctly placed
  # expect_true(all(!if_true_not_good))

})

test_that("sample ID malformed", {
  
  skip_cmdstan()

  counts_obj |>
  mutate(sample = if_else(sample %in% c("SCP424_pbmc1", "SCP424_pbmc2", "SCP345_860"), "SCP424_pbmc1", sample)) |>
    mutate(count = count |> as.integer()) |> 
  sccomp_estimate(
    formula_composition = ~ type ,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1, verbose = FALSE
  ) |>
  expect_warning("sccomp says: the input data frame does not have the same number") |>
    expect_error("sccomp says: You have duplicated")


  # Generate counts from negative binomial
  # From https://github.com/MangiolaLaboratory/sccomp/issues/178
  n_samples <- 100
  n_taxa <- 5
  dispersion = 0.3
  base_means <- rexp(n_taxa, rate = 1/100)
  counts <- matrix(
    rnbinom(
      n = n_samples * n_taxa,
      mu = rep(base_means, each = n_samples),
      size = 1 / 0.3
    ),
    nrow = n_samples
  )
  colnames(counts) <- paste0("taxon_", 1:n_taxa)
  
  counts_df <- as.data.frame(counts) %>%
    mutate(
      condition = rep(c("control", "treatment"), each = n_samples/2),
      sample = paste0("s", 1:n())
    ) %>%
    pivot_longer(
      -c("sample", condition),
      names_to = "taxon", values_to = "count")
  
  counts_df |> 
    mutate(count = count |> as.integer()) |> 
    sccomp_estimate(
      formula_composition = ~ condition,
      cell_group = "taxon",
      sample = "sample",
      abundance = "count", verbose = FALSE
    ) |> 
    expect_no_error()
  
  # For ~ 1 formula do not expect error
  model_without_association = 
    seurat_obj |>
    sccomp_estimate( 
      formula_composition = ~ 1, 
      sample = "sample", 
      cell_group = "cell_group", verbose = FALSE
    ) |> 
    expect_no_error()

})

test_that("LOO", {
  skip_cmdstan()
  
  library(loo)
  
  # Fit first model
  model_with_factor_association = 
    seurat_obj |>
    sccomp_estimate( 
      formula_composition = ~ type, 
      sample = "sample", 
      cell_group = "cell_group", 
      inference_method = "hmc",
      cores = 1, 
      enable_loo = TRUE, max_sampling_iterations = n_iterations, verbose = FALSE
    ) |> 
    expect_no_error()
  
  # Fit second model
  model_without_association = 
    seurat_obj |>
    sccomp_estimate( 
      formula_composition = ~ 1, 
      sample = "sample", 
      cell_group = "cell_group", 
      inference_method = "hmc",
      cores = 1, 
      enable_loo = TRUE, max_sampling_iterations = n_iterations, verbose = FALSE
    ) |> 
    expect_no_error()
  
  # Compare models
  loo_compare(
    attr(model_with_factor_association, "fit")$loo(),
    attr(model_without_association, "fit")$loo()
  ) |> 
    suppressWarnings() |> 
    expect_no_error()
  
})

test_that("use two methods", {

  skip_cmdstan()

    seurat_obj |>
    sccomp_estimate( 
      formula_composition = ~ 1, 
      sample = "sample", 
      cell_group = "cell_group", 
      inference_method = "hmc",
      cores = 1, max_sampling_iterations = n_iterations, verbose = FALSE
    ) |> 
    sccomp_remove_outliers(inference_method = "pathfinder", cores = 1, max_sampling_iterations = n_iterations) |> 
    expect_no_error()
  
})

test_that("get_design_matrix_with_na_handling properly handles NA values", {
  skip_cmdstan()
  # Create a test dataset with NA values in categorical variables
  test_data <- tibble(
    sample = paste0("sample_", 1:10),
    type = c("A", "B", "A", "B", "A", NA, NA, "C", "C", "B"),
    numeric_var = c(1.2, 0.8, 1.5, 0.6, 1.0, 0.9, 1.1, 1.3, 0.7, 0.5)
  )
  test_data$type <- as.factor(test_data$type)
  
  # Create a test dataset with NA at the beginning
  test_data_na_first <- tibble(
    sample = paste0("sample_", 1:10),
    type = c(NA, NA, "A", "B", "A", "B", "A", "C", "C", "B"),
    numeric_var = c(0.9, 1.1, 1.2, 0.8, 1.5, 0.6, 1.0, 1.3, 0.7, 0.5)
  )
  test_data_na_first$type <- as.factor(test_data_na_first$type)
  
  # Create the formula
  formula <- ~ type + numeric_var
  
  # Call the function for both datasets
  design_matrix <- sccomp:::get_design_matrix_with_na_handling(test_data, formula, sample)
  design_matrix_na_first <- sccomp:::get_design_matrix_with_na_handling(test_data_na_first, formula, sample)
  
  # Test colnames of design_matrix_na_first
  expect_equal(colnames(design_matrix_na_first), colnames(design_matrix))
  
  # Test dimensions
  expect_equal(nrow(design_matrix), 10)
  
  # Test that NA column is removed
  expect_false(any(grepl("typeNA", colnames(design_matrix))))
  
  # Test that rows with NA values get equal distribution across factor levels
  # Get row indices with NA in original data
  na_rows <- which(is.na(test_data$type))
  
  # Get type columns in design matrix
  type_cols <- grep("^type", colnames(design_matrix))
  
  # Number of visible levels in design matrix
  n_visible_levels <- length(type_cols)
  
  # Expected value accounting for reference level: 1/(n_visible_levels+1)
  # This matches current implementation which considers reference level
  expected_value <- 1/(n_visible_levels + 1)
  
  # For NA rows, each type column should have the expected value
  for (row in na_rows) {
    for (col in type_cols) {
      expect_equal(design_matrix[row, col], expected_value, 
                   tolerance = 1e-6,
                   info = paste("Row", row, "Column", col))
    }
  }
  
  # Test multiclass factor with 3 levels
  test_data2 <- tibble(
    sample = paste0("sample_", 1:10),
    type3 = c("X", "Y", "Z", "X", "Y", NA, NA, "Z", "X", "Y"),
    numeric_var = c(1.2, 0.8, 1.5, 0.6, 1.0, 0.9, 1.1, 1.3, 0.7, 0.5)
  )
  test_data2$type3 <- as.factor(test_data2$type3)
  
  formula2 <- ~ type3 + numeric_var
  
  design_matrix2 <- get_design_matrix_with_na_handling(test_data2, formula2, sample)
  
  # Test dimensions
  expect_equal(nrow(design_matrix2), 10)
  
  # Test that NA column is removed
  expect_false(any(grepl("type3NA", colnames(design_matrix2))))
  
  # Get NA rows
  na_rows2 <- which(is.na(test_data2$type3))
  
  # Get type columns in design matrix
  type_cols2 <- grep("^type3", colnames(design_matrix2))
  
  # Number of visible levels in design matrix
  n_visible_levels2 <- length(type_cols2)
  
  # For a 3-level factor with one reference level, each visible column should get 1/3
  expected_value2 <- 1/(n_visible_levels2 + 1)
  
  # For NA rows, each type column should have the expected value
  for (row in na_rows2) {
    for (col in type_cols2) {
      expect_equal(design_matrix2[row, col], expected_value2, 
                   tolerance = 1e-6,
                   info = paste("Row", row, "Column", col))
    }
  }
  
  # Test mix of categorical and numeric variables
  test_data3 <- tibble(
    sample = paste0("sample_", 1:10),
    type = c("A", "B", "A", "B", "A", NA, NA, "C", "C", "B"),
    numeric_var = c(1.2, 0.8, 1.5, 0.6, 1.0, NA, 1.1, 1.3, 0.7, 0.5)
  )
  test_data3$type <- as.factor(test_data3$type)
  
  formula3 <- ~ type + numeric_var
  
  # This should handle both NA in factors and numerics
  expect_error(design_matrix3 <- get_design_matrix_with_na_handling(test_data3, formula3, sample), NA)
})

test_that("renamed columns in Seurat input", {
  skip_cmdstan()
  
  # Create a copy of seurat_obj with renamed columns
  renamed_seurat = seurat_obj
  colnames(renamed_seurat[[]])[colnames(renamed_seurat[[]]) == "sample"] = "a"
  colnames(renamed_seurat[[]])[colnames(renamed_seurat[[]]) == "cell_group"] = "b"
  
  # Test that sccomp_estimate works with renamed columns
  renamed_estimate = 
    renamed_seurat |>
    sccomp_estimate(
      formula_composition = ~ type,
      formula_variability = ~ 1,
      sample = "a",
      cell_group = "b",
      cores = 1,
      mcmc_seed = 42,
      max_sampling_iterations = n_iterations,
      verbose = FALSE
    ) |>
    expect_no_error()
  
  # Verify the results contain the expected columns and data
  expect_true("b" %in% (renamed_estimate |> attr(".cell_group") |> as.character()))
  expect_true(all(c("a", "b") %in% colnames(renamed_seurat[[]])))
  
  # Test that we can still use the renamed columns for predictions
  renamed_estimate |>
    sccomp_predict(
      formula_composition = ~ type,
      new_data = renamed_seurat[[]] |> distinct(a, type),
      number_of_draws = 1
    ) |>
    expect_no_error()
})


