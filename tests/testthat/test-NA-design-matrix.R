n_iterations = 1000

test_that("test NA factors with interactions",{
  skip_cmdstan()
  
  library(dplyr)

  formula = ~ group2__ * type
  
  # Test that sccomp_estimate and sccomp_predict work with NA values and expect the ESS warning
    res =
      seurat_obj |>
      sccomp_estimate(
        formula_composition = formula,
        formula_variability = ~ 1,
        "sample", "cell_group",
        cores = 1,
        mcmc_seed = 42,   
        max_sampling_iterations = n_iterations,
        inference_method = "hmc", verbose=FALSE
      )
    
    # Create new data using one sample and setting continuous covariate to NA
    new_data = 
      res |> 
      attr("count_data") |> 
      sccomp:::.subset(sample)  |>
      filter(sample == "SI-GA-H1") |>
      mutate(type = NA)
    
    result = sccomp_predict(
      res,
      new_data = bind_rows(
        res |> attr("count_data") |> sccomp:::.subset(sample) |> filter(sample != "SI-GA-H1"),
        new_data
      ),
      formula_composition = formula
    )
  
  # Check that the result is not NULL
  expect_false(is.null(result))
})

test_that("test NA in both variables",{
  skip_cmdstan()
  
  library(dplyr)
  
  formula = ~ group2__ * type
  
  # Test that sccomp_estimate and sccomp_predict work with NA values and expect the ESS warning
    res =
      seurat_obj |>
      sccomp_estimate(
        formula_composition = formula,
        formula_variability = ~ 1,
        "sample", "cell_group",
        cores = 1,
        mcmc_seed = 42,   
        max_sampling_iterations = n_iterations,
        inference_method = "hmc", verbose=FALSE
      )
    
    # Create new data with NA in both variables
    new_data = 
      res |> 
      attr("count_data") |> 
      sccomp:::.subset(sample)  |>
      filter(sample == "SI-GA-H1") |>
      mutate(
        sample = "NEW_SAMPLE_1",
        type = NA,
        group2__ = NA
      )
    
    result = sccomp_predict(
      res,
      new_data = bind_rows(
        res |> attr("count_data") |> sccomp:::.subset(sample) |> filter(sample != "SI-GA-H1"),
        new_data
      ),
      formula_composition = formula
    )
  
  # Check that the result is not NULL
  expect_false(is.null(result))
})

test_that("test NA in three-way interactions",{
  skip_cmdstan()
  
  library(dplyr)
  
  formula = ~ group2__ * type * continuous_covariate
  
  # Test that sccomp_estimate and sccomp_predict work with NA values and expect the ESS warning
    res =
      seurat_obj |>
      sccomp_estimate(
        formula_composition = formula,
        formula_variability = ~ 1,
        "sample", "cell_group",
        cores = 1,
        mcmc_seed = 42,   
        max_sampling_iterations = n_iterations,
        inference_method = "hmc", verbose=FALSE
      )
    
    # Create new data using one sample and setting continuous covariate to NA
    new_data = 
      res |> 
      attr("count_data") |> 
      sccomp:::.subset(sample)  |>
      filter(sample == "SI-GA-H1") |>
      mutate(
        sample = "NEW_SAMPLE_2",  # Use unique sample name
        type = NA,
        continuous_covariate = NA
      )
    
    result = sccomp_predict(
      res,
      new_data = bind_rows(
        res |> attr("count_data") |> sccomp:::.subset(sample) |> filter(sample != "SI-GA-H1"),
        new_data
      ),
      formula_composition = formula
    )

  
  # Check that the result is not NULL
  expect_false(is.null(result))
})
