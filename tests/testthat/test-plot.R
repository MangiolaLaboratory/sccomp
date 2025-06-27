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
      max_sampling_iterations = n_iterations, verbose=FALSE
    )
  
  my_estimate_with_variance = 
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ type,
      formula_variability = ~ type,
      "sample", "cell_group",
      
      cores = 1, 
      inference_method = "pathfinder",
      max_sampling_iterations = n_iterations, verbose=FALSE
    )
}

# Test for plot_1d_intervals function
test_that("plot_1d_intervals function works correctly", {
  # skip_cmdstan()
  
  my_estimate |> 
    sccomp_test() |> 
    plot_1D_intervals(
      significance_threshold = 0.025
    ) |> 
    expect_s3_class("patchwork")
})

# Test for plot_2d_intervals function
test_that("plot_2d_intervals function works correctly", {
  # skip_cmdstan()
  
  my_estimate_with_variance |> 
    sccomp_test() |> 
    plot_2D_intervals(
      significance_threshold = 0.025
    ) |>
    expect_s3_class("ggplot")
})

# Test for show_fdr_message parameter in plot functions
test_that("show_fdr_message parameter works correctly in plot_1D_intervals", {
  # skip_cmdstan()
  
  # Test with show_fdr_message = TRUE (default)
  plot_with_message <- my_estimate |> 
    sccomp_test() |> 
    plot_1D_intervals(
      significance_threshold = 0.025,
      show_fdr_message = TRUE
    )
  
  expect_s3_class(plot_with_message, "patchwork")
  
  # Test with show_fdr_message = FALSE
  plot_without_message <- my_estimate |> 
    sccomp_test() |> 
    plot_1D_intervals(
      significance_threshold = 0.025,
      show_fdr_message = FALSE
    )
  
  expect_s3_class(plot_without_message, "patchwork")
  
  # Verify that both plots are created successfully (no errors)
  expect_no_error(plot_with_message)
  expect_no_error(plot_without_message)
})

test_that("show_fdr_message parameter works correctly in plot_2D_intervals", {
  # skip_cmdstan()
  
  # Test with show_fdr_message = TRUE (default)
  plot_with_message <- my_estimate_with_variance |> 
    sccomp_test() |> 
    plot_2D_intervals(
      significance_threshold = 0.025,
      show_fdr_message = TRUE
    )
  
  expect_s3_class(plot_with_message, "ggplot")
  
  # Test with show_fdr_message = FALSE
  plot_without_message <- my_estimate_with_variance |> 
    sccomp_test() |> 
    plot_2D_intervals(
      significance_threshold = 0.025,
      show_fdr_message = FALSE
    )
  
  expect_s3_class(plot_without_message, "ggplot")
  
  # Verify that both plots are created successfully (no errors)
  expect_no_error(plot_with_message)
  expect_no_error(plot_without_message)
})

test_that("show_fdr_message parameter accepts logical values", {
  # skip_cmdstan()
  
  # Test with TRUE
  expect_no_error(
    my_estimate |> 
      sccomp_test() |> 
      plot_1D_intervals(show_fdr_message = TRUE)
  )
  
  # Test with FALSE
  expect_no_error(
    my_estimate |> 
      sccomp_test() |> 
      plot_1D_intervals(show_fdr_message = FALSE)
  )
  
  # Test with TRUE for 2D plots
  expect_no_error(
    my_estimate_with_variance |> 
      sccomp_test() |> 
      plot_2D_intervals(show_fdr_message = TRUE)
  )
  
  # Test with FALSE for 2D plots
  expect_no_error(
    my_estimate_with_variance |> 
      sccomp_test() |> 
      plot_2D_intervals(show_fdr_message = FALSE)
  )
})

test_that("plot functions work with different significance thresholds", {
  # skip_cmdstan()
  
  # Test plot_1D_intervals with different thresholds
  expect_no_error(
    my_estimate |> 
      sccomp_test() |> 
      plot_1D_intervals(significance_threshold = 0.01)
  )
  
  expect_no_error(
    my_estimate |> 
      sccomp_test() |> 
      plot_1D_intervals(significance_threshold = 0.1)
  )
  
  # Test plot_2D_intervals with different thresholds
  expect_no_error(
    my_estimate_with_variance |> 
      sccomp_test() |> 
      plot_2D_intervals(significance_threshold = 0.01)
  )
  
  expect_no_error(
    my_estimate_with_variance |> 
      sccomp_test() |> 
      plot_2D_intervals(significance_threshold = 0.1)
  )
})

test_that("significance_statistic argument works for plot_1D_intervals", {
  # skip_cmdstan()
  expect_no_error(
    my_estimate |> 
      sccomp_test() |> 
      plot_1D_intervals(significance_statistic = "FDR")
  )
  expect_no_error(
    my_estimate |> 
      sccomp_test() |> 
      plot_1D_intervals(significance_statistic = "pH0")
  )
})

test_that("significance_statistic argument works for plot_2D_intervals", {
  # skip_cmdstan()
  expect_no_error(
    my_estimate_with_variance |> 
      sccomp_test() |> 
      plot_2D_intervals(significance_statistic = "FDR")
  )
  expect_no_error(
    my_estimate_with_variance |> 
      sccomp_test() |> 
      plot_2D_intervals(significance_statistic = "pH0")
  )
})

test_that("show_fdr_message argument works for plot_1D_intervals and plot_2D_intervals", {
  # skip_cmdstan()
  expect_no_error(
    my_estimate |> 
      sccomp_test() |> 
      plot_1D_intervals(significance_statistic = "FDR", show_fdr_message = TRUE)
  )
  expect_no_error(
    my_estimate |> 
      sccomp_test() |> 
      plot_1D_intervals(significance_statistic = "FDR", show_fdr_message = FALSE)
  )
  expect_no_error(
    my_estimate_with_variance |> 
      sccomp_test() |> 
      plot_2D_intervals(significance_statistic = "FDR", show_fdr_message = TRUE)
  )
  expect_no_error(
    my_estimate_with_variance |> 
      sccomp_test() |> 
      plot_2D_intervals(significance_statistic = "FDR", show_fdr_message = FALSE)
  )
})

test_that("significance_statistic and show_fdr_message work via plot() S3 method", {
  # skip_cmdstan()
  expect_no_error(
    plot(
      my_estimate |> sccomp_test(),
      significance_statistic = "FDR",
      show_fdr_message = TRUE
    )
  )
  expect_no_error(
    plot(
      my_estimate |> sccomp_test(),
      significance_statistic = "pH0",
      show_fdr_message = TRUE
    )
  )
  expect_no_error(
    plot(
      my_estimate |> sccomp_test(),
      significance_statistic = "FDR",
      show_fdr_message = FALSE
    )
  )
  expect_no_error(
    plot(
      my_estimate |> sccomp_test(),
      significance_statistic = "pH0",
      show_fdr_message = FALSE
    )
  )
}) 