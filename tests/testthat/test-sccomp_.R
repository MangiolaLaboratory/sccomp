library(dplyr)
library(sccomp)
data("seurat_obj")
data("sce_obj")
data("counts_obj")


# test_that("dirichlet multinomial outliers",{
#
#   library(dplyr)
#   library(sccomp)
#   library(digest)
#
#   res =
#     sccomp::counts_obj  |>
#     sccomp_estimate(
#       formula = ~ type,
#       sample, cell_type, count
#     )
#
#   expect_equal(
#     res |>
#       distinct(cell_type, significant) |>
#       pull(significant) |>
#       digest(algo="md5"),
#     "f6cef772af43198f586e15c96b2f1239"
#   )
#
# })
# test_that("counts dirichlet multinomial outlier VB",{
#
#   if(interactive()){
#     res =
#       counts_obj  |>
#       sccomp_estimate(
#         formula = ~ type,
#         sample, cell_group, count,
#         noise_model = "dirichlet_multinomial",
#         cores = 1
#       )
#   }
#
# })

test_that("Generate data",{


  seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ type ,
      formula_variability = ~ 1,
      sample, cell_group,
      approximate_posterior_inference = FALSE,
      cores = 1,
      mcmc_seed = 42,
      max_sampling_iterations = 1000
    ) |>

    sccomp_replicate() |>
    nrow() |>
    expect_equal(600)

  # With grouping
  seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ 0 + type + (type | group__),
      formula_variability = ~ 1,
      sample, cell_group,
      approximate_posterior_inference = FALSE,
      cores = 1,
      mcmc_seed = 42,     
      max_sampling_iterations = 1000
    ) |>


    sccomp_replicate(~ 0 + type) |>
    nrow() |>
    expect_equal(600)


})

test_that("outliers",{
  
  res =
    seurat_obj |>
    ## filter(cell_group %in% c("NK cycling", "B immature")) |>
    sccomp_estimate(
      formula_composition = ~ type,
      formula_variability = ~ 1,
      sample, cell_group,
      approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42, 
      verbose = TRUE, 
      max_sampling_iterations = 1000
    )
  
})

test_that("multilevel multi beta binomial from Seurat",{

  res =
    seurat_obj |>
    ## filter(cell_group %in% c("NK cycling", "B immature")) |>
    sccomp_estimate(
      formula_composition = ~ type + (1 | group__),
      formula_variability = ~ 1,
      sample, cell_group,
      #approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42,     
      max_sampling_iterations = 1000
    )

  # Check order
  res |>
    filter(parameter == "typehealthy") |>
    arrange(desc(abs(c_effect))) |>
    slice(1:3) |>
    pull(cell_group) |>
    sort() |> 
    expect_equal(c( "CD4 cm high cytokine", "CD4 ribosome" , "Mono NKG7 2"  ) |> sort())

  # Check convergence
  res |>
    filter(c_R_k_hat > 4) |>
    nrow() |>
    expect_equal(0)

  res =
    seurat_obj |>
    ## filter(cell_group %in% c("NK cycling", "B immature")) |>
    sccomp_estimate(
      formula_composition = ~ 0 + type + (type | group__),
      formula_variability = ~ 1,
      sample, cell_group,
      #approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42,    
      max_sampling_iterations = 1000
    )

  # res |>
  #   filter(parameter == "typecancer - typehealthy") |>
  #   arrange(desc(abs(c_effect))) |>
  #   slice(1:3) |>
  #   pull(cell_group) |>
  #   sort() |>
  #   expect_equal(c("B mem"  ,  "CD4 cm high cytokine" ,"CD4 ribosome"         ))

  # Check convergence
  res |>
    filter(c_R_k_hat > 4) |>
    nrow() |>
    expect_equal(0)

})

test_that("multilevel multi beta binomial from Seurat with intercept and continuous covariate",{


  res =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ continuous_covariate + (1 + continuous_covariate | group__),
      formula_variability = ~ 1,
      sample, cell_group,
      #approximate_posterior_inference = "all",
      cores = 1,
      #mcmc_seed = 42,   
      max_sampling_iterations = 1000
    )

    expect_in(
      "T gd2",
      res |>
        filter(parameter == "continuous_covariate") |>
        arrange(desc(abs(c_effect))) |>
        slice(1:3) |>
        pull(cell_group)
    )


  # Check convergence
  res |>
    filter(c_R_k_hat > 4) |>
    nrow() |>
    expect_equal(0)


})



# test_that("wrongly-set groups",{
#
#   # library(tidyseurat)
#   # seurat_obj =
#   #   seurat_obj |>
#   #   nest(data = -c(sample, type)) |>
#   #   mutate(group__wrong = c(1,1,1,1,1, 2,2,2,2,2, 1,1,1,1,1, 2,2,2,2,2) |> as.character()) |>
#   #   unnest(data)
#
#     expect_error(
#       object =
#         seurat_obj |>
#         ## filter(cell_group %in% c("NK cycling", "B immature")) |>
#         sccomp_estimate(
#           formula_composition = ~ 0 + type + (type | group__wrong),
#           formula_variability = ~ 1,
#           sample, cell_group,
#           approximate_posterior_inference = "all",
#           contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
#           cores = 1,
#           mcmc_seed = 42,       max_sampling_iterations = 1000
#         ) ,
#       regexp = "should not be shared"
#     )
#
# })

test_that("multi beta binomial from Seurat",{

  res =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~  type,
      formula_variability = ~ 1,
      sample, cell_group,
      approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42,    
      max_sampling_iterations = 1000
    )

  res |>
    filter(parameter == "typehealthy") |>
    arrange(desc(abs(c_effect))) |>
    slice(1) |>
    pull(cell_group) |>
    expect_equal(c("B mem"  ))

  # Check convergence
  res |>
    filter(c_R_k_hat > 4) |>
    nrow() |>
    expect_equal(0)

})


test_that("multi beta binomial contrasts from Seurat",{

  res = seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ 0 + type,
      formula_variability = ~ 1,
      sample, cell_group,
      approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42,     
      max_sampling_iterations = 1000
    ) |> 
    sccomp_test(contrasts = c("typecancer - typehealthy", "typehealthy - typecancer") )

  
  expect_equal(
    res[1,"c_effect"] |> as.numeric(),
    -res[2,"c_effect"] |> as.numeric()
  )

})

test_that("remove unwanted variation",{

  library(tidyseurat)

  data =
    seurat_obj |>

    # Add batch
    nest(data = -c(sample, type)) |>
    mutate(batch = rep(c(0,1), 10)) |>
    unnest(data)

  # Estimate
  estimate =
    data |>
    sccomp_estimate(
      formula_composition = ~ type + batch,
      formula_variability = ~ 1,
      sample, cell_group,
      approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42,    
      max_sampling_iterations = 1000
    )

  estimate |>
    sccomp_remove_unwanted_variation(~ type)

})

test_that("multi beta binomial from SCE",{

    res =
      sce_obj |>
    sccomp_estimate(
      formula_composition = ~ type,
      formula_variability = ~ 1,
      sample,
      cell_group,
      approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42,      
      max_sampling_iterations = 1000
    )

  res |>
    filter(parameter == "typehealthy") |>
    arrange(desc(abs(c_effect))) |>
    slice(1) |>
    pull(cell_group) |>
    expect_equal(c("B mem"  ))

  # Check convergence
  res |>
    filter(c_R_k_hat > 4) |>
    nrow() |>
    expect_equal(0)
})

res_composition =
  seurat_obj[[]] |>
  sccomp_estimate(
    formula_composition = ~ type,
    formula_variability = ~ 1,
    sample,
    cell_group,
    approximate_posterior_inference = "all",
    cores = 1,
    mcmc_seed = 42,   
    max_sampling_iterations = 1000
  )

res_composition_variability =
  seurat_obj[[]] |>
  sccomp_estimate(
    formula_composition = ~ type,
    formula_variability = ~ type,
    sample,
    cell_group,
    approximate_posterior_inference = "all",
    cores = 1,
    mcmc_seed = 42,    
    max_sampling_iterations = 1000
  )

test_that("multi beta binomial from metadata",{

  res_composition  |>
    filter(parameter == "typehealthy") |>
    arrange(desc(abs(c_effect))) |>
    slice(1) |>
    pull(cell_group) |>
    expect_equal(c("B mem"  ))

  # Check convergence
  res_composition |>
    filter(c_R_k_hat > 4) |>
    nrow() |>
    expect_equal(0)

})

test_that("plot test composition",{

  res_composition |> 
    sccomp_test() |> 
    plot_summary()


})

test_that("plot test variability",{

  res_composition_variability |> 
    sccomp_test() |> 
    plot_summary()


})

test_that("test constrasts",{


  estimate =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ type ,
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = FALSE,
      cores = 1,
      mcmc_seed = 42,       max_sampling_iterations = 1000
    )

  new_test =
    estimate |>
    sccomp_test() 

})





