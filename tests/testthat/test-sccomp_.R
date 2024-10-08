library(dplyr)
library(tidyr)
library(sccomp)
data("seurat_obj")
data("sce_obj")
data("counts_obj")

set.seed(42)

if (instantiate::stan_cmdstan_exists()){
  
  my_estimate = 
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ continuous_covariate * type ,
      formula_variability = ~ 1,
      sample, cell_group,
      
      cores = 1, 
      inference_method = "pathfinder",
      # mcmc_seed = 42,
      max_sampling_iterations = 1000, verbose=FALSE
    )
  
  my_estimate_full_bayes = 
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ continuous_covariate * type ,
      formula_variability = ~ 1,
      sample, cell_group,
      cores = 1, 
      inference_method = "hmc",
      mcmc_seed = 42,
      max_sampling_iterations = 1000, verbose=FALSE
    )
  
  my_estimate_no_intercept = 
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ 0 + type + continuous_covariate,
      formula_variability = ~ 1,
      sample, cell_group,
      cores = 1,
      mcmc_seed = 42,
      max_sampling_iterations = 1000, verbose=FALSE
    )
  
  my_estimate_random = 
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ type + (type | group__),
      
      formula_variability = ~ 1,
      sample, cell_group,
      cores = 1,
      mcmc_seed = 42,     
      max_sampling_iterations = 1000, verbose=FALSE
    )
  
  # my_estimate_random2 = 
  # 	seurat_obj |>
  # 	sccomp_estimate(
  # 		formula_composition = ~ 1 +  type + (1 + type | group__),
  # 		formula_variability = ~ 1,
  # 		sample, cell_group,
  # 		cores = 1,
  # 		mcmc_seed = 42,     
  # 		max_sampling_iterations = 1000
  # 	)
  # 
  # my_estimate_random3 = 
  # 	seurat_obj |>
  # 	sccomp_estimate(
  # 		formula_composition = ~  type + (1 | group__),
  # 		formula_variability = ~ 1,
  # 		sample, cell_group,
  # 		cores = 1,
  # 		mcmc_seed = 42,     
  # 		max_sampling_iterations = 1000
  # 	)
  
  
  
}


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
   
  new_data_tibble = new_data_seurat[[]] |> distinct(sample, type, continuous_covariate)
  
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
      new_data = new_data_seurat
    ) |>
    nrow() |>
    expect_equal(60)
  
})

test_that("outliers",{
  
  skip_cmdstan()
  
  my_estimate |>
    sccomp_remove_outliers(
      cores = 1, 
      max_sampling_iterations = 1000,
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
      sample, cell_group,
      cores = 1,
      mcmc_seed = 42,     
      max_sampling_iterations = 1000, verbose=FALSE
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
  #debugonce(sccomp:::fit_model)
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
      sample, cell_group,
      cores = 1,
      mcmc_seed = 42,     
      max_sampling_iterations = 1000
    )
  
  
})


test_that("multilevel multi beta binomial from Seurat with intercept and continuous covariate",{

  skip_cmdstan()
  
  library(sccomp)

  res =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ continuous_covariate + (1 + continuous_covariate | group__),
      formula_variability = ~ 1,
      sample, cell_group,
      cores = 1,
      mcmc_seed = 42,   
      max_sampling_iterations = 1000,
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
#           contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
#           cores = 1,
#           mcmc_seed = 42,       max_sampling_iterations = 1000
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

test_that("remove unwanted variation",{

  skip_cmdstan()
  
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
      cores = 1,
      mcmc_seed = 43,    
      max_sampling_iterations = 1000, verbose = FALSE
    )

  estimate |>
    sccomp_remove_unwanted_variation(formula_composition = ~ type)

})

test_that("multi beta binomial from SCE",{

  skip_cmdstan()
  
    res =
      sce_obj |>
    sccomp_estimate(
      formula_composition = ~ type,
      formula_variability = ~ 1,
      sample,
      cell_group,
      cores = 1,
      mcmc_seed = 42,      
      max_sampling_iterations = 1000, verbose = FALSE
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
    sample,
    cell_group,
    cores = 1,
    mcmc_seed = 42,   
    max_sampling_iterations = 1000, verbose = FALSE
  )

res_composition_variability =
  seurat_obj[[]] |>
  sccomp_estimate(
    formula_composition = ~ type,
    formula_variability = ~ type,
    sample,
    cell_group,
    cores = 1,
    mcmc_seed = 42,    
    max_sampling_iterations = 1000, verbose = FALSE
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
})

# Test for plot_1d_intervals function
test_that("plot_1d_intervals function works correctly", {
  
  skip_cmdstan()

    my_estimate |> 
    sccomp_test() |> 
    plot_1D_intervals(
    significance_threshold = 0.025
  ) |> 
  expect_s3_class("patchwork")
})

# Test for plot_2d_intervals function
test_that("plot_2d_intervals function works correctly", {
  
  skip_cmdstan()
  
  my_estimate |> 
    sccomp_test() |> 
    plot_2D_intervals(
      significance_threshold = 0.025
    ) |>
  expect_s3_class("ggplot")
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
    expect_warning("sccomp says: for columns which have special characters") 
  

})



test_that("sccomp_proportional_fold_change",{
  
  skip_cmdstan()
  
  
  my_estimate |> 
    sccomp_proportional_fold_change(formula_composition = ~  type, from =  "healthy", to = "cancer") |> 
    expect_no_error()
  
 
  
  
})



# fit = 
# 	seurat_obj |> 
# 	dplyr::count(sample, cell_group, continuous_covariate, type, group__) |> 
# 	with_groups(sample, ~ .x |> mutate(tot = sum(n))) |> 
# 	rename(group = group__) |> 
# 	filter(cell_group=="B mem") |> 
# 	brm(
# 		n | trials(tot) ~ type + (1 + type| group), 
# 		data = _, 
# 		family=beta_binomial(), 
# 		cores = 4
# 	)
