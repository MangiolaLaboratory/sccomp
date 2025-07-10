# compare_branches_job.R
#
# Launch two RStudio background jobs – one for each branch – that install the
# selected sccomp branch, run a lightweight model fit on the built-in
# `counts_obj` test data, and save the results to disk.  When both jobs are
# done you can load the two *.rds files and compare ESS / estimates.
#
# ‑- How to use -------------------------------------------------------------
# 1. Adjust `branches` or `repo` if needed.
# 2. Source this script (or run it line-by-line).  The jobs will appear in the
#    RStudio Jobs pane and run in the background so your console stays free.
# 3. When the jobs finish, run `compare_results()` (defined at the bottom) to
#    load and visualise the results.
# --------------------------------------------------------------------------

if (!requireNamespace("job", quietly = TRUE)) {
  install.packages("job")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

library(sccomp)
library(cellNexus)
library(dplyr)
library(stringr)
library(tidyr)
library(job)


# --------------------------------------------------------------------------
# Run master branch only ---------------------------------------------------
# --------------------------------------------------------------------------
run_branch <- function(my_df, branch = "master", 
                       repo = "stemangiola/sccomp",
                       model_formula = ~ sex + age_days_scaled,
                       output_samples = 400) {
  
  
  
  
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  if (!requireNamespace("posterior", quietly = TRUE)) install.packages("posterior")
  
  remove.packages("sccomp")  
  
  library(devtools)
  message("Installing branch ", branch)
  devtools::install_github(repo, ref = branch, force = TRUE)
  
  library(sccomp)
  library(posterior)
  
  clear_stan_model_cache()
  system("rm -r ~/.sccomp_models")
  
  # Check available columns
  
  
  
  # Run sccomp_estimate with actual test data
  start_time <- Sys.time()
  cat("Starting sccomp_estimate...\n")
  
  # debugonce(sccomp:::fit_model)
  
  
  
  estimate <- sccomp_estimate(
    my_df,
    formula_composition = ~
      sex*age_days_scaled + self_reported_ethnicity +
      (1 | dataset_id) +
      (sex*age_days_scaled | tissue_groups),
    sample = "sample_id",
    abundance = "n",
    cell_group = "cell_type_unified_ensemble",
    inference_method = "hmc",
    mcmc_seed = 123,
    verbose = TRUE, 
    refresh = 5, 
    cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
    
  )
  end_time = Sys.time()
  # estimate |> attr("fit") = NULL
  
  # return the estimate object and time taken
  list(estimate = estimate, time = end_time - start_time)
  
} 

# Load cellNexus data
full_df =
  cellNexus::get_metadata() |>
  filter(
    tissue_groups == "blood",
    str_detect(assay, "10x"),
    disease == "normal",
    !empty_droplet,
    alive
    # , scDblFinder.class == "singlet" # Uncomment if this column exists
  )
my_samples = full_df |> distinct(sample_id) |> head(300) |> pull(1)
my_df =
  full_df |>
  filter(sample_id %in% my_samples) |>
  count(sample_id, cell_type_unified_ensemble, sex, age_days, self_reported_ethnicity, dataset_id, tissue_groups) |>
  mutate(age_days_scaled = (age_days |> sqrt() - mean(age_days |> sqrt())) / sd(age_days |> sqrt())) |>
  mutate(n = as.integer(n)) |>
  collect() |>
  drop_na() |>
  droplevels()

job::job({ result_master = run_branch(my_df, "master")  })

job::job({ result_sum_to_zero = run_branch(my_df, "sum-to-zero-variable")  })

library(ggplot2)
result_master$estimate |> 
  select(cell_type_unified_ensemble, parameter, factor, c_ess_bulk_master = c_ess_bulk, c_rhat_master = c_rhat) |> 
  left_join(
    result_sum_to_zero$estimate |> select(cell_type_unified_ensemble, parameter, factor, c_ess_bulk_sum_to_zero = c_ess_bulk, c_rhat_sum_to_zero = c_rhat)
  ) |> 
  ggplot(aes(c_ess_bulk_master, c_ess_bulk_sum_to_zero)) +
  geom_point(aes(color = cell_type_unified_ensemble)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm", color = "red")

result_master$estimate |> 
  select(cell_type_unified_ensemble, parameter, factor, c_ess_bulk_master = c_ess_bulk, c_rhat_master = c_rhat, c_effect_master = c_effect) |> 
  left_join(
    result_sum_to_zero$estimate |> select(cell_type_unified_ensemble, parameter, factor, c_ess_bulk_sum_to_zero = c_ess_bulk, c_rhat_sum_to_zero = c_rhat, c_effect_sum_to_zero = c_effect)
  )|> arrange(c_ess_bulk_master) |> View()
