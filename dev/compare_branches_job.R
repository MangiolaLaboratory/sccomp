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

# Install required packages if not already installed
if (!requireNamespace("job", quietly = TRUE)) {
  install.packages("job")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Load required libraries
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
  
  # Install required packages for the function
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  if (!requireNamespace("posterior", quietly = TRUE)) install.packages("posterior")
  
  # Remove existing sccomp package to ensure clean installation
  remove.packages("sccomp")  
  
  library(devtools)
  message("Installing branch ", branch)
  # Install the specified branch from GitHub
  devtools::install_github(repo, ref = branch, force = TRUE)
  
  library(sccomp)
  library(posterior)
  
  # Clear Stan model cache to ensure fresh model compilation
  clear_stan_model_cache()
  system("rm -r ~/.sccomp_models")
  
  # Check available columns
  
  
  
  # Run sccomp_estimate with actual test data
  start_time <- Sys.time()
  cat("Starting sccomp_estimate...\n")
  
  # debugonce(sccomp:::fit_model)
  
  
  
  # Fit the compositional model using sccomp_estimate
  # This model includes:
  # - Fixed effects: sex, age_days_scaled, self_reported_ethnicity
  # - Interaction: sex*age_days_scaled
  # - Random effects: (1 | dataset_id) for dataset-level variation
  # - Random effects: (sex*age_days_scaled | tissue_groups) for tissue-specific effects
  estimate <- sccomp_estimate(
    my_df,
    formula_composition = ~
      sex*age_days_scaled + self_reported_ethnicity +
      (1 | dataset_id) +
      (sex*age_days_scaled | tissue_groups),
    sample = "sample_id",
    abundance = "n",
    cell_group = "cell_type_unified_ensemble",
    inference_method = "hmc",  # Hamiltonian Monte Carlo
    mcmc_seed = 123,  # Set seed for reproducibility
    verbose = TRUE, 
    refresh = 5,  # Print progress every 5 iterations
    cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))  # Use SLURM cores if available
    
  )
  end_time = Sys.time()
  # estimate |> attr("fit") = NULL
  
  # return the estimate object and time taken
  list(estimate = estimate, time = end_time - start_time)
  
} 

# Load cellNexus data and prepare for analysis
# Filter for blood tissue samples with 10x technology, normal disease state
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

# Select first 300 samples for analysis
my_samples = full_df |> distinct(sample_id) |> head(300) |> pull(1)

# Prepare the data frame for sccomp analysis
# - Count cells by sample and cell type
# - Include metadata: sex, age, ethnicity, dataset, tissue
# - Scale age using square root transformation and z-standardization
my_df =
  full_df |>
  filter(sample_id %in% my_samples) |>
  count(sample_id, cell_type_unified_ensemble, sex, age_days, self_reported_ethnicity, dataset_id, tissue_groups) |>
  mutate(age_days_scaled = (age_days |> sqrt() - mean(age_days |> sqrt())) / sd(age_days |> sqrt())) |>
  mutate(n = as.integer(n)) |>
  collect() |>
  drop_na() |>
  droplevels()

# Launch background job for master branch
job::job({ result_master = run_branch(my_df, "master")  })

# Launch background job for sum-to-zero-variable branch
job::job({ result_sum_to_zero = run_branch(my_df, "sum-to-zero-variable")  })

# Load ggplot2 for visualization
library(ggplot2)

# Compare ESS (Effective Sample Size) between branches
# Create scatter plot of ESS values from both branches
result_master$estimate |> 
  select(cell_type_unified_ensemble, parameter, factor, c_ess_bulk_master = c_ess_bulk, c_rhat_master = c_rhat) |> 
  left_join(
    result_sum_to_zero$estimate |> select(cell_type_unified_ensemble, parameter, factor, c_ess_bulk_sum_to_zero = c_ess_bulk, c_rhat_sum_to_zero = c_rhat)
  ) |> 
  ggplot(aes(c_ess_bulk_master, c_ess_bulk_sum_to_zero)) +
  geom_point(aes(color = cell_type_unified_ensemble)) +
  geom_abline(intercept = 0, slope = 1) +  # Perfect agreement line
  geom_smooth(method = "lm", color = "red")  # Linear trend line

# Create detailed comparison table with ESS, R-hat, and effect estimates
# Sort by master branch ESS to identify parameters with different convergence
result_master$estimate |> 
  select(cell_type_unified_ensemble, parameter, factor, c_ess_bulk_master = c_ess_bulk, c_rhat_master = c_rhat, c_effect_master = c_effect) |> 
  left_join(
    result_sum_to_zero$estimate |> select(cell_type_unified_ensemble, parameter, factor, c_ess_bulk_sum_to_zero = c_ess_bulk, c_rhat_sum_to_zero = c_rhat, c_effect_sum_to_zero = c_effect)
  )|> arrange(c_ess_bulk_master) |> View()
