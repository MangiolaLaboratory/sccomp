library(sccomp)
library(tidyverse)
library(magrittr)
library(tidybulk)



# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = args[2]
add_outliers = as.integer(args[3])
add_hyperpriors = as.integer(args[4])

print(add_hyperpriors)

# exposures = counts_obj %>% group_by(sample) %>% summarise(s=sum(count)) %>% pull(s) %>% sort %>% head(-1)
beta_0 = readRDS("dev/beta_0.rds")

if(add_hyperpriors==0){

  # uninformative
  prior_mean_variable_association = list(intercept = c(0, 2), slope = c(0,  1), standard_deviation = c(20, 40))

} else if(add_hyperpriors==1) {
  # PREVIOUS
  prior_mean_variable_association = list(intercept = c(4.9205842, 0.11952838), slope = c(-0.7608308,  0.09501801), standard_deviation = c(37.45106,76.65637))

  # s41587-020-0602-4_COVID_19
  # prior_mean_variable_association = list(intercept = c(3.8357711, 1.05), slope = c(-0.9338016,  0.10), standard_deviation = c(5.390190,8.746909))

} else if(add_hyperpriors==2) {
  # BRCA1_s41467-021-21783-3
  prior_mean_variable_association = list(intercept = c(5.8281775, 0.1457238), slope = c(-0.8904407,  0.1061705), standard_deviation = c(53.71133,66.01899))
} else if(add_hyperpriors==3) {
  # Wrong hyperprior
  prior_mean_variable_association = list(intercept = c(10, 0.15), slope = c(1,  0.10), standard_deviation = c(37.45106,76.65637))
}

print(prior_mean_variable_association)

# Iterate over runs
readRDS(input_file) %>%
  mutate( results_sccomp = map2(
    data, run,
    ~  .x %>%
      mutate(generated_counts = as.integer(generated_counts)) %>%
      sccomp_glm(
        formula_composition = ~type, .sample = sample, .cell_group = cell_type, .count = generated_counts,
        check_outliers = add_outliers==1,
        approximate_posterior_inference = "none",
        percent_false_positive = 0.001 * 100,
        mcmc_seed = .y * 2,
        prior_mean_variable_association = prior_mean_variable_association,
        cores = 4
      )
  )) %>%
  saveRDS(output_file)


