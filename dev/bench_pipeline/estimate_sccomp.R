library(sccomp)
library(tidyverse)
library(magrittr)
library(tidybulk)
library(tictoc)


# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = args[2]
add_outliers = as.integer(args[3])


# exposures = counts_obj %>% group_by(sample) %>% summarise(s=sum(count)) %>% pull(s) %>% sort %>% head(-1)
# beta_0 = readRDS("~/PostDoc/sccomp_dev/dev/beta_0.rds")


# Iterate over runs
readRDS(input_file) %>%
  mutate( results_sccomp = map2(
    data, run,
    ~ {

      tic()

      results =
        .x %>%
        mutate(generated_counts = as.integer(generated_counts)) %>%
        sccomp_glm(formula_composition = ~type, .sample = sample, .cell_group = cell_type, .count = generated_counts,
          check_outliers = add_outliers==1,
          approximate_posterior_inference = "none",
          percent_false_positive = 0.001 * 100,
          mcmc_seed = .y * 2,
          prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(5.06983, 8.549324))
        )

      time_arr = toc()

      results |> tidybulk:::add_attr(time_arr$toc - time_arr$tic, "execution_time_seconds")
    }
  )) %>%
  saveRDS(output_file)


