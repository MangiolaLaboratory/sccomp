library(sccomp)
library(tidyverse)
library(magrittr)
library(tidybulk)



# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = args[2]
add_outliers = as.integer(args[3])


# exposures = counts_obj %>% group_by(sample) %>% summarise(s=sum(count)) %>% pull(s) %>% sort %>% head(-1)
beta_0 = readRDS("dev/beta_0.rds")


# Iterate over runs
readRDS(input_file) %>%
  mutate( results_sccomp = map2(
    data, run,
    ~  .x %>%
      sccomp_glm(
        ~type,
        sample, cell_type, generated_counts,
        check_outliers = add_outliers==1,
        approximate_posterior_inference = FALSE,
        percent_false_positive = 0.001 * 100,
        seed = .y * 2,
        prior_mean_variable_association =  list(intercept = c(5.6260004, 0.1), slope = c(-0.6940178,  0.01), standard_deviation = c(0.816423129, 0.01))
      )
  )) %>%
  saveRDS(output_file)


