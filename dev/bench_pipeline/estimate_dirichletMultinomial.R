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
  mutate( results_DirichletMultinomial = map2(
    data, run,
    ~  .x %>%
      mutate(.value = as.integer(.value)) %>%
      sccomp_glm(
        ~type,
        sample, cell_type, .value,
        check_outliers = FALSE,
        approximate_posterior_inference = FALSE,
        percent_false_positive = 0.001 * 100,
        noise_model="dirichlet_multinomial",
        seed = .y * 2
      )
  )) %>%
  saveRDS(output_file)

