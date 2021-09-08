library(tidyverse)
library(magrittr)
library(tidybulk)
library(speckle)
library(edgeR)
library(limma)


# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = args[2]
add_outliers = as.integer(args[3])


# exposures = counts_obj %>% group_by(sample) %>% summarise(s=sum(count)) %>% pull(s) %>% sort %>% head(-1)
beta_0 = readRDS("dev/beta_0.rds")


# Iterate over runs
readRDS(input_file) %>%
  # edgeR
  mutate( results_edger = map(
    data,
    ~  .x %>%
      mutate(across(c(sample, cell_type), ~ as.character(.x))) %>%
      test_differential_abundance(
        ~ type,
        sample, cell_type, .value
      ) %>%
      pivot_transcript(cell_type)
  )) %>%

  # robust edgeR
  mutate( results_edgerRobust = map(
    data,
    ~  .x %>%
      mutate(across(c(sample, cell_type), ~ as.character(.x))) %>%
      test_differential_abundance(
        ~ type,
        sample, cell_type, .value,
        method = "edger_robust_likelihood_ratio"
      ) %>%
      pivot_transcript(cell_type)
  )) %>%
  # voom
  mutate( results_voom = map(
    data,
    ~  .x %>%
      mutate(across(c(sample, cell_type), ~ as.character(.x))) %>%
      # group_by(sample) %>%
      # mutate(proportion = (.value+1)/sum(.value+1)) %>%
      # ungroup(sample) %>%
      # mutate(rate = proportion %>% boot::logit()) %>%
      # mutate(rate = rate - min(rate)) %>%
      # mutate(counts = rate %>% exp()) %>%
      test_differential_abundance(
        ~ type,
        sample, cell_type, .value,
        method = "limma_voom",
      ) %>%
      pivot_transcript(cell_type)
  )) %>%

  # speckle
  mutate( results_speckle = map(
    data,
    ~  .x %>%
      mutate(cell_number = map(.value, ~ 1:.x)) %>%
      unnest(cell_number) %>%
      unite("cell", c(sample, cell_type, cell_number), remove = FALSE) %>%
      propeller(clusters =  .$cell_type, sample = .$sample, group = .$type) %>%
      as_tibble(rownames="cell_type")
  )) %>%
  saveRDS(output_file)

