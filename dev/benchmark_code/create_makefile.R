#~/unix3XX/third_party_sofware/cctools-7.1.5-x86_64-centos6/bin/makeflow -T torque  --do-not-save-failed-output dev/TCGA_makeflow_pipeline/makefile_ARMET_TCGA.makeflow

library(tidyverse)
library(glue)

root_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/sccomp/dev/"
code_directory = glue("{root_directory}benchmark_code/")
results_directory = glue("{root_directory}benchmark_results/")

tab="\t"

commands_df =
  expand_grid(
  slope =seq(0.1, 3, length.out = 10),
  n_samples = c(10, 20),
  n_cell_type = c(10, 20),
  max_cell_counts_per_sample = c(1000, 3000),
  add_outliers = c(0, 1)
) %>%

  # Add algorithm
  expand_grid(method = c("sccomp", "dirichletMultinomial", "rest")) %>%

  # Make commands
  unite("file_prefix" , c(slope, n_samples, n_cell_type, max_cell_counts_per_sample, add_outliers), remove = FALSE) %>%
  mutate(input_file = sprintf("input_%s.rds", file_prefix)) %>%
  mutate(output_file = sprintf("output_%s_%s.rds", file_prefix, method)) %>%
  mutate(create_input_command = glue("{results_directory}{input_file}:\n{tab}Rscript {code_directory}create_input.R {slope} {n_samples} {n_cell_type} {max_cell_counts_per_sample} {add_outliers} {results_directory}{input_file}")) %>%
  mutate(estimate_command = glue("{results_directory}{output_file}:\n{tab}Rscript {code_directory}estimate_{method}.R {results_directory}{input_file} {results_directory}{output_file} {add_outliers}"))


# Create input
"CATEGORY=create_input\nMEMORY=10024\nCORES=4\nWALL_TIME=1500" %>%
  c(
    commands_df %>% distinct(create_input_command) %>% pull(create_input_command)
  ) %>%

  # Estimate
  c("CATEGORY=estimate_sccomp\nMEMORY=30024\nCORES=2\nWALL_TIME=36000") %>%
  c(
    commands_df %>% distinct(estimate_command) %>% pull(estimate_command)
  ) %>%


  write_lines(glue("{code_directory}/makefile.makeflow"))

