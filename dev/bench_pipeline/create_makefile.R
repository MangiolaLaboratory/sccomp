#~/unix3XX/third_party_sofware/cctools-7.1.5-x86_64-centos6/bin/makeflow -T torque  --do-not-save-failed-output dev/TCGA_makeflow_pipeline/makefile_ARMET_TCGA.makeflow

library(tidyverse)
library(glue)

root_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/sccomp/dev/"
code_directory = glue("{root_directory}benchmark_code/")
results_directory = glue("{root_directory}benchmark_results/")

tab="\t"

commands_df =
  expand_grid(
  slope =seq(0.1, 2, length.out = 20),
  n_samples = c(2, 4, 6, 10, 15, 20),
  n_cell_type = c(5, 10, 20),
  max_cell_counts_per_sample = c(1000),
  add_outliers = c(0, 1)
) %>%

  # Add algorithm
  expand_grid(method = c("sccomp", "dirichletMultinomial", "rest")) %>%

  # Make commands
  unite("file_prefix" , c(slope, n_samples, n_cell_type, max_cell_counts_per_sample, add_outliers), remove = FALSE) %>%
  mutate(input_file = sprintf("input_%s.rds", file_prefix)) %>%
  mutate(output_file = sprintf("output_%s_%s.rds", file_prefix, method)) %>%
  mutate(parsed_file = glue("parsed_{file_prefix}.rds")) %>%
  mutate(auc_file = glue("auc_{file_prefix}.rds")) %>%
  mutate(create_input_command = glue("{results_directory}{input_file}:\n{tab}Rscript {code_directory}create_input.R {slope} {n_samples} {n_cell_type} {max_cell_counts_per_sample} {add_outliers} {results_directory}{input_file}")) %>%
  mutate(estimate_command = glue("{results_directory}{output_file}:{results_directory}{input_file}\n{tab}Rscript {code_directory}estimate_{method}.R {results_directory}{input_file} {results_directory}{output_file} {add_outliers}")) %>%

  # Filter 2 samples for other methods
  filter(!(method!="sccomp" & n_samples <3))

# Create input
"CATEGORY=create_input\nMEMORY=10024\nCORES=4\nWALL_TIME=1500" %>%
  c(
    commands_df %>% distinct(create_input_command) %>% pull(create_input_command)
  ) %>%

  # Estimate sccomp
  c("CATEGORY=estimate_sccomp\nMEMORY=30024\nCORES=12") %>%
  c(
    commands_df %>% filter(method %in% c("sccomp", "dirichletMultinomial")) %>%  distinct(estimate_command) %>% pull(estimate_command)
  ) %>%

  # Estimate sccomp
  c("CATEGORY=estimate_rest\nMEMORY=10024\nCORES=4\nWALL_TIME=1500") %>%
  c(
    commands_df %>% filter(method %in% c("rest")) %>% distinct(estimate_command) %>% pull(estimate_command)
  ) %>%

  # Calculate significance
  c("CATEGORY=significance\nMEMORY=10024\nCORES=4\nWALL_TIME=1500") %>%
  c(

      commands_df %>%

      # PArse separately because now I have
      mutate(output_file = glue("{results_directory}{output_file}")) %>%

      pivot_wider(names_from = method, values_from = c(input_file, output_file, estimate_command)) %>%
      tidyr::unite( "output_files", contains("output_file"), sep=" ", na.rm = TRUE) %>%

      mutate(parse_estimates_command = glue("{results_directory}{parsed_file}:{output_files}\n{tab}Rscript {code_directory}parse_estimates.R {slope} {n_samples} {n_cell_type} {max_cell_counts_per_sample} {add_outliers} {results_directory}{parsed_file} {output_files}")) %>%
      distinct(parse_estimates_command) %>%
      pull(parse_estimates_command)
  ) %>%

  # Calculate AUC
  c("CATEGORY=auc\nMEMORY=10024\nCORES=4\nWALL_TIME=1500") %>%
  c(

    commands_df %>%
      select(-c(input_file, output_file, estimate_command)) %>%
      distinct() %>%
      mutate(auc_command = glue("{results_directory}{auc_file}:{results_directory}{parsed_file}\n{tab}Rscript {code_directory}calculate_auc.R {results_directory}{parsed_file} {results_directory}{auc_file}")) %>%
      distinct(auc_command) %>%
      pull(auc_command)
  ) %>%


  write_lines(glue("{code_directory}/makefile.makeflow"))

