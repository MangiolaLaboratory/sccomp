library(tidyverse)
library(magrittr)
library(tidybulk)
library(speckle)
library(edgeR)
library(limma)
library(broom)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = args[2]
add_outliers = as.integer(args[3])


# exposures = counts_obj %>% group_by(sample) %>% summarise(s=sum(count)) %>% pull(s) %>% sort %>% head(-1)
beta_0 = readRDS("dev/beta_0.rds")
source("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/sccomp/dev/quasi_binomial/counts_to_quasi_binomial_estimate.R", echo=TRUE)
source("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/sccomp/dev/rlm_robust/rlm_estimate.R", echo=TRUE)

# Iterate over runs
readRDS(input_file) %>%

  # # edgeR
  # mutate( results_edger = map(
  #   data,
  #   ~  .x %>%
  #     mutate(across(c(sample, cell_type), ~ as.character(.x))) %>%
  #     test_differential_abundance(
  #       ~ type,
  #       sample, cell_type, generated_counts
  #     ) %>%
  #     pivot_transcript(cell_type)
  # )) %>%
  #
  # # deseq2
  # mutate( results_deseq2 = map(
  #   data,
  #   ~ {
  #     input =  .x %>%
  #       mutate(across(c(sample, cell_type), ~ as.character(.x))) %>%
  #       mutate(type = factor(type))
  #
  #     # Sometimes gives error for non inversability
  #     tryCatch({
  #       test_differential_abundance(
  #         input,
  #         ~ type,
  #         sample, cell_type, generated_counts,
  #         method = "deseq2"
  #       ) %>%
  #       pivot_transcript(cell_type)
  #     },
  #       error=function(cond) {
  #
  #         # Return empty dataset
  #         return(
  #           tibble(
  #             cell_type = character(),
  #             coefficients = list(),
  #             baseMean = numeric(),
  #             log2FoldChange = numeric(),
  #             lfcSE    = numeric(),
  #             stat  = numeric(),
  #             pvalue  = numeric(),
  #             padj = numeric()
  #           )
  #         )
  #     })
  #
  #   }
  # )) %>%
  #
  # # voom
  # mutate( results_voom = map(
  #   data,
  #   ~  .x %>%
  #     mutate(across(c(sample, cell_type), ~ as.character(.x))) %>%
  #     # group_by(sample) %>%
  #     # mutate(proportion = (generated_counts+1)/sum(generated_counts+1)) %>%
  #     # ungroup(sample) %>%
  #     # mutate(rate = proportion %>% boot::logit()) %>%
  #     # mutate(rate = rate - min(rate)) %>%
  #     # mutate(counts = rate %>% exp()) %>%
  #     test_differential_abundance(
  #       ~ type,
  #       sample, cell_type, generated_counts,
  #       method = "limma_voom",
  #     ) %>%
  #     pivot_transcript(cell_type)
  # )) %>%

  # logit linear
  mutate( results_logitLinear = map(
    data,
    ~  .x %>%
      mutate(across(c(sample, cell_type), ~ as.character(.x))) %>%
      group_by(sample) %>%
      mutate(proportion = (generated_counts+1)/sum(generated_counts+1)) %>%
      ungroup(sample) %>%
      mutate(rate = proportion %>% boot::logit()) %>%
      nest(data = -cell_type) %>%
      mutate(fit = map(
        data,
        ~ lm(rate ~ type, data = .x) %>%
          broom::tidy() %>%
          filter(term =="type"))) %>%
      unnest(fit) %>%
      unnest(data) %>%
      pivot_transcript(cell_type)
  )) %>%

  # T test on proportions
  mutate( results_ttest = map(
    data,
    ~  .x %>%
      mutate(across(c(sample, cell_type), ~ as.character(.x))) %>%
      group_by(sample) %>%
      mutate(proportion = (generated_counts+1)/sum(generated_counts+1)) %>%
      ungroup(sample) %>%
      nest(data = -cell_type) %>%
      mutate(fit = map(
        data,
        ~ t.test(
          filter(.x, type==1) %>% pull(proportion),
          filter(.x, type==0) %>% pull(proportion)
        ) %>%
          broom::tidy()
      )) %>%
      unnest(fit) %>%
      unnest(data) %>%
      pivot_transcript(cell_type)
  )) %>%

  # speckle
  mutate( results_speckle = map(
    data,
    ~  .x %>%
      mutate(cell_number = map(generated_counts, ~ 1:.x)) %>%
      unnest(cell_number) %>%
      unite("cell", c(sample, cell_type, cell_number), remove = FALSE) %>%
      propeller(clusters =  .$cell_type, sample = .$sample, group = .$type) %>%
      as_tibble(rownames="cell_type")
  )) %>%

  # Quasi biinomial
  mutate( results_quasiBinomial = map(
    data,
    ~   counts_to_quasi_binomial(.x, sample, cell_type, generated_counts, type)
  )) %>%

  # rlm
  mutate( results_rlm = map(
    data,
    ~   counts_to_rlm(.x, sample, cell_type, generated_counts, type)
  )) %>%


  saveRDS(output_file)

