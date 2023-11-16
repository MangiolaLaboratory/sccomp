library(tidyverse)
library(tidyseurat)
library(sccomp)
library(job)
library(patchwork)


prior_mean_overdispersion_association = list(intercept = c(5, 5), slope = c(0,  5), standard_deviation = c(5,5))

job({

  load("data/counts_obj.rda")
    counts_obj  |>
    mutate(is_benign = type=="benign") |>
      rename(cell_type = cell_group) %>%
    sccomp_glm(
      formula = ~ is_benign,
      sample, cell_type, count,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_overdispersion_association = prior_mean_overdispersion_association
    ) %>%
      saveRDS("dev/data_integration/estimate_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds")
})

job({
    readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
    rename(type = `Sample Type`) %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_overdispersion_association = prior_mean_overdispersion_association
    ) %>%
    saveRDS("dev/data_integration/estimate_GSE139829_uveal_melanoma.rds")

})

job({
    readRDS("dev/data_integration/SCP1288_renal_cell_carcinoma.rds")  |>
    tidyseurat::filter(!is.na(sample) & !is.na(cell_type) & !is.na(sex))  |>
    sccomp_glm(
      formula = ~ sex,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_overdispersion_association = prior_mean_overdispersion_association
    ) %>%
    saveRDS("dev/data_integration/estimate_SCP1288_renal_cell_carcinoma.rds")
})

job({
    readRDS("dev/data_integration/SCP1039_bc_cells.rds")  |>
    #mutate(type = subtype=="TNBC") %>%
    mutate(type = factor(subtype, levels = c("TNBC", "HER2+", "ER+"))) %>%
    mutate(cell_type = celltype_subset) %>%
    sccomp_glm(
      formula = ~ type,
      formula_variability = ~ type,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_overdispersion_association = prior_mean_overdispersion_association
    ) %>%
    saveRDS("dev/data_integration/estimate_SCP1039_bc_cells.rds")
})

job({
    readRDS("dev/data_integration/s41587-020-0602-4_COVID_19.rds")  |>
    mutate(is_critical = severity=="critical") %>%
    sccomp_glm(
      formula = ~ is_critical,
      formula_variability = ~is_critical,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_overdispersion_association = prior_mean_overdispersion_association
    ) %>%
    saveRDS("dev/data_integration/estimate_s41587-020-0602-4_COVID_19.rds")
})

job({
    readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
    sccomp_glm(
      formula = ~ time,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_overdispersion_association = prior_mean_overdispersion_association
    ) %>%
    saveRDS("dev/data_integration/estimate_GSE120575_melanoma.rds")
})

job({

  library(tidySingleCellExperiment)

  readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/sccomp/dev/data_integration/BRCA1_s41467-021-21783-3.rds") %>%
    filter(ptime %>% is.na() %>% `!`) %>%

    # Scale ptime
    mutate(ptime = scales::rescale(ptime)) %>%
    rename(cell_type = CellTypesFinal) %>%
    rename(sample = Sample) %>%
    sccomp_glm(
      formula = ~ ptime,
      sample, cell_type ,
      approximate_posterior_inference = FALSE,
      variance_association = FALSE,
      prior_mean_overdispersion_association = prior_mean_overdispersion_association
    ) %>%
    saveRDS("dev/data_integration/estimate_BRCA1_s41467-021-21783-3.rds")
})
