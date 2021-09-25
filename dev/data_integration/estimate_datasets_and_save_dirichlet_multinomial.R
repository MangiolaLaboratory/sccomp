library(tidyverse)
library(tidyseurat)
library(sccomp)
library(job)
library(patchwork)



job({

  load("data/counts_obj.rda")

  counts_obj %>%
    group_by(sample) %>%
    mutate(proportion = (count+1)/sum(count+1)) %>%
    ungroup(sample) %>%
    mutate(proportion_logit = boot::logit(proportion)) %>%
    group_by(cell_group) %>%
    summarise(mean = mean(proportion_logit), variance = sd(proportion_logit)) %>%
    ggplot(aes(mean, variance)) +
    geom_point() +
    geom_smooth(method="lm") +
    theme_bw()


    counts_obj  |>
    mutate(is_benign = type=="benign") |>
    sccomp_glm(
      formula = ~ is_benign,
      sample, cell_group, count,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2)),
     noise_model = "dirichlet_multinomial"
    ) %>%
      saveRDS("dev/data_integration/estimate_dirichlet_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds")
})

job({
  estimate_UVM =
    readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
    rename(type = `Sample Type`) %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2)),
      noise_model = "dirichlet_multinomial"
    ) %>%
    saveRDS("dev/data_integration/estimate_dirichlet_GSE139829_uveal_melanoma.rds")

})

job({
  estimate_renal_cell_carcinoma =
    readRDS("dev/data_integration/SCP1288_renal_cell_carcinoma.rds")  |>
    tidyseurat::filter(!is.na(sample) & !is.na(cell_type) & !is.na(sex))  |>
    sccomp_glm(
      formula = ~ sex,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2)),
      noise_model = "dirichlet_multinomial"
    ) %>%
    saveRDS("dev/data_integration/estimate_dirichlet_SCP1288_renal_cell_carcinoma.rds")
})

job({
  estimate_bc_cells =
    readRDS("dev/data_integration/SCP1039_bc_cells.rds")  |>
    mutate(type = subtype=="TNBC") %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2)),
      noise_model = "dirichlet_multinomial"
    ) %>%
    saveRDS("dev/data_integration/estimate_dirichlet_SCP1039_bc_cells.rds")
})

job({
  estimate_COVID =
    readRDS("dev/data_integration/s41587-020-0602-4_COVID_19.rds")  |>
    mutate(is_critical = severity=="critical") %>%
    sccomp_glm(
      formula = ~ is_critical,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2)),
      noise_model = "dirichlet_multinomial"
    ) %>%
    saveRDS("dev/data_integration/estimate_dirichlet_s41587-020-0602-4_COVID_19.rds")
})

job({
  estimate_melanoma =
    readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
    sccomp_glm(
      formula = ~ time,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2)),
      noise_model = "dirichlet_multinomial"
    ) %>%
    saveRDS("dev/data_integration/estimate_dirichlet_GSE120575_melanoma.rds")
})
