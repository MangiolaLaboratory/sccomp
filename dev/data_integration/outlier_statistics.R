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

  estimate_oligo =
  counts_obj  |>
  mutate(is_benign = type=="benign") |>
  sccomp_glm(
    formula = ~ is_benign,
    sample, cell_group, count,
    approximate_posterior_inference = FALSE,
    prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
  )
})
estimate_oligo %>% attr("mean_concentration_association")
# [1]  5.6260004 -0.6940178
# prec_sd  = 0.3312485

job({
  estimate_UVM =
  readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
    rename(type = `Sample Type`) %>%
  sccomp_glm(
    formula = ~ type,
    sample, cell_type,
    approximate_posterior_inference = FALSE,
    prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
  )
})
estimate_UVM %>% attr("mean_concentration_association")
# [1]  3.0423138 -0.6920534
# prec_sd  = 0.1987102

job({
  estimate_renal_cell_carcinoma =
    readRDS("dev/data_integration/SCP1288_renal_cell_carcinoma.rds")  |>
    tidyseurat::filter(!is.na(sample) & !is.na(cell_type) & !is.na(sex))  |>
    sccomp_glm(
      formula = ~ sex,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    )
})
estimate_renal_cell_carcinoma %>% attr("mean_concentration_association")
# [1]  4.1541595    -0.7367941
# prec_sd  =  0.5060364

job({
  estimate_bc_cells =
    readRDS("dev/data_integration/SCP1039_bc_cells.rds")  |>
    mutate(type = subtype=="TNBC") %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    )
})
estimate_bc_cells %>% attr("mean_concentration_association")
# [1]  3.2800052 -0.7575131
# prec_sd  = 1.0363162

job({
  estimate_COVID =
    readRDS("dev/data_integration/s41587-020-0602-4_COVID_19.rds")  |>
    mutate(is_critical = severity=="critical") %>%
    sccomp_glm(
      formula = ~ is_critical,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    )
})
estimate_COVID %>% attr("mean_concentration_association")
# [1]  3.8746815    -0.8179763
# prec_sd  = 0.6306725

job({
  estimate_melanoma =
    readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
    sccomp_glm(
      formula = ~ time,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    )
})
estimate_melanoma %>% attr("mean_concentration_association")
# [1]  2.2340179    -0.7466973
# prec_sd  = 0.6198534


# Plor trends
size_df =

    list(
      estimate_oligo = estimate_oligo %>% nrow(),
      estimate_UVM = estimate_UVM  %>% nrow(),
      estimate_renal_cell_carcinoma = estimate_renal_cell_carcinoma  %>% nrow(),
      estimate_bc_cells = estimate_bc_cells  %>% nrow(),
      estimate_COVID = estimate_COVID  %>% nrow(),
      estimate_melanoma = estimate_melanoma  %>% nrow()
    ) %>%
      enframe("dataset", "number_of_cell_types") %>%
      unnest(number_of_cell_types)

plot_outlier =
  list(
  estimate_oligo %>% mutate(dataset = "estimate_oligo"),
  estimate_UVM %>% mutate(dataset = "estimate_UVM"),
  estimate_renal_cell_carcinoma %>% mutate(dataset = "estimate_renal_cell_carcinoma"),
  estimate_bc_cells %>% mutate(dataset = "estimate_bc_cells"),
  estimate_COVID %>% mutate(dataset = "estimate_COVID"),
  estimate_melanoma %>% mutate(dataset = "estimate_melanoma")
) %>%
  reduce(bind_rows) %>%
  unnest(outliers) %>%
  filter(outlier) %>%
  count(dataset, significant) %>%
  with_groups(dataset, ~ mutate(.x, nn = sum(n))) %>%
  mutate(dataset = gsub("estimate_", "", dataset)) %>%

  ggplot() +
  geom_bar(aes(forcats::fct_reorder(dataset, nn, .desc = TRUE ),  n, fill=significant  ), stat = "identity", position = "stack")+
  geom_point(aes(dataset, number_of_cell_types), data = size_df %>% mutate(dataset = gsub("estimate_", "", dataset))  ) +
  scale_fill_manual(values = c("grey", "#e11f28")) +
  xlab("Dataset") +
  ylab("Count of outliers") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle=30, hjust = 1))


ggsave(
  "dev/outliers_across_datasets_WEHI_seminar_2.pdf",
  plot = plot_outlier,
  units = c("mm"),
  width = 80 ,
  height = 183 ,
  limitsize = FALSE
)
