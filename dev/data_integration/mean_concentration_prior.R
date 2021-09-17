library(tidyverse)
library(tidyseurat)
library(sccomp)
library(job)
library(patchwork)




# Fit
job({
  estimate_benign =
    readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/expanded_analyses_with_control/benign_adjusted_cell_type.rds") |>
  sccomp_glm(
    formula = ~ 1,
    .sample = sample, .cell_group = curated_cell_type_pretty,
    approximate_posterior_inference = F,
    prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
  )
})
# estimate_benign %>% attr("mean_concentration_association")
# [1]  4.6424601 -0.8061759
# prec_sd = 0.8282289

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
  filter(type != "benign") |>
  sccomp_glm(
    formula = ~ 1,
    sample, cell_group, count,
    approximate_posterior_inference = F,
    prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
  )
})
# estimate_oligo %>% attr("mean_concentration_association")
# [1]  5.6260004 -0.6940178
# prec_sd  = 0.3312485

job({
  estimate_UVM =
  readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
  sccomp_glm(
    formula = ~ 1,
    sample, cell_type,
    approximate_posterior_inference = F,
    prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
  )
})
# estimate_UVM %>% attr("mean_concentration_association")
# [1]  3.0423138 -0.6920534
# prec_sd  = 0.1987102

job({
  estimate_renal_cell_carcinoma =
    readRDS("dev/data_integration/SCP1288_renal_cell_carcinoma.rds")  |>
    tidyseurat::filter(!is.na(sample) & !is.na(cell_type))  |>
    sccomp_glm(
      formula = ~ 1,
      sample, cell_type,
      approximate_posterior_inference = F,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    )
})
# estimate_renal_cell_carcinoma %>% attr("mean_concentration_association")
# [1]  4.1541595    -0.7367941
# prec_sd  =  0.5060364

job({
  estimate_bc_cells =
    readRDS("dev/data_integration/SCP1039_bc_cells.rds")  |>
    sccomp_glm(
      formula = ~ 1,
      sample, cell_type,
      approximate_posterior_inference = F,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    )
})
# estimate_bc_cells %>% attr("mean_concentration_association")
# [1]  3.2800052 -0.7575131
# prec_sd  = 1.0363162

job({
  estimate_COVID =
    readRDS("dev/data_integration/s41587-020-0602-4_COVID_19.rds")  |>
    sccomp_glm(
      formula = ~ 1,
      sample, cell_type,
      approximate_posterior_inference = F,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    )
})
# estimate_COVID %>% attr("mean_concentration_association")
# [1]  3.8746815    -0.8179763
# prec_sd  = 0.6306725

job({
  estimate_melanoma =
    readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
    sccomp_glm(
      formula = ~ 1,
      sample, cell_type,
      approximate_posterior_inference = F,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    )
})
# estimate_melanoma %>% attr("mean_concentration_association")
# [1]  2.2340179    -0.7466973
# prec_sd  = 0.6198534


# intercept = c(4.6424601, 5.6260004, 3.0423138)
# slope = c(-0.8061759, -0.6940178, -0.6920534)
# standard_deviation = c(0.8282289, 0.3312485, 0.1987102)
#
# prior_mean_variable_association = list(
#   intercept = c(mean(intercept), sd(intercept)),
#   slope = c(mean(slope), sd(slope)),
#   standard_deviation = c(mean(standard_deviation), sd(standard_deviation))
# )

# save(prior_mean_variable_association, file="data/prior_mean_variable_association.rda", compress = "xz")

# Plor trends
slopes_df =

    list(
      estimate_benign = estimate_benign %>% attr("mean_concentration_association"),
      estimate_oligo = estimate_oligo %>% attr("mean_concentration_association"),
      estimate_UVM = estimate_UVM %>% attr("mean_concentration_association"),
      estimate_renal_cell_carcinoma = estimate_renal_cell_carcinoma %>% attr("mean_concentration_association"),
      estimate_bc_cells = estimate_bc_cells %>% attr("mean_concentration_association"),
      estimate_COVID = estimate_COVID %>% attr("mean_concentration_association"),
      estimate_melanoma = estimate_melanoma %>% attr("mean_concentration_association")
    ) %>%
      enframe("dataset") %>%
      mutate(value = map(
        value,
        ~ t(.x) %>%
          as.data.frame %>%
          setNames(c("intercept", "slope", "standard_deviation"))
      )) %>%
      unnest(value)

plot_association_all =
  list(
  estimate_benign %>% mutate(dataset = "estimate_benign"),
  estimate_oligo %>% mutate(dataset = "estimate_oligo"),
  estimate_UVM %>% mutate(dataset = "estimate_UVM"),
  estimate_renal_cell_carcinoma %>% mutate(dataset = "estimate_renal_cell_carcinoma"),
  estimate_bc_cells %>% mutate(dataset = "estimate_bc_cells"),
  estimate_COVID %>% mutate(dataset = "estimate_COVID"),
  estimate_melanoma %>% mutate(dataset = "estimate_melanoma")
) %>%
  reduce(bind_rows) %>%
  unnest(concentration) %>%
  ggplot(aes(`.median_(Intercept)`, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`, color=dataset),  alpha = 0.4) +
  geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`, color=dataset), alpha = 0.4) +
  geom_point(size=0.1) +
  geom_abline(aes(intercept =  intercept, slope = slope, color = dataset), linetype = "dashed", data = slopes_df) +
  scale_color_brewer(palette="Set1") +
  xlab("Category rate") +
  ylab("Category log-concentration") +
  theme_bw() +
  theme(legend.position = "none")

plot_association_one =
  list(
  estimate_benign %>% mutate(dataset = "estimate_benign")
) %>%
  reduce(bind_rows) %>%
  unnest(concentration) %>%
  ggplot(aes(`.median_(Intercept)`, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`, color=dataset),  alpha = 0.4) +
  geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`, color=dataset), alpha = 0.4) +
  geom_point(size=0.1) +
  geom_abline(aes(intercept =  intercept, slope = slope), linetype = "dashed", color="grey", data = slopes_df %>% filter(dataset == "estimate_benign")) +
  scale_color_brewer(palette="Set1") +
  guides(color='none') +
  xlab("Category rate") +
  ylab("Category log-concentration") +
  theme_bw()+
  theme(legend.position = "bottom")


# CyTOF
estimate_CyTOF =
  readRDS("dev/data_integration/Cytof_raw.rds") %>%
  as.data.frame() %>%
  setNames(c("cell_cluster", "sample", "count")) %>%
  as_tibble() %>%
  tidyr::extract(sample, c("factor_of_interest", "patient", "batch"), "([A-Z]+)([0-9]+)_B(.+)", remove = FALSE) %>%
  filter(batch=="1") %>%
  sccomp_glm(
    formula = ~ 1,
    sample, cell_cluster,count,
    approximate_posterior_inference = F,
    prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
  )
# estimate_CyTOF %>% attr("mean_concentration_association")
# [1]  4.2623500    -1.1306609
# prec_sd  = 0.4845913

estimate_CyTOF %>%
  unnest(concentration) %>%
  ggplot(aes(`.median_(Intercept)`, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`),  alpha = 0.4) +
  geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`), alpha = 0.4) +
  geom_point(size=0.1) +
  geom_abline(intercept =  4.2623500, slope =  -1.1306609, linetype = "dashed") +
  scale_color_brewer(palette="Set1") +
  xlab("Category rate") +
  ylab("Category log-concentration") +
  theme_bw() +
  theme(legend.position = "none")






(
  plot_association_one +  plot_association_all &
    plot_layout(guides = "collect") &
    theme(legend.position = 'bottom', legend.key.size = unit(0.2, 'cm'))
) %>%
  ggsave(
    "dev/data_integration/association_plot.pdf",
    plot = .,
    units = c("mm"),
    width = 183 ,
    height = 110 ,
    limitsize = FALSE
  )
