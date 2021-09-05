library(tidyverse)
library(tidyseurat)
library(sccomp)
library(job)

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


# datasets_benign  =
#   benign_alive %>%
#   nest(data = - sample) %>%
#   mutate(data = map(
#     data,
#     ~ .x %>%
#       tidysc::adjust_abundance(
#         ~ subsets_Mito_percent + mito_RPS
#       ) %>%
#       Seurat::RunPCA() %>%
#       tidysc::cluster_elements()
#   ))

# Fit
job({
  estimate_benign =
    readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/expanded_analyses_with_control/benign_adjusted_cell_type.rds") |>
  sccomp_glm(
    formula = ~ 1,
    .sample = sample, .cell_group = curated_cell_type_pretty,
    approximate_posterior_inference = F
  )
})

estimate_benign %>% attr("mean_concentration_association")
# [1]  4.6424601 -0.8061759
# prec_sd = 0.8282289

job({
  estimate_oligo =
  counts_obj  |>
  filter(type != "benign") |>
  sccomp_glm(
    formula = ~ 1,
    sample, cell_group, count,
    approximate_posterior_inference = F
  )
})

estimate_oligo %>% attr("mean_concentration_association")
# [1]  5.6260004 -0.6940178
# prec_sd  = 0.3312485

job({
  estimate_UVM =
  readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
  sccomp_glm(
    formula = ~ 1,
    sample, cell_type,
    approximate_posterior_inference = F
  )
})

estimate_UVM %>% attr("mean_concentration_association")
# [1]  3.0423138 -0.6920534
# prec_sd  = 0.1987102

intercept = c(4.6424601, 5.6260004, 3.0423138)
slope = c(-0.8061759, -0.6940178, -0.6920534)
standard_deviation = c(0.8282289, 0.3312485, 0.1987102)

prior_mean_variable_association = list(
  intercept = c(mean(intercept), sd(intercept)),
  slope = c(mean(slope), sd(slope)),
  standard_deviation = c(mean(standard_deviation), sd(standard_deviation))
)

save(prior_mean_variable_association, file="data/prior_mean_variable_association.rda", compress = "xz")
