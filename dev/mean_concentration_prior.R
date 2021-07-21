library(tidyverse)
library(tidyseurat)
library(sccomp)


load("data/counts_obj.rda")
library(tidyverse)

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

benign <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/expanded_analyses_with_control/benign_adjusted_cell_type.rds")

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
estimate =
  benign |>
  sccomp_glm(
    formula = ~ 1,
    .sample = sample, .cell_group = curated_cell_type_pretty,
    approximate_posterior_inference = F
  )

estimate %>% attr("mean_concentration_association")
# [1]  4.6424601 -0.8061759

estimate_oligo =
  counts_obj  |>
  filter(type != "benign") |>
  sccomp_glm(
    formula = ~ 1,
    sample, cell_group, count,
    approximate_posterior_inference = F
  )

estimate_oligo %>% attr("mean_concentration_association")
# [1]  5.6260004 -0.6940178
# prec_sd  = 0.816423129
