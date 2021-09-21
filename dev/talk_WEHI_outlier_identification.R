


library(dplyr)
library(sccomp)
library(ggplot2)
library(forcats)
library(tidyr)
library(job)
data("seurat_obj")
data("sce_obj")
data("counts_obj")

job({
  res =
  counts_obj %>%
  sccomp_glm(
    ~ type,
    sample, cell_group, count,
    approximate_posterior_inference = FALSE,
    percent_false_positive = 1
  )

})






data_for_plot =
  res %>%
  tidyr::unnest(count_data) %>%
  #left_join(counts_obj %>% distinct(sample, type), by = c("sample")) %>%
  group_by(sample) %>%
  mutate(proportion = (count+1)/sum(count+1)) %>%
  ungroup(sample)



plot_input_data =
  ggplot() +
  geom_boxplot(
    aes(type, proportion),
    outlier.shape = NA,
    data = data_for_plot %>% filter(!outlier)
  ) +
  geom_jitter(aes(type, proportion), size = 1, data = data_for_plot) +
  facet_wrap(~ forcats::fct_reorder(cell_group, desc(abs(.median_typecancer))), scale="free_y") +
  scale_y_continuous(trans="logit") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  theme_bw() +
  theme(
    strip.background =element_rect(fill="white", color="white"),
    legend.position = "bottom", axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

plot_outliers =
  ggplot() +
  geom_boxplot(
    aes(type, proportion),
    outlier.shape = NA,
    data = data_for_plot %>% filter(!outlier)
  ) +
  geom_jitter(aes(type, proportion, color=outlier), size = 1, data = data_for_plot) +
  facet_wrap(~ forcats::fct_reorder(cell_group, desc(abs(.median_typecancer))), scale="free_y") +
  scale_y_continuous(trans="logit") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  theme_bw() +
  theme(
    strip.background =element_rect(fill="white", color="white"),
    legend.position = "bottom", axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


plot_fitted_outlier =
  ggplot() +
  geom_boxplot(
    aes(type, proportion, fill=significant),
    outlier.shape = NA,
    data = data_for_plot %>% filter(!outlier)
  ) +
  geom_jitter(aes(type, proportion, color=outlier), size = 1, data = data_for_plot) +
  facet_wrap(~ forcats::fct_reorder(cell_group, desc(abs(.median_typecancer))), scale="free_y") +
  scale_y_continuous(trans="logit") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  theme_bw() +
  theme(
    strip.background =element_rect(fill="white", color="white"),
    legend.position = "bottom", axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )



# rng =  rstan::gqs(
#   sccomp:::stanmodels$glm_multi_beta_binomial_generate_date,
#   #rstan::stan_model("inst/stan/glm_multi_beta_binomial_generate_date.stan"),
#   draws =  as.matrix(fit_object$fit),
#   data = fit_object$data_for_model
# )

# plot_generated =
#   sccomp:::draws_to_tibble_x_y(rng, "counts", "N", "M") %>%
#   filter(.draw == 1) %>%
#   left_join(tibble(N = 1:20, type = factor(fit_object$data_for_model$X[,2]))) %>%
#   mutate(M = as.character(M)) %>%
#   left_join(res %>% select(.median_typecancer) %>% mutate(M = as.character(1:n())), by="M") %>%
#   rename(sample = N, count = .value) %>%
#   with_groups(sample, ~ mutate(.x, proportion = (count+1)/sum(count+1)) )


job({
  fit_object =
    counts_obj %>%
    sccomp:::sccomp_glm(
      ~ type,
      sample, cell_group, count,
      approximate_posterior_inference = FALSE, variance_association = TRUE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    )
})

fit_object %>%
  simulate_data(sample, cell_group) %>%
  unnest(generated_data) %>%

  # Add type and median
  left_join(counts_obj %>% distinct(sample, type)) %>%
  left_join(fit_object %>% distinct(cell_group , .median_typecancer)) %>%

  ggplot() +
  geom_boxplot(
    aes(type, generated_proportions ),
    outlier.shape = NA
  ) +
  geom_jitter(aes(type, generated_proportions ), size = 1) +
  facet_wrap(~ forcats::fct_reorder(cell_group , desc(abs(.median_typecancer))), scale="free_y") +
  scale_y_continuous(trans="logit") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  xlab("Biological condition") +
  ylab("Cell-group generated_proportions ") +
  theme_bw() +
  theme(
    strip.background =element_rect(fill="white", color="white"),
    legend.position = "bottom", axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

plot_theoretical_distribution =
  sccomp:::draws_to_tibble_x_y(rng, "counts", "N", "M") %>%
  filter(.draw < 10) %>%
  left_join(tibble(N = 1:20, type = factor(fit_object$data_for_model$X[,2]))) %>%
  mutate(M = as.character(M)) %>%
  left_join(res %>% select(.median_typecancer) %>% mutate(M = as.character(1:n())), by="M") %>%
  mutate(sample = as.character(N)) %>%
  rename( count = .value) %>%
  with_groups(c(sample, .draw), ~ mutate(.x, proportion = (count+1)/sum(count+1)) ) %>%
  left_join(
    res %>%
      mutate(M = as.character(1:n())) %>%
      unnest(outliers) %>%
      nest(data = -sample) %>%
      mutate(sample = as.character(1:n())) %>% unnest(data)  %>%
      select(count_outlier = count, M, sample, outlier) %>%
      filter(outlier)
  ) %>%


  ggplot(aes(sample, count + 1)) +
  stat_halfeye(size = 0) +
  geom_point(aes(sample, count_outlier + 1), color = "#e11f28") +
  facet_wrap(~ forcats::fct_reorder(M, desc(abs(.median_typecancer))), scale="free_y") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  scale_y_log10() +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  theme_bw() +
  theme(
    strip.background =element_rect(fill="white", color="white"),
    legend.position = "bottom", axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )



ggsave(
  "dev/plot_input_data_WEHI_talk.pdf",
  plot = plot_input_data,
  units = c("mm"),
  width = 183 ,
  height = 183 ,
  limitsize = FALSE
)


ggsave(
  "dev/plot_outliers_WEHI_talk.pdf",
  plot = plot_outliers,
  units = c("mm"),
  width = 183 ,
  height = 183 ,
  limitsize = FALSE
)

ggsave(
  "dev/plot_fitted_outlier_WEHI_talk.pdf",
  plot = plot_fitted_outlier,
  units = c("mm"),
  width = 183 ,
  height = 183 ,
  limitsize = FALSE
)

ggsave(
  "dev/plot_generated_WEHI_talk.pdf",
  plot = plot_generated,
  units = c("mm"),
  width = 183 ,
  height = 183,
  limitsize = FALSE
)


ggsave(
  "dev/plot_theoretical_distribution_WEHI_talk.pdf",
  plot = plot_theoretical_distribution,
  units = c("mm"),
  width = 183 ,
  height = 183,
  limitsize = FALSE
)










# Composition with healthy
benign_adjusted_cell_type <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/expanded_analyses_with_control/benign_adjusted_cell_type.rds")
cancer_counts_alive_seurat <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/expanded_analyses_with_control/cancer_only_analyses/cancer_counts_alive_seurat.rds")
integrated_counts_curated = readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/oligo_breast/expanded_analyses_with_control/cancer_only_analyses/integrated_counts_curated.rds")



glm_input_benign =
  cancer_counts_alive_seurat  %>%
  left_join(integrated_counts_curated %>% as_tibble() %>% select(cell, curated_cell_type_pretty )) %>%
  bind_rows(
    benign_adjusted_cell_type
  ) %>%
  mutate(
    curated_cell_type_pretty = case_when(
      curated_cell_type_pretty %in% c("CD8 em 1", "CD8 em 2", "CD8 em 3") ~ "CD8 em",
      grepl("CD4 cm", curated_cell_type_pretty) ~ "CD4 cm",
      grepl("Mono NKG7", curated_cell_type_pretty) ~ "Mono NKG7",
      grepl("Mac M1", curated_cell_type_pretty) ~ "Mac M1",
      curated_cell_type_pretty == "CD4 em high cytokine" ~ "CD4 em",

      TRUE ~ curated_cell_type_pretty
    )
  ) %>%
  mutate(curated_cell_type_pretty = if_else(is.na(curated_cell_type_pretty), "none", curated_cell_type_pretty)) %>%
  mutate(is_benign = type=="benign") %>%
  filter(curated_cell_type_pretty!="none")

rm(
  benign_adjusted_cell_type,
  cancer_counts_alive_seurat,
  integrated_counts_curated
)
gc()

job({
  estimate_benign_WEHI_talk =
    glm_input_benign  %>%
    sccomp_glm(~is_benign, sample, curated_cell_type_pretty,
               approximate_posterior_inference = FALSE,
               prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    )

})

job({
  fit_object =
    glm_input_benign %>%
    sccomp:::estimate_multi_beta_binomial_glm(
      ~is_benign, sample, curated_cell_type_pretty,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    )
})




data_for_plot =
  estimate_benign_WEHI_talk %>%
  #filter(curated_cell_type_pretty!="none") %>%
  arrange(desc(abs(.median_is_benignTRUE))) %>%
  tidyr::unnest(count_data ) %>%
  #left_join(glm_input_benign %>% distinct(sample, type, is_benign), by = c( "sample")) %>%
  with_groups(sample, ~ mutate(.x, proportion = (count+1)/sum(count+1)) )




plot_composition =
  ggplot() +
  geom_boxplot(aes(is_benign, generated_proportions ), color="blue", data = data_for_plot %>% unnest(generated_data ) ,  height=0, alpha=0.5) +
  geom_boxplot(
    aes(is_benign, proportion, fill=significant),
    outlier.shape = NA,
    data = data_for_plot %>% filter(!outlier)
  ) +
  geom_jitter(aes(is_benign, proportion, color=outlier), size = 1, data = data_for_plot) +
  facet_wrap(~ forcats::fct_reorder(curated_cell_type_pretty, desc(abs(.median_is_benignTRUE))), scale="free_y", ncol=4) +
  scale_y_continuous(trans="logit") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  theme_bw() +
  theme(
    strip.background =element_rect(fill="white", color="white"),
    legend.position = "bottom", axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


rng =  rstan::gqs(
  sccomp:::stanmodels$glm_multi_beta_binomial_generate_date,
  #rstan::stan_model("inst/stan/glm_multi_beta_binomial_generate_date.stan"),
  draws =  as.matrix(fit_object$fit),
  data = fit_object$data_for_model
)

both_data_sets =
  sccomp:::draws_to_tibble_x_y(rng, "counts", "N", "M") %>%
  filter(.draw == 1) %>%
  left_join(tibble(M = 1:24, type = data_for_plot$X[,2])) %>%
  rename(sample = N, count = .value) %>%
  with_groups(sample, ~ mutate(.x, proportion = (count+1)/sum(count+1)) ) %>%
  geom_boxplot(
    aes(is_benign, proportion, fill=significant),
    outlier.shape = NA
  ) +
  geom_jitter(aes(is_benign, proportion, color=outlier), size = 1) +
  facet_wrap(~ factor(M), scale="free_y", ncol=4) +
  scale_y_continuous(trans="logit") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  theme_bw() +
  theme(
    strip.background =element_rect(fill="white", color="white"),
    legend.position = "bottom", axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


data_for_plot_benign =
  estimate_benign %>%
  filter(curated_cell_type_pretty!="none") %>%
  arrange(desc(abs(.median_is_benignTRUE))) %>%
  #head(16) %>%
  tidyr::unnest(outliers) %>%
  left_join(glm_input_benign %>% count(sample, curated_cell_type_pretty, type, is_benign), by = c("curated_cell_type_pretty", "sample")) %>%
  filter(is_benign %>% is.na %>% `!`) %>%
  group_by(sample) %>%
  mutate(proportion = (n+1)/sum(n+1)) %>%
  ungroup(sample) %>%
  mutate(type = case_when(type=="OMBC" ~ "Stable", type=="MBC" ~ "Progressive", TRUE ~ type)) %>%
  mutate(type = forcats::fct_relevel(type, c("Stable", "Progressive"))) %>%
  mutate(type = forcats::fct_relevel(type, c("benign", "Stable", "Progressive"))) %>%
  mutate(is_benign = case_when(is_benign ~ "Healthy", TRUE ~ "Cancer")) %>%
  mutate(is_benign = is_benign %>% factor(levels = c("Healthy", "Cancer"))) %>%
  filter(curated_cell_type_pretty != "Stem") %>%
  filter(curated_cell_type_pretty != "CD4 ribosome")

# Plot for top three for cancer
composition_benign_for_cancer =
  ggplot() +
  geom_boxplot(
    aes(is_benign, proportion, fill=type),
    outlier.shape = NA,
    data = data_for_plot_benign %>% filter(curated_cell_type_pretty %in% c("CD4 em", "MAIT", "T gd2", "CD8 em")) %>% filter(!outlier)
  ) +
  geom_jitter(aes(is_benign, proportion, color=outlier), size = 1, data = data_for_plot_benign%>% filter(curated_cell_type_pretty %in% c("CD4 em", "MAIT", "T gd2", "CD8 em"))) +
  facet_wrap(~ forcats::fct_reorder(curated_cell_type_pretty, desc(abs(.median_is_benignTRUE))), scale="free_y", ncol=4) +
  scale_y_continuous(trans="logit") +
  scale_color_manual(values = c("black", "#e11f28")) +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  theme(
    strip.background =element_rect(fill="white", color="white"),
    legend.position = "bottom", axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


plot_pairs = fit_object$fit %>% pairs(pars = c("beta", "prec_coeff") )

ggsave(
  "dev/approximation_dirichlet_multinomial_pairs.pdf", plot = plot_pairs,
  units = c("mm"),
  width = 183/2 ,
  height = 183/2 ,
  limitsize = FALSE
)
