library(patchwork)
library(bayestestR)
library(yardstick)
library(forcats)
library(glue)
library(purrr)
#
# dir("dev/benchmark_results_f58cd58973261636d29f0462eb9b406d762f3895 345a2afa986458672dd95f5dff3f5d973bf1d994//", pattern = "auc", full.names = TRUE) %>%
#   map_dfr(~ .x %>% readRDS()) %>%
#   nest(data = -c(add_outliers, max_cell_counts_per_sample)) %>%
#   pull(data) %>% .[[2]] %>%
#
#
#   nest(data = -c(  slope, n_samples, n_cell_type )) %>%
#   mutate(min_auc = map_dbl(data, ~ min(.x$auc))) %>%
#   unnest(data) %>%
#
#   mutate(`Accuracy (AUC) relative to the smallest` = auc / min_auc) %>%
#
#
#   ggplot(aes(slope, `Accuracy (AUC) relative to the smallest`, color=name)) +
#   geom_line() +
#   facet_grid( n_cell_type ~ n_samples, scale="free_y") +
#   theme_bw()

result_directory = "dev/benchmark_results_e3498c0f9349ba9df28728253b3c1963fcda4b4c"
result_directory = "dev/benchmark_results"


plot_auc =
  dir(result_directory, pattern = "auc", full.names = TRUE) %>%
  map_dfr(~ .x %>% readRDS()) %>%
  nest(data = c(name, auc)) %>%
  mutate(random_auc = map_dbl(data, ~ .x %>% filter(name=="random") %>% pull(auc))) %>%
  unnest(data) %>%

  # filter(name !="edgerRobust" & name !="speckle") %>%
  filter(name !="edgerRobust" & name !="speckle") %>%
  filter(n_samples >2) %>%
  mutate(name = if_else(name=="random", "zrandom", name)) %>%
  # Choose one mode
  nest(data = -c(add_outliers, max_cell_counts_per_sample)) %>%
  mutate(plot = map(
    data,
    ~ .x %>%

      mutate(`Accuracy (AUC)` = auc ) %>%
      mutate(n_cell_type_ = glue("{n_cell_type} cell types"), n_samples_ = glue("{n_samples} samples")) %>%
      rename(Algorithm = name) %>%

      ggplot(aes(slope, `Accuracy (AUC)`, color=Algorithm)) +
      geom_line(size=0.3) +
      facet_grid( fct_reorder(n_cell_type_,n_cell_type)  ~ fct_reorder(n_samples_, n_samples), scale="free_y") +
      geom_hline(yintercept = 0.1, linetype="dashed", color="grey") +
      scale_color_brewer(palette="Paired") +
      theme_bw() +
      ylab("Accuracy (AUC) up to 0.1 false-positive-rate") +
      xlab("Slope") +
      theme(
        strip.background =element_rect(fill="white", color="white")
      ) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  )) %>%
  pull(plot) %>%
  patchwork::wrap_plots(ncol = 1)

ggsave(
  "dev/plot_auc_WEHI_seminar.pdf",
  plot = plot_auc,
  units = c("mm"),
  width = 183 ,
  height = 183 ,
  limitsize = FALSE
)

plot_ROC =
  readRDS(glue("{result_directory}/parsed_1.5_10_10_1000_1.rds")) %>%
  mutate(roc = map(
    df_for_roc,
    ~ .x %>%
      roc_curve(event_level = 'second', truth = significant, probability)
  )) %>%
  select(-df_for_roc) %>%
  unnest(roc) %>%

  filter(name !="edgerRobust" & name !="speckle") %>%
  rename(Algorithm = name) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity, color=Algorithm)) +
  geom_rect( mapping=aes(xmin=0, xmax=0.1, ymin=0, ymax=1), fill="#EBECF0", color=NA, alpha=0.1) +
  geom_path() +
  geom_abline(lty = 3) +
  scale_color_brewer(palette="Paired") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_fixed()+
  ylab("True positive rate") +
  xlab("False positive rate")

ggsave(
  "dev/plot_ROC_WEHI_seminar.pdf",
  plot = plot_ROC,
  units = c("mm"),
  width = 183 ,
  height = 183/2 ,
  limitsize = FALSE
)
# readRDS("dev/benchmark_results/output_0.1_15_10_1000_0_rest.rds")  %>%
#
#
#   mutate(df_for_roc = map2(
#     data, results_edger,
#     ~ left_join(
#       .x %>%
#         unnest(coefficients) %>%
#         distinct(cell_type, beta_1),
#
#       .y  %>%
#         select(cell_type, false_discovery_rate= PValue , estimate = logFC),
#       by = "cell_type"
#     )
#   )) %>%
#   select(df_for_roc) %>%
#
#   unnest(df_for_roc) %>%
#   mutate(false_discovery_rate = if_else(beta_1 * estimate < 0, 1, false_discovery_rate)) %>%
#   mutate(significant = (beta_1 != 0) %>% as.factor()) %>%
#   mutate(probability = 1-false_discovery_rate)  %>%
#   #mutate(probability = runif(n(), 0, 1)) %>%
#   roc_curve(event_level = 'second', truth = significant, probability) %>%
#   ggplot(
#     aes(
#       x = 1 - specificity,
#       y = sensitivity
#     )
#   ) + # plot with 2 ROC curves for each model
#   geom_line(size = 1.1) +
#   geom_abline(slope = 1, intercept = 0, size = 0.4) +
#   scale_color_manual(values = c("#48466D", "#3D84A8")) +
#   coord_fixed() +
#   theme_cowplot()
#
# readRDS("dev/benchmark_results/output_0.3_15_10_1000_0_sccomp.rds")  %>%
#
#
#   mutate(df_for_roc = map2(
#     data, results_sccomp,
#     ~ left_join(
#       .x %>%
#         unnest(coefficients) %>%
#         distinct(cell_type, beta_1),
#
#       .y  %>%
#         select(cell_type,  false_discovery_rate , estimate = .median_type),
#       by = "cell_type"
#     )
#   )) %>%
#   select(df_for_roc) %>%
#
#   unnest(df_for_roc) %>%
#   mutate(false_discovery_rate = if_else(beta_1 * estimate < 0, 1, false_discovery_rate)) %>%
#   mutate(significant = (beta_1 != 0) %>% as.factor()) %>%
#   mutate(probability = 1-false_discovery_rate)  %>%
#   #mutate(probability = runif(n(), 0, 1)) %>%
#   roc_curve(event_level = 'second', truth = significant, probability) %>%
#   ggplot(
#     aes(
#       x = 1 - specificity,
#       y = sensitivity
#     )
#   ) + # plot with 2 ROC curves for each model
#   geom_line(size = 1.1) +
#   geom_abline(slope = 1, intercept = 0, size = 0.4) +
#   scale_color_manual(values = c("#48466D", "#3D84A8")) +
#   coord_fixed() +
#   theme_cowplot()
