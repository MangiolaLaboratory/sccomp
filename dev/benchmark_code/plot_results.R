


dir("dev/benchmark_results_f58cd58973261636d29f0462eb9b406d762f3895 345a2afa986458672dd95f5dff3f5d973bf1d994//", pattern = "auc", full.names = TRUE) %>%
  map_dfr(~ .x %>% readRDS()) %>%
  nest(data = -c(add_outliers, max_cell_counts_per_sample)) %>%
  pull(data) %>% .[[2]] %>%


  nest(data = -c(  slope, n_samples, n_cell_type )) %>%
  mutate(min_auc = map_dbl(data, ~ min(.x$auc))) %>%
  unnest(data) %>%

  mutate(`Accuracy (AUC) relative to the smallest` = auc / min_auc) %>%


  ggplot(aes(slope, `Accuracy (AUC) relative to the smallest`, color=name)) +
  geom_line() +
  facet_grid( n_cell_type ~ n_samples, scale="free_y") +
  theme_bw()


dir("dev/benchmark_results", pattern = "auc", full.names = TRUE) %>%
  map_dfr(~ .x %>% readRDS()) %>%
  nest(data = -c(add_outliers, max_cell_counts_per_sample)) %>%
  pull(data) %>% .[[2]] %>%


  nest(data = -c(  slope, n_samples, n_cell_type )) %>%
  mutate(min_auc = map_dbl(data, ~ min(.x$auc))) %>%
  unnest(data) %>%

  mutate(`Accuracy (AUC) relative to the smallest` = auc / min_auc) %>%

  ggplot(aes(slope, auc, color=name)) +
  geom_line() +
  facet_grid( n_cell_type ~ n_samples, scale="free_y") +
  theme_bw()


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
