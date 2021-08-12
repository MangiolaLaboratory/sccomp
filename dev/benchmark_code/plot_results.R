


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
