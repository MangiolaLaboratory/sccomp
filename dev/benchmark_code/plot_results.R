


dir("dev/benchmark_results/", pattern = "auc", full.names = TRUE) %>%
  map_dfr(~ .x %>% readRDS()) %>%
  nest(data = -c(add_outliers, max_cell_counts_per_sample)) %>%
  pull(data) %>% .[[1]] %>%
  ggplot(aes(slope, auc, color=name)) +
  geom_line() +
  facet_grid(n_samples ~ n_cell_type)
