
res = readRDS("dev/data_integration/estimate_dirichlet_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds") 


data_proportion =
  res %>%
  distinct() %>%
  unnest(count_data) %>%
  select(cell_group, sample, outlier, count, type, significant) %>%
  with_groups(sample, ~ mutate(.x, proportion = (count)/sum(count)) ) 

simulated_proportion = 
  res %>%
  simulate_data(sample, cell_group, number_of_draws = 10) %>%
  unnest(generated_data) %>% 
  
  left_join(data_proportion %>% distinct(type, sample))



ggplot() +
  
  geom_boxplot(
    aes(type, generated_proportions),
    outlier.shape = NA, alpha=0.2,
    data = simulated_proportion, color="blue"
  ) + 
  geom_jitter(aes(type, generated_proportions), color="blue" ,alpha=0.2, size = 0.6, data = simulated_proportion) + 
  
  geom_boxplot(
    aes(type, proportion, fill=significant),
    outlier.shape = NA, 
    data = data_proportion |> filter(!outlier)
  ) + 
  geom_jitter(aes(type, proportion, color=outlier), size = 1, data = data_proportion) + 
  
  facet_wrap(~ interaction(cell_group), scale="free_y") +
  scale_y_continuous(trans="logit") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  xlab("Biological condition") + 
  ylab("Cell-group proportion") + 
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))
