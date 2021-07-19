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
