


library(rstan)
library(gtools)
library(tidyverse)
library(cellsig)
library(magrittr)
library(tidybayes)

# Real data
sccomp::cell_counts %>%
  nest(data = -sample) %>%
  mutate(sample_tot = map_int(data, ~ sum(.x$count))) %>%
  unnest(data) %>%
  mutate(frac_logit = (count/sample_tot) %>% boot::logit() ) %>%

  group_by(cell_type, type) %>%
  summarise(m=mean(frac_logit), s=sd(frac_logit)) %>%
  ggplot(aes(m, s, color=type)) +
  geom_point() +
  geom_smooth(method="lm")


input_df =
  sccomp::cell_counts %>%
  mutate(count = count + 1L) %>%
  nest(data = -sample) %>%
  mutate(sample_tot = map_int(data, ~ sum(.x$count))) %>%
  unnest(data) %>%
  select(-c(phenotype)) %>%
  spread(cell_type, count)

f4 = stan(file = "inst/stan/glm_multi_beta_binomial.stan",
          data = list(
            N = 20,
            M = 36,
            tot = input_df$sample_tot,
            y = input_df %>% select(-type, -sample_tot) %>% nanny::as_matrix(rownames = sample),
            X = input_df %>% select(sample, type) %>% model.matrix(~ type, data=.)
          ),
          cores = 4
)

f4 %>% plot(pars="beta")
f4 %>% plot(pars="precision")


f3 %>%
  spread_draws(precision[M], beta[C,M]) %>%
  filter(C==1) %>%
  mean_qi() %>%
  ggplot(aes(beta, precision)) +
  geom_errorbar(aes(ymin = precision.lower, ymax = precision.upper)) +
  geom_point() +
  geom_smooth(method="lm")
