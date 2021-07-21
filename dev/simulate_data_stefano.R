library(sccomp)
library(tidyverse)
library(magrittr)
library(tidybulk)

input_data =
  expand_grid(
    sample = 1:10, cell_type = 1:10
  ) %>%
  nest(d = -cell_type) %>%
  mutate(beta_0 = rnorm(n = n(), 0, 1)) %>%
  mutate(beta_1 = (cell_type == 1) %>% as.integer %>% multiply_by(0.5)) %>%
  nest(coefficients = starts_with("beta_")) %>%
  unnest(d) %>%
  nest(d = -sample) %>%
  mutate(type = sample(c(0,1), size = n(), replace = TRUE)) %>%
  mutate(tot_count = sample(600:10000, size = n(), replace = TRUE)) %>%
  unnest(d)


# debugonce(simulate_data)
sim_data =
  simulate_data(input_data,
                formula = ~ type ,
                sample,
                cell_type, tot_count, coefficients
  )

sim_data %>%
  group_by(sample) %>%
  mutate(proportion = (.value+1)/sum(.value+1)) %>%
  ungroup(sample) %>%
  ggplot(aes(factor(type), proportion)) +
  geom_boxplot() +
  facet_wrap(~ cell_type)


# dmbvs
source(file.path("dev/dmbvs-master/dmbvs-master/code", "wrapper.R"))
source(file.path("dev/dmbvs-master/dmbvs-master/code", "helper_functions.R"))
install.packages("dirmult")

simdata = simulate_dirichlet_multinomial_regression(n_obs = 100, n_vars = 50,
                                                    n_taxa = 50, n_relevant_vars = 5,
                                                    n_relevant_taxa = 5)
results = dmbvs(XX = simdata$XX[,-1][,1, drop=FALSE], YY = simdata$YY,
                intercept_variance = 10, slab_variance = 10,
                bb_alpha = 0.02, bb_beta = 1.98, GG = 1100L, thin = 10L, burn = 100L,
                exec = file.path(".", "dev/dmbvs-master/dmbvs-master/code", "dmbvs.x"), output_location = "dev/dmbvs-master")

# edgeR
sim_data %>%
  mutate(across(c(sample, cell_type), ~ as.character(.x))) %>%
  test_differential_abundance(
    ~ type,
    sample, cell_type, .value
  ) %>%
  pivot_transcript(cell_type) %>% filter(FDR<0.05)

# Voom Limma
sim_data %>%
  mutate(across(c(sample, cell_type), ~ as.character(.x))) %>%
  group_by(sample) %>%
  mutate(proportion = (.value+1)/sum(.value+1)) %>%
  ungroup(sample) %>%
  mutate(rate = proportion %>% boot::logit()) %>%
  mutate(rate = rate - min(rate)) %>%
  mutate(counts = rate %>% exp()) %>%
  test_differential_abundance(
    ~ type,
    sample, cell_type, counts,
    method = "limma_voom",
  ) %>%
  pivot_transcript(cell_type) %>%
  filter(adj.P.Val<0.05)

# Speckle
library(speckle)

sim_data %>%
  mutate(cell_number = map(.value, ~ 1:.x)) %>%
  unnest(cell_number) %>%
  unite("cell", c(sample, cell_type, cell_number), remove = FALSE) %>%
  propeller(clusters =  .$cell_type, sample = .$sample, group = .$type) %>%
  as_tibble(rownames="curated_cell_type")

# sccomp
sim_data %>%
  mutate(.value = as.integer(.value)) %>%
  sccomp_glm(
    ~type,
    sample, cell_type, .value,
    check_outliers = FALSE,
    approximate_posterior_inference = FALSE
  ) %>%
  dplyr::select(contains("type"))
