library(sccomp)
library(tidyverse)
library(magrittr)
library(tidybulk)
library(speckle)
library(edgeR)
library(limma)

input_data =
  expand_grid(
    sample = 1:12, cell_type = 1:12
  ) %>%
  nest(d = -cell_type) %>%
  mutate(beta_0 = rnorm(n = n(), 0, 1)) %>%
  mutate(  beta_1 = case_when(
      cell_type %in% 1:3 ~ 1,
      cell_type %in% 4:6 ~ -1,
      TRUE ~ 0
    )) %>%
  mutate(beta_1 = beta_1 %>% multiply_by(0.5)) %>%
  nest(coefficients = starts_with("beta_")) %>%
  unnest(d) %>%
  nest(d = -sample) %>%
  mutate(type = sample(c(0,1), size = n(), replace = TRUE)) %>%
  mutate(tot_count = sample(600:10000, size = n(), replace = TRUE)) %>%
  unnest(d)


# # debugonce(simulate_data)
# sim_data =
#   simulate_data(input_data,
#                 formula = ~ type ,
#                 sample,
#                 cell_type, tot_count, coefficients
#   )
#
# sim_data %>%
#   group_by(sample) %>%
#   mutate(proportion = (.value+1)/sum(.value+1)) %>%
#   ungroup(sample) %>%
#   ggplot(aes(factor(type), proportion)) +
#   geom_boxplot() +
#   facet_wrap(~ cell_type)



# dmbvs
# source(file.path("dev/dmbvs-master/dmbvs-master/code", "wrapper.R"))
# source(file.path("dev/dmbvs-master/dmbvs-master/code", "helper_functions.R"))
# install.packages("dirmult")
#
# simdata = simulate_dirichlet_multinomial_regression(n_obs = 100, n_vars = 50,
#                                                     n_taxa = 50, n_relevant_vars = 5,
#                                                     n_relevant_taxa = 5)
# results = dmbvs(XX = simdata$XX[,-1][,1, drop=FALSE], YY = simdata$YY,
#                 intercept_variance = 10, slab_variance = 10,
#                 bb_alpha = 0.02, bb_beta = 1.98, GG = 1100L, thin = 10L, burn = 100L,
#                 exec = file.path(".", "dev/dmbvs-master/dmbvs-master/code", "dmbvs.x"), output_location = "dev/dmbvs-master")

#probs = c(0, 0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
probs = seq(0, 0.2,length.out = 20)

# Iterate over runs
benchmark =
  tibble(run = 1:20) %>%
  mutate(data = map(
    run,
    ~ simulate_data(input_data,
                    formula = ~ type ,
                    sample,
                    cell_type, tot_count, coefficients,
                    seed = .x * 2
    )
  )) %>%

  # edgeR
  mutate( results_edger = map(
    data,
    ~  .x %>%
      mutate(across(c(sample, cell_type), ~ as.character(.x))) %>%
      test_differential_abundance(
        ~ type,
        sample, cell_type, .value
      ) %>%
      pivot_transcript(cell_type)
  )) %>%

  # voom
  mutate( results_voom = map(
    data,
    ~  .x %>%
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
      pivot_transcript(cell_type)
  )) %>%

  # speckle
  mutate( results_speckle = map(
    data,
    ~  .x %>%
      mutate(cell_number = map(.value, ~ 1:.x)) %>%
      unnest(cell_number) %>%
      unite("cell", c(sample, cell_type, cell_number), remove = FALSE) %>%
      propeller(clusters =  .$cell_type, sample = .$sample, group = .$type) %>%
      as_tibble(rownames="cell_type")
  )) %>%

  # sccomp
  mutate( results_sccomp = map2(
    data, run,
    ~  .x %>%
      mutate(.value = as.integer(.value)) %>%
      sccomp:::estimate_multi_beta_binomial_glm(
        ~type,
        sample, cell_type, .value,
        check_outliers = FALSE,
        approximate_posterior_inference = FALSE,
        percent_false_positive = 0.001 * 100,
        seed = .y * 2
      )
  ))

saveRDS(benchmark, "dev/benchmark.rds")



benchmark_hypothesis =
  benchmark %>%
  mutate(probs = map(run, ~ !!probs)) %>%
  unnest(probs) %>%
  mutate(hypothesis_edger = map2(results_edger, probs, ~ mutate(.x, positive = FDR<.y))) %>%
  mutate(hypothesis_voom = map2(results_voom, probs, ~ mutate(.x, positive = adj.P.Val<.y))) %>%
  mutate(hypothesis_speckle = map2(results_speckle, probs, ~ mutate(.x, positive = FDR<.y))) %>%
  mutate(hypothesis_sccomp = map2(
    results_sccomp, probs,
    ~ sccomp:::hypothesis_test_multi_beta_binomial_glm(
      sample, cell_type,
      .x$fit,
      .x$data_for_model, percent_false_positive = .y * 100,
      check_outliers = F,
      .x$truncation_df2
    )  %>%
      mutate(positive = (.lower_type * .upper_type) > 0)
  )) %>%
  dplyr::select(-contains("results"))



benchmark_hypothesis %>%
  dplyr::select(-contains("results")) %>%
  pivot_longer(contains("hypothesis"),names_prefix = "hypothesis_" ) %>%
  mutate(accuracy_df = map2(
    data, value             ,
    ~ left_join(

      .x %>% unnest(coefficients) %>% dplyr::select(cell_type, beta_1) %>% distinct %>% mutate(cell_type = as.character(cell_type)),
      .y %>% dplyr::select(cell_type, positive),
      by="cell_type"

    )

  )) %>%
  mutate(TP = map_int(accuracy_df, ~ .x %>% filter(positive & (beta_1 != 0)) %>% nrow())) %>%
  mutate(FP = map_int(accuracy_df, ~ .x %>% filter(positive & (beta_1 == 0)) %>% nrow())) %>%
  mutate(total_true_positive = 6, total_true_negative = 6) %>%
  mutate(FP_rate = FP/total_true_negative) %>%
  mutate(TP_rate = TP/total_true_positive) %>%
  group_by(name, probs) %>%
  summarise(mean_FP_rate = mean(FP_rate), mean_TP_rate = mean(TP_rate)) %>%
  arrange(mean_FP_rate) %>%
  ggplot(aes(mean_FP_rate, mean_TP_rate)) +
  geom_line(aes(color = name)) +
  xlim(0,0.1)






benchmark %>%
  mutate(accuracy_df = map2(
    data, results_sccomp             ,
    ~ left_join(

      .x %>% unnest(coefficients) %>% dplyr::select(cell_type, beta_1) %>% distinct %>% mutate(cell_type = as.character(cell_type)),
      .y %>% dplyr::select(cell_type, positive),
      by="cell_type"

    )

  )) %>%
  mutate(TP = map_int(accuracy_df, ~ .x %>% filter(positive & (beta_1 != 0)) %>% nrow())) %>%
  mutate(FP = map_int(accuracy_df, ~ .x %>% filter(positive & (beta_1 == 0)) %>% nrow()))


