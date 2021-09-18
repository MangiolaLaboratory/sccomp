library(sccomp)
library(tidyverse)
library(magrittr)
library(tidybulk)
library(speckle)
library(edgeR)
library(limma)
library(DirichletReg)
# exposures = counts_obj %>% group_by(sample) %>% summarise(s=sum(count)) %>% pull(s) %>% sort %>% head(-1)
beta_0 = readRDS("dev/beta_0.rds")




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
# probs = seq(0, 0.2,length.out = 20)

outliers =
  readRDS("dev/oligo_estimate_wth_benign.rds") %>%
  unnest(outliers) %>%
  filter(outlier) %>%
  mutate(ratio = (count+1)/.median)

outliers %>% filter(ratio>1) %>% pull(ratio) %>% mean
outliers %>% filter(ratio>1) %>% pull(ratio) %>% sd
outliers %>% filter(ratio<1) %>% pull(ratio) %>% mean
outliers %>% filter(ratio<1) %>% pull(ratio) %>% sd

# Iterate over runs
benchmark =
  tibble(run = 1:50) %>%
  mutate(data = map(
    run,
    ~ {
      input_data =
        expand_grid(
          sample = 1:10, cell_type = 1:10
        ) %>%
        nest(d = -cell_type) %>%
        mutate(beta_0 = sample(beta_0, size = n())) %>% # rnorm(n = n(), 0, 1)) %>%
        mutate(  beta_1 = case_when(
          cell_type %in% c(1, 3) ~ 1,
          cell_type %in% c(2, 4) ~ -1,
          TRUE ~ 0
        )) %>%
        mutate(beta_1 = beta_1 %>% multiply_by(0.8)) %>%
        nest(coefficients = starts_with("beta_")) %>%
        unnest(d) %>%
        nest(d = -sample) %>%
        mutate(type = sample(c(0,1), size = n(), replace = TRUE)) %>%
        mutate(tot_count = sample(600:1000, size = n(), replace = TRUE)) %>%
        unnest(d) %>%
        mutate(sample = as.character(sample), cell_type = as.character(cell_type))

      my_simulated_data =
        simulate_data(input_data,
                    formula = ~ type ,
                    sample,
                    cell_type, tot_count, coefficients,
                    seed = .x * 2
      )

      # Add multipliers
      ratio_changing =
        my_simulated_data %>%
        unnest(coefficients) %>%
        filter((beta_1<0 & type ==1) | (beta_1>0 & type ==0)) %>%
        group_by(cell_type) %>%
        sample_n(1) %>%
        mutate(ratio =  rnorm(n(), 7.6, 2.9)) %>%
        ungroup %>%
        select(sample, cell_type, ratio)

      ratio_NON_changing =
        my_simulated_data %>%
        unnest(coefficients) %>%
        filter(beta_1==0) %>%
        sample_frac(0.1) %>%
        mutate(ratio =  rnorm(n()/2, 7.6, 2.9) %>% c(rnorm(n()/2, 0.03572171, 0.01107917)) ) %>%
        ungroup %>%
        select(sample, cell_type, ratio)

      my_simulated_data %>%
        left_join(
          bind_rows(ratio_changing, ratio_NON_changing),
          by = c("sample", "cell_type")
        ) %>%
        replace_na(list(ratio = 1)) %>%
        mutate(.value = (.value * ratio) %>% as.integer()) %>%
        rowwise() %>%
        mutate(.value = max(0, .value)) %>%
        ungroup() %>%
        select(-ratio)


    }
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
      # group_by(sample) %>%
      # mutate(proportion = (.value+1)/sum(.value+1)) %>%
      # ungroup(sample) %>%
      # mutate(rate = proportion %>% boot::logit()) %>%
      # mutate(rate = rate - min(rate)) %>%
      # mutate(counts = rate %>% exp()) %>%
      test_differential_abundance(
        ~ type,
        sample, cell_type, .value,
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
      sccomp_glm(
        ~type,
        sample, cell_type, .value,
        check_outliers = TRUE,
        approximate_posterior_inference = FALSE,
        percent_false_positive = 0.001 * 100,
        seed = .y * 2
      )
  )) %>%
  mutate( results_DirichletMultinomial = map2(
    data, run,
    ~  .x %>%
      mutate(.value = as.integer(.value)) %>%
      sccomp_glm(
        ~type,
        sample, cell_type, .value,
        check_outliers = FALSE,
        approximate_posterior_inference = FALSE,
        percent_false_positive = 0.001 * 100,
        noise_model="dirichlet_multinomial",
        seed = .y * 2
      )
  ))


benchmark %>%
  pull(data) %>%
  .[[1]] %>%
  group_by(sample) %>%
  mutate(proportion = (.value+1)/sum(.value+1)) %>%
  ungroup(sample) %>%
  ggplot(aes(factor(type), proportion, fill = cell_type %in% 1:4)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.5) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  facet_wrap(~ as.integer(cell_type), ncol=5) +
  theme_bw() +
  theme(
    strip.background =element_rect(fill="white", color="white"),
    legend.position = "none",
  )

ggsave(
  "dev/example_boxplot_yes_outliers.pdf",
  units = c("mm"),
  width = 183/2 ,
  height = 183/2,
  limitsize = FALSE
)

probs = seq(0, 0.1,length.out = 50) %>% c(seq(0.1, 1,length.out = 50))

benchmark_hypothesis =
  benchmark %>%
  dplyr::mutate(probs  = map(run, ~ !!probs)) %>%
  unnest((probs) ) %>%
  mutate(hypothesis_edger = map2(results_edger, (probs) , ~ .x %>% arrange(FDR) %>% mutate(positive = FDR<.y) %>% mutate(trend = logFC))) %>%
  mutate(hypothesis_voom = map2(results_voom, (probs) , ~.x %>% arrange(adj.P.Val) %>% mutate(positive = adj.P.Val<.y) %>% mutate(trend = logFC))) %>%
  mutate(hypothesis_speckle = map2(results_speckle, (probs) , ~ .x %>% arrange(FDR) %>% mutate(positive = FDR<.y) %>% mutate(trend = -Tstatistic    ))) %>%
  mutate(hypothesis_sccomp = map2(
    results_sccomp, (probs) ,
    ~ .x  %>%
      arrange(false_discovery_rate) %>%
      mutate(positive = false_discovery_rate<.y) %>%
      mutate(trend = .median_type )
  )) %>%
  mutate(hypothesis_DirichletMultinomial  = map2(
    results_DirichletMultinomial , (probs) ,
    ~ .x  %>%
      arrange(false_discovery_rate) %>%
      mutate(positive = false_discovery_rate<(.y/10)) %>%
      mutate(trend = .median_type )
  )) %>%
  dplyr::select(-contains("results"))



benchmark_hypothesis %>%
  dplyr::select(-contains("results")) %>%
  pivot_longer(contains("hypothesis"),names_prefix = "hypothesis_" ) %>%
  mutate(accuracy_df = map2(
    data, value             ,
    ~ left_join(

      .x %>% unnest(coefficients) %>% dplyr::select(cell_type, beta_1) %>% distinct %>% mutate(cell_type = as.character(cell_type)),
      .y %>% dplyr::select(cell_type, positive, trend),
      by="cell_type"

    )

  )) %>%
  mutate(TP = map_int(accuracy_df, ~ .x %>% filter(positive & (beta_1 != 0) & (beta_1 * trend)>0) %>% nrow())) %>%
  mutate(FP = map_int(accuracy_df, ~ .x %>% filter(

    # Positive when not
    (positive & (beta_1 == 0)) |

    # Positive when yes but wrong direction
    (positive & (beta_1 != 0) & (beta_1 * trend)<0)
  ) %>% nrow())) %>%
  mutate(total_true_positive = 6, total_true_negative = 24-6) %>%
  mutate(FP_rate = FP/total_true_negative) %>%
  mutate(TP_rate = TP/total_true_positive) %>%
  group_by(name, probs) %>%
  summarise(mean_FP_rate = mean(FP_rate), mean_TP_rate = mean(TP_rate)) %>%
  ungroup() %>%
  arrange(mean_FP_rate) %>%
  rename(Algorithm = name) %>%
  distinct() %>%
  ggplot(aes(mean_FP_rate, mean_TP_rate)) +
  geom_line(aes(color = Algorithm)) +
  scale_x_continuous(limits =c(0,0.10)) +
  scale_color_brewer(palette="Set1") +
  ylab("Mean true-positive proportion") +
  xlab("Mean false-positive proportion") +
  theme_bw() +
  theme(legend.position = "bottom")

saveRDS(benchmark, "dev/benchmark_slope0.8_samples10_celltypes10_max1000exposure_50replicates_outliers.rds")


ggsave(
  "dev/roc_slope0.8_samples10_celltypes10_max1000exposure_50replicates_outliers.pdf",
  units = c("mm"),
  width = 183/2 ,
  height = 183/2,
  limitsize = FALSE
)


