library(sccomp)
#library(tidyverse)

library(dplyr)
library(tidyr)
library(purrr)


library(magrittr)
library(tidybulk)



# Read arguments
args = commandArgs(trailingOnly=TRUE)
slope = as.numeric(args[1])
n_samples = as.integer(args[2])
n_cell_type = as.integer(args[3])
max_cell_counts_per_sample = as.integer(args[4])
add_outliers = as.integer(args[5])
output_file = args[6]

outlier_probability = 0.1


# exposures = counts_obj %>% group_by(sample) %>% summarise(s=sum(count)) %>% pull(s) %>% sort %>% head(-1)
beta_0 = readRDS("dev/beta_0.rds")

is_odd = function(x) { ( seq_len(length(x)) %% 2 ) == 1 }

n_differentially_abundant = ceiling(n_cell_type*0.4)

set_factor_of_interest = function(n){
  c(rep(0, ceiling(n/2)), rep(1, ceiling(n/2))) %>%  sample(size = n) %>% .[1:n]
}

# Iterate over runs
tibble(run = 1:50) %>%
  mutate(data = map(
    run,
    ~ {
      input_data =
        expand_grid(
          sample = 1:n_samples, cell_type = 1:n_cell_type
        ) %>%
        nest(d = -cell_type) %>%
        mutate(beta_0 = 0 ) %>%  # sample(beta_0, size = n())) %>% # rnorm(n = n(), 0, 1)) %>%
        mutate(  beta_1 = case_when(
          cell_type %in% seq_len(n_differentially_abundant)[is_odd(seq_len(n_differentially_abundant))] ~ 1,
          cell_type %in% seq_len(n_differentially_abundant)[!is_odd(seq_len(n_differentially_abundant))]  ~ -1,
          TRUE ~ 0
        )) %>%
        mutate(beta_1 = beta_1 %>% multiply_by(slope)) %>%
        nest(coefficients = starts_with("beta_")) %>%
        unnest(d) %>%
        nest(d = -sample) %>%
        mutate(type = set_factor_of_interest(n())) %>%
        mutate(tot_count = sample(400:1000, size = n(), replace = TRUE)) %>%
        unnest(d) %>%
        mutate(sample = as.character(sample), cell_type = as.character(cell_type)) %>%
        unnest(coefficients)

      my_simulated_data =
        simulate_data(input_data,
                      readRDS("dev/oligo_breast_estimate.rds") %>%
                        tidybulk:::add_attr("noise_model", "logit_normal_multinomial"), # Use logit normal
                      formula = ~ type ,
                      .sample = sample,
                      .cell_group = cell_type,
                      .coefficients = c(beta_0, beta_1),
                      mcmc_seed = .x * 2
        )

      # Add outliers
      if(add_outliers==1){

        # Add multipliers
        ratio_changing =
          my_simulated_data %>%
          filter((beta_1<0 & type ==1) | (beta_1>0 & type ==0)) %>%

          # 0.2 because I am just taking half of the samples and I double the outlier probability
          sample_frac(outlier_probability*2) %>%
          mutate(ratio =  rnorm(n(), 7.6, 2.9)) %>%
          select(sample, cell_type, ratio)

        ratio_NON_changing =
          my_simulated_data %>%
          filter(beta_1==0) %>%
          sample_frac(outlier_probability) %>%
          mutate(ratio =  rnorm(ceiling(n()/2), 7.6, 2.9) %>% c(rnorm(floor(n()/2), 0.03572171, 0.01107917)) ) %>%
          ungroup %>%
          select(sample, cell_type, ratio)

        my_simulated_data %>%
          left_join(
            bind_rows(ratio_changing, ratio_NON_changing),
            by = c("sample", "cell_type")
          ) %>%
          replace_na(list(ratio = 1)) %>%
          mutate(generated_counts = (generated_counts * ratio) %>% as.integer()) %>%
          rowwise() %>%
          mutate(generated_counts = max(0, generated_counts) %>% as.integer()) %>%
          ungroup()
      }

      else {
        my_simulated_data
      }


    }
  ))  %>%

  saveRDS(output_file)



# Plot
#

# my_simulated_data2 |>
#   left_join(input_data |> distinct(type, sample)) %>%
#   ggplot(aes(factor(type), generated_proportions, fill=type)) +
#
#   geom_boxplot(
#     aes( fill=type),
#     outlier.shape = NA, alpha=0.2
#   ) +
#   geom_jitter(
#     aes( color=type) ,
#     alpha=0.2, size = 0.6,
#     data = simulated_proportion
#   ) +
#
#   facet_wrap(~ interaction(cell_type), scale="free_y") +
#   scale_y_continuous(trans="logit") +
#
#   xlab("Biological condition") +
#   ylab("Cell-group proportion") +
#   theme_bw() +
#   theme(strip.background =element_rect(fill="white"), legend.position = "bottom")
