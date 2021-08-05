library(sccomp)
library(tidyverse)
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

# exposures = counts_obj %>% group_by(sample) %>% summarise(s=sum(count)) %>% pull(s) %>% sort %>% head(-1)
beta_0 = readRDS("dev/beta_0.rds")

is_odd = function(x) { ( seq_len(length(x)) %% 2 ) == 1 }

n_differentially_abundant = ceiling(n_cell_type*0.4)

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
        mutate(beta_0 = sample(beta_0, size = n())) %>% # rnorm(n = n(), 0, 1)) %>%
        mutate(  beta_1 = case_when(
          cell_type %in% seq_len(n_differentially_abundant)[is_odd(seq_len(n_differentially_abundant))] ~ 1,
          cell_type %in% seq_len(n_differentially_abundant)[!is_odd(seq_len(n_differentially_abundant))]  ~ -1,
          TRUE ~ 0
        )) %>%
        mutate(beta_1 = beta_1 %>% multiply_by(slope)) %>%
        nest(coefficients = starts_with("beta_")) %>%
        unnest(d) %>%
        nest(d = -sample) %>%
        mutate(type = sample(c(0,1), size = n(), replace = TRUE)) %>%
        mutate(tot_count = sample(400:1000, size = n(), replace = TRUE)) %>%
        unnest(d) %>%
        mutate(sample = as.character(sample), cell_type = as.character(cell_type))

      my_simulated_data =
        simulate_data(input_data,
                      formula = ~ type ,
                      sample,
                      cell_type, tot_count, coefficients,
                      seed = .x * 2
        )

      # Add outliers
      if(add_outliers==1){

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

      else {
        my_simulated_data
      }


    }
  ))  %>%

  saveRDS(output_file)

