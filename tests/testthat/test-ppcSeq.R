context('ppcseq')

test_that("first test",{

  # library(furrr)
  # plan(multisession, workers=20)
  library(dplyr)
  library(sccomp)
  library(digest)
  #debugonce(sccomp_glm)

  res =
    sccomp::cell_counts %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type, count
    )

  expect_equal(
    res %>%
      distinct(cell_type, significant) %>%
      pull(significant) %>%
      digest(algo="md5"),
    "f6cef772af43198f586e15c96b2f1239"
  )

})





library(rstan)
library(gtools)
library(tidyverse)
library(cellsig)
library(magrittr)
library(tidybayes)
PFI_all_cancers = readRDS("~/PostDoc/supervision/yuhan/BLCA_IL2NK/dev/test_simulation_makeflow_pipeline/PFI_all_cancers.rds")


S = 100
slope = 1

# beta = cellsig:::get_alpha(slope = slope, 1, 1:10)
beta = matrix(c(c(0,0,0,0,0,0,1,1,1,1), c(1,0,0,0,0,0,0,0,0,1)), ncol=2)

X = cellsig:::get_survival_X(S, PFI_all_cancers %>% mutate(PFI.time.2 = sqrt(PFI.time.2)), center = TRUE)[,c("intercept", "real_days")] %>% as.matrix()

logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}
softmax <- function (x) {
  exp(x - logsumexp(x))
}

y =
  X %*% t(beta) %>%
  apply(1, softmax) %>%
  t %>%
  `*` (120) %>%
  as.data.frame() %>%
  apply(1, function(x) gtools::rdirichlet(1, x)) %>%
  t()

y %>%
  as_tibble() %>%
  rowid_to_column("n") %>%
  gather(m, p, -n) %>%
  left_join(X %>% as_tibble() %>% rowid_to_column("n")) %>%
  ggplot(aes(real_days, p, color=m)) +
  geom_line()

f = stan(file = "inst/stan/glm_multi_beta.stan",
     data = list(
       N = S,
       M = 10,
       y = y,
       X = X
     ),
     cores = 4
    )

f %>% plot(pars="beta")
f %>% plot(pars="precision")





# Generate quantities
q = stan(file = "inst/stan/glm_multi_beta_generate_data.stan",
         data = list(
           N = S,
           M = 10,
           C = 2,
           X = X,
           beta = beta %>% t(),
           precision = c(30, 30, 30, 30, 30, 30, 5, 5, 5, 5)
         ), algorithm = "Fixed_param")

y2 = q %>% rstan::extract("y") %$% y %>% .[1,,] %>% apply(1, function(x) x/sum(x)) %>% t()
y2 = q %>% gather_draws(y[N, M]) %>% filter(.draw < 100) %>% mean_qi() %>% select(N, M, .value) %>% spread(M, .value) %>% nanny::as_matrix(rownames = N)  %>% apply(1, function(x) x/sum(x)) %>% t()

y2 %>%
  as_tibble() %>%
  rowid_to_column("n") %>%
  gather(m, p, -n) %>%
  left_join(X %>% as_tibble() %>% rowid_to_column("n")) %>%
  ggplot(aes(real_days, p, color=m)) +
  geom_line()

f2 = stan(file = "inst/stan/glm_multi_beta.stan",
         data = list(
           N = S,
           M = 10,
           y = y2,
           X = X
         ),
         cores = 4
)

f2 %>% plot(pars="beta")
f2 %>% plot(pars="precision")

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
  mutate(frac = (count/sample_tot) ) %>%
  select(-c(count, sample_tot, phenotype)) %>%
  spread(cell_type, frac)

f3 = stan(file = "inst/stan/glm_multi_beta.stan",
          data = list(
            N = 20,
            M = 36,
            y = input_df %>% select(-type) %>% nanny::as_matrix(rownames = sample),
            X = input_df %>% select(sample, type) %>% model.matrix(~ type, data=.)
          ),
          cores = 4
)

f3 %>% plot(pars="beta")
f3 %>% plot(pars="precision")


f3 %>%
  spread_draws(precision[M], beta[C,M]) %>%
  filter(C==1) %>%
  mean_qi() %>%
  ggplot(aes(beta, precision)) +
  geom_errorbar(aes(ymin = precision.lower, ymax = precision.upper)) +
  geom_point() +
  geom_smooth(method="lm")
