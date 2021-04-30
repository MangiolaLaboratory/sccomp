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
PFI_all_cancers = readRDS("~/PostDoc/supervision/yuhan/BLCA_IL2NK/dev/test_simulation_makeflow_pipeline/PFI_all_cancers.rds")


S = 20
slope = 1

beta = cellsig:::get_alpha(slope = slope, 1, 1:10)
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
  `*` (80) %>%
  as.data.frame() %>%
  apply(1, function(x) gtools::rdirichlet(1, x)) %>%
  t()


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

y %>%
  as_tibble() %>%
  rowid_to_column("n") %>%
  gather(m, p, -n) %>%
  left_join(X %>% as_tibble() %>% rowid_to_column("n")) %>%
  ggplot(aes(real_days, p, color=m)) +
  geom_line()



# Generate quantities
q = stan(file = "inst/stan/glm_multi_beta_generate_data.stan",
         data = list(
           N = S,
           M = 10,
           C = 2,
           X = X,
           beta = beta %>% t(),
           precision = c(5, rep(30, 9))
         ), algorithm = "Fixed_param")

y2 = q %>% rstan::extract("y") %$% y %>% .[1,,] %>% apply(1, function(x) x/sum(x)) %>% t()

y2 %>%
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
           y = y2,
           X = X
         ),
         cores = 4
)

f %>% plot(pars="beta")
f %>% plot(pars="precision")
