library(tidyverse)
library(extraDistr)
library(TailRank)

prop_logit = 0
sd = 1
draws = 20000
p = 0.995

x =
  tibble(
  mu_logit = rnorm(draws, prop_logit, sd),
  precision_log = rnorm(draws, 4, 1)
) %>%
  mutate(
    mu = boot::inv.logit(mu_logit),
    precision = exp(precision_log)
  ) %>%
tibble(
  alpha = mu*precision,
  beta = (1-mu)*precision,
  N = 500
) %>%

  # Ground truth
  mutate(r = pmap_int(
    list(alpha, beta, N),
    ~ rbb(1,  N = ..3,u =  ..1 ,v = ..2 )

  )) %>%
  mutate(quantile(r, p)) %>%

# Quantile method
  mutate(q = pmap_int(
    list(alpha, beta, N),
    ~ qbb(p = p,  N =..3, u =  ..1 ,v = ..2 )

  )) %>%
  mutate(mean(q)) %>%

    # Hypersampling method
  mutate(r2 = pmap(
    list(alpha, beta, N),
    ~ rbb(10,  N = ..3, u =..1 , v = ..2 )

  )) %>%
  unnest(r2) %>%
  mutate(quantile(r2, p))


sampling = sample(1:nrow(x), 20000, replace = TRUE)
rbb(length(sampling),  N = 500, u =x$alpha[sampling] , v = x$beta[sampling] ) %>%
  quantile(p)
