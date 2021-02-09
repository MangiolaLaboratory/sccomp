sccomp - outlier-free compositional analysis
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidyseurat/workflows/R-CMD-check/badge.svg)](https://github.com/stemangiola/tidyseurat/actions/)
<!-- badges: end -->

``` r
library(furrr)
```

    ## Loading required package: future

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(sccomp)
library(tidyr)
library(ggplot2)
plan(multisession, workers=10)
```

``` r
res = 
  sccomp::cell_counts %>%
  sccomp_glm(
    formula = ~ type,
    sample, cell_type, count
  )
```

    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.06 seconds (Warm-up)
    ## Chain 1:                0.71 seconds (Sampling)
    ## Chain 1:                10.77 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.92 seconds (Warm-up)
    ## Chain 1:                0.55 seconds (Sampling)
    ## Chain 1:                10.47 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.04 seconds (Warm-up)
    ## Chain 1:                0.52 seconds (Sampling)
    ## Chain 1:                10.56 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.84 seconds (Warm-up)
    ## Chain 1:                0.7 seconds (Sampling)
    ## Chain 1:                10.54 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.01 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 100 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.34 seconds (Warm-up)
    ## Chain 1:                0.68 seconds (Sampling)
    ## Chain 1:                11.02 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.21 seconds (Warm-up)
    ## Chain 1:                0.59 seconds (Sampling)
    ## Chain 1:                10.8 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.7 seconds (Warm-up)
    ## Chain 1:                0.6 seconds (Sampling)
    ## Chain 1:                10.3 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.31 seconds (Warm-up)
    ## Chain 1:                0.58 seconds (Sampling)
    ## Chain 1:                10.89 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.01 seconds (Warm-up)
    ## Chain 1:                0.57 seconds (Sampling)
    ## Chain 1:                10.58 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.15 seconds (Warm-up)
    ## Chain 1:                0.69 seconds (Sampling)
    ## Chain 1:                10.84 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.74 seconds (Warm-up)
    ## Chain 1:                0.39 seconds (Sampling)
    ## Chain 1:                10.13 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.86 seconds (Warm-up)
    ## Chain 1:                0.71 seconds (Sampling)
    ## Chain 1:                10.57 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.84 seconds (Warm-up)
    ## Chain 1:                0.58 seconds (Sampling)
    ## Chain 1:                10.42 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.78 seconds (Warm-up)
    ## Chain 1:                0.42 seconds (Sampling)
    ## Chain 1:                10.2 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.46 seconds (Warm-up)
    ## Chain 1:                0.64 seconds (Sampling)
    ## Chain 1:                10.1 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.53 seconds (Warm-up)
    ## Chain 1:                0.66 seconds (Sampling)
    ## Chain 1:                10.19 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.77 seconds (Warm-up)
    ## Chain 1:                0.5 seconds (Sampling)
    ## Chain 1:                10.27 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.84 seconds (Warm-up)
    ## Chain 1:                0.33 seconds (Sampling)
    ## Chain 1:                10.17 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.88 seconds (Warm-up)
    ## Chain 1:                0.47 seconds (Sampling)
    ## Chain 1:                10.35 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.95 seconds (Warm-up)
    ## Chain 1:                0.56 seconds (Sampling)
    ## Chain 1:                10.51 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.16 seconds (Warm-up)
    ## Chain 1:                0.73 seconds (Sampling)
    ## Chain 1:                10.89 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.94 seconds (Warm-up)
    ## Chain 1:                0.57 seconds (Sampling)
    ## Chain 1:                10.51 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.16 seconds (Warm-up)
    ## Chain 1:                0.71 seconds (Sampling)
    ## Chain 1:                10.87 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.3 seconds (Warm-up)
    ## Chain 1:                0.73 seconds (Sampling)
    ## Chain 1:                11.03 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.99 seconds (Warm-up)
    ## Chain 1:                0.78 seconds (Sampling)
    ## Chain 1:                10.77 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.98 seconds (Warm-up)
    ## Chain 1:                0.56 seconds (Sampling)
    ## Chain 1:                10.54 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.66 seconds (Warm-up)
    ## Chain 1:                0.56 seconds (Sampling)
    ## Chain 1:                10.22 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.91 seconds (Warm-up)
    ## Chain 1:                0.57 seconds (Sampling)
    ## Chain 1:                10.48 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.53 seconds (Warm-up)
    ## Chain 1:                0.7 seconds (Sampling)
    ## Chain 1:                10.23 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.52 seconds (Warm-up)
    ## Chain 1:                0.75 seconds (Sampling)
    ## Chain 1:                10.27 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.24 seconds (Warm-up)
    ## Chain 1:                0.62 seconds (Sampling)
    ## Chain 1:                10.86 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.91 seconds (Warm-up)
    ## Chain 1:                0.62 seconds (Sampling)
    ## Chain 1:                11.53 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.68 seconds (Warm-up)
    ## Chain 1:                0.65 seconds (Sampling)
    ## Chain 1:                11.33 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.76 seconds (Warm-up)
    ## Chain 1:                0.63 seconds (Sampling)
    ## Chain 1:                11.39 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.29 seconds (Warm-up)
    ## Chain 1:                0.9 seconds (Sampling)
    ## Chain 1:                11.19 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.01 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 100 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.41 seconds (Warm-up)
    ## Chain 1:                0.72 seconds (Sampling)
    ## Chain 1:                11.13 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.68 seconds (Warm-up)
    ## Chain 1:                0.57 seconds (Sampling)
    ## Chain 1:                11.25 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.63 seconds (Warm-up)
    ## Chain 1:                0.8 seconds (Sampling)
    ## Chain 1:                11.43 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.47 seconds (Warm-up)
    ## Chain 1:                0.47 seconds (Sampling)
    ## Chain 1:                10.94 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.76 seconds (Warm-up)
    ## Chain 1:                0.69 seconds (Sampling)
    ## Chain 1:                11.45 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.61 seconds (Warm-up)
    ## Chain 1:                0.57 seconds (Sampling)
    ## Chain 1:                10.18 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.95 seconds (Warm-up)
    ## Chain 1:                0.67 seconds (Sampling)
    ## Chain 1:                10.62 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.19 seconds (Warm-up)
    ## Chain 1:                0.78 seconds (Sampling)
    ## Chain 1:                10.97 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.91 seconds (Warm-up)
    ## Chain 1:                0.57 seconds (Sampling)
    ## Chain 1:                10.48 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.28 seconds (Warm-up)
    ## Chain 1:                0.54 seconds (Sampling)
    ## Chain 1:                10.82 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.01 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 100 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.08 seconds (Warm-up)
    ## Chain 1:                0.46 seconds (Sampling)
    ## Chain 1:                10.54 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.72 seconds (Warm-up)
    ## Chain 1:                0.68 seconds (Sampling)
    ## Chain 1:                10.4 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.41 seconds (Warm-up)
    ## Chain 1:                0.61 seconds (Sampling)
    ## Chain 1:                11.02 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.04 seconds (Warm-up)
    ## Chain 1:                0.66 seconds (Sampling)
    ## Chain 1:                10.7 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.94 seconds (Warm-up)
    ## Chain 1:                0.7 seconds (Sampling)
    ## Chain 1:                10.64 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.45 seconds (Warm-up)
    ## Chain 1:                0.56 seconds (Sampling)
    ## Chain 1:                10.01 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.68 seconds (Warm-up)
    ## Chain 1:                0.5 seconds (Sampling)
    ## Chain 1:                10.18 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.24 seconds (Warm-up)
    ## Chain 1:                0.58 seconds (Sampling)
    ## Chain 1:                10.82 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.7 seconds (Warm-up)
    ## Chain 1:                0.72 seconds (Sampling)
    ## Chain 1:                10.42 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.95 seconds (Warm-up)
    ## Chain 1:                0.5 seconds (Sampling)
    ## Chain 1:                10.45 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.34 seconds (Warm-up)
    ## Chain 1:                0.58 seconds (Sampling)
    ## Chain 1:                10.92 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.74 seconds (Warm-up)
    ## Chain 1:                0.4 seconds (Sampling)
    ## Chain 1:                10.14 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.52 seconds (Warm-up)
    ## Chain 1:                0.35 seconds (Sampling)
    ## Chain 1:                10.87 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.01 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 100 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.09 seconds (Warm-up)
    ## Chain 1:                0.55 seconds (Sampling)
    ## Chain 1:                10.64 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.58 seconds (Warm-up)
    ## Chain 1:                0.42 seconds (Sampling)
    ## Chain 1:                10 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.31 seconds (Warm-up)
    ## Chain 1:                0.68 seconds (Sampling)
    ## Chain 1:                10.99 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.26 seconds (Warm-up)
    ## Chain 1:                0.38 seconds (Sampling)
    ## Chain 1:                10.64 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.29 seconds (Warm-up)
    ## Chain 1:                0.53 seconds (Sampling)
    ## Chain 1:                10.82 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.74 seconds (Warm-up)
    ## Chain 1:                0.65 seconds (Sampling)
    ## Chain 1:                11.39 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.66 seconds (Warm-up)
    ## Chain 1:                0.57 seconds (Sampling)
    ## Chain 1:                11.23 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.16 seconds (Warm-up)
    ## Chain 1:                0.63 seconds (Sampling)
    ## Chain 1:                10.79 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.07 seconds (Warm-up)
    ## Chain 1:                0.83 seconds (Sampling)
    ## Chain 1:                10.9 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.37 seconds (Warm-up)
    ## Chain 1:                0.58 seconds (Sampling)
    ## Chain 1:                10.95 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 11.01 seconds (Warm-up)
    ## Chain 1:                0.75 seconds (Sampling)
    ## Chain 1:                11.76 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.5 seconds (Warm-up)
    ## Chain 1:                0.66 seconds (Sampling)
    ## Chain 1:                11.16 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.83 seconds (Warm-up)
    ## Chain 1:                0.72 seconds (Sampling)
    ## Chain 1:                10.55 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.69 seconds (Warm-up)
    ## Chain 1:                0.59 seconds (Sampling)
    ## Chain 1:                11.28 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.67 seconds (Warm-up)
    ## Chain 1:                0.73 seconds (Sampling)
    ## Chain 1:                10.4 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.32 seconds (Warm-up)
    ## Chain 1:                0.65 seconds (Sampling)
    ## Chain 1:                10.97 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.36 seconds (Warm-up)
    ## Chain 1:                0.62 seconds (Sampling)
    ## Chain 1:                10.98 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.93 seconds (Warm-up)
    ## Chain 1:                0.43 seconds (Sampling)
    ## Chain 1:                10.36 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.07 seconds (Warm-up)
    ## Chain 1:                0.6 seconds (Sampling)
    ## Chain 1:                10.67 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.14 seconds (Warm-up)
    ## Chain 1:                0.6 seconds (Sampling)
    ## Chain 1:                10.74 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.78 seconds (Warm-up)
    ## Chain 1:                0.33 seconds (Sampling)
    ## Chain 1:                10.11 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.55 seconds (Warm-up)
    ## Chain 1:                0.55 seconds (Sampling)
    ## Chain 1:                10.1 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.01 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 100 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.88 seconds (Warm-up)
    ## Chain 1:                0.62 seconds (Sampling)
    ## Chain 1:                10.5 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.78 seconds (Warm-up)
    ## Chain 1:                0.73 seconds (Sampling)
    ## Chain 1:                10.51 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.62 seconds (Warm-up)
    ## Chain 1:                0.56 seconds (Sampling)
    ## Chain 1:                10.18 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.52 seconds (Warm-up)
    ## Chain 1:                0.64 seconds (Sampling)
    ## Chain 1:                10.16 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.76 seconds (Warm-up)
    ## Chain 1:                0.66 seconds (Sampling)
    ## Chain 1:                10.42 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.74 seconds (Warm-up)
    ## Chain 1:                0.6 seconds (Sampling)
    ## Chain 1:                10.34 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.08 seconds (Warm-up)
    ## Chain 1:                0.51 seconds (Sampling)
    ## Chain 1:                10.59 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.56 seconds (Warm-up)
    ## Chain 1:                0.59 seconds (Sampling)
    ## Chain 1:                10.15 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.94 seconds (Warm-up)
    ## Chain 1:                0.59 seconds (Sampling)
    ## Chain 1:                10.53 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.02 seconds (Warm-up)
    ## Chain 1:                0.65 seconds (Sampling)
    ## Chain 1:                10.67 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.86 seconds (Warm-up)
    ## Chain 1:                0.74 seconds (Sampling)
    ## Chain 1:                10.6 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.37 seconds (Warm-up)
    ## Chain 1:                0.68 seconds (Sampling)
    ## Chain 1:                11.05 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.9 seconds (Warm-up)
    ## Chain 1:                0.59 seconds (Sampling)
    ## Chain 1:                10.49 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.73 seconds (Warm-up)
    ## Chain 1:                0.79 seconds (Sampling)
    ## Chain 1:                10.52 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.03 seconds (Warm-up)
    ## Chain 1:                0.67 seconds (Sampling)
    ## Chain 1:                10.7 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.91 seconds (Warm-up)
    ## Chain 1:                0.61 seconds (Sampling)
    ## Chain 1:                10.52 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.94 seconds (Warm-up)
    ## Chain 1:                0.65 seconds (Sampling)
    ## Chain 1:                10.59 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.58 seconds (Warm-up)
    ## Chain 1:                0.58 seconds (Sampling)
    ## Chain 1:                11.16 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 9.54 seconds (Warm-up)
    ## Chain 1:                0.64 seconds (Sampling)
    ## Chain 1:                10.18 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_dirichlet_multinomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 1050 [  0%]  (Warmup)
    ## Chain 1: Iteration:  105 / 1050 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  210 / 1050 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  315 / 1050 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  420 / 1050 [ 40%]  (Warmup)
    ## Chain 1: Iteration:  525 / 1050 [ 50%]  (Warmup)
    ## Chain 1: Iteration:  630 / 1050 [ 60%]  (Warmup)
    ## Chain 1: Iteration:  735 / 1050 [ 70%]  (Warmup)
    ## Chain 1: Iteration:  840 / 1050 [ 80%]  (Warmup)
    ## Chain 1: Iteration:  945 / 1050 [ 90%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 1050 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 1050 / 1050 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 10.2 seconds (Warm-up)
    ## Chain 1:                0.58 seconds (Sampling)
    ## Chain 1:                10.78 seconds (Total)
    ## Chain 1:

These are the cell\_types for which inference was biased by the presence
of outliers

``` r
res %>%
  filter(significant != significant_pre_filtering) %>%
  distinct(cell_type, phenotype)
```

    ## # A tibble: 3 x 2
    ##   phenotype                                     cell_type
    ##   <fct>                                         <fct>    
    ## 1 Monocyte:CD14+_IL1B_CXCL8_inflammatory        M2       
    ## 2 T_cell:CD4+_central_memory_JUN_FOS_JUNB_DUSP1 CD4 3    
    ## 3 T_cell:CD8+_non_activated                     CD8 3

``` r
data_for_plot = 
  res %>%
  rename(lower = `2.5%`, upper = `97.5%`, median=`50%`) %>%
  select(cell_type, lower, upper, median, quantiles_pre_filtering, significant, significant_pre_filtering) %>%
  unnest(quantiles_pre_filtering) %>%
  distinct() 
```

``` r
data_for_plot %>%
  ggplot(aes(x=`50%`, y=cell_type)) +
  geom_vline(xintercept = 0, colour="grey") +
  geom_errorbar(aes(xmin=`2.5%`, xmax=`97.5%`, color=significant_pre_filtering)) +
  geom_point() +
  theme_bw() +
  xlab("Credible interval slope") +
  ggtitle("Before outlier filtering")
```

![](man/figures/unnamed-chunk-7-1.png)<!-- -->

``` r
data_for_plot %>%
  ggplot(aes(x=median, y=cell_type)) +
  geom_vline(xintercept = 0, colour="grey") +
  geom_errorbar(aes(xmin=lower, xmax=upper, color=significant)) +
  geom_point() +
  theme_bw() +
  xlab("Credible interval slope") +
  ggtitle("After outlier filtering")
```

![](man/figures/unnamed-chunk-8-1.png)<!-- -->
