# plot

This function plots a summary of the results of the model.

## Usage

``` r
# S3 method for class 'sccomp_tbl'
plot(
  x,
  significance_threshold = 0.05,
  test_composition_above_logit_fold_change = attr(.data,
    "test_composition_above_logit_fold_change"),
  significance_statistic = c("FDR", "pH0"),
  show_fdr_message = TRUE,
  ...
)
```

## Arguments

- x:

  A tibble including a cell_group name column \| sample name column \|
  read counts column \| factor columns \| Pvalue column \| a
  significance column

- significance_threshold:

  Numeric value specifying the significance threshold for highlighting
  differences. Default is 0.025.

- test_composition_above_logit_fold_change:

  A positive integer. It is the effect threshold used for the hypothesis
  test. A value of 0.2 correspond to a change in cell proportion of 10%
  for a cell type with baseline proportion of 50%. That is, a cell type
  goes from 45% to 50%. When the baseline proportion is closer to 0 or 1
  this effect thrshold has consistent value in the logit uncontrained
  scale.

- significance_statistic:

  Character vector indicating which statistic to highlight. Default is
  "FDR".

- show_fdr_message:

  Logical. Whether to show the Bayesian FDR interpretation message on
  the plot. Default is TRUE.

- ...:

  For internal use

## Value

A `ggplot`

## References

S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z.
Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust
differential composition and variability analysis for single-cell data,
Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120,
https://doi.org/10.1073/pnas.2203828120 (2023).

## Examples

``` r
print("cmdstanr is needed to run this example.")
#> [1] "cmdstanr is needed to run this example."
# Note: Before running the example, ensure that the 'cmdstanr' package is installed:
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))

# \donttest{
  if (instantiate::stan_cmdstan_exists()) {
    data("counts_obj")

    estimate = sccomp_estimate(
      counts_obj,
      ~ type, ~1, "sample", "cell_group", "count",
      cores = 1
    )

    # estimate |> plot()
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Precompiled model not found. Compiling the model...
#> Running make /tmp/Rtmputl3kG/model-38122eab5363 "STAN_THREADS=TRUE" \
#>   "STANCFLAGS += --include-paths=/tmp/Rtmputl3kG/temp_libpath3812d1f3955/sccomp/stan --name='glm_multi_beta_binomial_model'"
#> 
#> --- Translating Stan model to C++ code ---
#> bin/stanc --include-paths=/tmp/Rtmputl3kG/temp_libpath3812d1f3955/sccomp/stan --name='glm_multi_beta_binomial_model' --o=/tmp/Rtmputl3kG/model-38122eab5363.hpp /tmp/Rtmputl3kG/model-38122eab5363.stan
#> 
#> --- Compiling C++ code ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -I stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -I stan/lib/stan_math/lib/eigen_3.4.0 -I stan/lib/stan_math/lib/boost_1.87.0 -I stan/lib/stan_math/lib/sundials_6.1.1/include -I stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS          -c -Wno-ignored-attributes   -x c++ -o /tmp/Rtmputl3kG/model-38122eab5363.o /tmp/Rtmputl3kG/model-38122eab5363.hpp
#> 
#> --- Linking model ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -I stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -I stan/lib/stan_math/lib/eigen_3.4.0 -I stan/lib/stan_math/lib/boost_1.87.0 -I stan/lib/stan_math/lib/sundials_6.1.1/include -I stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS               -Wl,-L,"/home/runner/.cmdstan/cmdstan-2.38.0/stan/lib/stan_math/lib/tbb"   -Wl,-rpath,"/home/runner/.cmdstan/cmdstan-2.38.0/stan/lib/stan_math/lib/tbb"      /tmp/Rtmputl3kG/model-38122eab5363.o src/cmdstan/main_threads.o       -ltbb   stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_nvecserial.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_cvodes.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_idas.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_kinsol.a  stan/lib/stan_math/lib/tbb/libtbb.so.2 -o /tmp/Rtmputl3kG/model-38122eab5363
#> rm /tmp/Rtmputl3kG/model-38122eab5363.hpp /tmp/Rtmputl3kG/model-38122eab5363.o
#> Model compiled and saved to cache successfully.
#> Path [1] :Initial log joint density = -481591.693540 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.548e-02   3.349e-01    1.000e+00  1.000e+00      3299 -3.704e+03 -3.704e+03                   
#> Path [1] :Best Iter: [56] ELBO (-3703.871840) evaluations: (3299) 
#> Path [2] :Initial log joint density = -481748.738466 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.389e-03   2.295e-01    1.000e+00  1.000e+00      3065 -3.707e+03 -3.703e+03                   
#> Path [2] :Best Iter: [55] ELBO (-3702.729221) evaluations: (3065) 
#> Path [3] :Initial log joint density = -482597.964046 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.082e-02   3.605e-01    1.000e+00  1.000e+00      3261 -3.699e+03 -3.708e+03                   
#> Path [3] :Best Iter: [55] ELBO (-3698.692066) evaluations: (3261) 
#> Path [4] :Initial log joint density = -481746.681762 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.409e-03   2.139e-01    1.000e+00  1.000e+00      3189 -3.709e+03 -3.703e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3702.905641) evaluations: (3189) 
#> Path [5] :Initial log joint density = -481578.137437 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.587e-03   1.834e-01    9.227e-01  9.227e-01      2950 -3.705e+03 -3.719e+03                   
#> Path [5] :Best Iter: [47] ELBO (-3705.339247) evaluations: (2950) 
#> Path [6] :Initial log joint density = -483875.477861 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.098e-02   3.726e-01    1.000e+00  1.000e+00      3596 -3.698e+03 -3.707e+03                   
#> Path [6] :Best Iter: [55] ELBO (-3697.786841) evaluations: (3596) 
#> Path [7] :Initial log joint density = -481670.678360 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.853e-03   2.302e-01    6.237e-01  6.237e-01      3260 -3.705e+03 -3.709e+03                   
#> Path [7] :Best Iter: [56] ELBO (-3705.059536) evaluations: (3260) 
#> Path [8] :Initial log joint density = -481871.158637 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.185e-02   3.000e-01    1.000e+00  1.000e+00      3572 -3.699e+03 -3.704e+03                   
#> Path [8] :Best Iter: [60] ELBO (-3699.361472) evaluations: (3572) 
#> Path [9] :Initial log joint density = -481562.755362 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.109e-02   2.717e-01    5.129e-01  1.000e+00      3115 -3.707e+03 -3.716e+03                   
#> Path [9] :Best Iter: [39] ELBO (-3706.517543) evaluations: (3115) 
#> Path [10] :Initial log joint density = -481670.762014 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.158e-03   2.920e-01    1.000e+00  1.000e+00      2911 -3.708e+03 -3.716e+03                   
#> Path [10] :Best Iter: [48] ELBO (-3708.281679) evaluations: (2911) 
#> Path [11] :Initial log joint density = -485447.340233 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.084e-03   2.467e-01    1.000e+00  1.000e+00      3214 -3.709e+03 -3.707e+03                   
#> Path [11] :Best Iter: [55] ELBO (-3706.677840) evaluations: (3214) 
#> Path [12] :Initial log joint density = -481434.559358 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      4.558e-03   2.447e-01    1.000e+00  1.000e+00      3008 -3.708e+03 -3.724e+03                   
#> Path [12] :Best Iter: [46] ELBO (-3708.379004) evaluations: (3008) 
#> Path [13] :Initial log joint density = -484087.061783 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      5.234e-03   1.740e-01    1.000e+00  1.000e+00      3416 -3.704e+03 -3.708e+03                   
#> Path [13] :Best Iter: [55] ELBO (-3704.005635) evaluations: (3416) 
#> Path [14] :Initial log joint density = -481678.267564 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      5.423e-03   1.952e-01    1.000e+00  1.000e+00      2831 -3.708e+03 -3.711e+03                   
#> Path [14] :Best Iter: [43] ELBO (-3707.932600) evaluations: (2831) 
#> Path [15] :Initial log joint density = -482423.890422 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.255e-02   2.735e-01    1.000e+00  1.000e+00      3513 -3.696e+03 -3.704e+03                   
#> Path [15] :Best Iter: [57] ELBO (-3696.428362) evaluations: (3513) 
#> Path [16] :Initial log joint density = -481853.944358 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.136e-03   1.778e-01    1.000e+00  1.000e+00      3157 -3.707e+03 -3.715e+03                   
#> Path [16] :Best Iter: [44] ELBO (-3707.231036) evaluations: (3157) 
#> Path [17] :Initial log joint density = -483082.905936 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.471e-02   3.080e-01    1.000e+00  1.000e+00      3547 -3.697e+03 -3.701e+03                   
#> Path [17] :Best Iter: [57] ELBO (-3697.437522) evaluations: (3547) 
#> Path [18] :Initial log joint density = -481872.548760 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.222e-03   2.616e-01    9.841e-01  9.841e-01      3271 -3.697e+03 -3.713e+03                   
#> Path [18] :Best Iter: [56] ELBO (-3696.991463) evaluations: (3271) 
#> Path [19] :Initial log joint density = -481778.755615 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.116e-02   2.560e-01    1.000e+00  1.000e+00      3320 -3.705e+03 -3.704e+03                   
#> Path [19] :Best Iter: [58] ELBO (-3704.229455) evaluations: (3320) 
#> Path [20] :Initial log joint density = -481620.547631 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.046e-03   2.458e-01    1.000e+00  1.000e+00      3160 -3.705e+03 -3.707e+03                   
#> Path [20] :Best Iter: [48] ELBO (-3704.662678) evaluations: (3160) 
#> Path [21] :Initial log joint density = -481930.410947 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.867e-03   2.074e-01    1.000e+00  1.000e+00      3220 -3.704e+03 -3.705e+03                   
#> Path [21] :Best Iter: [55] ELBO (-3704.286265) evaluations: (3220) 
#> Path [22] :Initial log joint density = -481638.587112 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      7.921e-03   2.849e-01    8.459e-01  8.459e-01      2549 -3.707e+03 -3.720e+03                   
#> Path [22] :Best Iter: [39] ELBO (-3706.987514) evaluations: (2549) 
#> Path [23] :Initial log joint density = -482118.464888 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.959e-03   2.039e-01    1.000e+00  1.000e+00      3426 -3.702e+03 -3.702e+03                   
#> Path [23] :Best Iter: [59] ELBO (-3701.690432) evaluations: (3426) 
#> Path [24] :Initial log joint density = -481778.219369 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.134e-02   3.187e-01    1.000e+00  1.000e+00      2832 -3.708e+03 -3.714e+03                   
#> Path [24] :Best Iter: [50] ELBO (-3708.337530) evaluations: (2832) 
#> Path [25] :Initial log joint density = -481764.271730 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.174e-02   2.592e-01    1.000e+00  1.000e+00      3272 -3.706e+03 -3.705e+03                   
#> Path [25] :Best Iter: [55] ELBO (-3704.892597) evaluations: (3272) 
#> Path [26] :Initial log joint density = -481774.958556 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.584e-03   3.203e-01    1.000e+00  1.000e+00      3137 -3.705e+03 -3.709e+03                   
#> Path [26] :Best Iter: [48] ELBO (-3705.349923) evaluations: (3137) 
#> Path [27] :Initial log joint density = -482583.637151 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      5.042e-03   2.271e-01    1.000e+00  1.000e+00      3618 -3.702e+03 -3.710e+03                   
#> Path [27] :Best Iter: [57] ELBO (-3702.400771) evaluations: (3618) 
#> Path [28] :Initial log joint density = -481723.352534 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.456e-02   3.291e-01    1.000e+00  1.000e+00      3357 -3.698e+03 -3.704e+03                   
#> Path [28] :Best Iter: [58] ELBO (-3698.147026) evaluations: (3357) 
#> Path [29] :Initial log joint density = -481776.672020 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.046e-02   2.248e-01    1.000e+00  1.000e+00      3543 -3.700e+03 -3.699e+03                   
#> Path [29] :Best Iter: [60] ELBO (-3699.252206) evaluations: (3543) 
#> Path [30] :Initial log joint density = -481620.998102 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.000e-03   2.039e-01    1.000e+00  1.000e+00      3074 -3.708e+03 -3.702e+03                   
#> Path [30] :Best Iter: [55] ELBO (-3702.239281) evaluations: (3074) 
#> Path [31] :Initial log joint density = -481406.523787 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.004e-02   2.717e-01    7.881e-01  7.881e-01      2863 -3.707e+03 -3.718e+03                   
#> Path [31] :Best Iter: [41] ELBO (-3707.161354) evaluations: (2863) 
#> Path [32] :Initial log joint density = -482297.012882 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.025e-02   2.858e-01    1.000e+00  1.000e+00      3311 -3.704e+03 -3.706e+03                   
#> Path [32] :Best Iter: [56] ELBO (-3703.929082) evaluations: (3311) 
#> Path [33] :Initial log joint density = -481815.365950 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.372e-03   2.040e-01    1.000e+00  1.000e+00      3332 -3.703e+03 -3.706e+03                   
#> Path [33] :Best Iter: [55] ELBO (-3703.495450) evaluations: (3332) 
#> Path [34] :Initial log joint density = -481570.693824 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.505e-03   1.762e-01    1.000e+00  1.000e+00      3521 -3.698e+03 -3.709e+03                   
#> Path [34] :Best Iter: [56] ELBO (-3697.519428) evaluations: (3521) 
#> Path [35] :Initial log joint density = -482106.915971 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.189e-02   2.220e-01    1.000e+00  1.000e+00      3330 -3.702e+03 -3.704e+03                   
#> Path [35] :Best Iter: [57] ELBO (-3702.192787) evaluations: (3330) 
#> Path [36] :Initial log joint density = -481532.168548 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.144e-02   2.091e-01    1.000e+00  1.000e+00      2935 -3.709e+03 -3.711e+03                   
#> Path [36] :Best Iter: [48] ELBO (-3709.065981) evaluations: (2935) 
#> Path [37] :Initial log joint density = -481702.550342 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.515e-03   2.010e-01    1.000e+00  1.000e+00      2944 -3.710e+03 -3.712e+03                   
#> Path [37] :Best Iter: [47] ELBO (-3709.641366) evaluations: (2944) 
#> Path [38] :Initial log joint density = -482996.098210 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.678e-03   2.939e-01    1.000e+00  1.000e+00      3191 -3.707e+03 -3.708e+03                   
#> Path [38] :Best Iter: [49] ELBO (-3706.901778) evaluations: (3191) 
#> Path [39] :Initial log joint density = -483845.405007 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.128e-02   3.007e-01    1.000e+00  1.000e+00      3400 -3.701e+03 -3.706e+03                   
#> Path [39] :Best Iter: [55] ELBO (-3700.913097) evaluations: (3400) 
#> Path [40] :Initial log joint density = -481904.676665 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      3.892e-03   2.305e-01    7.568e-01  7.568e-01      3356 -3.702e+03 -3.713e+03                   
#> Path [40] :Best Iter: [57] ELBO (-3701.979184) evaluations: (3356) 
#> Path [41] :Initial log joint density = -481558.174944 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.474e-02   3.153e-01    1.000e+00  1.000e+00      3163 -3.702e+03 -3.707e+03                   
#> Path [41] :Best Iter: [55] ELBO (-3701.833966) evaluations: (3163) 
#> Path [42] :Initial log joint density = -481586.057638 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.110e-03   1.917e-01    1.000e+00  1.000e+00      3244 -3.704e+03 -3.706e+03                   
#> Path [42] :Best Iter: [55] ELBO (-3704.239994) evaluations: (3244) 
#> Path [43] :Initial log joint density = -482319.267209 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      4.611e-03   2.006e-01    1.000e+00  1.000e+00      3363 -3.703e+03 -3.711e+03                   
#> Path [43] :Best Iter: [56] ELBO (-3702.872951) evaluations: (3363) 
#> Path [44] :Initial log joint density = -481304.464573 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      3.950e-03   1.938e-01    7.304e-01  7.304e-01      2910 -3.707e+03 -3.721e+03                   
#> Path [44] :Best Iter: [43] ELBO (-3706.940704) evaluations: (2910) 
#> Path [45] :Initial log joint density = -481600.997964 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.271e-02   3.455e-01    1.000e+00  1.000e+00      2783 -3.709e+03 -3.713e+03                   
#> Path [45] :Best Iter: [50] ELBO (-3709.307636) evaluations: (2783) 
#> Path [46] :Initial log joint density = -481489.480441 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      4.707e-03   1.093e-01    8.153e-01  8.153e-01      3358 -3.698e+03 -3.717e+03                   
#> Path [46] :Best Iter: [57] ELBO (-3698.254954) evaluations: (3358) 
#> Path [47] :Initial log joint density = -481949.776350 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.207e-02   2.406e-01    1.000e+00  1.000e+00      3238 -3.699e+03 -3.699e+03                   
#> Path [47] :Best Iter: [56] ELBO (-3698.606069) evaluations: (3238) 
#> Path [48] :Initial log joint density = -481385.989372 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      6.272e-03   2.248e-01    4.590e-01  1.000e+00      2869 -3.709e+03 -3.717e+03                   
#> Path [48] :Best Iter: [46] ELBO (-3709.277588) evaluations: (2869) 
#> Path [49] :Initial log joint density = -482299.926994 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.703e-02   3.929e-01    1.000e+00  1.000e+00      3209 -3.701e+03 -3.708e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3700.686816) evaluations: (3209) 
#> Path [50] :Initial log joint density = -481389.624713 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.207e-02   1.597e-01    1.000e+00  1.000e+00      3495 -3.700e+03 -3.699e+03                   
#> Path [50] :Best Iter: [60] ELBO (-3698.793939) evaluations: (3495) 
#> Finished in  13.4 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
