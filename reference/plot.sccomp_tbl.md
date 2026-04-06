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
  significance_statistic = c("pH0", "FDR"),
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
  "pH0".

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
#> Running make /tmp/RtmpT5zwuv/model-3b264ed0c124 "STAN_THREADS=TRUE" \
#>   "STANCFLAGS += --include-paths=/tmp/RtmpT5zwuv/temp_libpath3b26437ccd76/sccomp/stan --name='glm_multi_beta_binomial_model'"
#> 
#> --- Translating Stan model to C++ code ---
#> bin/stanc --include-paths=/tmp/RtmpT5zwuv/temp_libpath3b26437ccd76/sccomp/stan --name='glm_multi_beta_binomial_model' --o=/tmp/RtmpT5zwuv/model-3b264ed0c124.hpp /tmp/RtmpT5zwuv/model-3b264ed0c124.stan
#> 
#> --- Compiling C++ code ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -I stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -I stan/lib/stan_math/lib/eigen_3.4.0 -I stan/lib/stan_math/lib/boost_1.87.0 -I stan/lib/stan_math/lib/sundials_6.1.1/include -I stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS          -c -Wno-ignored-attributes   -x c++ -o /tmp/RtmpT5zwuv/model-3b264ed0c124.o /tmp/RtmpT5zwuv/model-3b264ed0c124.hpp
#> 
#> --- Linking model ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -I stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -I stan/lib/stan_math/lib/eigen_3.4.0 -I stan/lib/stan_math/lib/boost_1.87.0 -I stan/lib/stan_math/lib/sundials_6.1.1/include -I stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS               -Wl,-L,"/home/runner/.cmdstan/cmdstan-2.38.0/stan/lib/stan_math/lib/tbb"   -Wl,-rpath,"/home/runner/.cmdstan/cmdstan-2.38.0/stan/lib/stan_math/lib/tbb"      /tmp/RtmpT5zwuv/model-3b264ed0c124.o src/cmdstan/main_threads.o       -ltbb   stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_nvecserial.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_cvodes.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_idas.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_kinsol.a  stan/lib/stan_math/lib/tbb/libtbb.so.2 -o /tmp/RtmpT5zwuv/model-3b264ed0c124
#> rm /tmp/RtmpT5zwuv/model-3b264ed0c124.o /tmp/RtmpT5zwuv/model-3b264ed0c124.hpp
#> Model compiled and saved to cache successfully.
#> Path [1] :Initial log joint density = -482696.034405 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.100e-03   2.290e-01    1.000e+00  1.000e+00      3293 -3.685e+03 -3.684e+03                   
#> Path [1] :Best Iter: [57] ELBO (-3684.264959) evaluations: (3293) 
#> Path [2] :Initial log joint density = -482575.629154 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.346e-02   3.095e-01    9.455e-01  9.455e-01      2909 -3.692e+03 -3.713e+03                   
#> Path [2] :Best Iter: [42] ELBO (-3692.063021) evaluations: (2909) 
#> Path [3] :Initial log joint density = -481768.326965 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.788e+05      1.081e-02   3.007e-01    1.000e+00  1.000e+00      3950 -3.681e+03 -3.688e+03                   
#> Path [3] :Best Iter: [64] ELBO (-3681.380507) evaluations: (3950) 
#> Path [4] :Initial log joint density = -483335.763416 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.257e-03   3.415e-01    9.756e-01  9.756e-01      3277 -3.683e+03 -3.696e+03                   
#> Path [4] :Best Iter: [57] ELBO (-3682.657146) evaluations: (3277) 
#> Path [5] :Initial log joint density = -482409.383290 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.023e-02   2.490e-01    1.000e+00  1.000e+00      3025 -3.690e+03 -3.688e+03                   
#> Path [5] :Best Iter: [55] ELBO (-3688.194889) evaluations: (3025) 
#> Path [6] :Initial log joint density = -483034.637676 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.291e-02   2.509e-01    1.000e+00  1.000e+00      3192 -3.682e+03 -3.683e+03                   
#> Path [6] :Best Iter: [56] ELBO (-3682.307131) evaluations: (3192) 
#> Path [7] :Initial log joint density = -482041.683904 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      3.626e-03   2.176e-01    7.938e-01  7.938e-01      3449 -3.685e+03 -3.697e+03                   
#> Path [7] :Best Iter: [56] ELBO (-3684.717885) evaluations: (3449) 
#> Path [8] :Initial log joint density = -482390.384478 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.647e-03   2.162e-01    1.000e+00  1.000e+00      3348 -3.687e+03 -3.684e+03                   
#> Path [8] :Best Iter: [57] ELBO (-3684.445957) evaluations: (3348) 
#> Path [9] :Initial log joint density = -482173.820708 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.374e-02   2.931e-01    9.618e-01  9.618e-01      3145 -3.686e+03 -3.692e+03                   
#> Path [9] :Best Iter: [55] ELBO (-3685.834287) evaluations: (3145) 
#> Path [10] :Initial log joint density = -481441.169881 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.334e-02   2.824e-01    1.000e+00  1.000e+00      3037 -3.692e+03 -3.696e+03                   
#> Path [10] :Best Iter: [51] ELBO (-3691.837471) evaluations: (3037) 
#> Path [11] :Initial log joint density = -481548.264723 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      2.061e-03   2.833e-01    5.955e-01  5.955e-01      3498 -3.682e+03 -3.698e+03                   
#> Path [11] :Best Iter: [59] ELBO (-3682.330804) evaluations: (3498) 
#> Path [12] :Initial log joint density = -481524.784047 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.039e-02   2.557e-01    1.000e+00  1.000e+00      2783 -3.690e+03 -3.689e+03                   
#> Path [12] :Best Iter: [51] ELBO (-3689.127762) evaluations: (2783) 
#> Path [13] :Initial log joint density = -481755.600981 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.614e-02   2.741e-01    1.000e+00  1.000e+00      3006 -3.689e+03 -3.698e+03                   
#> Path [13] :Best Iter: [38] ELBO (-3688.890253) evaluations: (3006) 
#> Path [14] :Initial log joint density = -481518.977061 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      1.309e-02   3.621e-01    1.000e+00  1.000e+00      3863 -3.681e+03 -3.688e+03                   
#> Path [14] :Best Iter: [62] ELBO (-3681.316662) evaluations: (3863) 
#> Path [15] :Initial log joint density = -481838.545746 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      7.195e-03   1.784e-01    9.482e-01  9.482e-01      3621 -3.682e+03 -3.693e+03                   
#> Path [15] :Best Iter: [59] ELBO (-3682.303562) evaluations: (3621) 
#> Path [16] :Initial log joint density = -481807.480647 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.166e-02   2.890e-01    1.000e+00  1.000e+00      3109 -3.681e+03 -3.685e+03                   
#> Path [16] :Best Iter: [55] ELBO (-3681.156938) evaluations: (3109) 
#> Path [17] :Initial log joint density = -481855.963124 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      5.883e-03   1.505e-01    1.000e+00  1.000e+00      2832 -3.689e+03 -3.694e+03                   
#> Path [17] :Best Iter: [49] ELBO (-3688.606235) evaluations: (2832) 
#> Path [18] :Initial log joint density = -481552.613951 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      8.159e-03   2.736e-01    8.565e-01  8.565e-01      2762 -3.691e+03 -3.701e+03                   
#> Path [18] :Best Iter: [50] ELBO (-3690.900149) evaluations: (2762) 
#> Path [19] :Initial log joint density = -481348.770091 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.339e-02   3.058e-01    1.000e+00  1.000e+00      3595 -3.685e+03 -3.685e+03                   
#> Path [19] :Best Iter: [58] ELBO (-3684.957376) evaluations: (3595) 
#> Path [20] :Initial log joint density = -483393.624223 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.056e-03   2.242e-01    3.985e-01  1.000e+00      3278 -3.682e+03 -3.689e+03                   
#> Path [20] :Best Iter: [57] ELBO (-3682.365390) evaluations: (3278) 
#> Path [21] :Initial log joint density = -481642.359415 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.065e-02   2.092e-01    1.000e+00  1.000e+00      3470 -3.685e+03 -3.687e+03                   
#> Path [21] :Best Iter: [59] ELBO (-3685.159872) evaluations: (3470) 
#> Path [22] :Initial log joint density = -481852.343994 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      8.277e-03   2.264e-01    8.262e-01  8.262e-01      3565 -3.683e+03 -3.693e+03                   
#> Path [22] :Best Iter: [55] ELBO (-3683.295230) evaluations: (3565) 
#> Path [23] :Initial log joint density = -481533.756819 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      2.880e-03   2.915e-01    6.307e-01  6.307e-01      3256 -3.692e+03 -3.697e+03                   
#> Path [23] :Best Iter: [55] ELBO (-3691.769636) evaluations: (3256) 
#> Path [24] :Initial log joint density = -481794.794276 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      3.867e-03   2.724e-01    6.208e-01  6.208e-01      2957 -3.692e+03 -3.718e+03                   
#> Path [24] :Best Iter: [52] ELBO (-3691.882932) evaluations: (2957) 
#> Path [25] :Initial log joint density = -481609.573590 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      3.110e-03   2.848e-01    6.856e-01  6.856e-01      3316 -3.683e+03 -3.699e+03                   
#> Path [25] :Best Iter: [56] ELBO (-3683.001738) evaluations: (3316) 
#> Path [26] :Initial log joint density = -481555.661178 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.914e-02   2.974e-01    1.000e+00  1.000e+00      3419 -3.682e+03 -3.687e+03                   
#> Path [26] :Best Iter: [55] ELBO (-3682.492770) evaluations: (3419) 
#> Path [27] :Initial log joint density = -482142.853546 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      7.747e-03   2.396e-01    9.451e-01  9.451e-01      3546 -3.685e+03 -3.692e+03                   
#> Path [27] :Best Iter: [58] ELBO (-3684.834771) evaluations: (3546) 
#> Path [28] :Initial log joint density = -483596.432907 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.226e-03   2.002e-01    1.000e+00  1.000e+00      3346 -3.685e+03 -3.683e+03                   
#> Path [28] :Best Iter: [57] ELBO (-3682.927344) evaluations: (3346) 
#> Path [29] :Initial log joint density = -481538.161655 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.630e-03   2.330e-01    1.000e+00  1.000e+00      3134 -3.692e+03 -3.692e+03                   
#> Path [29] :Best Iter: [55] ELBO (-3691.514265) evaluations: (3134) 
#> Path [30] :Initial log joint density = -485181.467539 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.481e-03   2.349e-01    1.000e+00  1.000e+00      3354 -3.690e+03 -3.684e+03                   
#> Path [30] :Best Iter: [57] ELBO (-3683.865788) evaluations: (3354) 
#> Path [31] :Initial log joint density = -482130.334454 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.349e-03   2.137e-01    1.000e+00  1.000e+00      3363 -3.682e+03 -3.681e+03                   
#> Path [31] :Best Iter: [59] ELBO (-3681.394210) evaluations: (3363) 
#> Path [32] :Initial log joint density = -481732.952881 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.503e-02   3.932e-01    1.000e+00  1.000e+00      3119 -3.690e+03 -3.689e+03                   
#> Path [32] :Best Iter: [55] ELBO (-3688.856491) evaluations: (3119) 
#> Path [33] :Initial log joint density = -482727.895332 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.286e-02   2.347e-01    1.000e+00  1.000e+00      3554 -3.685e+03 -3.687e+03                   
#> Path [33] :Best Iter: [59] ELBO (-3684.711455) evaluations: (3554) 
#> Path [34] :Initial log joint density = -481397.656790 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.747e-02   3.465e-01    1.000e+00  1.000e+00      3234 -3.683e+03 -3.689e+03                   
#> Path [34] :Best Iter: [55] ELBO (-3683.479865) evaluations: (3234) 
#> Path [35] :Initial log joint density = -481710.669373 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      4.242e-03   1.928e-01    8.071e-01  8.071e-01      2995 -3.691e+03 -3.708e+03                   
#> Path [35] :Best Iter: [43] ELBO (-3691.257188) evaluations: (2995) 
#> Path [36] :Initial log joint density = -481616.593440 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.148e-02   1.962e-01    1.000e+00  1.000e+00      2974 -3.693e+03 -3.684e+03                   
#> Path [36] :Best Iter: [55] ELBO (-3684.171631) evaluations: (2974) 
#> Path [37] :Initial log joint density = -484285.716139 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.788e+05      1.286e-02   1.327e-01    1.000e+00  1.000e+00      4011 -3.681e+03 -3.681e+03                   
#> Path [37] :Best Iter: [62] ELBO (-3681.242130) evaluations: (4011) 
#> Path [38] :Initial log joint density = -481727.117331 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.142e-03   1.787e-01    1.000e+00  1.000e+00      3045 -3.689e+03 -3.695e+03                   
#> Path [38] :Best Iter: [50] ELBO (-3688.936326) evaluations: (3045) 
#> Path [39] :Initial log joint density = -483834.429459 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.780e-02   2.433e-01    1.000e+00  1.000e+00      3872 -3.681e+03 -3.685e+03                   
#> Path [39] :Best Iter: [61] ELBO (-3681.206106) evaluations: (3872) 
#> Path [40] :Initial log joint density = -481820.175988 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.732e-03   1.626e-01    1.000e+00  1.000e+00      3399 -3.689e+03 -3.696e+03                   
#> Path [40] :Best Iter: [45] ELBO (-3689.280284) evaluations: (3399) 
#> Path [41] :Initial log joint density = -481411.052556 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.146e-03   2.733e-01    9.700e-01  9.700e-01      2934 -3.688e+03 -3.703e+03                   
#> Path [41] :Best Iter: [51] ELBO (-3687.914719) evaluations: (2934) 
#> Path [42] :Initial log joint density = -481793.895469 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.422e-02   2.643e-01    1.000e+00  1.000e+00      3356 -3.682e+03 -3.684e+03                   
#> Path [42] :Best Iter: [55] ELBO (-3682.071874) evaluations: (3356) 
#> Path [43] :Initial log joint density = -481561.109533 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      4.362e-03   2.874e-01    6.251e-01  6.251e-01      3325 -3.684e+03 -3.692e+03                   
#> Path [43] :Best Iter: [55] ELBO (-3683.725104) evaluations: (3325) 
#> Path [44] :Initial log joint density = -481410.120520 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.284e-02   1.953e-01    1.000e+00  1.000e+00      3385 -3.680e+03 -3.680e+03                   
#> Path [44] :Best Iter: [57] ELBO (-3679.923481) evaluations: (3385) 
#> Path [45] :Initial log joint density = -482106.783428 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      3.404e-03   2.365e-01    6.931e-01  6.931e-01      3277 -3.683e+03 -3.696e+03                   
#> Path [45] :Best Iter: [55] ELBO (-3683.041037) evaluations: (3277) 
#> Path [46] :Initial log joint density = -481673.153145 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.541e-02   3.889e-01    9.763e-01  9.763e-01      2832 -3.690e+03 -3.699e+03                   
#> Path [46] :Best Iter: [51] ELBO (-3689.829585) evaluations: (2832) 
#> Path [47] :Initial log joint density = -481598.770257 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      8.729e-03   2.193e-01    1.000e+00  1.000e+00      2704 -3.691e+03 -3.690e+03                   
#> Path [47] :Best Iter: [51] ELBO (-3690.185746) evaluations: (2704) 
#> Path [48] :Initial log joint density = -481464.098525 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.083e-02   2.050e-01    1.000e+00  1.000e+00      3153 -3.690e+03 -3.683e+03                   
#> Path [48] :Best Iter: [56] ELBO (-3683.082964) evaluations: (3153) 
#> Path [49] :Initial log joint density = -481668.302441 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.900e-03   2.179e-01    1.000e+00  1.000e+00      3002 -3.691e+03 -3.692e+03                   
#> Path [49] :Best Iter: [38] ELBO (-3691.221904) evaluations: (3002) 
#> Path [50] :Initial log joint density = -481662.774713 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.155e-03   2.013e-01    1.000e+00  1.000e+00      3220 -3.689e+03 -3.690e+03                   
#> Path [50] :Best Iter: [53] ELBO (-3688.555504) evaluations: (3220) 
#> Finished in  13.5 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
