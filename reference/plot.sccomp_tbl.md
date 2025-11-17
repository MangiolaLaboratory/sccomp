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
#> Running make /tmp/RtmpFf0i1u/model-3ca024f58da9 "STAN_THREADS=TRUE" \
#>   "STANCFLAGS += --include-paths=/tmp/RtmpFf0i1u/temp_libpath3ca03c17963f/sccomp/stan --name='glm_multi_beta_binomial_model'"
#> 
#> --- Translating Stan model to C++ code ---
#> bin/stanc --include-paths=/tmp/RtmpFf0i1u/temp_libpath3ca03c17963f/sccomp/stan --name='glm_multi_beta_binomial_model' --o=/tmp/RtmpFf0i1u/model-3ca024f58da9.hpp /tmp/RtmpFf0i1u/model-3ca024f58da9.stan
#> 
#> --- Compiling C++ code ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -I stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -I stan/lib/stan_math/lib/eigen_3.4.0 -I stan/lib/stan_math/lib/boost_1.87.0 -I stan/lib/stan_math/lib/sundials_6.1.1/include -I stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS          -c -Wno-ignored-attributes   -x c++ -o /tmp/RtmpFf0i1u/model-3ca024f58da9.o /tmp/RtmpFf0i1u/model-3ca024f58da9.hpp
#> 
#> --- Linking model ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -I stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -I stan/lib/stan_math/lib/eigen_3.4.0 -I stan/lib/stan_math/lib/boost_1.87.0 -I stan/lib/stan_math/lib/sundials_6.1.1/include -I stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS               -Wl,-L,"/home/runner/.cmdstan/cmdstan-2.37.0/stan/lib/stan_math/lib/tbb"   -Wl,-rpath,"/home/runner/.cmdstan/cmdstan-2.37.0/stan/lib/stan_math/lib/tbb"      /tmp/RtmpFf0i1u/model-3ca024f58da9.o src/cmdstan/main_threads.o       -ltbb   stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_nvecserial.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_cvodes.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_idas.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_kinsol.a  stan/lib/stan_math/lib/tbb/libtbb.so.2 -o /tmp/RtmpFf0i1u/model-3ca024f58da9
#> rm /tmp/RtmpFf0i1u/model-3ca024f58da9.hpp /tmp/RtmpFf0i1u/model-3ca024f58da9.o
#> Model compiled and saved to cache successfully.
#> Path [1] :Initial log joint density = -481862.988576 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.531e-02   2.088e-01    1.000e+00  1.000e+00      3593 -3.701e+03 -3.700e+03                   
#> Path [1] :Best Iter: [59] ELBO (-3700.231028) evaluations: (3593) 
#> Path [2] :Initial log joint density = -482290.585323 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.589e-03   2.126e-01    1.000e+00  1.000e+00      3072 -3.708e+03 -3.713e+03                   
#> Path [2] :Best Iter: [48] ELBO (-3707.987035) evaluations: (3072) 
#> Path [3] :Initial log joint density = -481684.913644 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      5.032e-03   2.097e-01    9.523e-01  9.523e-01      2798 -3.705e+03 -3.717e+03                   
#> Path [3] :Best Iter: [48] ELBO (-3705.289391) evaluations: (2798) 
#> Path [4] :Initial log joint density = -481700.614186 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.037e-02   2.289e-01    1.000e+00  1.000e+00      3293 -3.702e+03 -3.709e+03                   
#> Path [4] :Best Iter: [56] ELBO (-3701.509975) evaluations: (3293) 
#> Path [5] :Initial log joint density = -481696.455772 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.517e-03   1.964e-01    1.000e+00  1.000e+00      3245 -3.702e+03 -3.698e+03                   
#> Path [5] :Best Iter: [57] ELBO (-3697.824691) evaluations: (3245) 
#> Path [6] :Initial log joint density = -481317.790661 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.652e-02   3.997e-01    1.000e+00  1.000e+00      3082 -3.706e+03 -3.717e+03                   
#> Path [6] :Best Iter: [43] ELBO (-3706.376535) evaluations: (3082) 
#> Path [7] :Initial log joint density = -481456.506182 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.445e-02   3.905e-01    1.000e+00  1.000e+00      3149 -3.698e+03 -3.704e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3698.257585) evaluations: (3149) 
#> Path [8] :Initial log joint density = -481503.105373 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.569e-02   3.697e-01    1.000e+00  1.000e+00      3136 -3.700e+03 -3.708e+03                   
#> Path [8] :Best Iter: [56] ELBO (-3699.858984) evaluations: (3136) 
#> Path [9] :Initial log joint density = -481319.923982 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.860e-02   4.522e-01    1.000e+00  1.000e+00      3177 -3.708e+03 -3.706e+03                   
#> Path [9] :Best Iter: [55] ELBO (-3705.645267) evaluations: (3177) 
#> Path [10] :Initial log joint density = -481742.182993 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.243e-03   1.800e-01    1.000e+00  1.000e+00      3080 -3.704e+03 -3.716e+03                   
#> Path [10] :Best Iter: [42] ELBO (-3703.921100) evaluations: (3080) 
#> Path [11] :Initial log joint density = -483411.221409 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.616e-03   2.562e-01    1.000e+00  1.000e+00      3344 -3.706e+03 -3.710e+03                   
#> Path [11] :Best Iter: [45] ELBO (-3705.960071) evaluations: (3344) 
#> Path [12] :Initial log joint density = -481627.486259 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      5.817e-03   2.465e-01    8.348e-01  8.348e-01      3382 -3.702e+03 -3.715e+03                   
#> Path [12] :Best Iter: [57] ELBO (-3701.983691) evaluations: (3382) 
#> Path [13] :Initial log joint density = -482097.177815 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.426e-02   2.661e-01    1.000e+00  1.000e+00      3480 -3.698e+03 -3.699e+03                   
#> Path [13] :Best Iter: [59] ELBO (-3697.720469) evaluations: (3480) 
#> Path [14] :Initial log joint density = -482761.842072 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.270e-02   3.223e-01    1.000e+00  1.000e+00      3503 -3.701e+03 -3.706e+03                   
#> Path [14] :Best Iter: [59] ELBO (-3700.502409) evaluations: (3503) 
#> Path [15] :Initial log joint density = -481410.820999 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.646e-03   2.095e-01    4.505e-01  1.000e+00      2906 -3.709e+03 -3.720e+03                   
#> Path [15] :Best Iter: [47] ELBO (-3708.806655) evaluations: (2906) 
#> Path [16] :Initial log joint density = -483251.937308 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      5.046e-03   2.137e-01    8.471e-01  8.471e-01      3544 -3.704e+03 -3.708e+03                   
#> Path [16] :Best Iter: [58] ELBO (-3703.780707) evaluations: (3544) 
#> Path [17] :Initial log joint density = -482259.803342 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.613e-03   2.576e-01    1.000e+00  1.000e+00      2995 -3.709e+03 -3.718e+03                   
#> Path [17] :Best Iter: [47] ELBO (-3709.439874) evaluations: (2995) 
#> Path [18] :Initial log joint density = -482480.541580 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.615e-03   2.669e-01    1.000e+00  1.000e+00      3373 -3.701e+03 -3.711e+03                   
#> Path [18] :Best Iter: [56] ELBO (-3701.263421) evaluations: (3373) 
#> Path [19] :Initial log joint density = -481757.440327 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.705e-03   2.535e-01    1.000e+00  1.000e+00      3357 -3.701e+03 -3.700e+03                   
#> Path [19] :Best Iter: [59] ELBO (-3700.260288) evaluations: (3357) 
#> Path [20] :Initial log joint density = -481673.755792 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.370e-02   3.480e-01    1.000e+00  1.000e+00      3495 -3.700e+03 -3.702e+03                   
#> Path [20] :Best Iter: [55] ELBO (-3700.221937) evaluations: (3495) 
#> Path [21] :Initial log joint density = -482673.895844 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.760e-03   1.553e-01    7.828e-01  7.828e-01      3162 -3.702e+03 -3.710e+03                   
#> Path [21] :Best Iter: [55] ELBO (-3701.737683) evaluations: (3162) 
#> Path [22] :Initial log joint density = -481798.629188 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      5.909e-03   2.428e-01    7.363e-01  7.363e-01      3130 -3.708e+03 -3.723e+03                   
#> Path [22] :Best Iter: [40] ELBO (-3707.584316) evaluations: (3130) 
#> Path [23] :Initial log joint density = -481768.942323 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.012e-03   2.389e-01    9.176e-01  9.176e-01      2914 -3.706e+03 -3.714e+03                   
#> Path [23] :Best Iter: [39] ELBO (-3705.554292) evaluations: (2914) 
#> Path [24] :Initial log joint density = -481611.961978 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.520e-02   3.540e-01    1.000e+00  1.000e+00      3179 -3.708e+03 -3.716e+03                   
#> Path [24] :Best Iter: [34] ELBO (-3708.214485) evaluations: (3179) 
#> Path [25] :Initial log joint density = -483188.315722 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.651e-02   2.759e-01    1.000e+00  1.000e+00      3555 -3.699e+03 -3.703e+03                   
#> Path [25] :Best Iter: [59] ELBO (-3698.767478) evaluations: (3555) 
#> Path [26] :Initial log joint density = -482594.925406 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.085e-03   1.665e-01    1.000e+00  1.000e+00      3271 -3.706e+03 -3.701e+03                   
#> Path [26] :Best Iter: [57] ELBO (-3700.573582) evaluations: (3271) 
#> Path [27] :Initial log joint density = -481746.267152 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.546e-03   2.785e-01    9.161e-01  9.161e-01      3017 -3.710e+03 -3.711e+03                   
#> Path [27] :Best Iter: [40] ELBO (-3710.097717) evaluations: (3017) 
#> Path [28] :Initial log joint density = -482838.234823 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.023e-02   2.678e-01    1.000e+00  1.000e+00      3428 -3.703e+03 -3.703e+03                   
#> Path [28] :Best Iter: [59] ELBO (-3703.398838) evaluations: (3428) 
#> Path [29] :Initial log joint density = -485360.105071 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.472e-03   2.295e-01    1.000e+00  1.000e+00      3223 -3.708e+03 -3.705e+03                   
#> Path [29] :Best Iter: [55] ELBO (-3705.013329) evaluations: (3223) 
#> Path [30] :Initial log joint density = -481451.816285 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.781e-03   1.466e-01    1.000e+00  1.000e+00      2918 -3.707e+03 -3.710e+03                   
#> Path [30] :Best Iter: [52] ELBO (-3706.920935) evaluations: (2918) 
#> Path [31] :Initial log joint density = -481904.295606 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.087e-02   2.142e-01    1.000e+00  1.000e+00      3317 -3.704e+03 -3.709e+03                   
#> Path [31] :Best Iter: [55] ELBO (-3704.421405) evaluations: (3317) 
#> Path [32] :Initial log joint density = -481764.569795 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      2.846e-03   2.455e-01    5.932e-01  5.932e-01      3445 -3.702e+03 -3.713e+03                   
#> Path [32] :Best Iter: [58] ELBO (-3701.737994) evaluations: (3445) 
#> Path [33] :Initial log joint density = -481864.539955 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.209e-02   3.722e-01    1.000e+00  1.000e+00      2632 -3.708e+03 -3.723e+03                   
#> Path [33] :Best Iter: [47] ELBO (-3707.904145) evaluations: (2632) 
#> Path [34] :Initial log joint density = -481760.691268 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.542e-02   2.821e-01    1.000e+00  1.000e+00      3161 -3.702e+03 -3.708e+03                   
#> Path [34] :Best Iter: [55] ELBO (-3702.218239) evaluations: (3161) 
#> Path [35] :Initial log joint density = -481218.809659 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.471e-02   2.909e-01    1.000e+00  1.000e+00      3106 -3.707e+03 -3.702e+03                   
#> Path [35] :Best Iter: [56] ELBO (-3702.162515) evaluations: (3106) 
#> Path [36] :Initial log joint density = -481595.662328 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.296e-03   2.326e-01    1.000e+00  1.000e+00      2937 -3.709e+03 -3.716e+03                   
#> Path [36] :Best Iter: [51] ELBO (-3708.867823) evaluations: (2937) 
#> Path [37] :Initial log joint density = -481671.595051 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      3.304e-03   2.103e-01    8.298e-01  8.298e-01      3391 -3.702e+03 -3.714e+03                   
#> Path [37] :Best Iter: [58] ELBO (-3701.524678) evaluations: (3391) 
#> Path [38] :Initial log joint density = -481589.764069 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      9.696e-03   2.576e-01    1.000e+00  1.000e+00      3053 -3.701e+03 -3.700e+03                   
#> Path [38] :Best Iter: [56] ELBO (-3700.451759) evaluations: (3053) 
#> Path [39] :Initial log joint density = -481836.297725 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.167e-02   4.165e-01    9.490e-01  9.490e-01      3576 -3.700e+03 -3.710e+03                   
#> Path [39] :Best Iter: [58] ELBO (-3700.066422) evaluations: (3576) 
#> Path [40] :Initial log joint density = -481544.800496 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.980e-03   2.153e-01    1.000e+00  1.000e+00      2909 -3.708e+03 -3.712e+03                   
#> Path [40] :Best Iter: [49] ELBO (-3708.077640) evaluations: (2909) 
#> Path [41] :Initial log joint density = -481585.005091 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      3.188e-03   2.584e-01    6.529e-01  6.529e-01      3449 -3.700e+03 -3.712e+03                   
#> Path [41] :Best Iter: [57] ELBO (-3700.151413) evaluations: (3449) 
#> Path [42] :Initial log joint density = -481936.638697 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.143e-02   2.201e-01    9.419e-01  9.419e-01      3313 -3.701e+03 -3.712e+03                   
#> Path [42] :Best Iter: [56] ELBO (-3700.910754) evaluations: (3313) 
#> Path [43] :Initial log joint density = -481642.464190 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.181e-02   2.687e-01    7.419e-01  7.419e-01      3318 -3.703e+03 -3.712e+03                   
#> Path [43] :Best Iter: [56] ELBO (-3703.238736) evaluations: (3318) 
#> Path [44] :Initial log joint density = -481785.839566 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.603e-02   1.972e-01    1.000e+00  1.000e+00      2915 -3.706e+03 -3.712e+03                   
#> Path [44] :Best Iter: [52] ELBO (-3705.757674) evaluations: (2915) 
#> Path [45] :Initial log joint density = -481993.445979 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      2.383e-03   2.872e-01    5.858e-01  5.858e-01      3084 -3.707e+03 -3.721e+03                   
#> Path [45] :Best Iter: [46] ELBO (-3707.224495) evaluations: (3084) 
#> Path [46] :Initial log joint density = -481915.419186 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      3.342e-03   2.768e-01    6.118e-01  6.118e-01      3356 -3.701e+03 -3.716e+03                   
#> Path [46] :Best Iter: [57] ELBO (-3701.083257) evaluations: (3356) 
#> Path [47] :Initial log joint density = -481764.274299 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      3.836e-03   3.522e-01    6.411e-01  6.411e-01      3257 -3.709e+03 -3.716e+03                   
#> Path [47] :Best Iter: [55] ELBO (-3708.791195) evaluations: (3257) 
#> Path [48] :Initial log joint density = -481797.449839 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      8.975e-03   1.846e-01    1.000e+00  1.000e+00      3847 -3.701e+03 -3.706e+03                   
#> Path [48] :Best Iter: [58] ELBO (-3701.167893) evaluations: (3847) 
#> Path [49] :Initial log joint density = -481668.280115 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.037e-02   1.832e-01    1.000e+00  1.000e+00      3045 -3.709e+03 -3.711e+03                   
#> Path [49] :Best Iter: [40] ELBO (-3708.607085) evaluations: (3045) 
#> Path [50] :Initial log joint density = -482256.963741 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.239e-03   2.154e-01    1.000e+00  1.000e+00      3156 -3.704e+03 -3.704e+03                   
#> Path [50] :Best Iter: [56] ELBO (-3704.449832) evaluations: (3156) 
#> Finished in  13.4 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
# }
```
