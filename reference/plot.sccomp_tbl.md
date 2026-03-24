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
#> Running make /tmp/RtmpJhaTpp/model-3c7c247fbd28 "STAN_THREADS=TRUE" \
#>   "STANCFLAGS += --include-paths=/tmp/RtmpJhaTpp/temp_libpath3c7c6ec8aa3a/sccomp/stan --name='glm_multi_beta_binomial_model'"
#> 
#> --- Translating Stan model to C++ code ---
#> bin/stanc --include-paths=/tmp/RtmpJhaTpp/temp_libpath3c7c6ec8aa3a/sccomp/stan --name='glm_multi_beta_binomial_model' --o=/tmp/RtmpJhaTpp/model-3c7c247fbd28.hpp /tmp/RtmpJhaTpp/model-3c7c247fbd28.stan
#> 
#> --- Compiling C++ code ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -I stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -I stan/lib/stan_math/lib/eigen_3.4.0 -I stan/lib/stan_math/lib/boost_1.87.0 -I stan/lib/stan_math/lib/sundials_6.1.1/include -I stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS          -c -Wno-ignored-attributes   -x c++ -o /tmp/RtmpJhaTpp/model-3c7c247fbd28.o /tmp/RtmpJhaTpp/model-3c7c247fbd28.hpp
#> 
#> --- Linking model ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -I stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -I stan/lib/stan_math/lib/eigen_3.4.0 -I stan/lib/stan_math/lib/boost_1.87.0 -I stan/lib/stan_math/lib/sundials_6.1.1/include -I stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS               -Wl,-L,"/home/runner/.cmdstan/cmdstan-2.38.0/stan/lib/stan_math/lib/tbb"   -Wl,-rpath,"/home/runner/.cmdstan/cmdstan-2.38.0/stan/lib/stan_math/lib/tbb"      /tmp/RtmpJhaTpp/model-3c7c247fbd28.o src/cmdstan/main_threads.o       -ltbb   stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_nvecserial.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_cvodes.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_idas.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_kinsol.a  stan/lib/stan_math/lib/tbb/libtbb.so.2 -o /tmp/RtmpJhaTpp/model-3c7c247fbd28
#> rm /tmp/RtmpJhaTpp/model-3c7c247fbd28.o /tmp/RtmpJhaTpp/model-3c7c247fbd28.hpp
#> Model compiled and saved to cache successfully.
#> Path [1] :Initial log joint density = -482697.585689 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.103e-03   2.227e-01    4.960e-01  1.000e+00      3529 -3.702e+03 -3.709e+03                   
#> Path [1] :Best Iter: [57] ELBO (-3701.910821) evaluations: (3529) 
#> Path [2] :Initial log joint density = -482578.349549 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.385e-03   2.405e-01    8.054e-01  8.054e-01      3074 -3.708e+03 -3.715e+03                   
#> Path [2] :Best Iter: [41] ELBO (-3707.507444) evaluations: (3074) 
#> Path [3] :Initial log joint density = -481769.787432 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      9.547e-03   2.167e-01    1.000e+00  1.000e+00      3692 -3.701e+03 -3.704e+03                   
#> Path [3] :Best Iter: [56] ELBO (-3700.652232) evaluations: (3692) 
#> Path [4] :Initial log joint density = -483337.744884 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.175e-03   1.136e-01    7.654e-01  7.654e-01      3192 -3.703e+03 -3.713e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3702.971816) evaluations: (3192) 
#> Path [5] :Initial log joint density = -482409.776555 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.119e-03   1.939e-01    1.000e+00  1.000e+00      2943 -3.706e+03 -3.707e+03                   
#> Path [5] :Best Iter: [48] ELBO (-3705.666406) evaluations: (2943) 
#> Path [6] :Initial log joint density = -483035.512893 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.842e-03   2.670e-01    4.429e-01  1.000e+00      3193 -3.703e+03 -3.712e+03                   
#> Path [6] :Best Iter: [56] ELBO (-3703.146646) evaluations: (3193) 
#> Path [7] :Initial log joint density = -482043.155979 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.169e-03   2.218e-01    7.620e-01  7.620e-01      3236 -3.708e+03 -3.714e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3707.588260) evaluations: (3236) 
#> Path [8] :Initial log joint density = -482391.731935 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      4.748e-03   1.988e-01    8.978e-01  8.978e-01      3162 -3.703e+03 -3.710e+03                   
#> Path [8] :Best Iter: [55] ELBO (-3703.014301) evaluations: (3162) 
#> Path [9] :Initial log joint density = -482175.842822 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.771e-03   2.692e-01    1.000e+00  1.000e+00      2890 -3.707e+03 -3.714e+03                   
#> Path [9] :Best Iter: [41] ELBO (-3707.229435) evaluations: (2890) 
#> Path [10] :Initial log joint density = -481441.540729 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      2.388e-03   2.476e-01    6.062e-01  6.062e-01      3201 -3.707e+03 -3.716e+03                   
#> Path [10] :Best Iter: [48] ELBO (-3706.904551) evaluations: (3201) 
#> Path [11] :Initial log joint density = -481548.580400 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.516e-03   2.396e-01    1.000e+00  1.000e+00      3136 -3.706e+03 -3.710e+03                   
#> Path [11] :Best Iter: [55] ELBO (-3706.198227) evaluations: (3136) 
#> Path [12] :Initial log joint density = -481525.496964 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      7.158e-03   2.684e-01    4.475e-01  1.000e+00      2951 -3.706e+03 -3.716e+03                   
#> Path [12] :Best Iter: [51] ELBO (-3705.844345) evaluations: (2951) 
#> Path [13] :Initial log joint density = -481756.269322 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.423e-03   2.226e-01    1.000e+00  1.000e+00      3386 -3.700e+03 -3.703e+03                   
#> Path [13] :Best Iter: [55] ELBO (-3699.580277) evaluations: (3386) 
#> Path [14] :Initial log joint density = -481520.824720 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      6.054e-03   2.129e-01    1.000e+00  1.000e+00      3650 -3.700e+03 -3.706e+03                   
#> Path [14] :Best Iter: [57] ELBO (-3700.076453) evaluations: (3650) 
#> Path [15] :Initial log joint density = -481839.345681 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      5.362e-03   2.311e-01    7.190e-01  7.190e-01      3531 -3.701e+03 -3.715e+03                   
#> Path [15] :Best Iter: [59] ELBO (-3700.963115) evaluations: (3531) 
#> Path [16] :Initial log joint density = -481807.784454 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.206e-03   1.956e-01    1.000e+00  1.000e+00      3049 -3.707e+03 -3.701e+03                   
#> Path [16] :Best Iter: [55] ELBO (-3701.110068) evaluations: (3049) 
#> Path [17] :Initial log joint density = -481857.067887 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      9.417e-03   2.441e-01    1.000e+00  1.000e+00      2673 -3.703e+03 -3.716e+03                   
#> Path [17] :Best Iter: [49] ELBO (-3703.381826) evaluations: (2673) 
#> Path [18] :Initial log joint density = -481553.304735 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.068e-02   2.818e-01    1.000e+00  1.000e+00      2673 -3.708e+03 -3.711e+03                   
#> Path [18] :Best Iter: [46] ELBO (-3708.039432) evaluations: (2673) 
#> Path [19] :Initial log joint density = -481351.567710 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.360e-02   3.392e-01    1.000e+00  1.000e+00      3540 -3.702e+03 -3.705e+03                   
#> Path [19] :Best Iter: [58] ELBO (-3702.192222) evaluations: (3540) 
#> Path [20] :Initial log joint density = -483394.378138 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.467e-02   2.439e-01    1.000e+00  1.000e+00      3280 -3.702e+03 -3.703e+03                   
#> Path [20] :Best Iter: [55] ELBO (-3702.298097) evaluations: (3280) 
#> Path [21] :Initial log joint density = -481642.596882 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.113e-02   2.261e-01    1.000e+00  1.000e+00      3138 -3.704e+03 -3.701e+03                   
#> Path [21] :Best Iter: [56] ELBO (-3701.081915) evaluations: (3138) 
#> Path [22] :Initial log joint density = -481853.375685 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      8.183e-03   2.714e-01    7.809e-01  7.809e-01      3391 -3.703e+03 -3.714e+03                   
#> Path [22] :Best Iter: [59] ELBO (-3702.915954) evaluations: (3391) 
#> Path [23] :Initial log joint density = -481535.968081 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      2.031e-03   3.144e-01    5.680e-01  5.680e-01      3456 -3.703e+03 -3.712e+03                   
#> Path [23] :Best Iter: [55] ELBO (-3702.696872) evaluations: (3456) 
#> Path [24] :Initial log joint density = -481796.467538 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.057e-02   2.284e-01    1.000e+00  1.000e+00      3009 -3.705e+03 -3.725e+03                   
#> Path [24] :Best Iter: [52] ELBO (-3705.432476) evaluations: (3009) 
#> Path [25] :Initial log joint density = -481611.099513 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.006e-02   2.277e-01    1.000e+00  1.000e+00      3447 -3.705e+03 -3.699e+03                   
#> Path [25] :Best Iter: [59] ELBO (-3699.452570) evaluations: (3447) 
#> Path [26] :Initial log joint density = -481556.642906 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.228e-03   2.212e-01    9.979e-01  9.979e-01      3176 -3.701e+03 -3.709e+03                   
#> Path [26] :Best Iter: [55] ELBO (-3701.455597) evaluations: (3176) 
#> Path [27] :Initial log joint density = -482144.883365 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.353e-03   2.355e-01    1.000e+00  1.000e+00      3242 -3.700e+03 -3.703e+03                   
#> Path [27] :Best Iter: [56] ELBO (-3700.341129) evaluations: (3242) 
#> Path [28] :Initial log joint density = -483598.112459 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      4.775e-03   1.761e-01    1.000e+00  1.000e+00      3565 -3.700e+03 -3.705e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3700.289128) evaluations: (3565) 
#> Path [29] :Initial log joint density = -481539.146412 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.485e-03   2.523e-01    1.000e+00  1.000e+00      2966 -3.706e+03 -3.720e+03                   
#> Path [29] :Best Iter: [45] ELBO (-3706.428949) evaluations: (2966) 
#> Path [30] :Initial log joint density = -485183.306530 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      8.768e-03   2.493e-01    1.000e+00  1.000e+00      3621 -3.699e+03 -3.709e+03                   
#> Path [30] :Best Iter: [58] ELBO (-3699.303220) evaluations: (3621) 
#> Path [31] :Initial log joint density = -482131.585338 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.454e-02   3.162e-01    1.000e+00  1.000e+00      3312 -3.699e+03 -3.704e+03                   
#> Path [31] :Best Iter: [57] ELBO (-3698.536073) evaluations: (3312) 
#> Path [32] :Initial log joint density = -481733.535936 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.732e-03   1.812e-01    9.037e-01  9.037e-01      3290 -3.700e+03 -3.709e+03                   
#> Path [32] :Best Iter: [55] ELBO (-3700.294572) evaluations: (3290) 
#> Path [33] :Initial log joint density = -482730.340083 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.804e-03   1.931e-01    8.771e-01  8.771e-01      3501 -3.704e+03 -3.708e+03                   
#> Path [33] :Best Iter: [57] ELBO (-3704.039835) evaluations: (3501) 
#> Path [34] :Initial log joint density = -481398.129490 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.949e-02   2.299e-01    1.000e+00  1.000e+00      3284 -3.702e+03 -3.701e+03                   
#> Path [34] :Best Iter: [56] ELBO (-3701.347990) evaluations: (3284) 
#> Path [35] :Initial log joint density = -481710.750034 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.338e-03   2.426e-01    9.911e-01  9.911e-01      2913 -3.709e+03 -3.718e+03                   
#> Path [35] :Best Iter: [51] ELBO (-3708.757650) evaluations: (2913) 
#> Path [36] :Initial log joint density = -481616.618758 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.227e-02   2.699e-01    1.000e+00  1.000e+00      2810 -3.709e+03 -3.712e+03                   
#> Path [36] :Best Iter: [46] ELBO (-3708.680265) evaluations: (2810) 
#> Path [37] :Initial log joint density = -484287.006341 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.646e-03   2.249e-01    1.000e+00  1.000e+00      3284 -3.702e+03 -3.704e+03                   
#> Path [37] :Best Iter: [55] ELBO (-3702.470523) evaluations: (3284) 
#> Path [38] :Initial log joint density = -481727.719183 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.269e-02   2.297e-01    1.000e+00  1.000e+00      2874 -3.706e+03 -3.712e+03                   
#> Path [38] :Best Iter: [50] ELBO (-3706.261498) evaluations: (2874) 
#> Path [39] :Initial log joint density = -483835.992828 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.605e-02   3.743e-01    1.000e+00  1.000e+00      3287 -3.699e+03 -3.710e+03                   
#> Path [39] :Best Iter: [56] ELBO (-3698.964021) evaluations: (3287) 
#> Path [40] :Initial log joint density = -481820.968994 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.416e-02   3.489e-01    9.705e-01  9.705e-01      3309 -3.702e+03 -3.712e+03                   
#> Path [40] :Best Iter: [55] ELBO (-3701.517605) evaluations: (3309) 
#> Path [41] :Initial log joint density = -481412.421257 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      5.429e-03   1.881e-01    1.000e+00  1.000e+00      2919 -3.705e+03 -3.717e+03                   
#> Path [41] :Best Iter: [49] ELBO (-3705.274412) evaluations: (2919) 
#> Path [42] :Initial log joint density = -481796.091583 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.337e-02   2.803e-01    1.000e+00  1.000e+00      3270 -3.702e+03 -3.703e+03                   
#> Path [42] :Best Iter: [57] ELBO (-3702.436661) evaluations: (3270) 
#> Path [43] :Initial log joint density = -481562.831985 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      3.634e-03   2.076e-01    7.568e-01  7.568e-01      3190 -3.704e+03 -3.712e+03                   
#> Path [43] :Best Iter: [55] ELBO (-3704.313541) evaluations: (3190) 
#> Path [44] :Initial log joint density = -481410.981700 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.351e-03   2.256e-01    1.000e+00  1.000e+00      3497 -3.698e+03 -3.703e+03                   
#> Path [44] :Best Iter: [56] ELBO (-3697.780960) evaluations: (3497) 
#> Path [45] :Initial log joint density = -482107.669716 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      3.063e-03   2.150e-01    6.963e-01  6.963e-01      3363 -3.700e+03 -3.713e+03                   
#> Path [45] :Best Iter: [57] ELBO (-3700.489172) evaluations: (3363) 
#> Path [46] :Initial log joint density = -481675.959709 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.168e-02   1.734e-01    1.000e+00  1.000e+00      3047 -3.706e+03 -3.705e+03                   
#> Path [46] :Best Iter: [54] ELBO (-3705.372286) evaluations: (3047) 
#> Path [47] :Initial log joint density = -481600.687927 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      6.621e-03   2.813e-01    7.506e-01  7.506e-01      2626 -3.706e+03 -3.724e+03                   
#> Path [47] :Best Iter: [44] ELBO (-3705.501604) evaluations: (2626) 
#> Path [48] :Initial log joint density = -481466.997242 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.566e-03   2.711e-01    1.000e+00  1.000e+00      2986 -3.709e+03 -3.720e+03                   
#> Path [48] :Best Iter: [46] ELBO (-3708.598329) evaluations: (2986) 
#> Path [49] :Initial log joint density = -481669.059102 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      2.435e-03   2.421e-01    7.462e-01  7.462e-01      3071 -3.706e+03 -3.713e+03                   
#> Path [49] :Best Iter: [39] ELBO (-3706.350767) evaluations: (3071) 
#> Path [50] :Initial log joint density = -481664.771845 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.238e-02   2.845e-01    1.000e+00  1.000e+00      3207 -3.707e+03 -3.705e+03                   
#> Path [50] :Best Iter: [56] ELBO (-3705.071724) evaluations: (3207) 
#> Finished in  13.4 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
