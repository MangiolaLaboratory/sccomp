# plot

This function plots a summary of the results of the model.

## Usage

``` r
# S3 method for class 'sccomp_tbl'
plot(
  x,
  significance_threshold = 0.05,
  test_composition_above_logit_fold_change = attr(x,
    "test_composition_above_logit_fold_change"),
  significance_statistic = c("pH0", "FDR"),
  show_fdr_message = TRUE,
  add_marginal_density = TRUE,
  omit_ci = FALSE,
  sort_by = c("none", "effect", "significance", "alphabetical"),
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
  differences. Default is 0.05.

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

- add_marginal_density:

  Logical. Whether to add marginal density plots on adjusted panels in
  2D intervals. Default is TRUE.

- omit_ci:

  Logical. Whether to omit credible interval error bars from 2D interval
  plots. Default is FALSE.

- sort_by:

  Character vector indicating how to sort taxa. Options are "none"
  (default), "effect" (by effect size), "significance" (by FDR/pH0), or
  "alphabetical".

- ...:

  For internal use

## Value

A list containing ggplot objects

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
    ) |>
    sccomp_test()

    plots = estimate |> plot()
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Precompiled model not found. Compiling the model...
#> Running make /tmp/Rtmp9NpqYX/model-3c5f450f6f9f "STAN_THREADS=TRUE" \
#>   "STANCFLAGS += --include-paths=/tmp/Rtmp9NpqYX/temp_libpath3c5f51bcedf7/sccomp/stan --name='glm_multi_beta_binomial_model'"
#> 
#> --- Translating Stan model to C++ code ---
#> bin/stanc --include-paths=/tmp/Rtmp9NpqYX/temp_libpath3c5f51bcedf7/sccomp/stan --name='glm_multi_beta_binomial_model' --o=/tmp/Rtmp9NpqYX/model-3c5f450f6f9f.hpp /tmp/Rtmp9NpqYX/model-3c5f450f6f9f.stan
#> 
#> --- Compiling C++ code ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -I stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -I stan/lib/stan_math/lib/eigen_3.4.0 -I stan/lib/stan_math/lib/boost_1.87.0 -I stan/lib/stan_math/lib/sundials_6.1.1/include -I stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS          -c -Wno-ignored-attributes   -x c++ -o /tmp/Rtmp9NpqYX/model-3c5f450f6f9f.o /tmp/Rtmp9NpqYX/model-3c5f450f6f9f.hpp
#> 
#> --- Linking model ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -I stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -I stan/lib/stan_math/lib/eigen_3.4.0 -I stan/lib/stan_math/lib/boost_1.87.0 -I stan/lib/stan_math/lib/sundials_6.1.1/include -I stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS               -Wl,-L,"/home/runner/.cmdstan/cmdstan-2.39.0/stan/lib/stan_math/lib/tbb"   -Wl,-rpath,"/home/runner/.cmdstan/cmdstan-2.39.0/stan/lib/stan_math/lib/tbb"      /tmp/Rtmp9NpqYX/model-3c5f450f6f9f.o src/cmdstan/main_threads.o       -ltbb   stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_nvecserial.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_cvodes.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_idas.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_kinsol.a  stan/lib/stan_math/lib/tbb/libtbb.so.2 -o /tmp/Rtmp9NpqYX/model-3c5f450f6f9f
#> rm /tmp/Rtmp9NpqYX/model-3c5f450f6f9f.hpp /tmp/Rtmp9NpqYX/model-3c5f450f6f9f.o
#> Model compiled and saved to cache successfully.
#> Path [1] :Initial log joint density = -482684.608336 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      4.975e-03   1.654e-01    1.000e+00  1.000e+00      3647 -3.686e+03 -3.700e+03                   
#> Path [1] :Best Iter: [57] ELBO (-3685.549632) evaluations: (3647) 
#> Path [2] :Initial log joint density = -482561.457083 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      9.095e-03   2.185e-01    1.000e+00  1.000e+00      3607 -3.688e+03 -3.690e+03                   
#> Path [2] :Best Iter: [57] ELBO (-3688.007264) evaluations: (3607) 
#> Path [3] :Initial log joint density = -481749.189963 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      9.061e-03   2.076e-01    1.000e+00  1.000e+00      4488 -3.689e+03 -3.690e+03                   
#> Path [3] :Best Iter: [68] ELBO (-3689.253919) evaluations: (4488) 
#> Path [4] :Initial log joint density = -483320.748753 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      7.086e-03   2.663e-01    1.000e+00  1.000e+00      3527 -3.690e+03 -3.699e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3690.190567) evaluations: (3527) 
#> Path [5] :Initial log joint density = -482399.422218 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      6.119e-03   2.000e-01    1.000e+00  1.000e+00      3930 -3.689e+03 -3.696e+03                   
#> Path [5] :Best Iter: [62] ELBO (-3688.734423) evaluations: (3930) 
#> Path [6] :Initial log joint density = -483019.809595 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      8.618e-03   3.440e-01    4.259e-01  1.000e+00      3505 -3.686e+03 -3.699e+03                   
#> Path [6] :Best Iter: [58] ELBO (-3685.765876) evaluations: (3505) 
#> Path [7] :Initial log joint density = -482018.264036 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      3.885e-03   1.400e-01    1.000e+00  1.000e+00      3906 -3.689e+03 -3.700e+03                   
#> Path [7] :Best Iter: [56] ELBO (-3688.832008) evaluations: (3906) 
#> Path [8] :Initial log joint density = -482374.997006 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      6.576e-03   2.588e-01    1.000e+00  1.000e+00      3290 -3.692e+03 -3.689e+03                   
#> Path [8] :Best Iter: [57] ELBO (-3688.887758) evaluations: (3290) 
#> Path [9] :Initial log joint density = -482059.982656 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      8.454e-03   2.236e-01    1.000e+00  1.000e+00      4693 -3.682e+03 -3.692e+03                   
#> Path [9] :Best Iter: [69] ELBO (-3682.237067) evaluations: (4693) 
#> Path [10] :Initial log joint density = -481409.068382 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      1.637e-02   2.779e-01    1.000e+00  1.000e+00      3335 -3.694e+03 -3.691e+03                   
#> Path [10] :Best Iter: [57] ELBO (-3691.140906) evaluations: (3335) 
#> Path [11] :Initial log joint density = -481523.925214 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.510e-02   2.129e-01    1.000e+00  1.000e+00      4890 -3.682e+03 -3.687e+03                   
#> Path [11] :Best Iter: [67] ELBO (-3681.987219) evaluations: (4890) 
#> Path [12] :Initial log joint density = -481457.463724 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      3.743e-02   2.417e-01    1.000e+00  1.000e+00      3830 -3.687e+03 -3.691e+03                   
#> Path [12] :Best Iter: [60] ELBO (-3686.748651) evaluations: (3830) 
#> Path [13] :Initial log joint density = -481707.229458 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.474e-02   2.203e-01    1.000e+00  1.000e+00      3927 -3.684e+03 -3.693e+03                   
#> Path [13] :Best Iter: [61] ELBO (-3683.676035) evaluations: (3927) 
#> Path [14] :Initial log joint density = -481496.498206 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      7.280e-03   2.011e-01    1.000e+00  1.000e+00      4537 -3.689e+03 -3.693e+03                   
#> Path [14] :Best Iter: [63] ELBO (-3688.720073) evaluations: (4537) 
#> Path [15] :Initial log joint density = -481822.682912 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      7.476e-03   2.028e-01    8.173e-01  8.173e-01      4509 -3.685e+03 -3.697e+03                   
#> Path [15] :Best Iter: [69] ELBO (-3684.847877) evaluations: (4509) 
#> Path [16] :Initial log joint density = -481734.465675 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      3.635e-02   2.471e-01    9.530e-01  9.530e-01      4578 -3.685e+03 -3.691e+03                   
#> Path [16] :Best Iter: [68] ELBO (-3684.928170) evaluations: (4578) 
#> Path [17] :Initial log joint density = -481815.576551 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.879e-02   3.479e-01    1.000e+00  1.000e+00      4054 -3.682e+03 -3.691e+03                   
#> Path [17] :Best Iter: [64] ELBO (-3682.072810) evaluations: (4054) 
#> Path [18] :Initial log joint density = -481527.793631 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.787e+05      1.534e-02   2.973e-01    7.559e-01  7.559e-01      2952 -3.695e+03 -3.710e+03                   
#> Path [18] :Best Iter: [47] ELBO (-3694.852286) evaluations: (2952) 
#> Path [19] :Initial log joint density = -481324.374543 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.581e-02   2.709e-01    1.000e+00  1.000e+00      4268 -3.688e+03 -3.695e+03                   
#> Path [19] :Best Iter: [68] ELBO (-3688.478955) evaluations: (4268) 
#> Path [20] :Initial log joint density = -483384.861001 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      1.258e-02   3.231e-01    1.000e+00  1.000e+00      3277 -3.687e+03 -3.692e+03                   
#> Path [20] :Best Iter: [57] ELBO (-3687.093954) evaluations: (3277) 
#> Path [21] :Initial log joint density = -481627.944841 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.729e-02   2.518e-01    1.000e+00  1.000e+00      4274 -3.687e+03 -3.693e+03                   
#> Path [21] :Best Iter: [59] ELBO (-3686.580638) evaluations: (4274) 
#> Path [22] :Initial log joint density = -481840.378006 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.129e-02   1.835e-01    1.000e+00  1.000e+00      4602 -3.689e+03 -3.688e+03                   
#> Path [22] :Best Iter: [72] ELBO (-3688.191868) evaluations: (4602) 
#> Path [23] :Initial log joint density = -481474.490297 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      5.312e-03   2.114e-01    1.000e+00  1.000e+00      3792 -3.686e+03 -3.698e+03                   
#> Path [23] :Best Iter: [59] ELBO (-3686.361993) evaluations: (3792) 
#> Path [24] :Initial log joint density = -481685.694033 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      9.011e-03   2.096e-01    1.000e+00  1.000e+00      4229 -3.685e+03 -3.688e+03                   
#> Path [24] :Best Iter: [58] ELBO (-3684.895413) evaluations: (4229) 
#> Path [25] :Initial log joint density = -481593.036311 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.059e-02   2.588e-01    1.000e+00  1.000e+00      4584 -3.686e+03 -3.694e+03                   
#> Path [25] :Best Iter: [59] ELBO (-3686.084079) evaluations: (4584) 
#> Path [26] :Initial log joint density = -481468.123512 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.397e-02   1.991e-01    1.000e+00  1.000e+00      4320 -3.687e+03 -3.692e+03                   
#> Path [26] :Best Iter: [66] ELBO (-3686.863068) evaluations: (4320) 
#> Path [27] :Initial log joint density = -482124.310358 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.468e-02   2.620e-01    1.000e+00  1.000e+00      3964 -3.686e+03 -3.691e+03                   
#> Path [27] :Best Iter: [61] ELBO (-3686.328287) evaluations: (3964) 
#> Path [28] :Initial log joint density = -483580.845892 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      6.859e-03   2.130e-01    9.962e-01  9.962e-01      3577 -3.688e+03 -3.697e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3688.369555) evaluations: (3577) 
#> Path [29] :Initial log joint density = -481516.987702 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.463e-02   2.083e-01    9.644e-01  9.644e-01      4279 -3.682e+03 -3.690e+03                   
#> Path [29] :Best Iter: [63] ELBO (-3682.254477) evaluations: (4279) 
#> Path [30] :Initial log joint density = -485163.608998 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      6.165e-03   1.612e-01    9.780e-01  9.780e-01      4321 -3.684e+03 -3.695e+03                   
#> Path [30] :Best Iter: [66] ELBO (-3684.379081) evaluations: (4321) 
#> Path [31] :Initial log joint density = -482114.883887 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.067e-02   2.287e-01    1.000e+00  1.000e+00      4570 -3.687e+03 -3.694e+03                   
#> Path [31] :Best Iter: [69] ELBO (-3687.291662) evaluations: (4570) 
#> Path [32] :Initial log joint density = -481698.523713 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      1.012e-02   2.423e-01    1.000e+00  1.000e+00      3281 -3.687e+03 -3.693e+03                   
#> Path [32] :Best Iter: [57] ELBO (-3687.229066) evaluations: (3281) 
#> Path [33] :Initial log joint density = -482712.657755 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.160e-02   2.371e-01    1.000e+00  1.000e+00      3756 -3.687e+03 -3.692e+03                   
#> Path [33] :Best Iter: [60] ELBO (-3687.481926) evaluations: (3756) 
#> Path [34] :Initial log joint density = -481370.061339 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      1.273e-02   2.371e-01    1.000e+00  1.000e+00      3026 -3.694e+03 -3.703e+03                   
#> Path [34] :Best Iter: [46] ELBO (-3693.858742) evaluations: (3026) 
#> Path [35] :Initial log joint density = -481683.980603 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.296e-02   2.418e-01    1.000e+00  1.000e+00      4013 -3.684e+03 -3.688e+03                   
#> Path [35] :Best Iter: [63] ELBO (-3684.061613) evaluations: (4013) 
#> Path [36] :Initial log joint density = -481603.277248 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      9.016e-03   2.232e-01    1.000e+00  1.000e+00      3270 -3.688e+03 -3.690e+03                   
#> Path [36] :Best Iter: [55] ELBO (-3688.325077) evaluations: (3270) 
#> Path [37] :Initial log joint density = -484269.036972 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.170e-02   1.713e-01    8.561e-01  8.561e-01      4475 -3.686e+03 -3.697e+03                   
#> Path [37] :Best Iter: [63] ELBO (-3685.606470) evaluations: (4475) 
#> Path [38] :Initial log joint density = -481638.597791 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      1.509e-02   3.169e-01    1.000e+00  1.000e+00      3388 -3.688e+03 -3.691e+03                   
#> Path [38] :Best Iter: [55] ELBO (-3687.971376) evaluations: (3388) 
#> Path [39] :Initial log joint density = -483817.787059 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      8.223e-03   2.175e-01    1.000e+00  1.000e+00      4274 -3.684e+03 -3.696e+03                   
#> Path [39] :Best Iter: [64] ELBO (-3684.293437) evaluations: (4274) 
#> Path [40] :Initial log joint density = -481782.150644 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.520e-02   1.923e-01    8.023e-01  8.023e-01      4649 -3.684e+03 -3.696e+03                   
#> Path [40] :Best Iter: [70] ELBO (-3684.375738) evaluations: (4649) 
#> Path [41] :Initial log joint density = -481364.515649 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      9.273e-03   2.429e-01    1.000e+00  1.000e+00      3005 -3.693e+03 -3.704e+03                   
#> Path [41] :Best Iter: [51] ELBO (-3692.596145) evaluations: (3005) 
#> Path [42] :Initial log joint density = -481772.790055 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.485e-02   2.041e-01    1.000e+00  1.000e+00      3596 -3.690e+03 -3.694e+03                   
#> Path [42] :Best Iter: [59] ELBO (-3689.985023) evaluations: (3596) 
#> Path [43] :Initial log joint density = -481528.626810 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      9.730e-03   2.784e-01    7.574e-01  7.574e-01      3251 -3.691e+03 -3.706e+03                   
#> Path [43] :Best Iter: [55] ELBO (-3691.404418) evaluations: (3251) 
#> Path [44] :Initial log joint density = -481357.142392 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      7.798e-03   1.619e-01    1.000e+00  1.000e+00      4486 -3.684e+03 -3.689e+03                   
#> Path [44] :Best Iter: [62] ELBO (-3684.429712) evaluations: (4486) 
#> Path [45] :Initial log joint density = -482092.476441 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      6.198e-03   2.051e-01    1.000e+00  1.000e+00      4611 -3.688e+03 -3.694e+03                   
#> Path [45] :Best Iter: [64] ELBO (-3687.757854) evaluations: (4611) 
#> Path [46] :Initial log joint density = -481644.339052 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.802e-02   2.480e-01    9.323e-01  9.323e-01      4183 -3.684e+03 -3.697e+03                   
#> Path [46] :Best Iter: [63] ELBO (-3684.389588) evaluations: (4183) 
#> Path [47] :Initial log joint density = -481506.371585 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.493e-02   2.401e-01    1.000e+00  1.000e+00      3965 -3.686e+03 -3.685e+03                   
#> Path [47] :Best Iter: [65] ELBO (-3685.498397) evaluations: (3965) 
#> Path [48] :Initial log joint density = -481438.352186 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      2.003e-02   1.814e-01    1.000e+00  1.000e+00      4739 -3.683e+03 -3.686e+03                   
#> Path [48] :Best Iter: [69] ELBO (-3683.477908) evaluations: (4739) 
#> Path [49] :Initial log joint density = -481644.254142 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.149e-02   2.204e-01    8.773e-01  8.773e-01      3823 -3.686e+03 -3.699e+03                   
#> Path [49] :Best Iter: [62] ELBO (-3686.345555) evaluations: (3823) 
#> Path [50] :Initial log joint density = -481632.790800 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.231e-02   2.136e-01    1.000e+00  1.000e+00      5085 -3.682e+03 -3.692e+03                   
#> Path [50] :Best Iter: [69] ELBO (-3682.357712) evaluations: (5085) 
#> Finished in  15.6 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Joining with `by = join_by(cell_group, M, parameter)`
#> sccomp says: When visualising proportions, especially for complex models, consider setting `remove_unwanted_effects=TRUE`. This will adjust the proportions, preserving only the observed effect.
#> sccomp says: from version 2.1.25, the default `significance_statistic` for boxplots is `pH0` (previously `FDR`). Set `significance_statistic = "FDR"` to use the previous default.
#> Precompiled model not found. Compiling the model...
#> Running make /tmp/Rtmp9NpqYX/model-3c5f27ebe699 "STAN_THREADS=TRUE" \
#>   "STANCFLAGS += --include-paths=/tmp/Rtmp9NpqYX/temp_libpath3c5f51bcedf7/sccomp/stan --name='glm_multi_beta_binomial_generate_data_model'"
#> 
#> --- Translating Stan model to C++ code ---
#> bin/stanc --include-paths=/tmp/Rtmp9NpqYX/temp_libpath3c5f51bcedf7/sccomp/stan --name='glm_multi_beta_binomial_generate_data_model' --o=/tmp/Rtmp9NpqYX/model-3c5f27ebe699.hpp /tmp/Rtmp9NpqYX/model-3c5f27ebe699.stan
#> 
#> --- Compiling C++ code ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -I stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -I stan/lib/stan_math/lib/eigen_3.4.0 -I stan/lib/stan_math/lib/boost_1.87.0 -I stan/lib/stan_math/lib/sundials_6.1.1/include -I stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS          -c -Wno-ignored-attributes   -x c++ -o /tmp/Rtmp9NpqYX/model-3c5f27ebe699.o /tmp/Rtmp9NpqYX/model-3c5f27ebe699.hpp
#> 
#> --- Linking model ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -I stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -I stan/lib/stan_math/lib/eigen_3.4.0 -I stan/lib/stan_math/lib/boost_1.87.0 -I stan/lib/stan_math/lib/sundials_6.1.1/include -I stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS               -Wl,-L,"/home/runner/.cmdstan/cmdstan-2.39.0/stan/lib/stan_math/lib/tbb"   -Wl,-rpath,"/home/runner/.cmdstan/cmdstan-2.39.0/stan/lib/stan_math/lib/tbb"      /tmp/Rtmp9NpqYX/model-3c5f27ebe699.o src/cmdstan/main_threads.o       -ltbb   stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_nvecserial.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_cvodes.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_idas.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_kinsol.a  stan/lib/stan_math/lib/tbb/libtbb.so.2 -o /tmp/Rtmp9NpqYX/model-3c5f27ebe699
#> rm /tmp/Rtmp9NpqYX/model-3c5f27ebe699.o /tmp/Rtmp9NpqYX/model-3c5f27ebe699.hpp
#> Model compiled and saved to cache successfully.
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.102 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> Joining with `by = join_by(cell_group, sample)`
#> Joining with `by = join_by(cell_group, type)`
#> Warning: Ignoring unknown parameters: `median.linewidth`
#> === Single Model Parameters ===
#> 
#> (Intercept):
#>   v = -(5.613 + -0.563 × c)
#> 
# }
```
