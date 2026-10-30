# plot

This function plots a summary of the results of the model.

## Usage

``` r
# S3 method for class 'sccomp_tbl'
plot(
  x,
  significance_threshold = 0.05,
  test_composition_above_CLR_fold_change = attr(x,
    "test_composition_above_CLR_fold_change"),
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

- test_composition_above_CLR_fold_change:

  A positive integer. It is the effect threshold used for the hypothesis
  test. A value of 0.2 correspond to a change in cell proportion of 10%
  for a cell type with baseline proportion of 50%. That is, a cell type
  goes from 45% to 50%. When the baseline proportion is closer to 0 or 1
  this effect thrshold has consistent value in the centered log-ratio
  uncontrained scale.

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
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), typecancer
#>   variability Xa: 1 parameter: (Intercept)
#>   random effects: none
#> Precompiled model not found. Compiling the model...
#> Warning: The "threads" argument of stan_package_compile() was deprecated on 2026-04-07 (instantiate version 0.2.3.9003). Alternative: Set threads in cpp_options, e.g. cpp_options = list(stan_threads = TRUE).
#> Running make /tmp/Rtmpc7h9Sq/model-397212f38c9f "STAN_THREADS=TRUE" \
#>   "STANCFLAGS += --include-paths=/tmp/Rtmpc7h9Sq/temp_libpath39727d588b56/sccomp/stan --name='glm_multi_beta_binomial_model'"
#> 
#> --- Translating Stan model to C++ code ---
#> bin/stanc --include-paths=/tmp/Rtmpc7h9Sq/temp_libpath39727d588b56/sccomp/stan --name='glm_multi_beta_binomial_model' --o=/tmp/Rtmpc7h9Sq/model-397212f38c9f.hpp /tmp/Rtmpc7h9Sq/model-397212f38c9f.stan
#> 
#> --- Compiling C++ code ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -isystem stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -isystem stan/lib/stan_math/lib/eigen_5.0.1 -isystem stan/lib/stan_math/lib/boost_1.87.0 -isystem stan/lib/stan_math/lib/sundials_6.1.1/include -isystem stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS          -c -Wno-ignored-attributes   -x c++ -o /tmp/Rtmpc7h9Sq/model-397212f38c9f.o /tmp/Rtmpc7h9Sq/model-397212f38c9f.hpp
#> 
#> --- Linking model ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -isystem stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -isystem stan/lib/stan_math/lib/eigen_5.0.1 -isystem stan/lib/stan_math/lib/boost_1.87.0 -isystem stan/lib/stan_math/lib/sundials_6.1.1/include -isystem stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS               -Wl,-L,"/home/runner/.cmdstan/cmdstan-2.40.0/stan/lib/stan_math/lib/tbb"   -Wl,-rpath,"/home/runner/.cmdstan/cmdstan-2.40.0/stan/lib/stan_math/lib/tbb"      /tmp/Rtmpc7h9Sq/model-397212f38c9f.o src/cmdstan/main_threads.o       -ltbb   stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_nvecserial.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_cvodes.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_idas.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_kinsol.a  stan/lib/stan_math/lib/tbb/libtbb.so.2 -o /tmp/Rtmpc7h9Sq/model-397212f38c9f
#> rm /tmp/Rtmpc7h9Sq/model-397212f38c9f.hpp /tmp/Rtmpc7h9Sq/model-397212f38c9f.o
#> Model compiled and saved to cache successfully.
#> Path [1] :Initial log joint density = -482668.843522 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      6.010e-03   2.875e-01    1.000e+00  1.000e+00      3610 -3.701e+03 -3.707e+03                   
#> Path [1] :Best Iter: [57] ELBO (-3701.380487) evaluations: (3610) 
#> Path [2] :Initial log joint density = -482547.994646 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      2.572e-03   2.547e-01    5.602e-01  5.602e-01      4093 -3.702e+03 -3.718e+03                   
#> Path [2] :Best Iter: [63] ELBO (-3701.544580) evaluations: (4093) 
#> Path [3] :Initial log joint density = -481732.205120 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.330e-02   1.954e-01    1.000e+00  1.000e+00      4981 -3.704e+03 -3.702e+03                   
#> Path [3] :Best Iter: [74] ELBO (-3702.224466) evaluations: (4981) 
#> Path [4] :Initial log joint density = -483306.273052 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.082e-02   2.366e-01    9.924e-01  9.924e-01      3657 -3.707e+03 -3.711e+03                   
#> Path [4] :Best Iter: [60] ELBO (-3706.924556) evaluations: (3657) 
#> Path [5] :Initial log joint density = -482381.178273 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      2.419e-02   2.862e-01    8.895e-01  8.895e-01      4473 -3.701e+03 -3.710e+03                   
#> Path [5] :Best Iter: [65] ELBO (-3701.397849) evaluations: (4473) 
#> Path [6] :Initial log joint density = -483008.152855 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      2.688e-03   3.222e-01    5.936e-01  5.936e-01      3680 -3.707e+03 -3.718e+03                   
#> Path [6] :Best Iter: [60] ELBO (-3707.382233) evaluations: (3680) 
#> Path [7] :Initial log joint density = -482003.286811 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.231e-02   3.996e-01    9.935e-01  9.935e-01      3624 -3.704e+03 -3.715e+03                   
#> Path [7] :Best Iter: [61] ELBO (-3704.066766) evaluations: (3624) 
#> Path [8] :Initial log joint density = -482354.641520 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      8.920e-03   2.749e-01    9.117e-01  9.117e-01      3940 -3.702e+03 -3.711e+03                   
#> Path [8] :Best Iter: [62] ELBO (-3702.287527) evaluations: (3940) 
#> Path [9] :Initial log joint density = -482041.641793 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.322e-02   2.877e-01    1.000e+00  1.000e+00      4532 -3.702e+03 -3.706e+03                   
#> Path [9] :Best Iter: [66] ELBO (-3701.719866) evaluations: (4532) 
#> Path [10] :Initial log joint density = -481392.119848 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      5.900e-03   1.729e-01    1.000e+00  1.000e+00      4057 -3.702e+03 -3.707e+03                   
#> Path [10] :Best Iter: [62] ELBO (-3701.643886) evaluations: (4057) 
#> Path [11] :Initial log joint density = -481509.841516 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      8.429e-03   2.615e-01    9.339e-01  9.339e-01      4446 -3.699e+03 -3.710e+03                   
#> Path [11] :Best Iter: [68] ELBO (-3699.230260) evaluations: (4446) 
#> Path [12] :Initial log joint density = -481443.807429 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      5.226e-03   2.081e-01    6.664e-01  6.664e-01      3724 -3.702e+03 -3.714e+03                   
#> Path [12] :Best Iter: [56] ELBO (-3702.195471) evaluations: (3724) 
#> Path [13] :Initial log joint density = -481694.310495 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      6.889e-03   1.886e-01    7.281e-01  7.281e-01      4715 -3.697e+03 -3.705e+03                   
#> Path [13] :Best Iter: [70] ELBO (-3696.803686) evaluations: (4715) 
#> Path [14] :Initial log joint density = -481484.750162 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.214e-02   2.494e-01    9.036e-01  9.036e-01      5115 -3.701e+03 -3.714e+03                   
#> Path [14] :Best Iter: [73] ELBO (-3700.804616) evaluations: (5115) 
#> Path [15] :Initial log joint density = -481806.542922 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      6.593e-03   2.889e-01    6.430e-01  6.430e-01      5172 -3.700e+03 -3.711e+03                   
#> Path [15] :Best Iter: [75] ELBO (-3699.675591) evaluations: (5172) 
#> Path [16] :Initial log joint density = -481719.014292 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.223e-02   1.991e-01    1.000e+00  1.000e+00      4773 -3.700e+03 -3.703e+03                   
#> Path [16] :Best Iter: [68] ELBO (-3700.395502) evaluations: (4773) 
#> Path [17] :Initial log joint density = -481801.086303 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      2.549e-02   3.440e-01    1.000e+00  1.000e+00      4407 -3.701e+03 -3.711e+03                   
#> Path [17] :Best Iter: [68] ELBO (-3700.562508) evaluations: (4407) 
#> Path [18] :Initial log joint density = -481514.577420 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      2.969e-02   3.301e-01    1.000e+00  1.000e+00      3420 -3.701e+03 -3.711e+03                   
#> Path [18] :Best Iter: [57] ELBO (-3700.931084) evaluations: (3420) 
#> Path [19] :Initial log joint density = -481314.105809 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      2.520e-02   2.757e-01    1.000e+00  1.000e+00      4341 -3.701e+03 -3.709e+03                   
#> Path [19] :Best Iter: [68] ELBO (-3701.336799) evaluations: (4341) 
#> Path [20] :Initial log joint density = -483374.272613 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      1.170e-02   3.205e-01    1.000e+00  1.000e+00      3600 -3.704e+03 -3.712e+03                   
#> Path [20] :Best Iter: [58] ELBO (-3704.011218) evaluations: (3600) 
#> Path [21] :Initial log joint density = -481609.516037 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      2.066e-02   2.559e-01    8.647e-01  8.647e-01      4503 -3.703e+03 -3.710e+03                   
#> Path [21] :Best Iter: [70] ELBO (-3702.707299) evaluations: (4503) 
#> Path [22] :Initial log joint density = -481825.545707 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      5.883e-03   1.719e-01    8.946e-01  8.946e-01      4360 -3.703e+03 -3.714e+03                   
#> Path [22] :Best Iter: [67] ELBO (-3702.594009) evaluations: (4360) 
#> Path [23] :Initial log joint density = -481456.946909 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      9.019e-03   1.484e-01    1.000e+00  1.000e+00      4149 -3.703e+03 -3.712e+03                   
#> Path [23] :Best Iter: [55] ELBO (-3702.724070) evaluations: (4149) 
#> Path [24] :Initial log joint density = -481668.965553 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.443e-02   2.032e-01    1.000e+00  1.000e+00      4250 -3.702e+03 -3.706e+03                   
#> Path [24] :Best Iter: [58] ELBO (-3701.767142) evaluations: (4250) 
#> Path [25] :Initial log joint density = -481582.024939 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      2.174e-02   2.911e-01    4.685e-01  1.000e+00      4900 -3.701e+03 -3.711e+03                   
#> Path [25] :Best Iter: [69] ELBO (-3700.956670) evaluations: (4900) 
#> Path [26] :Initial log joint density = -481454.388897 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.350e-02   2.070e-01    1.000e+00  1.000e+00      4132 -3.705e+03 -3.703e+03                   
#> Path [26] :Best Iter: [66] ELBO (-3703.254958) evaluations: (4132) 
#> Path [27] :Initial log joint density = -482107.232775 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      2.424e-02   1.361e-01    1.000e+00  1.000e+00      4840 -3.698e+03 -3.702e+03                   
#> Path [27] :Best Iter: [68] ELBO (-3698.269651) evaluations: (4840) 
#> Path [28] :Initial log joint density = -483561.558517 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      7.591e-03   2.774e-01    1.000e+00  1.000e+00      3400 -3.705e+03 -3.709e+03                   
#> Path [28] :Best Iter: [56] ELBO (-3705.453213) evaluations: (3400) 
#> Path [29] :Initial log joint density = -481503.180672 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.978e-02   2.993e-01    1.000e+00  1.000e+00      4208 -3.701e+03 -3.704e+03                   
#> Path [29] :Best Iter: [65] ELBO (-3700.576809) evaluations: (4208) 
#> Path [30] :Initial log joint density = -485149.723514 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.675e-02   3.801e-01    1.000e+00  1.000e+00      4187 -3.700e+03 -3.709e+03                   
#> Path [30] :Best Iter: [65] ELBO (-3700.439192) evaluations: (4187) 
#> Path [31] :Initial log joint density = -482096.710937 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.571e-02   1.967e-01    1.000e+00  1.000e+00      5376 -3.699e+03 -3.704e+03                   
#> Path [31] :Best Iter: [72] ELBO (-3699.303645) evaluations: (5376) 
#> Path [32] :Initial log joint density = -481686.382658 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      6.368e-03   2.039e-01    7.848e-01  7.848e-01      4014 -3.701e+03 -3.714e+03                   
#> Path [32] :Best Iter: [65] ELBO (-3700.549630) evaluations: (4014) 
#> Path [33] :Initial log joint density = -482692.920844 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.341e-02   3.482e-01    1.000e+00  1.000e+00      3725 -3.701e+03 -3.711e+03                   
#> Path [33] :Best Iter: [60] ELBO (-3701.341628) evaluations: (3725) 
#> Path [34] :Initial log joint density = -481356.308585 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      6.500e-03   1.654e-01    1.000e+00  1.000e+00      3718 -3.699e+03 -3.704e+03                   
#> Path [34] :Best Iter: [59] ELBO (-3698.988471) evaluations: (3718) 
#> Path [35] :Initial log joint density = -481670.690009 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      8.126e-03   2.652e-01    5.834e-01  5.834e-01      4017 -3.703e+03 -3.709e+03                   
#> Path [35] :Best Iter: [58] ELBO (-3702.582492) evaluations: (4017) 
#> Path [36] :Initial log joint density = -481584.074715 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.171e-02   1.778e-01    1.000e+00  1.000e+00      3710 -3.701e+03 -3.704e+03                   
#> Path [36] :Best Iter: [59] ELBO (-3700.895095) evaluations: (3710) 
#> Path [37] :Initial log joint density = -484252.083022 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      6.734e-03   2.599e-01    7.538e-01  7.538e-01      4182 -3.701e+03 -3.714e+03                   
#> Path [37] :Best Iter: [65] ELBO (-3701.046708) evaluations: (4182) 
#> Path [38] :Initial log joint density = -481623.199422 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.512e-02   1.626e-01    1.000e+00  1.000e+00      3686 -3.701e+03 -3.700e+03                   
#> Path [38] :Best Iter: [62] ELBO (-3700.348715) evaluations: (3686) 
#> Path [39] :Initial log joint density = -483797.822835 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      2.077e-02   2.785e-01    4.101e-01  1.000e+00      4057 -3.700e+03 -3.712e+03                   
#> Path [39] :Best Iter: [64] ELBO (-3699.946789) evaluations: (4057) 
#> Path [40] :Initial log joint density = -481767.664046 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      5.246e-03   1.961e-01    8.505e-01  8.505e-01      4183 -3.698e+03 -3.710e+03                   
#> Path [40] :Best Iter: [66] ELBO (-3698.169261) evaluations: (4183) 
#> Path [41] :Initial log joint density = -481350.348589 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      6.043e-03   2.014e-01    1.000e+00  1.000e+00      3040 -3.710e+03 -3.715e+03                   
#> Path [41] :Best Iter: [51] ELBO (-3710.244115) evaluations: (3040) 
#> Path [42] :Initial log joint density = -481754.141679 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      6.695e-03   2.297e-01    7.351e-01  7.351e-01      4151 -3.702e+03 -3.711e+03                   
#> Path [42] :Best Iter: [64] ELBO (-3702.047903) evaluations: (4151) 
#> Path [43] :Initial log joint density = -481513.376213 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      1.018e-02   2.278e-01    9.570e-01  9.570e-01      3479 -3.706e+03 -3.715e+03                   
#> Path [43] :Best Iter: [57] ELBO (-3706.203835) evaluations: (3479) 
#> Path [44] :Initial log joint density = -481344.652153 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      8.439e-03   2.331e-01    1.000e+00  1.000e+00      4505 -3.700e+03 -3.703e+03                   
#> Path [44] :Best Iter: [66] ELBO (-3700.315194) evaluations: (4505) 
#> Path [45] :Initial log joint density = -482074.552173 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.953e-02   2.533e-01    1.000e+00  1.000e+00      4978 -3.699e+03 -3.708e+03                   
#> Path [45] :Best Iter: [74] ELBO (-3698.980654) evaluations: (4978) 
#> Path [46] :Initial log joint density = -481630.810617 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      6.555e-03   1.567e-01    1.000e+00  1.000e+00      3962 -3.700e+03 -3.711e+03                   
#> Path [46] :Best Iter: [62] ELBO (-3700.405705) evaluations: (3962) 
#> Path [47] :Initial log joint density = -481492.900042 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.488e-02   1.513e-01    1.000e+00  1.000e+00      4091 -3.702e+03 -3.705e+03                   
#> Path [47] :Best Iter: [64] ELBO (-3702.415377) evaluations: (4091) 
#> Path [48] :Initial log joint density = -481424.382881 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      8.963e-03   2.082e-01    9.387e-01  9.387e-01      4535 -3.699e+03 -3.707e+03                   
#> Path [48] :Best Iter: [69] ELBO (-3699.049711) evaluations: (4535) 
#> Path [49] :Initial log joint density = -481625.839119 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      5.962e-03   2.055e-01    8.480e-01  8.480e-01      3687 -3.703e+03 -3.712e+03                   
#> Path [49] :Best Iter: [60] ELBO (-3702.624824) evaluations: (3687) 
#> Path [50] :Initial log joint density = -481612.295489 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.124e-02   1.818e-01    1.000e+00  1.000e+00      5119 -3.701e+03 -3.703e+03                   
#> Path [50] :Best Iter: [75] ELBO (-3701.065805) evaluations: (5119) 
#> Finished in  14.5 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Joining with `by = join_by(cell_group, M, parameter)`
#> sccomp says: When visualising proportions, especially for complex models, consider setting `remove_unwanted_effects=TRUE`. This will adjust the proportions, preserving only the observed effect.
#> sccomp says: from version 2.1.25, the default `significance_statistic` for boxplots is `pH0` (previously `FDR`). Set `significance_statistic = "FDR"` to use the previous default.
#> Precompiled model not found. Compiling the model...
#> Warning: The "threads" argument of stan_package_compile() was deprecated on 2026-04-07 (instantiate version 0.2.3.9003). Alternative: Set threads in cpp_options, e.g. cpp_options = list(stan_threads = TRUE).
#> Running make /tmp/Rtmpc7h9Sq/model-397254ee2871 "STAN_THREADS=TRUE" \
#>   "STANCFLAGS += --include-paths=/tmp/Rtmpc7h9Sq/temp_libpath39727d588b56/sccomp/stan --name='glm_multi_beta_binomial_generate_data_model'"
#> 
#> --- Translating Stan model to C++ code ---
#> bin/stanc --include-paths=/tmp/Rtmpc7h9Sq/temp_libpath39727d588b56/sccomp/stan --name='glm_multi_beta_binomial_generate_data_model' --o=/tmp/Rtmpc7h9Sq/model-397254ee2871.hpp /tmp/Rtmpc7h9Sq/model-397254ee2871.stan
#> 
#> --- Compiling C++ code ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -isystem stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -isystem stan/lib/stan_math/lib/eigen_5.0.1 -isystem stan/lib/stan_math/lib/boost_1.87.0 -isystem stan/lib/stan_math/lib/sundials_6.1.1/include -isystem stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS          -c -Wno-ignored-attributes   -x c++ -o /tmp/Rtmpc7h9Sq/model-397254ee2871.o /tmp/Rtmpc7h9Sq/model-397254ee2871.hpp
#> 
#> --- Linking model ---
#> g++ -Wno-deprecated-declarations -std=c++17 -pthread -D_REENTRANT -Wno-sign-compare -Wno-ignored-attributes -Wno-class-memaccess     -DSTAN_THREADS -isystem stan/lib/stan_math/lib/tbb_2020.3/include    -O3 -I src -I stan/src -I stan/lib/rapidjson_1.1.0/ -I lib/CLI11-1.9.1/ -I stan/lib/stan_math/ -isystem stan/lib/stan_math/lib/eigen_5.0.1 -isystem stan/lib/stan_math/lib/boost_1.87.0 -isystem stan/lib/stan_math/lib/sundials_6.1.1/include -isystem stan/lib/stan_math/lib/sundials_6.1.1/src/sundials    -DBOOST_DISABLE_ASSERTS               -Wl,-L,"/home/runner/.cmdstan/cmdstan-2.40.0/stan/lib/stan_math/lib/tbb"   -Wl,-rpath,"/home/runner/.cmdstan/cmdstan-2.40.0/stan/lib/stan_math/lib/tbb"      /tmp/Rtmpc7h9Sq/model-397254ee2871.o src/cmdstan/main_threads.o       -ltbb   stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_nvecserial.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_cvodes.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_idas.a stan/lib/stan_math/lib/sundials_6.1.1/lib/libsundials_kinsol.a  stan/lib/stan_math/lib/tbb/libtbb.so.2 -o /tmp/Rtmpc7h9Sq/model-397254ee2871
#> rm /tmp/Rtmpc7h9Sq/model-397254ee2871.hpp /tmp/Rtmpc7h9Sq/model-397254ee2871.o
#> Model compiled and saved to cache successfully.
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.103 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> Joining with `by = join_by(cell_group, sample)`
#> Joining with `by = join_by(cell_group, type)`
#> Warning: Ignoring unknown parameters: `median.linewidth`
#> === Single Model Parameters ===
#> 
#> (Intercept):
#>   v = -(5.632 + -0.568 × c)
#> 
# }
```
