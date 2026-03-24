# sccomp_predict

This function replicates counts from a real-world dataset.

## Usage

``` r
sccomp_predict(
  fit,
  formula_composition = NULL,
  new_data = NULL,
  number_of_draws = 500,
  mcmc_seed = sample_seed(),
  summary_instead_of_draws = TRUE,
  robust = FALSE
)
```

## Arguments

- fit:

  The result of sccomp_estimate.

- formula_composition:

  A formula. The formula describing the model for differential
  abundance, for example ~treatment. This formula can be a sub-formula
  of your estimated model; in this case all other factor will be
  factored out.

- new_data:

  A sample-wise data frame including the column that represent the
  factors in your formula. If you want to predict proportions for 10
  samples, there should be 10 rows. T

- number_of_draws:

  An integer. How may copies of the data you want to draw from the model
  joint posterior distribution.

- mcmc_seed:

  An integer. Used for Markov-chain Monte Carlo reproducibility. By
  default a random number is sampled from 1 to 999999. This itself can
  be controlled by set.seed()

- summary_instead_of_draws:

  Return the summary values (i.e. mean and quantiles) of the predicted
  proportions, or return single draws. Single draws can be helful to
  better analyse the uncertainty of the prediction.

- robust:

  A logical. If TRUE, use robust statistics (median and median absolute
  deviation) instead of classical statistics (mean and standard
  deviation) for the summary calculations.

## Value

A tibble (`tbl`) with the following columns:

- **cell_group** - A character column representing the cell group being
  tested.

- **sample** - A factor column representing the sample name for which
  the predictions are made.

- **proportion_mean** - A numeric column representing the predicted mean
  (or median when robust=TRUE) proportions from the model.

- **proportion_lower** - A numeric column representing the lower bound
  (2.5%) of the 95% credible interval for the predicted proportions.

- **proportion_upper** - A numeric column representing the upper bound
  (97.5%) of the 95% credible interval for the predicted proportions.

- **unconstrained_mean** - A numeric column representing the mean
  unconstrained predictors (before softmax transformation).

- **unconstrained_lower** - A numeric column representing the lower
  bound (2.5%) of the 95% credible interval for the unconstrained
  predictors.

- **unconstrained_upper** - A numeric column representing the upper
  bound (97.5%) of the 95% credible interval for the unconstrained
  predictors.

- **unconstrained** - A numeric column (when
  summary_instead_of_draws=FALSE) representing individual draws of the
  unconstrained predictors.

- **proportion** - A numeric column (when
  summary_instead_of_draws=FALSE) representing individual draws of the
  predicted proportions.

- **.draw** - An integer column (when summary_instead_of_draws=FALSE)
  representing the draw index.

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
  if (instantiate::stan_cmdstan_exists() && .Platform$OS.type == "unix") {
    data("counts_obj")

    sccomp_estimate(
      counts_obj,
      ~ type, ~1, "sample", "cell_group", "count",
      cores = 1
    ) |>
    sccomp_predict()
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481736.465292 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.245e-03   1.980e-01    1.000e+00  1.000e+00      3193 -3.705e+03 -3.703e+03                   
#> Path [1] :Best Iter: [57] ELBO (-3702.565331) evaluations: (3193) 
#> Path [2] :Initial log joint density = -481860.688422 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      8.305e-03   2.726e-01    8.709e-01  8.709e-01      2724 -3.711e+03 -3.721e+03                   
#> Path [2] :Best Iter: [45] ELBO (-3710.584652) evaluations: (2724) 
#> Path [3] :Initial log joint density = -481646.542410 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.282e-03   2.254e-01    1.000e+00  1.000e+00      3396 -3.701e+03 -3.708e+03                   
#> Path [3] :Best Iter: [58] ELBO (-3701.003387) evaluations: (3396) 
#> Path [4] :Initial log joint density = -483351.891450 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.737e-03   1.862e-01    9.297e-01  9.297e-01      3277 -3.700e+03 -3.710e+03                   
#> Path [4] :Best Iter: [56] ELBO (-3700.166928) evaluations: (3277) 
#> Path [5] :Initial log joint density = -481491.369531 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.122e-02   3.609e-01    4.949e-01  1.000e+00      3124 -3.709e+03 -3.710e+03                   
#> Path [5] :Best Iter: [43] ELBO (-3708.636080) evaluations: (3124) 
#> Path [6] :Initial log joint density = -481661.919655 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.076e-02   2.913e-01    4.798e-01  1.000e+00      3123 -3.705e+03 -3.707e+03                   
#> Path [6] :Best Iter: [43] ELBO (-3705.205422) evaluations: (3123) 
#> Path [7] :Initial log joint density = -481587.336505 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.808e-03   3.538e-01    3.729e-01  1.000e+00      3127 -3.705e+03 -3.707e+03                   
#> Path [7] :Best Iter: [43] ELBO (-3705.254014) evaluations: (3127) 
#> Path [8] :Initial log joint density = -481535.379693 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.966e-03   2.348e-01    1.000e+00  1.000e+00      3187 -3.706e+03 -3.703e+03                   
#> Path [8] :Best Iter: [57] ELBO (-3702.952267) evaluations: (3187) 
#> Path [9] :Initial log joint density = -481825.565125 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.776e-02   3.184e-01    8.911e-01  8.911e-01      2983 -3.708e+03 -3.719e+03                   
#> Path [9] :Best Iter: [48] ELBO (-3708.159174) evaluations: (2983) 
#> Path [10] :Initial log joint density = -481501.000775 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      7.657e-03   2.326e-01    1.000e+00  1.000e+00      2672 -3.707e+03 -3.709e+03                   
#> Path [10] :Best Iter: [42] ELBO (-3707.317868) evaluations: (2672) 
#> Path [11] :Initial log joint density = -481524.238457 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      8.147e-03   2.613e-01    1.000e+00  1.000e+00      2731 -3.708e+03 -3.719e+03                   
#> Path [11] :Best Iter: [45] ELBO (-3708.016193) evaluations: (2731) 
#> Path [12] :Initial log joint density = -483138.185320 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.330e-03   2.200e-01    1.000e+00  1.000e+00      3155 -3.706e+03 -3.703e+03                   
#> Path [12] :Best Iter: [56] ELBO (-3702.944946) evaluations: (3155) 
#> Path [13] :Initial log joint density = -481576.111766 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      6.565e-03   1.443e-01    1.000e+00  1.000e+00      3478 -3.699e+03 -3.713e+03                   
#> Path [13] :Best Iter: [58] ELBO (-3699.258470) evaluations: (3478) 
#> Path [14] :Initial log joint density = -481764.989011 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.149e-03   2.340e-01    1.000e+00  1.000e+00      3030 -3.709e+03 -3.712e+03                   
#> Path [14] :Best Iter: [43] ELBO (-3708.899677) evaluations: (3030) 
#> Path [15] :Initial log joint density = -481387.741241 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.577e-03   2.368e-01    8.248e-01  8.248e-01      3010 -3.709e+03 -3.710e+03                   
#> Path [15] :Best Iter: [54] ELBO (-3708.517140) evaluations: (3010) 
#> Path [16] :Initial log joint density = -484897.723958 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      6.879e-03   1.901e-01    8.535e-01  8.535e-01      3527 -3.700e+03 -3.713e+03                   
#> Path [16] :Best Iter: [55] ELBO (-3699.954256) evaluations: (3527) 
#> Path [17] :Initial log joint density = -481749.069759 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      3.157e-03   3.016e-01    6.092e-01  6.092e-01      3151 -3.705e+03 -3.716e+03                   
#> Path [17] :Best Iter: [50] ELBO (-3705.420546) evaluations: (3151) 
#> Path [18] :Initial log joint density = -481923.580134 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.555e-02   2.769e-01    1.000e+00  1.000e+00      3273 -3.703e+03 -3.706e+03                   
#> Path [18] :Best Iter: [57] ELBO (-3702.926294) evaluations: (3273) 
#> Path [19] :Initial log joint density = -482255.116395 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.393e-02   2.714e-01    1.000e+00  1.000e+00      3327 -3.701e+03 -3.699e+03                   
#> Path [19] :Best Iter: [58] ELBO (-3699.435785) evaluations: (3327) 
#> Path [20] :Initial log joint density = -481546.240853 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.100e-03   2.068e-01    1.000e+00  1.000e+00      2944 -3.708e+03 -3.711e+03                   
#> Path [20] :Best Iter: [52] ELBO (-3708.200608) evaluations: (2944) 
#> Path [21] :Initial log joint density = -481984.349458 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.085e-02   3.039e-01    9.335e-01  9.335e-01      3220 -3.701e+03 -3.711e+03                   
#> Path [21] :Best Iter: [57] ELBO (-3700.894361) evaluations: (3220) 
#> Path [22] :Initial log joint density = -483060.662799 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.423e-03   2.285e-01    9.188e-01  9.188e-01      2989 -3.705e+03 -3.720e+03                   
#> Path [22] :Best Iter: [52] ELBO (-3705.414592) evaluations: (2989) 
#> Path [23] :Initial log joint density = -481516.159447 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.730e-03   2.545e-01    8.761e-01  8.761e-01      3080 -3.706e+03 -3.712e+03                   
#> Path [23] :Best Iter: [54] ELBO (-3705.728833) evaluations: (3080) 
#> Path [24] :Initial log joint density = -481385.707339 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.805e-03   1.690e-01    1.000e+00  1.000e+00      2863 -3.706e+03 -3.710e+03                   
#> Path [24] :Best Iter: [41] ELBO (-3705.904961) evaluations: (2863) 
#> Path [25] :Initial log joint density = -481626.134189 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      7.901e-03   2.339e-01    1.000e+00  1.000e+00      3478 -3.700e+03 -3.703e+03                   
#> Path [25] :Best Iter: [58] ELBO (-3700.311628) evaluations: (3478) 
#> Path [26] :Initial log joint density = -481832.934655 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      2.559e-02   3.115e-01    1.000e+00  1.000e+00      3653 -3.699e+03 -3.707e+03                   
#> Path [26] :Best Iter: [58] ELBO (-3699.288833) evaluations: (3653) 
#> Path [27] :Initial log joint density = -481618.829359 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.574e-03   2.248e-01    1.000e+00  1.000e+00      2993 -3.709e+03 -3.710e+03                   
#> Path [27] :Best Iter: [37] ELBO (-3708.527766) evaluations: (2993) 
#> Path [28] :Initial log joint density = -481547.034706 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.349e-02   4.568e-01    1.000e+00  1.000e+00      3074 -3.708e+03 -3.712e+03                   
#> Path [28] :Best Iter: [54] ELBO (-3708.430839) evaluations: (3074) 
#> Path [29] :Initial log joint density = -481485.073153 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.396e-02   3.879e-01    9.610e-01  9.610e-01      3101 -3.701e+03 -3.714e+03                   
#> Path [29] :Best Iter: [55] ELBO (-3700.851565) evaluations: (3101) 
#> Path [30] :Initial log joint density = -481709.001740 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.010e-03   2.297e-01    1.000e+00  1.000e+00      3026 -3.709e+03 -3.712e+03                   
#> Path [30] :Best Iter: [44] ELBO (-3708.672494) evaluations: (3026) 
#> Path [31] :Initial log joint density = -483175.094337 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.036e-02   3.089e-01    1.000e+00  1.000e+00      3158 -3.703e+03 -3.705e+03                   
#> Path [31] :Best Iter: [55] ELBO (-3702.585400) evaluations: (3158) 
#> Path [32] :Initial log joint density = -481389.434125 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.246e-03   1.964e-01    1.000e+00  1.000e+00      3026 -3.705e+03 -3.707e+03                   
#> Path [32] :Best Iter: [41] ELBO (-3705.329082) evaluations: (3026) 
#> Path [33] :Initial log joint density = -481637.819583 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      9.740e-03   1.906e-01    1.000e+00  1.000e+00      2886 -3.709e+03 -3.716e+03                   
#> Path [33] :Best Iter: [44] ELBO (-3709.471560) evaluations: (2886) 
#> Path [34] :Initial log joint density = -481601.725867 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.239e-02   2.980e-01    1.000e+00  1.000e+00      3134 -3.702e+03 -3.711e+03                   
#> Path [34] :Best Iter: [55] ELBO (-3701.944585) evaluations: (3134) 
#> Path [35] :Initial log joint density = -482984.485174 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.229e-03   2.188e-01    9.279e-01  9.279e-01      3268 -3.702e+03 -3.709e+03                   
#> Path [35] :Best Iter: [55] ELBO (-3701.606776) evaluations: (3268) 
#> Path [36] :Initial log joint density = -482364.704955 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.632e-02   3.706e-01    1.000e+00  1.000e+00      3437 -3.698e+03 -3.703e+03                   
#> Path [36] :Best Iter: [56] ELBO (-3697.836833) evaluations: (3437) 
#> Path [37] :Initial log joint density = -481709.538224 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.180e-02   3.242e-01    8.385e-01  8.385e-01      3163 -3.709e+03 -3.712e+03                   
#> Path [37] :Best Iter: [41] ELBO (-3708.864064) evaluations: (3163) 
#> Path [38] :Initial log joint density = -483548.081022 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.682e-02   2.059e-01    1.000e+00  1.000e+00      3510 -3.702e+03 -3.702e+03                   
#> Path [38] :Best Iter: [60] ELBO (-3701.718762) evaluations: (3510) 
#> Path [39] :Initial log joint density = -481657.934013 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.809e-03   2.065e-01    1.000e+00  1.000e+00      3227 -3.706e+03 -3.700e+03                   
#> Path [39] :Best Iter: [56] ELBO (-3700.297599) evaluations: (3227) 
#> Path [40] :Initial log joint density = -481567.525112 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      9.791e-03   1.822e-01    1.000e+00  1.000e+00      2866 -3.706e+03 -3.715e+03                   
#> Path [40] :Best Iter: [50] ELBO (-3706.201049) evaluations: (2866) 
#> Path [41] :Initial log joint density = -481676.575940 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      5.983e-03   2.137e-01    1.000e+00  1.000e+00      2827 -3.707e+03 -3.716e+03                   
#> Path [41] :Best Iter: [49] ELBO (-3706.601060) evaluations: (2827) 
#> Path [42] :Initial log joint density = -481967.920483 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.102e-02   2.692e-01    1.000e+00  1.000e+00      3518 -3.704e+03 -3.704e+03                   
#> Path [42] :Best Iter: [57] ELBO (-3703.559800) evaluations: (3518) 
#> Path [43] :Initial log joint density = -481440.955609 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      9.330e-03   2.678e-01    9.031e-01  9.031e-01      3215 -3.700e+03 -3.714e+03                   
#> Path [43] :Best Iter: [55] ELBO (-3700.102331) evaluations: (3215) 
#> Path [44] :Initial log joint density = -481860.382020 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      2.649e-02   2.111e-01    1.000e+00  1.000e+00      3334 -3.700e+03 -3.700e+03                   
#> Path [44] :Best Iter: [58] ELBO (-3699.608257) evaluations: (3334) 
#> Path [45] :Initial log joint density = -481561.106265 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.059e-02   3.338e-01    1.000e+00  1.000e+00      3155 -3.699e+03 -3.706e+03                   
#> Path [45] :Best Iter: [55] ELBO (-3698.590171) evaluations: (3155) 
#> Path [46] :Initial log joint density = -481680.333856 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.257e-03   1.803e-01    1.000e+00  1.000e+00      3053 -3.705e+03 -3.701e+03                   
#> Path [46] :Best Iter: [56] ELBO (-3701.137391) evaluations: (3053) 
#> Path [47] :Initial log joint density = -481543.598434 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.154e-02   3.887e-01    1.000e+00  1.000e+00      2944 -3.708e+03 -3.722e+03                   
#> Path [47] :Best Iter: [48] ELBO (-3708.422815) evaluations: (2944) 
#> Path [48] :Initial log joint density = -481753.026806 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.090e-03   2.074e-01    1.000e+00  1.000e+00      3161 -3.700e+03 -3.705e+03                   
#> Path [48] :Best Iter: [55] ELBO (-3700.169205) evaluations: (3161) 
#> Path [49] :Initial log joint density = -481559.321930 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.469e-03   2.472e-01    1.000e+00  1.000e+00      2921 -3.708e+03 -3.709e+03                   
#> Path [49] :Best Iter: [50] ELBO (-3708.362261) evaluations: (2921) 
#> Path [50] :Initial log joint density = -481570.241650 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.306e-03   3.328e-01    1.000e+00  1.000e+00      2913 -3.713e+03 -3.726e+03                   
#> Path [50] :Best Iter: [45] ELBO (-3712.981570) evaluations: (2913) 
#> Finished in  13.3 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 720 × 9
#>    sample type   cell_group proportion_mean proportion_lower proportion_upper
#>    <fct>  <fct>  <chr>                <dbl>            <dbl>            <dbl>
#>  1 10x_6K benign B1                 0.0593           0.0458           0.0752 
#>  2 10x_6K benign B2                 0.0378           0.0289           0.0472 
#>  3 10x_6K benign B3                 0.0125           0.00951          0.0159 
#>  4 10x_6K benign BM                 0.00676          0.00509          0.00881
#>  5 10x_6K benign CD4 1              0.0255           0.0214           0.0302 
#>  6 10x_6K benign CD4 2              0.0509           0.0420           0.0620 
#>  7 10x_6K benign CD4 3              0.0818           0.0614           0.103  
#>  8 10x_6K benign CD4 4              0.00169          0.00115          0.00243
#>  9 10x_6K benign CD4 5              0.0305           0.0231           0.0380 
#> 10 10x_6K benign CD8 1              0.112            0.0963           0.128  
#> # ℹ 710 more rows
#> # ℹ 3 more variables: unconstrained_mean <dbl>, unconstrained_lower <dbl>,
#> #   unconstrained_upper <dbl>
# }

```
