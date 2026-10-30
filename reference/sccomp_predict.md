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
  samples, there should be 10 rows. Continuous covariates, including
  those inside `s()` / `t2()`, should be in the original units of the
  fitted data. They are z-scored once from those fitted samples before
  the design is built.

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
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), typecancer
#>   variability Xa: 1 parameter: (Intercept)
#>   random effects: none
#> Loading model from cache...
#> Path [1] :Initial log joint density = -482914.268162 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      8.051e-03   2.683e-01    1.000e+00  1.000e+00      3327 -3.704e+03 -3.714e+03                   
#> Path [1] :Best Iter: [55] ELBO (-3704.320911) evaluations: (3327) 
#> Path [2] :Initial log joint density = -482036.045745 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      3.224e-02   2.326e-01    1.000e+00  1.000e+00      4594 -3.700e+03 -3.702e+03                   
#> Path [2] :Best Iter: [67] ELBO (-3700.497265) evaluations: (4594) 
#> Path [3] :Initial log joint density = -481485.892437 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.668e-02   2.574e-01    1.000e+00  1.000e+00      4165 -3.701e+03 -3.708e+03                   
#> Path [3] :Best Iter: [64] ELBO (-3701.282933) evaluations: (4165) 
#> Path [4] :Initial log joint density = -481640.932278 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      5.678e-03   1.328e-01    1.000e+00  1.000e+00      4607 -3.699e+03 -3.710e+03                   
#> Path [4] :Best Iter: [67] ELBO (-3699.378226) evaluations: (4607) 
#> Path [5] :Initial log joint density = -481411.773949 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      8.920e-03   1.372e-01    8.651e-01  8.651e-01      3925 -3.702e+03 -3.714e+03                   
#> Path [5] :Best Iter: [63] ELBO (-3701.944425) evaluations: (3925) 
#> Path [6] :Initial log joint density = -481609.490912 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      3.395e-03   2.794e-01    6.225e-01  6.225e-01      3987 -3.703e+03 -3.720e+03                   
#> Path [6] :Best Iter: [61] ELBO (-3702.746382) evaluations: (3987) 
#> Path [7] :Initial log joint density = -481973.953077 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.787e+05      6.123e-03   1.892e-01    9.743e-01  9.743e-01      5835 -3.697e+03 -3.705e+03                   
#> Path [7] :Best Iter: [79] ELBO (-3697.179480) evaluations: (5835) 
#> Path [8] :Initial log joint density = -482137.651357 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.932e-02   1.663e-01    1.000e+00  1.000e+00      5291 -3.702e+03 -3.703e+03                   
#> Path [8] :Best Iter: [74] ELBO (-3701.506360) evaluations: (5291) 
#> Path [9] :Initial log joint density = -481746.173174 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.668e-02   2.216e-01    1.000e+00  1.000e+00      4583 -3.700e+03 -3.709e+03                   
#> Path [9] :Best Iter: [69] ELBO (-3700.311426) evaluations: (4583) 
#> Path [10] :Initial log joint density = -482376.089014 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      7.017e-03   2.160e-01    1.000e+00  1.000e+00      4438 -3.700e+03 -3.708e+03                   
#> Path [10] :Best Iter: [67] ELBO (-3700.241492) evaluations: (4438) 
#> Path [11] :Initial log joint density = -481574.565221 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      7.150e-03   2.108e-01    1.000e+00  1.000e+00      3651 -3.702e+03 -3.714e+03                   
#> Path [11] :Best Iter: [58] ELBO (-3702.458344) evaluations: (3651) 
#> Path [12] :Initial log joint density = -481451.320036 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      8.972e-03   2.263e-01    1.000e+00  1.000e+00      3962 -3.700e+03 -3.708e+03                   
#> Path [12] :Best Iter: [63] ELBO (-3700.305119) evaluations: (3962) 
#> Path [13] :Initial log joint density = -481580.551034 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      8.162e-03   1.253e-01    4.828e-01  9.616e-01      4274 -3.703e+03 -3.717e+03                   
#> Path [13] :Best Iter: [65] ELBO (-3703.097805) evaluations: (4274) 
#> Path [14] :Initial log joint density = -481597.673989 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      4.097e-03   2.480e-01    6.435e-01  6.435e-01      4132 -3.700e+03 -3.710e+03                   
#> Path [14] :Best Iter: [63] ELBO (-3699.541445) evaluations: (4132) 
#> Path [15] :Initial log joint density = -481460.010815 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      5.093e-03   2.328e-01    4.028e-01  4.028e-01      4302 -3.705e+03 -3.705e+03                   
#> Path [15] :Best Iter: [68] ELBO (-3704.642207) evaluations: (4302) 
#> Path [16] :Initial log joint density = -481503.941402 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      6.775e-03   2.243e-01    4.406e-01  1.000e+00      3780 -3.704e+03 -3.715e+03                   
#> Path [16] :Best Iter: [61] ELBO (-3704.216355) evaluations: (3780) 
#> Path [17] :Initial log joint density = -481865.445071 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      9.982e-03   2.066e-01    9.571e-01  9.571e-01      3572 -3.703e+03 -3.712e+03                   
#> Path [17] :Best Iter: [59] ELBO (-3702.946444) evaluations: (3572) 
#> Path [18] :Initial log joint density = -482089.329452 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      1.197e-02   1.369e-01    8.843e-01  8.843e-01      5579 -3.698e+03 -3.712e+03                   
#> Path [18] :Best Iter: [73] ELBO (-3698.435248) evaluations: (5579) 
#> Path [19] :Initial log joint density = -481619.633161 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      6.755e-03   2.262e-01    1.000e+00  1.000e+00      3109 -3.709e+03 -3.712e+03                   
#> Path [19] :Best Iter: [48] ELBO (-3708.796664) evaluations: (3109) 
#> Path [20] :Initial log joint density = -481535.898138 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      8.281e-03   1.414e-01    1.000e+00  1.000e+00      3708 -3.703e+03 -3.707e+03                   
#> Path [20] :Best Iter: [59] ELBO (-3703.344334) evaluations: (3708) 
#> Path [21] :Initial log joint density = -481597.115542 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      8.120e-03   2.216e-01    1.000e+00  1.000e+00      3193 -3.710e+03 -3.717e+03                   
#> Path [21] :Best Iter: [48] ELBO (-3710.493497) evaluations: (3193) 
#> Path [22] :Initial log joint density = -481786.630182 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      1.534e-02   1.815e-01    1.000e+00  1.000e+00      5293 -3.697e+03 -3.701e+03                   
#> Path [22] :Best Iter: [76] ELBO (-3697.385394) evaluations: (5293) 
#> Path [23] :Initial log joint density = -481733.446996 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.341e-02   3.118e-01    4.899e-01  1.000e+00      4152 -3.703e+03 -3.712e+03                   
#> Path [23] :Best Iter: [60] ELBO (-3702.516488) evaluations: (4152) 
#> Path [24] :Initial log joint density = -482348.068153 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.645e-02   2.775e-01    1.000e+00  1.000e+00      4096 -3.701e+03 -3.706e+03                   
#> Path [24] :Best Iter: [62] ELBO (-3700.836153) evaluations: (4096) 
#> Path [25] :Initial log joint density = -481519.838599 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      4.380e-03   2.275e-01    8.799e-01  8.799e-01      4502 -3.702e+03 -3.706e+03                   
#> Path [25] :Best Iter: [68] ELBO (-3701.788949) evaluations: (4502) 
#> Path [26] :Initial log joint density = -481564.978581 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      2.075e-02   3.076e-01    1.000e+00  1.000e+00      4149 -3.701e+03 -3.703e+03                   
#> Path [26] :Best Iter: [57] ELBO (-3701.043943) evaluations: (4149) 
#> Path [27] :Initial log joint density = -481669.518808 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      7.542e-03   1.908e-01    1.000e+00  1.000e+00      4999 -3.699e+03 -3.703e+03                   
#> Path [27] :Best Iter: [71] ELBO (-3699.054872) evaluations: (4999) 
#> Path [28] :Initial log joint density = -484078.325531 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      7.507e-03   1.597e-01    1.000e+00  1.000e+00      4436 -3.700e+03 -3.705e+03                   
#> Path [28] :Best Iter: [68] ELBO (-3699.614413) evaluations: (4436) 
#> Path [29] :Initial log joint density = -485858.904301 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.343e-02   1.854e-01    1.000e+00  1.000e+00      4754 -3.701e+03 -3.701e+03                   
#> Path [29] :Best Iter: [71] ELBO (-3700.541382) evaluations: (4754) 
#> Path [30] :Initial log joint density = -481523.163870 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      6.664e-03   1.810e-01    7.366e-01  7.366e-01      4601 -3.700e+03 -3.716e+03                   
#> Path [30] :Best Iter: [66] ELBO (-3699.632828) evaluations: (4601) 
#> Path [31] :Initial log joint density = -481606.770362 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      8.062e-03   2.043e-01    1.000e+00  1.000e+00      4565 -3.703e+03 -3.707e+03                   
#> Path [31] :Best Iter: [67] ELBO (-3703.166413) evaluations: (4565) 
#> Path [32] :Initial log joint density = -481647.166033 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      6.777e-03   1.577e-01    1.000e+00  1.000e+00      3193 -3.712e+03 -3.720e+03                   
#> Path [32] :Best Iter: [53] ELBO (-3712.016842) evaluations: (3193) 
#> Path [33] :Initial log joint density = -485058.619486 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.195e-02   1.910e-01    1.000e+00  1.000e+00      4411 -3.700e+03 -3.704e+03                   
#> Path [33] :Best Iter: [66] ELBO (-3699.776634) evaluations: (4411) 
#> Path [34] :Initial log joint density = -481819.353020 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.787e+05      8.647e-03   1.482e-01    1.000e+00  1.000e+00      5609 -3.700e+03 -3.711e+03                   
#> Path [34] :Best Iter: [77] ELBO (-3700.295638) evaluations: (5609) 
#> Path [35] :Initial log joint density = -481427.318114 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      4.694e-03   1.853e-01    1.000e+00  1.000e+00      3334 -3.709e+03 -3.715e+03                   
#> Path [35] :Best Iter: [54] ELBO (-3708.840431) evaluations: (3334) 
#> Path [36] :Initial log joint density = -481556.565910 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.608e-02   1.679e-01    1.000e+00  1.000e+00      3994 -3.698e+03 -3.702e+03                   
#> Path [36] :Best Iter: [60] ELBO (-3697.950524) evaluations: (3994) 
#> Path [37] :Initial log joint density = -481833.202035 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      8.908e-03   1.812e-01    1.000e+00  1.000e+00      4884 -3.699e+03 -3.702e+03                   
#> Path [37] :Best Iter: [70] ELBO (-3698.898605) evaluations: (4884) 
#> Path [38] :Initial log joint density = -482001.821527 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      2.026e-02   3.265e-01    1.000e+00  1.000e+00      5114 -3.698e+03 -3.705e+03                   
#> Path [38] :Best Iter: [76] ELBO (-3698.200776) evaluations: (5114) 
#> Path [39] :Initial log joint density = -481789.721256 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      6.288e-03   1.257e-01    1.000e+00  1.000e+00      3927 -3.703e+03 -3.712e+03                   
#> Path [39] :Best Iter: [60] ELBO (-3702.588510) evaluations: (3927) 
#> Path [40] :Initial log joint density = -483002.567159 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      7.654e-03   2.242e-01    1.000e+00  1.000e+00      3604 -3.707e+03 -3.711e+03                   
#> Path [40] :Best Iter: [56] ELBO (-3706.563418) evaluations: (3604) 
#> Path [41] :Initial log joint density = -482652.991102 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      6.986e-03   1.474e-01    9.273e-01  9.273e-01      4767 -3.701e+03 -3.712e+03                   
#> Path [41] :Best Iter: [65] ELBO (-3701.087981) evaluations: (4767) 
#> Path [42] :Initial log joint density = -481584.817571 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.753e-02   2.192e-01    1.000e+00  1.000e+00      4784 -3.700e+03 -3.704e+03                   
#> Path [42] :Best Iter: [72] ELBO (-3700.403606) evaluations: (4784) 
#> Path [43] :Initial log joint density = -481753.084642 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      3.324e-02   2.568e-01    1.000e+00  1.000e+00      4216 -3.703e+03 -3.705e+03                   
#> Path [43] :Best Iter: [65] ELBO (-3703.064215) evaluations: (4216) 
#> Path [44] :Initial log joint density = -482981.464867 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      7.036e-03   1.846e-01    9.679e-01  9.679e-01      4782 -3.701e+03 -3.708e+03                   
#> Path [44] :Best Iter: [70] ELBO (-3701.358948) evaluations: (4782) 
#> Path [45] :Initial log joint density = -481724.996540 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      9.004e-03   1.961e-01    9.470e-01  9.470e-01      3418 -3.701e+03 -3.716e+03                   
#> Path [45] :Best Iter: [57] ELBO (-3701.497518) evaluations: (3418) 
#> Path [46] :Initial log joint density = -481763.228681 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      1.265e-02   1.644e-01    1.000e+00  1.000e+00      4801 -3.701e+03 -3.700e+03                   
#> Path [46] :Best Iter: [73] ELBO (-3700.032682) evaluations: (4801) 
#> Path [47] :Initial log joint density = -481524.329010 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.536e-02   2.372e-01    1.000e+00  1.000e+00      4327 -3.699e+03 -3.702e+03                   
#> Path [47] :Best Iter: [65] ELBO (-3699.404106) evaluations: (4327) 
#> Path [48] :Initial log joint density = -481971.501253 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83      -4.787e+05      7.535e-03   1.754e-01    7.497e-01  7.497e-01      5869 -3.700e+03 -3.713e+03                   
#> Path [48] :Best Iter: [81] ELBO (-3700.142768) evaluations: (5869) 
#> Path [49] :Initial log joint density = -481792.282091 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.062e-02   1.957e-01    7.749e-01  7.749e-01      4170 -3.699e+03 -3.714e+03                   
#> Path [49] :Best Iter: [57] ELBO (-3699.331011) evaluations: (4170) 
#> Path [50] :Initial log joint density = -481846.776134 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.453e-02   2.714e-01    9.621e-01  9.621e-01      4155 -3.700e+03 -3.708e+03                   
#> Path [50] :Best Iter: [64] ELBO (-3699.740226) evaluations: (4155) 
#> Finished in  14.6 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.518 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 720 × 9
#>    sample type   cell_group proportion_mean proportion_lower proportion_upper
#>    <fct>  <fct>  <chr>                <dbl>            <dbl>            <dbl>
#>  1 10x_6K benign B1                 0.0587           0.0462           0.0734 
#>  2 10x_6K benign B2                 0.0380           0.0287           0.0489 
#>  3 10x_6K benign B3                 0.0127           0.00954          0.0167 
#>  4 10x_6K benign BM                 0.00661          0.00502          0.00840
#>  5 10x_6K benign CD4 1              0.0253           0.0214           0.0295 
#>  6 10x_6K benign CD4 2              0.0512           0.0420           0.0620 
#>  7 10x_6K benign CD4 3              0.0816           0.0631           0.105  
#>  8 10x_6K benign CD4 4              0.00153          0.00104          0.00213
#>  9 10x_6K benign CD4 5              0.0310           0.0242           0.0398 
#> 10 10x_6K benign CD8 1              0.112            0.0939           0.132  
#> # ℹ 710 more rows
#> # ℹ 3 more variables: unconstrained_mean <dbl>, unconstrained_lower <dbl>,
#> #   unconstrained_upper <dbl>
# }

```
