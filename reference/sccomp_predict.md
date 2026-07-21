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
#> Path [1] :Initial log joint density = -482932.383905 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      9.747e-03   2.924e-01    1.000e+00  1.000e+00      3417 -3.689e+03 -3.694e+03                   
#> Path [1] :Best Iter: [55] ELBO (-3688.734569) evaluations: (3417) 
#> Path [2] :Initial log joint density = -482054.265586 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.361e-02   2.173e-01    1.000e+00  1.000e+00      4179 -3.686e+03 -3.689e+03                   
#> Path [2] :Best Iter: [65] ELBO (-3685.571776) evaluations: (4179) 
#> Path [3] :Initial log joint density = -481503.040027 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      8.378e-03   2.425e-01    1.000e+00  1.000e+00      4281 -3.685e+03 -3.687e+03                   
#> Path [3] :Best Iter: [64] ELBO (-3685.326098) evaluations: (4281) 
#> Path [4] :Initial log joint density = -481659.166283 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.135e-02   1.592e-01    8.119e-01  8.119e-01      4464 -3.685e+03 -3.697e+03                   
#> Path [4] :Best Iter: [67] ELBO (-3684.828249) evaluations: (4464) 
#> Path [5] :Initial log joint density = -481425.302402 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      1.459e-02   3.260e-01    1.000e+00  1.000e+00      3419 -3.691e+03 -3.701e+03                   
#> Path [5] :Best Iter: [58] ELBO (-3691.481763) evaluations: (3419) 
#> Path [6] :Initial log joint density = -481628.923079 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      7.793e-03   1.632e-01    1.000e+00  1.000e+00      3271 -3.690e+03 -3.689e+03                   
#> Path [6] :Best Iter: [56] ELBO (-3688.703673) evaluations: (3271) 
#> Path [7] :Initial log joint density = -481989.483880 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      3.041e-03   1.879e-01    7.104e-01  7.104e-01      4997 -3.685e+03 -3.698e+03                   
#> Path [7] :Best Iter: [69] ELBO (-3685.095710) evaluations: (4997) 
#> Path [8] :Initial log joint density = -482153.232193 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      5.868e-03   2.368e-01    7.788e-01  7.788e-01      4434 -3.686e+03 -3.696e+03                   
#> Path [8] :Best Iter: [69] ELBO (-3685.578236) evaluations: (4434) 
#> Path [9] :Initial log joint density = -481761.807464 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      9.786e-03   3.036e-01    1.000e+00  1.000e+00      4114 -3.688e+03 -3.697e+03                   
#> Path [9] :Best Iter: [64] ELBO (-3688.178688) evaluations: (4114) 
#> Path [10] :Initial log joint density = -482392.816714 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.874e-02   2.744e-01    9.947e-01  9.947e-01      4254 -3.684e+03 -3.692e+03                   
#> Path [10] :Best Iter: [65] ELBO (-3683.913586) evaluations: (4254) 
#> Path [11] :Initial log joint density = -481587.051888 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      5.122e-03   1.624e-01    1.000e+00  1.000e+00      3465 -3.690e+03 -3.695e+03                   
#> Path [11] :Best Iter: [55] ELBO (-3690.266291) evaluations: (3465) 
#> Path [12] :Initial log joint density = -481465.327929 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      5.800e-03   1.711e-01    1.000e+00  1.000e+00      4151 -3.684e+03 -3.691e+03                   
#> Path [12] :Best Iter: [63] ELBO (-3684.351069) evaluations: (4151) 
#> Path [13] :Initial log joint density = -481600.267960 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      4.922e-03   2.146e-01    8.522e-01  8.522e-01      3186 -3.697e+03 -3.699e+03                   
#> Path [13] :Best Iter: [45] ELBO (-3696.976711) evaluations: (3186) 
#> Path [14] :Initial log joint density = -481618.302860 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.341e-02   1.895e-01    1.000e+00  1.000e+00      4088 -3.686e+03 -3.690e+03                   
#> Path [14] :Best Iter: [64] ELBO (-3686.208341) evaluations: (4088) 
#> Path [15] :Initial log joint density = -481476.229753 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.002e-02   2.152e-01    1.000e+00  1.000e+00      3872 -3.687e+03 -3.691e+03                   
#> Path [15] :Best Iter: [55] ELBO (-3687.341433) evaluations: (3872) 
#> Path [16] :Initial log joint density = -481518.858871 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      8.033e-03   2.641e-01    1.000e+00  1.000e+00      3566 -3.688e+03 -3.694e+03                   
#> Path [16] :Best Iter: [59] ELBO (-3688.141450) evaluations: (3566) 
#> Path [17] :Initial log joint density = -481883.939822 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      1.228e-02   3.162e-01    9.840e-01  9.840e-01      3220 -3.692e+03 -3.701e+03                   
#> Path [17] :Best Iter: [55] ELBO (-3692.002099) evaluations: (3220) 
#> Path [18] :Initial log joint density = -482103.719604 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      6.510e-03   2.296e-01    7.939e-01  7.939e-01      4300 -3.687e+03 -3.699e+03                   
#> Path [18] :Best Iter: [66] ELBO (-3687.008337) evaluations: (4300) 
#> Path [19] :Initial log joint density = -481638.416075 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      9.151e-03   2.212e-01    1.000e+00  1.000e+00      3238 -3.693e+03 -3.688e+03                   
#> Path [19] :Best Iter: [57] ELBO (-3687.655541) evaluations: (3238) 
#> Path [20] :Initial log joint density = -481550.779889 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      1.251e-02   2.091e-01    1.000e+00  1.000e+00      3141 -3.692e+03 -3.691e+03                   
#> Path [20] :Best Iter: [55] ELBO (-3691.412394) evaluations: (3141) 
#> Path [21] :Initial log joint density = -481612.959251 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      6.369e-03   1.727e-01    1.000e+00  1.000e+00      3026 -3.696e+03 -3.698e+03                   
#> Path [21] :Best Iter: [47] ELBO (-3695.838859) evaluations: (3026) 
#> Path [22] :Initial log joint density = -481806.730128 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      2.369e-02   1.021e-01    1.000e+00  1.000e+00      5001 -3.683e+03 -3.686e+03                   
#> Path [22] :Best Iter: [73] ELBO (-3682.538157) evaluations: (5001) 
#> Path [23] :Initial log joint density = -481750.292965 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      7.803e-03   2.546e-01    7.415e-01  7.415e-01      3972 -3.684e+03 -3.697e+03                   
#> Path [23] :Best Iter: [62] ELBO (-3683.728571) evaluations: (3972) 
#> Path [24] :Initial log joint density = -482359.995264 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.935e-03   2.569e-01    5.615e-01  5.615e-01      4571 -3.683e+03 -3.696e+03                   
#> Path [24] :Best Iter: [69] ELBO (-3683.288229) evaluations: (4571) 
#> Path [25] :Initial log joint density = -481534.624613 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      5.806e-03   1.905e-01    1.000e+00  1.000e+00      3655 -3.688e+03 -3.701e+03                   
#> Path [25] :Best Iter: [55] ELBO (-3687.660407) evaluations: (3655) 
#> Path [26] :Initial log joint density = -481581.122858 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.562e-02   2.236e-01    1.000e+00  1.000e+00      4321 -3.689e+03 -3.694e+03                   
#> Path [26] :Best Iter: [65] ELBO (-3688.985127) evaluations: (4321) 
#> Path [27] :Initial log joint density = -481688.626338 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.528e-02   3.102e-01    1.000e+00  1.000e+00      4409 -3.689e+03 -3.692e+03                   
#> Path [27] :Best Iter: [63] ELBO (-3689.376300) evaluations: (4409) 
#> Path [28] :Initial log joint density = -484096.201495 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      2.872e-02   2.306e-01    1.000e+00  1.000e+00      4345 -3.684e+03 -3.685e+03                   
#> Path [28] :Best Iter: [65] ELBO (-3683.760704) evaluations: (4345) 
#> Path [29] :Initial log joint density = -485877.907487 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.532e-02   1.727e-01    1.000e+00  1.000e+00      4824 -3.684e+03 -3.687e+03                   
#> Path [29] :Best Iter: [65] ELBO (-3683.709301) evaluations: (4824) 
#> Path [30] :Initial log joint density = -481538.659713 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      6.377e-03   1.505e-01    1.000e+00  1.000e+00      4209 -3.690e+03 -3.687e+03                   
#> Path [30] :Best Iter: [66] ELBO (-3687.093941) evaluations: (4209) 
#> Path [31] :Initial log joint density = -481624.011986 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      6.628e-03   1.815e-01    8.322e-01  8.322e-01      3979 -3.686e+03 -3.699e+03                   
#> Path [31] :Best Iter: [59] ELBO (-3686.408727) evaluations: (3979) 
#> Path [32] :Initial log joint density = -481661.325022 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      1.226e-02   3.099e-01    1.000e+00  1.000e+00      3257 -3.689e+03 -3.697e+03                   
#> Path [32] :Best Iter: [56] ELBO (-3689.392225) evaluations: (3257) 
#> Path [33] :Initial log joint density = -485074.841298 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.208e-02   2.434e-01    1.000e+00  1.000e+00      4607 -3.686e+03 -3.686e+03                   
#> Path [33] :Best Iter: [68] ELBO (-3686.169805) evaluations: (4607) 
#> Path [34] :Initial log joint density = -481829.565353 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      7.924e-03   2.067e-01    1.000e+00  1.000e+00      4811 -3.685e+03 -3.691e+03                   
#> Path [34] :Best Iter: [67] ELBO (-3684.631209) evaluations: (4811) 
#> Path [35] :Initial log joint density = -481440.513069 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      1.193e-02   2.561e-01    1.000e+00  1.000e+00      3246 -3.690e+03 -3.692e+03                   
#> Path [35] :Best Iter: [55] ELBO (-3689.629779) evaluations: (3246) 
#> Path [36] :Initial log joint density = -481570.867358 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.326e-02   2.273e-01    8.892e-01  8.892e-01      3718 -3.685e+03 -3.694e+03                   
#> Path [36] :Best Iter: [60] ELBO (-3685.143988) evaluations: (3718) 
#> Path [37] :Initial log joint density = -481850.003401 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      9.136e-03   2.351e-01    4.322e-01  1.000e+00      4246 -3.683e+03 -3.693e+03                   
#> Path [37] :Best Iter: [67] ELBO (-3682.938444) evaluations: (4246) 
#> Path [38] :Initial log joint density = -482021.593965 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      3.630e-03   2.139e-01    7.163e-01  7.163e-01      4353 -3.685e+03 -3.699e+03                   
#> Path [38] :Best Iter: [67] ELBO (-3684.539506) evaluations: (4353) 
#> Path [39] :Initial log joint density = -481809.642552 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.787e+05      1.197e-02   3.336e-01    5.193e-01  1.000e+00      2878 -3.695e+03 -3.702e+03                   
#> Path [39] :Best Iter: [50] ELBO (-3695.427142) evaluations: (2878) 
#> Path [40] :Initial log joint density = -483016.240236 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      7.291e-03   1.799e-01    1.000e+00  1.000e+00      3523 -3.690e+03 -3.697e+03                   
#> Path [40] :Best Iter: [56] ELBO (-3689.631032) evaluations: (3523) 
#> Path [41] :Initial log joint density = -482663.361163 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.709e-02   1.188e-01    1.000e+00  1.000e+00      4539 -3.684e+03 -3.691e+03                   
#> Path [41] :Best Iter: [65] ELBO (-3683.875471) evaluations: (4539) 
#> Path [42] :Initial log joint density = -481602.926230 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.105e-02   3.305e-01    1.000e+00  1.000e+00      4056 -3.689e+03 -3.694e+03                   
#> Path [42] :Best Iter: [63] ELBO (-3688.640830) evaluations: (4056) 
#> Path [43] :Initial log joint density = -481769.746016 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.309e-02   1.865e-01    1.000e+00  1.000e+00      4466 -3.684e+03 -3.686e+03                   
#> Path [43] :Best Iter: [67] ELBO (-3684.028060) evaluations: (4466) 
#> Path [44] :Initial log joint density = -482997.383986 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      7.990e-03   2.702e-01    6.126e-01  6.126e-01      4629 -3.685e+03 -3.694e+03                   
#> Path [44] :Best Iter: [59] ELBO (-3684.610178) evaluations: (4629) 
#> Path [45] :Initial log joint density = -481740.610222 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      5.997e-03   1.740e-01    1.000e+00  1.000e+00      3331 -3.689e+03 -3.699e+03                   
#> Path [45] :Best Iter: [55] ELBO (-3688.942477) evaluations: (3331) 
#> Path [46] :Initial log joint density = -481780.003190 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      8.136e-03   2.390e-01    8.971e-01  8.971e-01      4016 -3.687e+03 -3.698e+03                   
#> Path [46] :Best Iter: [64] ELBO (-3687.089959) evaluations: (4016) 
#> Path [47] :Initial log joint density = -481539.541052 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.112e-02   2.495e-01    1.000e+00  1.000e+00      4511 -3.683e+03 -3.689e+03                   
#> Path [47] :Best Iter: [66] ELBO (-3683.058982) evaluations: (4511) 
#> Path [48] :Initial log joint density = -481987.938149 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      8.755e-03   2.245e-01    9.061e-01  9.061e-01      4468 -3.686e+03 -3.695e+03                   
#> Path [48] :Best Iter: [68] ELBO (-3686.244230) evaluations: (4468) 
#> Path [49] :Initial log joint density = -481810.491537 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      1.193e-02   1.514e-01    1.000e+00  1.000e+00      4894 -3.683e+03 -3.683e+03                   
#> Path [49] :Best Iter: [73] ELBO (-3682.868178) evaluations: (4894) 
#> Path [50] :Initial log joint density = -481867.653485 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      8.152e-03   2.171e-01    1.000e+00  1.000e+00      3966 -3.684e+03 -3.691e+03                   
#> Path [50] :Best Iter: [63] ELBO (-3683.778625) evaluations: (3966) 
#> Finished in  15.4 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.511 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 720 × 9
#>    sample type   cell_group proportion_mean proportion_lower proportion_upper
#>    <fct>  <fct>  <chr>                <dbl>            <dbl>            <dbl>
#>  1 10x_6K benign B1                 0.0583           0.0463           0.0734 
#>  2 10x_6K benign B2                 0.0381           0.0291           0.0487 
#>  3 10x_6K benign B3                 0.0127           0.00981          0.0160 
#>  4 10x_6K benign BM                 0.00675          0.00519          0.00867
#>  5 10x_6K benign CD4 1              0.0255           0.0215           0.0299 
#>  6 10x_6K benign CD4 2              0.0507           0.0412           0.0623 
#>  7 10x_6K benign CD4 3              0.0817           0.0628           0.104  
#>  8 10x_6K benign CD4 4              0.00165          0.00107          0.00239
#>  9 10x_6K benign CD4 5              0.0305           0.0236           0.0390 
#> 10 10x_6K benign CD8 1              0.111            0.0943           0.130  
#> # ℹ 710 more rows
#> # ℹ 3 more variables: unconstrained_mean <dbl>, unconstrained_lower <dbl>,
#> #   unconstrained_upper <dbl>
# }

```
