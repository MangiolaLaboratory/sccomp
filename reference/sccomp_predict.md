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
#> Path [1] :Initial log joint density = -481734.557481 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.782e-02   2.137e-01    1.000e+00  1.000e+00      3364 -3.683e+03 -3.683e+03                   
#> Path [1] :Best Iter: [59] ELBO (-3683.123617) evaluations: (3364) 
#> Path [2] :Initial log joint density = -481859.638338 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      7.881e-03   2.590e-01    8.462e-01  8.462e-01      2869 -3.691e+03 -3.706e+03                   
#> Path [2] :Best Iter: [51] ELBO (-3691.423923) evaluations: (2869) 
#> Path [3] :Initial log joint density = -481645.361562 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      7.029e-03   2.060e-01    8.714e-01  8.714e-01      3396 -3.686e+03 -3.700e+03                   
#> Path [3] :Best Iter: [58] ELBO (-3685.870491) evaluations: (3396) 
#> Path [4] :Initial log joint density = -483351.423837 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.088e-03   3.069e-01    1.000e+00  1.000e+00      3288 -3.686e+03 -3.698e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3686.302448) evaluations: (3288) 
#> Path [5] :Initial log joint density = -481488.987259 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.680e-03   3.606e-01    3.731e-01  1.000e+00      3158 -3.687e+03 -3.696e+03                   
#> Path [5] :Best Iter: [55] ELBO (-3687.070593) evaluations: (3158) 
#> Path [6] :Initial log joint density = -481659.676563 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.835e-03   2.918e-01    1.000e+00  1.000e+00      3207 -3.682e+03 -3.686e+03                   
#> Path [6] :Best Iter: [55] ELBO (-3682.089991) evaluations: (3207) 
#> Path [7] :Initial log joint density = -481583.848041 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      5.671e-03   2.458e-01    9.133e-01  9.133e-01      2833 -3.689e+03 -3.706e+03                   
#> Path [7] :Best Iter: [46] ELBO (-3688.856477) evaluations: (2833) 
#> Path [8] :Initial log joint density = -481535.021836 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.437e-02   3.491e-01    1.000e+00  1.000e+00      3272 -3.682e+03 -3.688e+03                   
#> Path [8] :Best Iter: [57] ELBO (-3681.522784) evaluations: (3272) 
#> Path [9] :Initial log joint density = -481825.435626 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.695e-03   2.858e-01    8.049e-01  8.049e-01      3119 -3.692e+03 -3.695e+03                   
#> Path [9] :Best Iter: [43] ELBO (-3691.756853) evaluations: (3119) 
#> Path [10] :Initial log joint density = -481498.967518 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.319e-03   2.079e-01    1.000e+00  1.000e+00      2959 -3.688e+03 -3.694e+03                   
#> Path [10] :Best Iter: [44] ELBO (-3688.222063) evaluations: (2959) 
#> Path [11] :Initial log joint density = -481521.112640 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      8.299e-03   2.774e-01    1.000e+00  1.000e+00      2731 -3.690e+03 -3.702e+03                   
#> Path [11] :Best Iter: [45] ELBO (-3690.410182) evaluations: (2731) 
#> Path [12] :Initial log joint density = -483135.835088 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.577e-03   2.328e-01    8.810e-01  8.810e-01      3155 -3.685e+03 -3.698e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3684.818956) evaluations: (3155) 
#> Path [13] :Initial log joint density = -481575.739498 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.113e-02   1.530e-01    1.000e+00  1.000e+00      3573 -3.682e+03 -3.684e+03                   
#> Path [13] :Best Iter: [58] ELBO (-3682.031083) evaluations: (3573) 
#> Path [14] :Initial log joint density = -481763.960910 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.415e-03   2.261e-01    1.000e+00  1.000e+00      2827 -3.691e+03 -3.694e+03                   
#> Path [14] :Best Iter: [46] ELBO (-3691.227838) evaluations: (2827) 
#> Path [15] :Initial log joint density = -481386.069627 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.000e-03   1.793e-01    8.857e-01  8.857e-01      3220 -3.685e+03 -3.695e+03                   
#> Path [15] :Best Iter: [57] ELBO (-3685.041133) evaluations: (3220) 
#> Path [16] :Initial log joint density = -484896.156105 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.226e-02   2.906e-01    1.000e+00  1.000e+00      3706 -3.682e+03 -3.692e+03                   
#> Path [16] :Best Iter: [62] ELBO (-3682.341422) evaluations: (3706) 
#> Path [17] :Initial log joint density = -481747.086581 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.214e-03   2.821e-01    7.805e-01  7.805e-01      3242 -3.687e+03 -3.698e+03                   
#> Path [17] :Best Iter: [56] ELBO (-3687.187680) evaluations: (3242) 
#> Path [18] :Initial log joint density = -481923.194865 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.743e-03   2.463e-01    8.445e-01  8.445e-01      3273 -3.686e+03 -3.698e+03                   
#> Path [18] :Best Iter: [55] ELBO (-3685.852024) evaluations: (3273) 
#> Path [19] :Initial log joint density = -482255.068876 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.292e-03   2.040e-01    8.844e-01  8.844e-01      3355 -3.683e+03 -3.692e+03                   
#> Path [19] :Best Iter: [56] ELBO (-3683.410975) evaluations: (3355) 
#> Path [20] :Initial log joint density = -481546.149047 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.032e-02   2.820e-01    9.455e-01  9.455e-01      2783 -3.690e+03 -3.698e+03                   
#> Path [20] :Best Iter: [49] ELBO (-3689.865635) evaluations: (2783) 
#> Path [21] :Initial log joint density = -481983.863610 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.217e-02   3.877e-01    1.000e+00  1.000e+00      3434 -3.684e+03 -3.692e+03                   
#> Path [21] :Best Iter: [59] ELBO (-3683.726125) evaluations: (3434) 
#> Path [22] :Initial log joint density = -483059.373740 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      2.324e-03   2.824e-01    6.305e-01  6.305e-01      3156 -3.690e+03 -3.698e+03                   
#> Path [22] :Best Iter: [55] ELBO (-3690.129439) evaluations: (3156) 
#> Path [23] :Initial log joint density = -481514.229393 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.468e-02   2.124e-01    1.000e+00  1.000e+00      3077 -3.692e+03 -3.686e+03                   
#> Path [23] :Best Iter: [55] ELBO (-3686.129104) evaluations: (3077) 
#> Path [24] :Initial log joint density = -481384.137221 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.921e-03   1.745e-01    8.792e-01  8.792e-01      2863 -3.690e+03 -3.701e+03                   
#> Path [24] :Best Iter: [41] ELBO (-3690.389443) evaluations: (2863) 
#> Path [25] :Initial log joint density = -481625.698424 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      1.485e-02   3.660e-01    1.000e+00  1.000e+00      3780 -3.682e+03 -3.684e+03                   
#> Path [25] :Best Iter: [63] ELBO (-3681.740819) evaluations: (3780) 
#> Path [26] :Initial log joint density = -481832.418436 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      8.098e-03   1.697e-01    1.000e+00  1.000e+00      3778 -3.684e+03 -3.686e+03                   
#> Path [26] :Best Iter: [59] ELBO (-3684.003621) evaluations: (3778) 
#> Path [27] :Initial log joint density = -481617.906730 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.139e-02   2.545e-01    9.137e-01  9.137e-01      3418 -3.685e+03 -3.693e+03                   
#> Path [27] :Best Iter: [57] ELBO (-3685.384374) evaluations: (3418) 
#> Path [28] :Initial log joint density = -481545.979805 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.246e-02   3.044e-01    1.000e+00  1.000e+00      2909 -3.691e+03 -3.697e+03                   
#> Path [28] :Best Iter: [52] ELBO (-3690.778931) evaluations: (2909) 
#> Path [29] :Initial log joint density = -481483.372353 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.058e-02   2.486e-01    1.000e+00  1.000e+00      3107 -3.685e+03 -3.683e+03                   
#> Path [29] :Best Iter: [56] ELBO (-3683.392447) evaluations: (3107) 
#> Path [30] :Initial log joint density = -481708.749372 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.824e-03   2.266e-01    7.577e-01  7.577e-01      3109 -3.684e+03 -3.692e+03                   
#> Path [30] :Best Iter: [55] ELBO (-3684.146137) evaluations: (3109) 
#> Path [31] :Initial log joint density = -483171.772948 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.666e-02   2.196e-01    1.000e+00  1.000e+00      3238 -3.685e+03 -3.685e+03                   
#> Path [31] :Best Iter: [57] ELBO (-3684.538992) evaluations: (3238) 
#> Path [32] :Initial log joint density = -481388.419279 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.422e-02   1.990e-01    1.000e+00  1.000e+00      3451 -3.681e+03 -3.682e+03                   
#> Path [32] :Best Iter: [59] ELBO (-3680.937473) evaluations: (3451) 
#> Path [33] :Initial log joint density = -481635.814908 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      4.940e-03   1.698e-01    1.000e+00  1.000e+00      3049 -3.689e+03 -3.702e+03                   
#> Path [33] :Best Iter: [46] ELBO (-3688.516342) evaluations: (3049) 
#> Path [34] :Initial log joint density = -481601.435952 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.258e-03   2.100e-01    1.000e+00  1.000e+00      3210 -3.690e+03 -3.693e+03                   
#> Path [34] :Best Iter: [55] ELBO (-3689.835050) evaluations: (3210) 
#> Path [35] :Initial log joint density = -482984.307984 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.348e-02   2.403e-01    1.000e+00  1.000e+00      3185 -3.694e+03 -3.686e+03                   
#> Path [35] :Best Iter: [55] ELBO (-3686.312367) evaluations: (3185) 
#> Path [36] :Initial log joint density = -482363.470665 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.358e-02   3.654e-01    1.000e+00  1.000e+00      3442 -3.681e+03 -3.688e+03                   
#> Path [36] :Best Iter: [56] ELBO (-3681.247203) evaluations: (3442) 
#> Path [37] :Initial log joint density = -481707.239198 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      7.260e-03   2.653e-01    1.000e+00  1.000e+00      2783 -3.692e+03 -3.711e+03                   
#> Path [37] :Best Iter: [48] ELBO (-3692.095405) evaluations: (2783) 
#> Path [38] :Initial log joint density = -483547.559866 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.501e-02   2.324e-01    1.000e+00  1.000e+00      3463 -3.685e+03 -3.685e+03                   
#> Path [38] :Best Iter: [57] ELBO (-3684.555250) evaluations: (3463) 
#> Path [39] :Initial log joint density = -481657.864204 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.519e-02   2.425e-01    1.000e+00  1.000e+00      3332 -3.688e+03 -3.687e+03                   
#> Path [39] :Best Iter: [57] ELBO (-3686.763529) evaluations: (3332) 
#> Path [40] :Initial log joint density = -481565.748249 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.752e-03   2.632e-01    9.388e-01  9.388e-01      3071 -3.691e+03 -3.696e+03                   
#> Path [40] :Best Iter: [44] ELBO (-3691.310433) evaluations: (3071) 
#> Path [41] :Initial log joint density = -481675.690961 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      4.037e-03   2.465e-01    7.902e-01  7.902e-01      2873 -3.690e+03 -3.703e+03                   
#> Path [41] :Best Iter: [46] ELBO (-3690.292354) evaluations: (2873) 
#> Path [42] :Initial log joint density = -481966.652117 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.832e-03   2.851e-01    9.318e-01  9.318e-01      3245 -3.685e+03 -3.692e+03                   
#> Path [42] :Best Iter: [55] ELBO (-3685.061706) evaluations: (3245) 
#> Path [43] :Initial log joint density = -481439.939422 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      7.776e-03   1.594e-01    8.314e-01  8.314e-01      3665 -3.680e+03 -3.689e+03                   
#> Path [43] :Best Iter: [55] ELBO (-3679.897141) evaluations: (3665) 
#> Path [44] :Initial log joint density = -481860.121007 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.681e-02   1.982e-01    1.000e+00  1.000e+00      3596 -3.684e+03 -3.682e+03                   
#> Path [44] :Best Iter: [60] ELBO (-3682.166176) evaluations: (3596) 
#> Path [45] :Initial log joint density = -481559.611205 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.962e-03   2.782e-01    1.000e+00  1.000e+00      3114 -3.690e+03 -3.688e+03                   
#> Path [45] :Best Iter: [55] ELBO (-3688.316597) evaluations: (3114) 
#> Path [46] :Initial log joint density = -481678.419036 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.969e-03   1.633e-01    1.000e+00  1.000e+00      3053 -3.691e+03 -3.694e+03                   
#> Path [46] :Best Iter: [53] ELBO (-3691.186058) evaluations: (3053) 
#> Path [47] :Initial log joint density = -481542.796923 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      7.045e-03   2.807e-01    9.563e-01  9.563e-01      2863 -3.691e+03 -3.705e+03                   
#> Path [47] :Best Iter: [40] ELBO (-3691.256401) evaluations: (2863) 
#> Path [48] :Initial log joint density = -481751.664310 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      4.822e-03   2.665e-01    6.480e-01  6.480e-01      3161 -3.688e+03 -3.697e+03                   
#> Path [48] :Best Iter: [55] ELBO (-3687.790586) evaluations: (3161) 
#> Path [49] :Initial log joint density = -481557.343307 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.328e-03   1.710e-01    8.144e-01  8.144e-01      3133 -3.688e+03 -3.687e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3687.108673) evaluations: (3133) 
#> Path [50] :Initial log joint density = -481568.993317 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.417e-02   1.858e-01    1.000e+00  1.000e+00      2957 -3.692e+03 -3.695e+03                   
#> Path [50] :Best Iter: [46] ELBO (-3692.173606) evaluations: (2957) 
#> Finished in  13.7 seconds.
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
#>  1 10x_6K benign B1                 0.0591           0.0453           0.0753 
#>  2 10x_6K benign B2                 0.0384           0.0294           0.0482 
#>  3 10x_6K benign B3                 0.0126           0.00947          0.0163 
#>  4 10x_6K benign BM                 0.00673          0.00511          0.00867
#>  5 10x_6K benign CD4 1              0.0254           0.0216           0.0300 
#>  6 10x_6K benign CD4 2              0.0509           0.0421           0.0623 
#>  7 10x_6K benign CD4 3              0.0823           0.0642           0.104  
#>  8 10x_6K benign CD4 4              0.00167          0.00116          0.00251
#>  9 10x_6K benign CD4 5              0.0305           0.0240           0.0385 
#> 10 10x_6K benign CD8 1              0.112            0.0963           0.128  
#> # ℹ 710 more rows
#> # ℹ 3 more variables: unconstrained_mean <dbl>, unconstrained_lower <dbl>,
#> #   unconstrained_upper <dbl>
# }

```
