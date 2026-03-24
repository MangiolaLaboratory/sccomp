# sccomp_test

This function test contrasts from a sccomp result.

## Usage

``` r
sccomp_test(
  .data,
  contrasts = NULL,
  percent_false_positive = 5,
  test_composition_above_logit_fold_change = 0.1,
  pass_fit = TRUE
)
```

## Arguments

- .data:

  A tibble. The result of sccomp_estimate.

- contrasts:

  A vector of character strings. For example if your formula is
  `~ 0 + treatment` and the factor treatment has values `yes` and `no`,
  your contrast could be "constrasts = c(treatmentyes - treatmentno)".

- percent_false_positive:

  A real between 0 and 100 non included. This used to identify outliers
  with a specific false positive rate.

- test_composition_above_logit_fold_change:

  A positive integer. It is the effect threshold used for the hypothesis
  test. A value of 0.2 correspond to a change in cell proportion of 10%
  for a cell type with baseline proportion of 50%. That is, a cell type
  goes from 45% to 50%. When the baseline proportion is closer to 0 or 1
  this effect thrshold has consistent value in the logit uncontrained
  scale.

- pass_fit:

  A boolean. Whether to pass the Stan fit as attribute in the output.
  Because the Stan fit can be very large, setting this to FALSE can be
  used to lower the memory imprint to save the output.

## Value

A tibble (`tbl`), with the following columns:

- cell_group - The cell groups being tested.

- parameter - The parameter being estimated from the design matrix
  described by the input formula_composition and formula_variability.

- factor - The covariate factor in the formula, if applicable (e.g., not
  present for Intercept or contrasts).

- c_lower - Lower (2.5%) quantile of the posterior distribution for a
  composition (c) parameter.

- c_effect - Mean of the posterior distribution for a composition (c)
  parameter.

- c_upper - Upper (97.5%) quantile of the posterior distribution for a
  composition (c) parameter.

- c_pH0 - Probability of the c_effect being smaller or bigger than the
  `test_composition_above_logit_fold_change` argument.

- c_FDR - False-discovery rate of the null hypothesis (no difference)
  for a composition (c).

- c_n_eff - Effective sample size - the number of independent draws in
  the sample, the higher the better.

- c_R_k_hat - R statistic, a measure of chain equilibrium, should be
  within 0.05 of 1.0.

- v_lower - Lower (2.5%) quantile of the posterior distribution for a
  variability (v) parameter.

- v_effect - Mean of the posterior distribution for a variability (v)
  parameter.

- v_upper - Upper (97.5%) quantile of the posterior distribution for a
  variability (v) parameter.

- v_pH0 - Probability of the null hypothesis (no difference) for a
  variability (v).

- v_FDR - False-discovery rate of the null hypothesis (no difference)
  for a variability (v).

- v_n_eff - Effective sample size for a variability (v) parameter.

- v_R_k_hat - R statistic for a variability (v) parameter.

- count_data - Nested input count data.

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

    estimates = sccomp_estimate(
      counts_obj,
      ~ 0 + type, ~1, "sample", "cell_group", "count",
      cores = 1
    ) |>
    sccomp_test("typecancer - typebenign")
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: typebenign, typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481404.331963 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              41      -4.788e+05      1.028e-02   1.860e-01    1.000e+00  1.000e+00      2087 -3.695e+03 -3.695e+03                   
#> Path [1] :Best Iter: [37] ELBO (-3694.529033) evaluations: (2087) 
#> Path [2] :Initial log joint density = -481776.333318 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      9.981e-03   1.823e-01    1.000e+00  1.000e+00      2763 -3.693e+03 -3.700e+03                   
#> Path [2] :Best Iter: [39] ELBO (-3693.301391) evaluations: (2763) 
#> Path [3] :Initial log joint density = -481200.892690 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      7.085e-03   1.438e-01    1.000e+00  1.000e+00      2305 -3.694e+03 -3.696e+03                   
#> Path [3] :Best Iter: [37] ELBO (-3693.571225) evaluations: (2305) 
#> Path [4] :Initial log joint density = -481595.725794 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.788e-03   2.046e-01    1.000e+00  1.000e+00      3005 -3.693e+03 -3.694e+03                   
#> Path [4] :Best Iter: [44] ELBO (-3692.646383) evaluations: (3005) 
#> Path [5] :Initial log joint density = -481327.650112 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.491e-02   2.149e-01    9.213e-01  9.213e-01      2483 -3.693e+03 -3.698e+03                   
#> Path [5] :Best Iter: [43] ELBO (-3693.096479) evaluations: (2483) 
#> Path [6] :Initial log joint density = -481716.895525 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      9.388e-03   2.149e-01    1.000e+00  1.000e+00      2688 -3.694e+03 -3.700e+03                   
#> Path [6] :Best Iter: [41] ELBO (-3694.292127) evaluations: (2688) 
#> Path [7] :Initial log joint density = -481422.880786 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      9.068e-03   1.714e-01    1.000e+00  1.000e+00      2842 -3.692e+03 -3.700e+03                   
#> Path [7] :Best Iter: [43] ELBO (-3691.605817) evaluations: (2842) 
#> Path [8] :Initial log joint density = -481782.021904 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.406e-02   2.216e-01    1.000e+00  1.000e+00      2858 -3.695e+03 -3.693e+03                   
#> Path [8] :Best Iter: [53] ELBO (-3692.511884) evaluations: (2858) 
#> Path [9] :Initial log joint density = -481593.365695 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.487e-03   1.792e-01    1.000e+00  1.000e+00      3099 -3.695e+03 -3.692e+03                   
#> Path [9] :Best Iter: [55] ELBO (-3691.853181) evaluations: (3099) 
#> Path [10] :Initial log joint density = -487900.064772 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.495e-03   1.267e-01    1.000e+00  1.000e+00      2981 -3.694e+03 -3.702e+03                   
#> Path [10] :Best Iter: [38] ELBO (-3693.797646) evaluations: (2981) 
#> Path [11] :Initial log joint density = -483379.187916 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      8.175e-03   2.151e-01    8.312e-01  8.312e-01      2504 -3.694e+03 -3.697e+03                   
#> Path [11] :Best Iter: [34] ELBO (-3694.208753) evaluations: (2504) 
#> Path [12] :Initial log joint density = -482566.965073 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.361e-02   1.933e-01    1.000e+00  1.000e+00      2723 -3.694e+03 -3.694e+03                   
#> Path [12] :Best Iter: [45] ELBO (-3694.166533) evaluations: (2723) 
#> Path [13] :Initial log joint density = -481229.693237 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      4.178e-03   2.554e-01    6.408e-01  6.408e-01      2267 -3.693e+03 -3.697e+03                   
#> Path [13] :Best Iter: [40] ELBO (-3692.854779) evaluations: (2267) 
#> Path [14] :Initial log joint density = -481516.039366 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      7.984e-03   1.507e-01    1.000e+00  1.000e+00      2626 -3.692e+03 -3.697e+03                   
#> Path [14] :Best Iter: [48] ELBO (-3692.392509) evaluations: (2626) 
#> Path [15] :Initial log joint density = -481180.203998 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              43      -4.788e+05      7.637e-03   2.720e-01    1.000e+00  1.000e+00      2190 -3.694e+03 -3.695e+03                   
#> Path [15] :Best Iter: [32] ELBO (-3693.617245) evaluations: (2190) 
#> Path [16] :Initial log joint density = -484040.125924 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.364e-02   2.593e-01    1.000e+00  1.000e+00      2965 -3.694e+03 -3.693e+03                   
#> Path [16] :Best Iter: [53] ELBO (-3693.393754) evaluations: (2965) 
#> Path [17] :Initial log joint density = -481391.376675 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.362e-02   1.886e-01    1.000e+00  1.000e+00      2789 -3.695e+03 -3.697e+03                   
#> Path [17] :Best Iter: [34] ELBO (-3694.519902) evaluations: (2789) 
#> Path [18] :Initial log joint density = -481671.108853 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      6.461e-03   1.942e-01    1.000e+00  1.000e+00      2778 -3.692e+03 -3.696e+03                   
#> Path [18] :Best Iter: [47] ELBO (-3691.836279) evaluations: (2778) 
#> Path [19] :Initial log joint density = -481521.546883 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              43      -4.788e+05      6.058e-03   2.829e-01    9.313e-01  9.313e-01      2339 -3.694e+03 -3.696e+03                   
#> Path [19] :Best Iter: [41] ELBO (-3693.808474) evaluations: (2339) 
#> Path [20] :Initial log joint density = -481685.638809 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.339e-03   1.851e-01    1.000e+00  1.000e+00      2820 -3.695e+03 -3.699e+03                   
#> Path [20] :Best Iter: [40] ELBO (-3694.730527) evaluations: (2820) 
#> Path [21] :Initial log joint density = -481985.437862 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.054e-02   1.479e-01    8.614e-01  8.614e-01      2596 -3.691e+03 -3.694e+03                   
#> Path [21] :Best Iter: [43] ELBO (-3691.178121) evaluations: (2596) 
#> Path [22] :Initial log joint density = -481627.233885 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.298e-02   2.079e-01    1.000e+00  1.000e+00      2744 -3.695e+03 -3.700e+03                   
#> Path [22] :Best Iter: [50] ELBO (-3694.648074) evaluations: (2744) 
#> Path [23] :Initial log joint density = -483891.446689 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.382e-02   2.202e-01    1.000e+00  1.000e+00      2772 -3.693e+03 -3.696e+03                   
#> Path [23] :Best Iter: [45] ELBO (-3693.261498) evaluations: (2772) 
#> Path [24] :Initial log joint density = -481636.912732 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.323e-02   2.606e-01    1.000e+00  1.000e+00      2918 -3.693e+03 -3.698e+03                   
#> Path [24] :Best Iter: [42] ELBO (-3692.971015) evaluations: (2918) 
#> Path [25] :Initial log joint density = -481562.340383 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      6.973e-03   1.929e-01    1.000e+00  1.000e+00      2544 -3.696e+03 -3.695e+03                   
#> Path [25] :Best Iter: [49] ELBO (-3695.270335) evaluations: (2544) 
#> Path [26] :Initial log joint density = -481417.247496 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.244e-02   2.109e-01    1.000e+00  1.000e+00      2587 -3.692e+03 -3.693e+03                   
#> Path [26] :Best Iter: [42] ELBO (-3692.470492) evaluations: (2587) 
#> Path [27] :Initial log joint density = -481295.040203 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      4.694e-03   1.940e-01    5.216e-01  5.216e-01      2378 -3.694e+03 -3.698e+03                   
#> Path [27] :Best Iter: [34] ELBO (-3693.866857) evaluations: (2378) 
#> Path [28] :Initial log joint density = -481328.863485 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      9.279e-03   2.055e-01    1.000e+00  1.000e+00      2371 -3.694e+03 -3.700e+03                   
#> Path [28] :Best Iter: [30] ELBO (-3694.165078) evaluations: (2371) 
#> Path [29] :Initial log joint density = -481395.141209 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      8.481e-03   2.686e-01    4.742e-01  9.982e-01      2517 -3.695e+03 -3.698e+03                   
#> Path [29] :Best Iter: [40] ELBO (-3694.910549) evaluations: (2517) 
#> Path [30] :Initial log joint density = -481575.487042 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.147e-02   2.060e-01    1.000e+00  1.000e+00      2750 -3.695e+03 -3.695e+03                   
#> Path [30] :Best Iter: [51] ELBO (-3694.664866) evaluations: (2750) 
#> Path [31] :Initial log joint density = -481517.656183 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.262e-02   2.217e-01    9.147e-01  9.147e-01      2798 -3.694e+03 -3.695e+03                   
#> Path [31] :Best Iter: [41] ELBO (-3694.037963) evaluations: (2798) 
#> Path [32] :Initial log joint density = -481247.310018 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      9.890e-03   2.451e-01    1.000e+00  1.000e+00      2623 -3.692e+03 -3.697e+03                   
#> Path [32] :Best Iter: [48] ELBO (-3691.637134) evaluations: (2623) 
#> Path [33] :Initial log joint density = -482380.856333 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      2.661e-02   1.939e-01    1.000e+00  1.000e+00      2376 -3.693e+03 -3.695e+03                   
#> Path [33] :Best Iter: [36] ELBO (-3693.217973) evaluations: (2376) 
#> Path [34] :Initial log joint density = -481318.177364 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      3.520e-03   1.586e-01    6.741e-01  6.741e-01      2651 -3.694e+03 -3.697e+03                   
#> Path [34] :Best Iter: [46] ELBO (-3693.545088) evaluations: (2651) 
#> Path [35] :Initial log joint density = -481597.970683 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      6.754e-03   1.699e-01    1.000e+00  1.000e+00      2530 -3.693e+03 -3.706e+03                   
#> Path [35] :Best Iter: [42] ELBO (-3693.263519) evaluations: (2530) 
#> Path [36] :Initial log joint density = -481545.764156 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.556e-03   1.403e-01    7.888e-01  7.888e-01      3101 -3.693e+03 -3.699e+03                   
#> Path [36] :Best Iter: [46] ELBO (-3693.221688) evaluations: (3101) 
#> Path [37] :Initial log joint density = -481306.600540 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      9.660e-03   1.838e-01    1.000e+00  1.000e+00      2626 -3.694e+03 -3.694e+03                   
#> Path [37] :Best Iter: [47] ELBO (-3693.727709) evaluations: (2626) 
#> Path [38] :Initial log joint density = -481467.312899 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      1.072e-02   2.098e-01    1.000e+00  1.000e+00      2500 -3.695e+03 -3.694e+03                   
#> Path [38] :Best Iter: [48] ELBO (-3693.588503) evaluations: (2500) 
#> Path [39] :Initial log joint density = -481431.548824 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      1.293e-02   2.037e-01    1.000e+00  1.000e+00      2612 -3.694e+03 -3.694e+03                   
#> Path [39] :Best Iter: [48] ELBO (-3694.105363) evaluations: (2612) 
#> Path [40] :Initial log joint density = -481444.690173 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      1.064e-02   1.498e-01    1.000e+00  1.000e+00      2561 -3.694e+03 -3.712e+03                   
#> Path [40] :Best Iter: [41] ELBO (-3694.377563) evaluations: (2561) 
#> Path [41] :Initial log joint density = -481435.188540 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      4.982e-03   1.480e-01    1.000e+00  1.000e+00      2337 -3.694e+03 -3.698e+03                   
#> Path [41] :Best Iter: [40] ELBO (-3694.465017) evaluations: (2337) 
#> Path [42] :Initial log joint density = -481814.819464 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      1.259e-02   1.694e-01    1.000e+00  1.000e+00      2563 -3.693e+03 -3.707e+03                   
#> Path [42] :Best Iter: [42] ELBO (-3692.628978) evaluations: (2563) 
#> Path [43] :Initial log joint density = -481327.491731 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              43      -4.788e+05      7.950e-03   1.773e-01    1.000e+00  1.000e+00      2144 -3.693e+03 -3.696e+03                   
#> Path [43] :Best Iter: [41] ELBO (-3693.443929) evaluations: (2144) 
#> Path [44] :Initial log joint density = -481315.051002 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      8.895e-03   1.611e-01    9.485e-01  9.485e-01      2369 -3.694e+03 -3.696e+03                   
#> Path [44] :Best Iter: [45] ELBO (-3694.031960) evaluations: (2369) 
#> Path [45] :Initial log joint density = -481363.608918 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.797e-02   2.253e-01    1.000e+00  1.000e+00      2591 -3.695e+03 -3.698e+03                   
#> Path [45] :Best Iter: [43] ELBO (-3695.340198) evaluations: (2591) 
#> Path [46] :Initial log joint density = -481841.946522 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.498e-03   1.752e-01    1.000e+00  1.000e+00      3053 -3.694e+03 -3.691e+03                   
#> Path [46] :Best Iter: [55] ELBO (-3690.775735) evaluations: (3053) 
#> Path [47] :Initial log joint density = -481467.543356 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.060e-02   2.908e-01    1.000e+00  1.000e+00      2965 -3.694e+03 -3.695e+03                   
#> Path [47] :Best Iter: [43] ELBO (-3694.050678) evaluations: (2965) 
#> Path [48] :Initial log joint density = -482154.643952 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              42      -4.788e+05      1.217e-02   2.914e-01    9.000e-01  9.000e-01      2153 -3.695e+03 -3.712e+03                   
#> Path [48] :Best Iter: [41] ELBO (-3694.583052) evaluations: (2153) 
#> Path [49] :Initial log joint density = -481966.336842 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.182e-03   1.562e-01    9.155e-01  9.155e-01      2986 -3.693e+03 -3.696e+03                   
#> Path [49] :Best Iter: [50] ELBO (-3693.417376) evaluations: (2986) 
#> Path [50] :Initial log joint density = -481925.113244 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      6.884e-03   1.222e-01    9.916e-01  9.916e-01      2603 -3.693e+03 -3.696e+03                   
#> Path [50] :Best Iter: [34] ELBO (-3692.582867) evaluations: (2603) 
#> Finished in  11.5 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
