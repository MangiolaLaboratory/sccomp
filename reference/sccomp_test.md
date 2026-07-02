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
#> Path [1] :Initial log joint density = -481550.964703 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.215e-02   1.306e-01    1.000e+00  1.000e+00      2762 -3.691e+03 -3.695e+03                   
#> Path [1] :Best Iter: [45] ELBO (-3691.119746) evaluations: (2762) 
#> Path [2] :Initial log joint density = -481385.598781 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      1.021e-02   1.936e-01    1.000e+00  1.000e+00      2291 -3.693e+03 -3.695e+03                   
#> Path [2] :Best Iter: [41] ELBO (-3692.670728) evaluations: (2291) 
#> Path [3] :Initial log joint density = -481689.252762 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.261e-02   3.038e-01    1.000e+00  1.000e+00      2718 -3.692e+03 -3.695e+03                   
#> Path [3] :Best Iter: [47] ELBO (-3691.694903) evaluations: (2718) 
#> Path [4] :Initial log joint density = -481221.961336 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      8.607e-03   1.559e-01    1.000e+00  1.000e+00      2848 -3.691e+03 -3.693e+03                   
#> Path [4] :Best Iter: [46] ELBO (-3691.459986) evaluations: (2848) 
#> Path [5] :Initial log joint density = -481897.820161 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.353e-02   1.394e-01    1.000e+00  1.000e+00      2797 -3.694e+03 -3.693e+03                   
#> Path [5] :Best Iter: [51] ELBO (-3693.196989) evaluations: (2797) 
#> Path [6] :Initial log joint density = -481366.164051 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      6.526e-03   2.972e-01    1.000e+00  1.000e+00      2446 -3.696e+03 -3.702e+03                   
#> Path [6] :Best Iter: [43] ELBO (-3695.937998) evaluations: (2446) 
#> Path [7] :Initial log joint density = -481795.852666 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.734e-03   2.199e-01    1.000e+00  1.000e+00      3003 -3.693e+03 -3.694e+03                   
#> Path [7] :Best Iter: [51] ELBO (-3692.726421) evaluations: (3003) 
#> Path [8] :Initial log joint density = -482341.651942 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.358e-02   2.140e-01    1.000e+00  1.000e+00      3139 -3.695e+03 -3.694e+03                   
#> Path [8] :Best Iter: [54] ELBO (-3693.895168) evaluations: (3139) 
#> Path [9] :Initial log joint density = -481405.459548 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      5.029e-03   2.065e-01    8.702e-01  8.702e-01      2480 -3.693e+03 -3.699e+03                   
#> Path [9] :Best Iter: [46] ELBO (-3692.539320) evaluations: (2480) 
#> Path [10] :Initial log joint density = -481670.147487 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      1.712e-02   2.847e-01    8.940e-01  8.940e-01      2302 -3.695e+03 -3.697e+03                   
#> Path [10] :Best Iter: [33] ELBO (-3694.804714) evaluations: (2302) 
#> Path [11] :Initial log joint density = -481408.450723 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.186e-02   2.386e-01    1.000e+00  1.000e+00      2588 -3.694e+03 -3.693e+03                   
#> Path [11] :Best Iter: [49] ELBO (-3693.151525) evaluations: (2588) 
#> Path [12] :Initial log joint density = -481818.869260 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.076e-03   2.425e-01    6.815e-01  6.815e-01      3150 -3.690e+03 -3.693e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3689.792577) evaluations: (3150) 
#> Path [13] :Initial log joint density = -481543.594764 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      5.154e-03   1.521e-01    9.321e-01  9.321e-01      2263 -3.692e+03 -3.695e+03                   
#> Path [13] :Best Iter: [32] ELBO (-3692.438874) evaluations: (2263) 
#> Path [14] :Initial log joint density = -481360.343731 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      9.323e-03   1.606e-01    1.000e+00  1.000e+00      2768 -3.692e+03 -3.694e+03                   
#> Path [14] :Best Iter: [36] ELBO (-3692.038022) evaluations: (2768) 
#> Path [15] :Initial log joint density = -483380.100470 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.747e-02   2.280e-01    1.000e+00  1.000e+00      2567 -3.695e+03 -3.694e+03                   
#> Path [15] :Best Iter: [47] ELBO (-3694.340136) evaluations: (2567) 
#> Path [16] :Initial log joint density = -481478.943130 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.214e-02   2.068e-01    1.000e+00  1.000e+00      2650 -3.693e+03 -3.692e+03                   
#> Path [16] :Best Iter: [50] ELBO (-3692.326944) evaluations: (2650) 
#> Path [17] :Initial log joint density = -481343.287031 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      1.522e-02   2.575e-01    1.000e+00  1.000e+00      2455 -3.693e+03 -3.694e+03                   
#> Path [17] :Best Iter: [41] ELBO (-3693.308165) evaluations: (2455) 
#> Path [18] :Initial log joint density = -481715.683558 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.018e-02   1.754e-01    1.000e+00  1.000e+00      2576 -3.694e+03 -3.694e+03                   
#> Path [18] :Best Iter: [35] ELBO (-3693.568238) evaluations: (2576) 
#> Path [19] :Initial log joint density = -481500.748351 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      6.984e-03   3.116e-01    1.000e+00  1.000e+00      2685 -3.693e+03 -3.695e+03                   
#> Path [19] :Best Iter: [47] ELBO (-3693.242205) evaluations: (2685) 
#> Path [20] :Initial log joint density = -481387.817261 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      6.235e-03   1.752e-01    9.612e-01  9.612e-01      2285 -3.692e+03 -3.693e+03                   
#> Path [20] :Best Iter: [42] ELBO (-3692.184865) evaluations: (2285) 
#> Path [21] :Initial log joint density = -481308.293194 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      7.550e-03   1.952e-01    1.000e+00  1.000e+00      2917 -3.690e+03 -3.696e+03                   
#> Path [21] :Best Iter: [50] ELBO (-3690.482060) evaluations: (2917) 
#> Path [22] :Initial log joint density = -481936.421709 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.779e-03   1.725e-01    1.000e+00  1.000e+00      3082 -3.691e+03 -3.689e+03                   
#> Path [22] :Best Iter: [56] ELBO (-3689.210072) evaluations: (3082) 
#> Path [23] :Initial log joint density = -481455.156733 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.276e-02   2.684e-01    1.000e+00  1.000e+00      2482 -3.695e+03 -3.696e+03                   
#> Path [23] :Best Iter: [31] ELBO (-3694.508723) evaluations: (2482) 
#> Path [24] :Initial log joint density = -481473.013563 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      7.937e-03   2.054e-01    1.000e+00  1.000e+00      2473 -3.696e+03 -3.699e+03                   
#> Path [24] :Best Iter: [40] ELBO (-3696.078103) evaluations: (2473) 
#> Path [25] :Initial log joint density = -485188.059219 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.324e-02   2.050e-01    1.000e+00  1.000e+00      2970 -3.692e+03 -3.697e+03                   
#> Path [25] :Best Iter: [45] ELBO (-3692.479850) evaluations: (2970) 
#> Path [26] :Initial log joint density = -481411.006264 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      1.003e-02   1.705e-01    1.000e+00  1.000e+00      2555 -3.694e+03 -3.695e+03                   
#> Path [26] :Best Iter: [42] ELBO (-3694.061621) evaluations: (2555) 
#> Path [27] :Initial log joint density = -481371.878982 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      9.217e-03   1.917e-01    1.000e+00  1.000e+00      2521 -3.695e+03 -3.696e+03                   
#> Path [27] :Best Iter: [33] ELBO (-3695.001587) evaluations: (2521) 
#> Path [28] :Initial log joint density = -481649.160075 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.178e-02   2.899e-01    1.000e+00  1.000e+00      3132 -3.690e+03 -3.689e+03                   
#> Path [28] :Best Iter: [56] ELBO (-3689.163167) evaluations: (3132) 
#> Path [29] :Initial log joint density = -481634.453742 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      8.673e-03   2.581e-01    1.000e+00  1.000e+00      2505 -3.694e+03 -3.695e+03                   
#> Path [29] :Best Iter: [33] ELBO (-3694.214914) evaluations: (2505) 
#> Path [30] :Initial log joint density = -481724.688391 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.657e-03   1.877e-01    8.231e-01  8.231e-01      3188 -3.691e+03 -3.693e+03                   
#> Path [30] :Best Iter: [56] ELBO (-3690.940497) evaluations: (3188) 
#> Path [31] :Initial log joint density = -481651.881676 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.112e-02   2.672e-01    1.000e+00  1.000e+00      3265 -3.690e+03 -3.690e+03                   
#> Path [31] :Best Iter: [57] ELBO (-3689.603117) evaluations: (3265) 
#> Path [32] :Initial log joint density = -481550.053877 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.843e-02   2.390e-01    1.000e+00  1.000e+00      2920 -3.693e+03 -3.697e+03                   
#> Path [32] :Best Iter: [43] ELBO (-3693.094125) evaluations: (2920) 
#> Path [33] :Initial log joint density = -481670.115613 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      9.061e-03   1.850e-01    1.000e+00  1.000e+00      2549 -3.694e+03 -3.693e+03                   
#> Path [33] :Best Iter: [49] ELBO (-3693.288278) evaluations: (2549) 
#> Path [34] :Initial log joint density = -481372.278767 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      3.658e-03   2.760e-01    5.446e-01  5.446e-01      2840 -3.693e+03 -3.698e+03                   
#> Path [34] :Best Iter: [37] ELBO (-3693.357880) evaluations: (2840) 
#> Path [35] :Initial log joint density = -482720.058701 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      8.720e-03   1.659e-01    1.000e+00  1.000e+00      2481 -3.693e+03 -3.699e+03                   
#> Path [35] :Best Iter: [36] ELBO (-3692.700852) evaluations: (2481) 
#> Path [36] :Initial log joint density = -481126.289521 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      6.883e-03   2.479e-01    8.787e-01  8.787e-01      2300 -3.695e+03 -3.699e+03                   
#> Path [36] :Best Iter: [37] ELBO (-3694.781392) evaluations: (2300) 
#> Path [37] :Initial log joint density = -481694.490191 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      3.138e-03   1.850e-01    5.745e-01  5.745e-01      3011 -3.693e+03 -3.696e+03                   
#> Path [37] :Best Iter: [46] ELBO (-3692.855677) evaluations: (3011) 
#> Path [38] :Initial log joint density = -481403.112427 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.089e-02   1.902e-01    1.000e+00  1.000e+00      2939 -3.694e+03 -3.694e+03                   
#> Path [38] :Best Iter: [50] ELBO (-3693.533715) evaluations: (2939) 
#> Path [39] :Initial log joint density = -484099.953729 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      2.870e-03   2.935e-01    6.515e-01  6.515e-01      2742 -3.695e+03 -3.706e+03                   
#> Path [39] :Best Iter: [37] ELBO (-3695.221615) evaluations: (2742) 
#> Path [40] :Initial log joint density = -481458.793692 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      9.392e-03   1.360e-01    1.000e+00  1.000e+00      2709 -3.693e+03 -3.696e+03                   
#> Path [40] :Best Iter: [39] ELBO (-3693.056862) evaluations: (2709) 
#> Path [41] :Initial log joint density = -481320.420513 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      6.512e-03   1.673e-01    1.000e+00  1.000e+00      2360 -3.693e+03 -3.694e+03                   
#> Path [41] :Best Iter: [43] ELBO (-3693.444026) evaluations: (2360) 
#> Path [42] :Initial log joint density = -481625.590365 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.202e-02   1.461e-01    1.000e+00  1.000e+00      2506 -3.694e+03 -3.703e+03                   
#> Path [42] :Best Iter: [45] ELBO (-3694.193957) evaluations: (2506) 
#> Path [43] :Initial log joint density = -483071.365093 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.132e-03   2.042e-01    1.000e+00  1.000e+00      2995 -3.692e+03 -3.695e+03                   
#> Path [43] :Best Iter: [47] ELBO (-3691.670208) evaluations: (2995) 
#> Path [44] :Initial log joint density = -481514.797675 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.266e-02   1.757e-01    1.000e+00  1.000e+00      3063 -3.692e+03 -3.691e+03                   
#> Path [44] :Best Iter: [55] ELBO (-3691.187085) evaluations: (3063) 
#> Path [45] :Initial log joint density = -481155.625838 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              42      -4.788e+05      7.973e-03   1.340e-01    8.380e-01  8.380e-01      2112 -3.693e+03 -3.695e+03                   
#> Path [45] :Best Iter: [33] ELBO (-3692.574433) evaluations: (2112) 
#> Path [46] :Initial log joint density = -481540.073453 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      8.527e-03   2.505e-01    7.942e-01  7.942e-01      2293 -3.693e+03 -3.693e+03                   
#> Path [46] :Best Iter: [44] ELBO (-3692.509145) evaluations: (2293) 
#> Path [47] :Initial log joint density = -481430.676854 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.628e-02   2.287e-01    8.818e-01  8.818e-01      2670 -3.693e+03 -3.702e+03                   
#> Path [47] :Best Iter: [46] ELBO (-3693.254037) evaluations: (2670) 
#> Path [48] :Initial log joint density = -481399.880373 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      2.900e-03   2.628e-01    5.587e-01  5.587e-01      2369 -3.694e+03 -3.700e+03                   
#> Path [48] :Best Iter: [29] ELBO (-3693.744810) evaluations: (2369) 
#> Path [49] :Initial log joint density = -482494.823197 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.747e-03   1.387e-01    1.000e+00  1.000e+00      3131 -3.693e+03 -3.697e+03                   
#> Path [49] :Best Iter: [52] ELBO (-3693.091660) evaluations: (3131) 
#> Path [50] :Initial log joint density = -481516.237413 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.747e-03   1.889e-01    3.934e-01  1.000e+00      3003 -3.695e+03 -3.696e+03                   
#> Path [50] :Best Iter: [52] ELBO (-3695.080001) evaluations: (3003) 
#> Finished in  11.9 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
