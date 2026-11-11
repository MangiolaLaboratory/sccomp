# Plot 1D Intervals for Cell-group Effects

This function creates a series of 1D interval plots for cell-group
effects, highlighting significant differences based on a given
significance threshold.

## Usage

``` r
sccomp_plot_intervals_1D(
  .data,
  factor = NULL,
  significance_threshold = 0.05,
  test_composition_above_CLR_fold_change = attr(.data,
    "test_composition_above_CLR_fold_change"),
  show_fdr_message = TRUE,
  significance_statistic = c("pH0", "FDR"),
  sort_by = c("none", "effect", "significance", "alphabetical")
)
```

## Arguments

- .data:

  Data frame containing the main data.

- factor:

  Optional character string selecting one model factor to plot. If
  provided, plots are restricted to that factor plus `(Intercept)`.

- significance_threshold:

  Numeric value specifying the significance threshold for highlighting
  differences.

- test_composition_above_CLR_fold_change:

  A positive integer. It is the effect threshold used for the hypothesis
  test. A value of 0.2 correspond to a change in cell proportion of 10%
  for a cell type with baseline proportion of 50%. That is, a cell type
  goes from 45% to 50%. When the baseline proportion is closer to 0 or 1
  this effect thrshold has consistent value in the centered log-ratio
  uncontrained scale.

- show_fdr_message:

  Logical. Whether to show the Bayesian FDR interpretation message on
  the plot. Default is TRUE.

- significance_statistic:

  Character vector indicating which statistic to highlight. Default is
  "pH0".

- sort_by:

  Character vector indicating how to sort taxa. Options are "none"
  (default), "effect" (by effect size), "significance" (by FDR/pH0), or
  "alphabetical".

## Value

A combined plot of 1D interval plots.

## Examples

``` r

print("cmdstanr is needed to run this example.")
#> [1] "cmdstanr is needed to run this example."

# \donttest{
  if (instantiate::stan_cmdstan_exists()) {
    data("counts_obj")

    estimate <- sccomp_estimate(
      counts_obj,
      ~ type,
      ~1,
      "sample",
      "cell_group",
      "count",
      cores = 1
    ) |>
    sccomp_test()

  # Example usage:
  my_plot = sccomp_plot_intervals_1D(estimate, sort_by = "effect")

  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), typecancer
#>   variability Xa: 1 parameter: (Intercept)
#>   random effects: none
#> Loading model from cache...
#> Path [1] :Initial log joint density = -482146.403357 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      9.691e-03   1.892e-01    1.000e+00  1.000e+00      4977 -3.701e+03 -3.706e+03                   
#> Path [1] :Best Iter: [67] ELBO (-3700.696601) evaluations: (4977) 
#> Path [2] :Initial log joint density = -481866.504334 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.129e-02   1.691e-01    1.000e+00  1.000e+00      4569 -3.700e+03 -3.706e+03                   
#> Path [2] :Best Iter: [70] ELBO (-3699.967536) evaluations: (4569) 
#> Path [3] :Initial log joint density = -482124.084995 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      1.262e-02   2.371e-01    1.000e+00  1.000e+00      5219 -3.700e+03 -3.703e+03                   
#> Path [3] :Best Iter: [74] ELBO (-3700.495716) evaluations: (5219) 
#> Path [4] :Initial log joint density = -481675.544247 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.419e-02   2.010e-01    1.000e+00  1.000e+00      3963 -3.700e+03 -3.702e+03                   
#> Path [4] :Best Iter: [57] ELBO (-3700.331463) evaluations: (3963) 
#> Path [5] :Initial log joint density = -481894.538747 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      2.358e-02   1.984e-01    1.000e+00  1.000e+00      5363 -3.699e+03 -3.702e+03                   
#> Path [5] :Best Iter: [77] ELBO (-3699.491863) evaluations: (5363) 
#> Path [6] :Initial log joint density = -483225.638432 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      9.203e-03   2.388e-01    1.000e+00  1.000e+00      5439 -3.698e+03 -3.699e+03                   
#> Path [6] :Best Iter: [76] ELBO (-3697.505833) evaluations: (5439) 
#> Path [7] :Initial log joint density = -484794.774409 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      6.677e-03   2.907e-01    1.000e+00  1.000e+00      3525 -3.706e+03 -3.714e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3706.394414) evaluations: (3525) 
#> Path [8] :Initial log joint density = -481563.172394 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      1.254e-02   2.709e-01    1.000e+00  1.000e+00      3251 -3.704e+03 -3.711e+03                   
#> Path [8] :Best Iter: [55] ELBO (-3704.340916) evaluations: (3251) 
#> Path [9] :Initial log joint density = -481585.032472 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      2.085e-02   2.053e-01    1.000e+00  1.000e+00      4534 -3.701e+03 -3.702e+03                   
#> Path [9] :Best Iter: [63] ELBO (-3701.095297) evaluations: (4534) 
#> Path [10] :Initial log joint density = -481846.005479 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      1.403e-02   2.481e-01    4.333e-01  1.000e+00      5460 -3.699e+03 -3.708e+03                   
#> Path [10] :Best Iter: [74] ELBO (-3698.854078) evaluations: (5460) 
#> Path [11] :Initial log joint density = -481565.177636 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      1.307e-02   3.399e-01    1.000e+00  1.000e+00      3207 -3.706e+03 -3.710e+03                   
#> Path [11] :Best Iter: [56] ELBO (-3706.081680) evaluations: (3207) 
#> Path [12] :Initial log joint density = -481552.880215 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.006e-02   2.581e-01    1.000e+00  1.000e+00      3960 -3.702e+03 -3.708e+03                   
#> Path [12] :Best Iter: [62] ELBO (-3702.381627) evaluations: (3960) 
#> Path [13] :Initial log joint density = -482716.875207 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.217e-02   1.862e-01    7.844e-01  7.844e-01      4767 -3.703e+03 -3.714e+03                   
#> Path [13] :Best Iter: [64] ELBO (-3703.246935) evaluations: (4767) 
#> Path [14] :Initial log joint density = -481581.758785 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.211e-02   1.655e-01    1.000e+00  1.000e+00      4285 -3.702e+03 -3.701e+03                   
#> Path [14] :Best Iter: [67] ELBO (-3701.075917) evaluations: (4285) 
#> Path [15] :Initial log joint density = -481443.042527 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.145e-02   1.933e-01    1.000e+00  1.000e+00      4131 -3.700e+03 -3.708e+03                   
#> Path [15] :Best Iter: [64] ELBO (-3700.016429) evaluations: (4131) 
#> Path [16] :Initial log joint density = -481637.295129 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.378e-02   1.380e-01    1.000e+00  1.000e+00      4084 -3.701e+03 -3.703e+03                   
#> Path [16] :Best Iter: [63] ELBO (-3701.186846) evaluations: (4084) 
#> Path [17] :Initial log joint density = -481759.251367 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.259e-02   3.110e-01    4.644e-01  1.000e+00      4229 -3.703e+03 -3.714e+03                   
#> Path [17] :Best Iter: [61] ELBO (-3702.733801) evaluations: (4229) 
#> Path [18] :Initial log joint density = -481722.201657 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      2.787e-02   2.308e-01    1.000e+00  1.000e+00      4251 -3.701e+03 -3.706e+03                   
#> Path [18] :Best Iter: [63] ELBO (-3701.291775) evaluations: (4251) 
#> Path [19] :Initial log joint density = -481695.684859 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      3.720e-02   2.091e-01    1.000e+00  1.000e+00      3744 -3.705e+03 -3.703e+03                   
#> Path [19] :Best Iter: [60] ELBO (-3702.851746) evaluations: (3744) 
#> Path [20] :Initial log joint density = -481934.783375 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              84      -4.787e+05      1.441e-02   2.482e-01    4.651e-01  1.000e+00      5987 -3.697e+03 -3.708e+03                   
#> Path [20] :Best Iter: [83] ELBO (-3696.530982) evaluations: (5987) 
#> Path [21] :Initial log joint density = -481467.835944 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      1.009e-02   2.912e-01    9.802e-01  9.802e-01      3563 -3.705e+03 -3.715e+03                   
#> Path [21] :Best Iter: [58] ELBO (-3705.366101) evaluations: (3563) 
#> Path [22] :Initial log joint density = -481451.751232 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.751e-02   2.428e-01    1.000e+00  1.000e+00      4037 -3.702e+03 -3.701e+03                   
#> Path [22] :Best Iter: [65] ELBO (-3700.778437) evaluations: (4037) 
#> Path [23] :Initial log joint density = -481551.510467 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.224e-02   1.777e-01    1.000e+00  1.000e+00      3959 -3.702e+03 -3.703e+03                   
#> Path [23] :Best Iter: [61] ELBO (-3701.652813) evaluations: (3959) 
#> Path [24] :Initial log joint density = -481380.911610 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.303e-02   2.657e-01    9.574e-01  9.574e-01      3927 -3.701e+03 -3.714e+03                   
#> Path [24] :Best Iter: [59] ELBO (-3701.462255) evaluations: (3927) 
#> Path [25] :Initial log joint density = -481688.682470 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      3.579e-03   2.112e-01    6.564e-01  6.564e-01      5314 -3.701e+03 -3.711e+03                   
#> Path [25] :Best Iter: [77] ELBO (-3700.506069) evaluations: (5314) 
#> Path [26] :Initial log joint density = -481493.071739 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      8.121e-03   2.108e-01    8.385e-01  8.385e-01      4993 -3.702e+03 -3.709e+03                   
#> Path [26] :Best Iter: [66] ELBO (-3702.336593) evaluations: (4993) 
#> Path [27] :Initial log joint density = -481444.361808 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      9.012e-03   2.069e-01    1.000e+00  1.000e+00      3445 -3.706e+03 -3.710e+03                   
#> Path [27] :Best Iter: [56] ELBO (-3705.507732) evaluations: (3445) 
#> Path [28] :Initial log joint density = -483780.310544 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      5.772e-03   1.679e-01    7.686e-01  7.686e-01      4410 -3.701e+03 -3.713e+03                   
#> Path [28] :Best Iter: [68] ELBO (-3701.309502) evaluations: (4410) 
#> Path [29] :Initial log joint density = -481549.009445 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      7.862e-03   1.912e-01    1.000e+00  1.000e+00      3718 -3.701e+03 -3.702e+03                   
#> Path [29] :Best Iter: [60] ELBO (-3700.852080) evaluations: (3718) 
#> Path [30] :Initial log joint density = -481615.336246 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      5.371e-03   3.217e-01    5.641e-01  5.641e-01      5053 -3.700e+03 -3.713e+03                   
#> Path [30] :Best Iter: [74] ELBO (-3700.109339) evaluations: (5053) 
#> Path [31] :Initial log joint density = -481273.418814 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.410e-02   1.426e-01    1.000e+00  1.000e+00      3943 -3.700e+03 -3.703e+03                   
#> Path [31] :Best Iter: [62] ELBO (-3700.327872) evaluations: (3943) 
#> Path [32] :Initial log joint density = -482398.026250 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.648e-02   2.665e-01    1.000e+00  1.000e+00      5348 -3.697e+03 -3.701e+03                   
#> Path [32] :Best Iter: [77] ELBO (-3696.804905) evaluations: (5348) 
#> Path [33] :Initial log joint density = -482655.642591 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      5.838e-03   2.311e-01    1.000e+00  1.000e+00      4042 -3.703e+03 -3.703e+03                   
#> Path [33] :Best Iter: [62] ELBO (-3702.639083) evaluations: (4042) 
#> Path [34] :Initial log joint density = -481186.803742 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.308e-02   1.720e-01    1.000e+00  1.000e+00      4805 -3.702e+03 -3.702e+03                   
#> Path [34] :Best Iter: [71] ELBO (-3701.961720) evaluations: (4805) 
#> Path [35] :Initial log joint density = -483884.088698 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      4.413e-03   1.417e-01    7.692e-01  7.692e-01      4443 -3.703e+03 -3.714e+03                   
#> Path [35] :Best Iter: [66] ELBO (-3702.529908) evaluations: (4443) 
#> Path [36] :Initial log joint density = -482074.659767 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      5.794e-03   2.205e-01    1.000e+00  1.000e+00      3676 -3.703e+03 -3.715e+03                   
#> Path [36] :Best Iter: [59] ELBO (-3703.181379) evaluations: (3676) 
#> Path [37] :Initial log joint density = -481559.813607 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.079e-02   2.101e-01    1.000e+00  1.000e+00      4056 -3.699e+03 -3.708e+03                   
#> Path [37] :Best Iter: [63] ELBO (-3698.881742) evaluations: (4056) 
#> Path [38] :Initial log joint density = -481480.578086 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.093e-02   3.099e-01    9.694e-01  9.694e-01      3784 -3.704e+03 -3.713e+03                   
#> Path [38] :Best Iter: [56] ELBO (-3704.310126) evaluations: (3784) 
#> Path [39] :Initial log joint density = -481380.110400 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.564e-02   2.485e-01    1.000e+00  1.000e+00      4206 -3.700e+03 -3.706e+03                   
#> Path [39] :Best Iter: [65] ELBO (-3700.120803) evaluations: (4206) 
#> Path [40] :Initial log joint density = -481716.426327 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      4.335e-03   1.914e-01    7.623e-01  7.623e-01      4056 -3.700e+03 -3.713e+03                   
#> Path [40] :Best Iter: [64] ELBO (-3699.734057) evaluations: (4056) 
#> Path [41] :Initial log joint density = -481525.540773 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      4.468e-03   1.919e-01    1.000e+00  1.000e+00      4566 -3.701e+03 -3.713e+03                   
#> Path [41] :Best Iter: [68] ELBO (-3700.692281) evaluations: (4566) 
#> Path [42] :Initial log joint density = -481697.935042 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.271e-02   2.681e-01    1.000e+00  1.000e+00      4706 -3.703e+03 -3.716e+03                   
#> Path [42] :Best Iter: [68] ELBO (-3702.799029) evaluations: (4706) 
#> Path [43] :Initial log joint density = -481636.137245 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      7.532e-03   2.506e-01    1.000e+00  1.000e+00      4359 -3.704e+03 -3.709e+03                   
#> Path [43] :Best Iter: [65] ELBO (-3703.904879) evaluations: (4359) 
#> Path [44] :Initial log joint density = -481721.783427 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      2.657e-03   2.683e-01    6.849e-01  6.849e-01      4482 -3.703e+03 -3.711e+03                   
#> Path [44] :Best Iter: [63] ELBO (-3702.738588) evaluations: (4482) 
#> Path [45] :Initial log joint density = -481425.824065 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.056e-02   2.215e-01    1.000e+00  1.000e+00      5014 -3.700e+03 -3.710e+03                   
#> Path [45] :Best Iter: [72] ELBO (-3700.419875) evaluations: (5014) 
#> Path [46] :Initial log joint density = -481365.926459 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      6.024e-03   2.588e-01    6.864e-01  6.864e-01      3932 -3.704e+03 -3.718e+03                   
#> Path [46] :Best Iter: [58] ELBO (-3703.515462) evaluations: (3932) 
#> Path [47] :Initial log joint density = -482299.527145 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      1.594e-02   3.137e-01    9.666e-01  9.666e-01      3413 -3.703e+03 -3.712e+03                   
#> Path [47] :Best Iter: [58] ELBO (-3703.448655) evaluations: (3413) 
#> Path [48] :Initial log joint density = -481711.548521 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.456e-02   2.650e-01    9.142e-01  9.142e-01      4937 -3.701e+03 -3.715e+03                   
#> Path [48] :Best Iter: [73] ELBO (-3701.236540) evaluations: (4937) 
#> Path [49] :Initial log joint density = -481555.468201 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      7.231e-03   2.422e-01    5.612e-01  5.612e-01      4264 -3.699e+03 -3.709e+03                   
#> Path [49] :Best Iter: [59] ELBO (-3699.232025) evaluations: (4264) 
#> Path [50] :Initial log joint density = -481505.644700 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.003e-02   2.027e-01    1.000e+00  1.000e+00      3794 -3.703e+03 -3.710e+03                   
#> Path [50] :Best Iter: [60] ELBO (-3703.161544) evaluations: (3794) 
#> Finished in  14.6 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Joining with `by = join_by(cell_group, M, parameter)`
# }

```
