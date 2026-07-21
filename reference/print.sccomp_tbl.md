# Print method for sccomp objects

Print method for sccomp objects.

The print method for sccomp objects provides a summary of the model
specifications, data dimensions, and convergence diagnostics.

The output is formatted to be easy to read and understand.

## Usage

``` r
# S3 method for class 'sccomp_tbl'
print(x, ...)
```

## Arguments

- x:

  A sccomp object

- ...:

  Additional arguments passed to print

## Value

The printed object

## Examples

``` r

# Note: Before running the example, ensure that the 'cmdstanr' package is installed:
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))

# \donttest{
if (instantiate::stan_cmdstan_exists()) {
  # Create a sccomp object
  data("counts_obj") 
  estimate <- sccomp_estimate(
    counts_obj,
    ~ type,
    ~1,
    "sample",
    "cell_group",
    "count",
    cores = 1
  )
  print(estimate)
}
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -482847.072075 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      4.375e-03   2.303e-01    8.201e-01  8.201e-01      3590 -3.687e+03 -3.699e+03                   
#> Path [1] :Best Iter: [59] ELBO (-3686.863874) evaluations: (3590) 
#> Path [2] :Initial log joint density = -481742.130766 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      7.542e-03   2.300e-01    1.000e+00  1.000e+00      4684 -3.687e+03 -3.687e+03                   
#> Path [2] :Best Iter: [65] ELBO (-3686.626645) evaluations: (4684) 
#> Path [3] :Initial log joint density = -481801.923933 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      7.783e-03   2.309e-01    9.268e-01  9.268e-01      4248 -3.687e+03 -3.697e+03                   
#> Path [3] :Best Iter: [66] ELBO (-3686.601190) evaluations: (4248) 
#> Path [4] :Initial log joint density = -481582.054843 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.574e-02   2.637e-01    1.000e+00  1.000e+00      4101 -3.688e+03 -3.689e+03                   
#> Path [4] :Best Iter: [57] ELBO (-3687.608146) evaluations: (4101) 
#> Path [5] :Initial log joint density = -481946.972454 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.335e-02   3.595e-01    1.000e+00  1.000e+00      4049 -3.687e+03 -3.694e+03                   
#> Path [5] :Best Iter: [65] ELBO (-3687.376276) evaluations: (4049) 
#> Path [6] :Initial log joint density = -481466.085570 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.551e-02   2.564e-01    1.000e+00  1.000e+00      4086 -3.688e+03 -3.687e+03                   
#> Path [6] :Best Iter: [66] ELBO (-3686.778373) evaluations: (4086) 
#> Path [7] :Initial log joint density = -482622.616017 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.479e-02   2.239e-01    1.000e+00  1.000e+00      3901 -3.689e+03 -3.689e+03                   
#> Path [7] :Best Iter: [59] ELBO (-3689.284358) evaluations: (3901) 
#> Path [8] :Initial log joint density = -481503.675225 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      2.447e-02   2.055e-01    1.000e+00  1.000e+00      5120 -3.681e+03 -3.688e+03                   
#> Path [8] :Best Iter: [70] ELBO (-3681.452815) evaluations: (5120) 
#> Path [9] :Initial log joint density = -481460.084686 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      6.991e-03   1.684e-01    1.000e+00  1.000e+00      4082 -3.684e+03 -3.691e+03                   
#> Path [9] :Best Iter: [64] ELBO (-3683.943133) evaluations: (4082) 
#> Path [10] :Initial log joint density = -481669.192074 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      9.140e-03   2.498e-01    7.438e-01  7.438e-01      5231 -3.680e+03 -3.697e+03                   
#> Path [10] :Best Iter: [75] ELBO (-3680.057837) evaluations: (5231) 
#> Path [11] :Initial log joint density = -481592.917387 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.472e-02   2.239e-01    7.951e-01  7.951e-01      4339 -3.684e+03 -3.695e+03                   
#> Path [11] :Best Iter: [61] ELBO (-3684.056021) evaluations: (4339) 
#> Path [12] :Initial log joint density = -482359.339713 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.319e-02   4.053e-01    1.000e+00  1.000e+00      3890 -3.687e+03 -3.699e+03                   
#> Path [12] :Best Iter: [63] ELBO (-3686.926325) evaluations: (3890) 
#> Path [13] :Initial log joint density = -483604.196560 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      4.416e-03   2.670e-01    6.799e-01  6.799e-01      3552 -3.685e+03 -3.701e+03                   
#> Path [13] :Best Iter: [58] ELBO (-3685.422928) evaluations: (3552) 
#> Path [14] :Initial log joint density = -481599.890423 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.658e-02   1.672e-01    1.000e+00  1.000e+00      4543 -3.685e+03 -3.682e+03                   
#> Path [14] :Best Iter: [71] ELBO (-3682.479400) evaluations: (4543) 
#> Path [15] :Initial log joint density = -482078.041258 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      7.805e-03   1.822e-01    1.000e+00  1.000e+00      4865 -3.686e+03 -3.695e+03                   
#> Path [15] :Best Iter: [71] ELBO (-3686.098796) evaluations: (4865) 
#> Path [16] :Initial log joint density = -481708.106043 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      7.663e-03   1.692e-01    1.000e+00  1.000e+00      3958 -3.687e+03 -3.695e+03                   
#> Path [16] :Best Iter: [64] ELBO (-3687.175336) evaluations: (3958) 
#> Path [17] :Initial log joint density = -481730.453668 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      8.534e-03   1.457e-01    1.000e+00  1.000e+00      4477 -3.684e+03 -3.691e+03                   
#> Path [17] :Best Iter: [62] ELBO (-3684.271148) evaluations: (4477) 
#> Path [18] :Initial log joint density = -481643.702460 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      7.301e-03   2.019e-01    1.000e+00  1.000e+00      3561 -3.687e+03 -3.687e+03                   
#> Path [18] :Best Iter: [60] ELBO (-3686.547510) evaluations: (3561) 
#> Path [19] :Initial log joint density = -481911.869291 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      1.149e-02   4.069e-01    9.398e-01  9.398e-01      3115 -3.687e+03 -3.699e+03                   
#> Path [19] :Best Iter: [55] ELBO (-3687.422928) evaluations: (3115) 
#> Path [20] :Initial log joint density = -481671.294913 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      1.464e-02   3.437e-01    1.000e+00  1.000e+00      3561 -3.685e+03 -3.696e+03                   
#> Path [20] :Best Iter: [58] ELBO (-3685.104025) evaluations: (3561) 
#> Path [21] :Initial log joint density = -482997.829748 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.325e-02   1.826e-01    1.000e+00  1.000e+00      4713 -3.684e+03 -3.689e+03                   
#> Path [21] :Best Iter: [70] ELBO (-3684.180587) evaluations: (4713) 
#> Path [22] :Initial log joint density = -481574.937501 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.787e+05      1.421e-02   2.976e-01    1.000e+00  1.000e+00      2754 -3.694e+03 -3.696e+03                   
#> Path [22] :Best Iter: [43] ELBO (-3694.142656) evaluations: (2754) 
#> Path [23] :Initial log joint density = -482162.372781 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.364e-02   2.510e-01    1.000e+00  1.000e+00      4876 -3.688e+03 -3.685e+03                   
#> Path [23] :Best Iter: [74] ELBO (-3685.370594) evaluations: (4876) 
#> Path [24] :Initial log joint density = -481417.348088 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      5.904e-03   1.657e-01    1.000e+00  1.000e+00      3212 -3.693e+03 -3.697e+03                   
#> Path [24] :Best Iter: [53] ELBO (-3693.022752) evaluations: (3212) 
#> Path [25] :Initial log joint density = -481791.898455 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      6.746e-03   1.980e-01    1.000e+00  1.000e+00      3925 -3.685e+03 -3.693e+03                   
#> Path [25] :Best Iter: [61] ELBO (-3685.462868) evaluations: (3925) 
#> Path [26] :Initial log joint density = -482442.929457 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      7.196e-03   2.197e-01    1.000e+00  1.000e+00      4157 -3.685e+03 -3.693e+03                   
#> Path [26] :Best Iter: [65] ELBO (-3684.678224) evaluations: (4157) 
#> Path [27] :Initial log joint density = -481810.893923 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      4.102e-02   2.801e-01    1.000e+00  1.000e+00      4846 -3.685e+03 -3.685e+03                   
#> Path [27] :Best Iter: [73] ELBO (-3684.722626) evaluations: (4846) 
#> Path [28] :Initial log joint density = -481761.959439 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      4.323e-03   1.444e-01    7.536e-01  7.536e-01      4660 -3.686e+03 -3.700e+03                   
#> Path [28] :Best Iter: [65] ELBO (-3685.864259) evaluations: (4660) 
#> Path [29] :Initial log joint density = -481541.582379 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.445e-02   1.801e-01    8.759e-01  8.759e-01      4303 -3.687e+03 -3.696e+03                   
#> Path [29] :Best Iter: [66] ELBO (-3687.039090) evaluations: (4303) 
#> Path [30] :Initial log joint density = -481905.484390 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      9.323e-03   2.322e-01    1.000e+00  1.000e+00      3384 -3.684e+03 -3.693e+03                   
#> Path [30] :Best Iter: [55] ELBO (-3684.394518) evaluations: (3384) 
#> Path [31] :Initial log joint density = -485021.918054 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      2.274e-03   2.345e-01    6.496e-01  6.496e-01      3776 -3.688e+03 -3.704e+03                   
#> Path [31] :Best Iter: [60] ELBO (-3687.873343) evaluations: (3776) 
#> Path [32] :Initial log joint density = -481734.420693 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.157e-02   2.238e-01    9.235e-01  9.235e-01      3675 -3.685e+03 -3.699e+03                   
#> Path [32] :Best Iter: [61] ELBO (-3685.388489) evaluations: (3675) 
#> Path [33] :Initial log joint density = -481489.588890 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      7.937e-03   1.725e-01    1.000e+00  1.000e+00      4567 -3.684e+03 -3.692e+03                   
#> Path [33] :Best Iter: [65] ELBO (-3684.342377) evaluations: (4567) 
#> Path [34] :Initial log joint density = -481336.211151 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      2.421e-02   1.803e-01    1.000e+00  1.000e+00      4251 -3.685e+03 -3.686e+03                   
#> Path [34] :Best Iter: [64] ELBO (-3685.353995) evaluations: (4251) 
#> Path [35] :Initial log joint density = -481810.294815 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.086e-02   1.492e-01    1.000e+00  1.000e+00      4867 -3.683e+03 -3.690e+03                   
#> Path [35] :Best Iter: [73] ELBO (-3683.321191) evaluations: (4867) 
#> Path [36] :Initial log joint density = -481520.248434 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      2.309e-02   2.736e-01    9.979e-01  9.979e-01      4505 -3.687e+03 -3.695e+03                   
#> Path [36] :Best Iter: [69] ELBO (-3687.382106) evaluations: (4505) 
#> Path [37] :Initial log joint density = -481859.824757 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      5.421e-03   1.636e-01    1.000e+00  1.000e+00      4238 -3.682e+03 -3.693e+03                   
#> Path [37] :Best Iter: [64] ELBO (-3681.799882) evaluations: (4238) 
#> Path [38] :Initial log joint density = -481796.287059 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      5.776e-03   2.297e-01    7.428e-01  7.428e-01      3260 -3.693e+03 -3.704e+03                   
#> Path [38] :Best Iter: [54] ELBO (-3692.678181) evaluations: (3260) 
#> Path [39] :Initial log joint density = -481714.110615 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      1.116e-02   1.887e-01    1.000e+00  1.000e+00      3227 -3.696e+03 -3.694e+03                   
#> Path [39] :Best Iter: [57] ELBO (-3694.035297) evaluations: (3227) 
#> Path [40] :Initial log joint density = -481589.630846 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.274e-02   3.126e-01    1.000e+00  1.000e+00      5134 -3.686e+03 -3.689e+03                   
#> Path [40] :Best Iter: [68] ELBO (-3686.409034) evaluations: (5134) 
#> Path [41] :Initial log joint density = -483433.697132 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      3.163e-03   3.223e-01    6.214e-01  6.214e-01      3343 -3.693e+03 -3.702e+03                   
#> Path [41] :Best Iter: [45] ELBO (-3692.544335) evaluations: (3343) 
#> Path [42] :Initial log joint density = -482104.910410 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.418e-02   2.849e-01    1.000e+00  1.000e+00      4365 -3.690e+03 -3.696e+03                   
#> Path [42] :Best Iter: [67] ELBO (-3689.712807) evaluations: (4365) 
#> Path [43] :Initial log joint density = -483688.076518 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      1.492e-02   2.492e-01    1.000e+00  1.000e+00      4768 -3.685e+03 -3.689e+03                   
#> Path [43] :Best Iter: [70] ELBO (-3685.216618) evaluations: (4768) 
#> Path [44] :Initial log joint density = -482045.543592 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.047e-02   2.154e-01    1.000e+00  1.000e+00      4209 -3.688e+03 -3.692e+03                   
#> Path [44] :Best Iter: [65] ELBO (-3687.634913) evaluations: (4209) 
#> Path [45] :Initial log joint density = -481414.227026 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      6.335e-03   1.827e-01    1.000e+00  1.000e+00      4407 -3.686e+03 -3.693e+03                   
#> Path [45] :Best Iter: [67] ELBO (-3685.892970) evaluations: (4407) 
#> Path [46] :Initial log joint density = -481854.536376 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.932e-02   2.525e-01    8.464e-01  8.464e-01      4600 -3.688e+03 -3.695e+03                   
#> Path [46] :Best Iter: [61] ELBO (-3687.888134) evaluations: (4600) 
#> Path [47] :Initial log joint density = -481988.086869 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.102e-02   1.820e-01    1.000e+00  1.000e+00      3590 -3.694e+03 -3.694e+03                   
#> Path [47] :Best Iter: [61] ELBO (-3694.397190) evaluations: (3590) 
#> Path [48] :Initial log joint density = -481659.296770 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.254e-02   3.009e-01    1.000e+00  1.000e+00      3478 -3.689e+03 -3.693e+03                   
#> Path [48] :Best Iter: [56] ELBO (-3688.629380) evaluations: (3478) 
#> Path [49] :Initial log joint density = -481676.493026 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      7.334e-03   1.934e-01    1.000e+00  1.000e+00      4789 -3.686e+03 -3.691e+03                   
#> Path [49] :Best Iter: [72] ELBO (-3686.490950) evaluations: (4789) 
#> Path [50] :Initial log joint density = -482183.554982 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      4.070e-03   2.915e-01    6.785e-01  6.785e-01      4742 -3.686e+03 -3.700e+03                   
#> Path [50] :Best Iter: [70] ELBO (-3686.189672) evaluations: (4742) 
#> Finished in  15.9 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> sccomp model
#> ============
#> 
#> Model specifications:
#>   Family: multi_beta_binomial 
#>   Composition formula: ~type 
#>   Variability formula: ~1 
#>   Inference method: pathfinder 
#> 
#> Data: Samples: 20   Cell groups: 36 
#> 
#> Column prefixes: c_ -> composition parameters  v_ -> variability parameters
#> 
#> Convergence diagnostics:
#>   For each parameter, n_eff is the effective sample size and R_k_hat is the potential
#>   scale reduction factor on split chains (at convergence, R_k_hat = 1).
#> 
#> # A tibble: 72 × 15
#>    cell_group parameter   factor c_lower c_effect  c_upper c_rhat c_ess_bulk
#>    <chr>      <chr>       <chr>    <dbl>    <dbl>    <dbl>  <dbl>      <dbl>
#>  1 B1         (Intercept) NA      0.947     1.20   1.45     1.00       3775.
#>  2 B1         typecancer  type   -0.933    -0.621 -0.303    1.00       3868.
#>  3 B2         (Intercept) NA      0.495     0.764  1.03     1.01       3905.
#>  4 B2         typecancer  type   -0.995    -0.665 -0.335    1.00       3691.
#>  5 B3         (Intercept) NA     -0.592    -0.336 -0.0674   1.00       3700.
#>  6 B3         typecancer  type   -0.603    -0.286  0.0311   1.00       3722.
#>  7 BM         (Intercept) NA     -1.23     -0.975 -0.723    1.00       3770.
#>  8 BM         typecancer  type   -0.614    -0.307  0.00820  1.000      4030.
#>  9 CD4 1      (Intercept) NA      0.209     0.368  0.533    1.00       4082.
#> 10 CD4 1      typecancer  type   -0.0254    0.191  0.413    1.00       3980.
#> # ℹ 62 more rows
#> # ℹ 7 more variables: c_ess_tail <dbl>, v_upper <dbl>, v_effect <dbl>,
#> #   v_lower <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>
# }
```
