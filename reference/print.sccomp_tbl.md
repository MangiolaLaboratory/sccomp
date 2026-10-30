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
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), typecancer
#>   variability Xa: 1 parameter: (Intercept)
#>   random effects: none
#> Loading model from cache...
#> Path [1] :Initial log joint density = -482833.232298 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.671e-02   2.798e-01    9.152e-01  9.152e-01      4206 -3.703e+03 -3.711e+03                   
#> Path [1] :Best Iter: [65] ELBO (-3702.596166) evaluations: (4206) 
#> Path [2] :Initial log joint density = -481728.691881 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.459e-02   2.342e-01    1.000e+00  1.000e+00      5257 -3.698e+03 -3.703e+03                   
#> Path [2] :Best Iter: [75] ELBO (-3698.021671) evaluations: (5257) 
#> Path [3] :Initial log joint density = -481783.552065 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      7.408e-03   2.098e-01    1.000e+00  1.000e+00      4150 -3.701e+03 -3.706e+03                   
#> Path [3] :Best Iter: [63] ELBO (-3701.161871) evaluations: (4150) 
#> Path [4] :Initial log joint density = -481566.613410 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      9.883e-03   2.136e-01    4.060e-01  1.000e+00      4057 -3.700e+03 -3.712e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3699.870620) evaluations: (4057) 
#> Path [5] :Initial log joint density = -481928.590011 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.665e-02   1.738e-01    1.000e+00  1.000e+00      4863 -3.700e+03 -3.699e+03                   
#> Path [5] :Best Iter: [74] ELBO (-3699.194032) evaluations: (4863) 
#> Path [6] :Initial log joint density = -481450.961872 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      6.634e-03   2.035e-01    1.000e+00  1.000e+00      4089 -3.699e+03 -3.706e+03                   
#> Path [6] :Best Iter: [62] ELBO (-3699.484471) evaluations: (4089) 
#> Path [7] :Initial log joint density = -482606.469559 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      3.672e-03   1.731e-01    7.404e-01  7.404e-01      4584 -3.699e+03 -3.716e+03                   
#> Path [7] :Best Iter: [68] ELBO (-3699.299438) evaluations: (4584) 
#> Path [8] :Initial log joint density = -481486.269206 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      1.946e-02   2.202e-01    1.000e+00  1.000e+00      5117 -3.696e+03 -3.702e+03                   
#> Path [8] :Best Iter: [75] ELBO (-3695.871642) evaluations: (5117) 
#> Path [9] :Initial log joint density = -481448.010455 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.219e-02   2.055e-01    4.445e-01  1.000e+00      4026 -3.699e+03 -3.708e+03                   
#> Path [9] :Best Iter: [65] ELBO (-3699.110828) evaluations: (4026) 
#> Path [10] :Initial log joint density = -481651.079473 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      2.147e-02   1.867e-01    1.000e+00  1.000e+00      4603 -3.701e+03 -3.702e+03                   
#> Path [10] :Best Iter: [67] ELBO (-3701.004639) evaluations: (4603) 
#> Path [11] :Initial log joint density = -481574.046218 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              94      -4.787e+05      2.686e-02   1.158e-01    1.000e+00  1.000e+00      7102 -3.694e+03 -3.699e+03                   
#> Path [11] :Best Iter: [90] ELBO (-3694.199546) evaluations: (7102) 
#> Path [12] :Initial log joint density = -482343.798688 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      6.694e-03   1.584e-01    1.000e+00  1.000e+00      4339 -3.701e+03 -3.713e+03                   
#> Path [12] :Best Iter: [66] ELBO (-3700.966067) evaluations: (4339) 
#> Path [13] :Initial log joint density = -483588.295898 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.312e-02   2.856e-01    8.768e-01  8.768e-01      4015 -3.702e+03 -3.710e+03                   
#> Path [13] :Best Iter: [63] ELBO (-3701.539995) evaluations: (4015) 
#> Path [14] :Initial log joint density = -481588.407827 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      6.536e-03   1.980e-01    1.000e+00  1.000e+00      4349 -3.701e+03 -3.713e+03                   
#> Path [14] :Best Iter: [61] ELBO (-3700.755526) evaluations: (4349) 
#> Path [15] :Initial log joint density = -482059.461207 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.787e+05      1.059e-02   1.438e-01    9.324e-01  9.324e-01      5616 -3.700e+03 -3.709e+03                   
#> Path [15] :Best Iter: [75] ELBO (-3700.124748) evaluations: (5616) 
#> Path [16] :Initial log joint density = -481689.102179 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.896e-02   1.640e-01    1.000e+00  1.000e+00      4900 -3.701e+03 -3.705e+03                   
#> Path [16] :Best Iter: [71] ELBO (-3701.307903) evaluations: (4900) 
#> Path [17] :Initial log joint density = -481712.404158 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      3.519e-03   2.220e-01    6.882e-01  6.882e-01      4487 -3.699e+03 -3.711e+03                   
#> Path [17] :Best Iter: [68] ELBO (-3699.289525) evaluations: (4487) 
#> Path [18] :Initial log joint density = -481630.636680 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      5.809e-03   2.328e-01    7.309e-01  7.309e-01      3927 -3.701e+03 -3.715e+03                   
#> Path [18] :Best Iter: [61] ELBO (-3701.258279) evaluations: (3927) 
#> Path [19] :Initial log joint density = -481892.515968 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      3.532e-03   2.551e-01    6.546e-01  6.546e-01      4008 -3.702e+03 -3.714e+03                   
#> Path [19] :Best Iter: [61] ELBO (-3702.368510) evaluations: (4008) 
#> Path [20] :Initial log joint density = -481653.336160 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      5.379e-03   1.769e-01    1.000e+00  1.000e+00      3310 -3.710e+03 -3.717e+03                   
#> Path [20] :Best Iter: [53] ELBO (-3709.724613) evaluations: (3310) 
#> Path [21] :Initial log joint density = -482985.675370 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.836e-02   2.807e-01    1.000e+00  1.000e+00      4475 -3.698e+03 -3.703e+03                   
#> Path [21] :Best Iter: [68] ELBO (-3698.435131) evaluations: (4475) 
#> Path [22] :Initial log joint density = -481559.739762 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      1.269e-02   2.018e-01    1.000e+00  1.000e+00      3533 -3.704e+03 -3.702e+03                   
#> Path [22] :Best Iter: [60] ELBO (-3701.787776) evaluations: (3533) 
#> Path [23] :Initial log joint density = -482149.261450 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      1.358e-02   2.081e-01    1.000e+00  1.000e+00      5511 -3.699e+03 -3.704e+03                   
#> Path [23] :Best Iter: [79] ELBO (-3699.399504) evaluations: (5511) 
#> Path [24] :Initial log joint density = -481404.604802 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      5.369e-03   1.684e-01    7.965e-01  7.965e-01      4026 -3.702e+03 -3.715e+03                   
#> Path [24] :Best Iter: [63] ELBO (-3701.685173) evaluations: (4026) 
#> Path [25] :Initial log joint density = -481776.066697 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      6.990e-03   1.864e-01    7.370e-01  7.370e-01      3830 -3.702e+03 -3.717e+03                   
#> Path [25] :Best Iter: [61] ELBO (-3701.901877) evaluations: (3830) 
#> Path [26] :Initial log joint density = -482429.250169 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      9.491e-03   2.119e-01    1.000e+00  1.000e+00      4169 -3.702e+03 -3.710e+03                   
#> Path [26] :Best Iter: [62] ELBO (-3701.919970) evaluations: (4169) 
#> Path [27] :Initial log joint density = -481794.562328 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      2.437e-02   3.490e-01    1.000e+00  1.000e+00      4053 -3.700e+03 -3.704e+03                   
#> Path [27] :Best Iter: [65] ELBO (-3700.497339) evaluations: (4053) 
#> Path [28] :Initial log joint density = -481747.549333 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.414e-02   3.585e-01    1.000e+00  1.000e+00      4589 -3.701e+03 -3.708e+03                   
#> Path [28] :Best Iter: [70] ELBO (-3700.575229) evaluations: (4589) 
#> Path [29] :Initial log joint density = -481524.650415 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      3.026e-02   2.837e-01    8.329e-01  8.329e-01      3878 -3.702e+03 -3.715e+03                   
#> Path [29] :Best Iter: [61] ELBO (-3702.430074) evaluations: (3878) 
#> Path [30] :Initial log joint density = -481892.241967 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.118e-02   2.475e-01    1.000e+00  1.000e+00      3908 -3.703e+03 -3.704e+03                   
#> Path [30] :Best Iter: [63] ELBO (-3703.211877) evaluations: (3908) 
#> Path [31] :Initial log joint density = -485006.283802 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      5.836e-03   2.522e-01    1.000e+00  1.000e+00      3552 -3.703e+03 -3.707e+03                   
#> Path [31] :Best Iter: [57] ELBO (-3703.247390) evaluations: (3552) 
#> Path [32] :Initial log joint density = -481716.192960 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.787e+05      9.612e-03   2.077e-01    1.000e+00  1.000e+00      5872 -3.698e+03 -3.701e+03                   
#> Path [32] :Best Iter: [84] ELBO (-3698.175982) evaluations: (5872) 
#> Path [33] :Initial log joint density = -481472.295653 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.416e-02   2.538e-01    1.000e+00  1.000e+00      4567 -3.700e+03 -3.705e+03                   
#> Path [33] :Best Iter: [65] ELBO (-3700.109896) evaluations: (4567) 
#> Path [34] :Initial log joint density = -481319.747264 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      2.469e-02   1.759e-01    1.000e+00  1.000e+00      4176 -3.703e+03 -3.703e+03                   
#> Path [34] :Best Iter: [62] ELBO (-3702.771824) evaluations: (4176) 
#> Path [35] :Initial log joint density = -481798.613940 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      6.647e-03   1.210e-01    1.000e+00  1.000e+00      5550 -3.700e+03 -3.712e+03                   
#> Path [35] :Best Iter: [75] ELBO (-3700.356883) evaluations: (5550) 
#> Path [36] :Initial log joint density = -481507.343900 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.109e-02   2.037e-01    1.000e+00  1.000e+00      4451 -3.701e+03 -3.704e+03                   
#> Path [36] :Best Iter: [63] ELBO (-3701.366187) evaluations: (4451) 
#> Path [37] :Initial log joint density = -481844.037160 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      8.972e-03   2.305e-01    1.000e+00  1.000e+00      3978 -3.701e+03 -3.702e+03                   
#> Path [37] :Best Iter: [60] ELBO (-3701.451445) evaluations: (3978) 
#> Path [38] :Initial log joint density = -481779.791486 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      2.003e-02   1.757e-01    1.000e+00  1.000e+00      3691 -3.707e+03 -3.706e+03                   
#> Path [38] :Best Iter: [61] ELBO (-3706.314288) evaluations: (3691) 
#> Path [39] :Initial log joint density = -481697.421788 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      9.759e-03   2.366e-01    1.000e+00  1.000e+00      4136 -3.705e+03 -3.706e+03                   
#> Path [39] :Best Iter: [64] ELBO (-3704.874202) evaluations: (4136) 
#> Path [40] :Initial log joint density = -481572.455630 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.787e+05      1.353e-02   2.262e-01    1.000e+00  1.000e+00      5621 -3.698e+03 -3.701e+03                   
#> Path [40] :Best Iter: [80] ELBO (-3697.756357) evaluations: (5621) 
#> Path [41] :Initial log joint density = -483422.909587 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      2.621e-03   2.886e-01    6.839e-01  6.839e-01      3642 -3.700e+03 -3.716e+03                   
#> Path [41] :Best Iter: [58] ELBO (-3699.673306) evaluations: (3642) 
#> Path [42] :Initial log joint density = -482084.515250 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.214e-02   2.252e-01    1.000e+00  1.000e+00      5069 -3.701e+03 -3.705e+03                   
#> Path [42] :Best Iter: [74] ELBO (-3701.097322) evaluations: (5069) 
#> Path [43] :Initial log joint density = -483670.624247 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      3.199e-03   2.555e-01    6.940e-01  6.940e-01      4594 -3.701e+03 -3.713e+03                   
#> Path [43] :Best Iter: [69] ELBO (-3700.557472) evaluations: (4594) 
#> Path [44] :Initial log joint density = -482028.198463 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      1.194e-02   2.269e-01    1.000e+00  1.000e+00      4806 -3.701e+03 -3.702e+03                   
#> Path [44] :Best Iter: [69] ELBO (-3701.425641) evaluations: (4806) 
#> Path [45] :Initial log joint density = -481398.450451 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      2.203e-02   3.010e-01    4.867e-01  1.000e+00      4570 -3.700e+03 -3.711e+03                   
#> Path [45] :Best Iter: [71] ELBO (-3700.040571) evaluations: (4570) 
#> Path [46] :Initial log joint density = -481837.490441 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      4.189e-02   2.188e-01    1.000e+00  1.000e+00      4170 -3.703e+03 -3.706e+03                   
#> Path [46] :Best Iter: [62] ELBO (-3703.148405) evaluations: (4170) 
#> Path [47] :Initial log joint density = -481972.297060 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.575e-02   1.524e-01    1.000e+00  1.000e+00      4341 -3.700e+03 -3.699e+03                   
#> Path [47] :Best Iter: [69] ELBO (-3699.310841) evaluations: (4341) 
#> Path [48] :Initial log joint density = -481644.126224 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      2.951e-03   2.342e-01    6.637e-01  6.637e-01      3513 -3.705e+03 -3.718e+03                   
#> Path [48] :Best Iter: [59] ELBO (-3705.389756) evaluations: (3513) 
#> Path [49] :Initial log joint density = -481658.894977 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.787e+05      1.097e-02   2.488e-01    1.000e+00  1.000e+00      5586 -3.701e+03 -3.709e+03                   
#> Path [49] :Best Iter: [80] ELBO (-3700.652592) evaluations: (5586) 
#> Path [50] :Initial log joint density = -482168.127638 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      1.666e-02   2.046e-01    1.000e+00  1.000e+00      5473 -3.699e+03 -3.698e+03                   
#> Path [50] :Best Iter: [79] ELBO (-3697.674844) evaluations: (5473) 
#> Finished in  15.4 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
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
#>    cell_group parameter   factor  c_lower c_effect c_upper c_rhat c_ess_bulk
#>    <chr>      <chr>       <chr>     <dbl>    <dbl>   <dbl>  <dbl>      <dbl>
#>  1 B1         (Intercept) NA      0.978      1.22   1.46    1.000      3808.
#>  2 B1         typecancer  type   -0.911     -0.588 -0.256   1.00       3749.
#>  3 B2         (Intercept) NA      0.513      0.786  1.06    1.00       3790.
#>  4 B2         typecancer  type   -0.986     -0.639 -0.299   1.000      3863.
#>  5 B3         (Intercept) NA     -0.588     -0.325 -0.0673  1.00       3756.
#>  6 B3         typecancer  type   -0.555     -0.250  0.0684  1.000      3614.
#>  7 BM         (Intercept) NA     -1.23      -0.973 -0.708   1.00       3844.
#>  8 BM         typecancer  type   -0.598     -0.259  0.0674  1.00       4007.
#>  9 CD4 1      (Intercept) NA      0.221      0.384  0.549   1.00       4000.
#> 10 CD4 1      typecancer  type   -0.00719    0.222  0.445   1.00       3989.
#> # ℹ 62 more rows
#> # ℹ 7 more variables: c_ess_tail <dbl>, v_upper <dbl>, v_effect <dbl>,
#> #   v_lower <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>
# }
```
