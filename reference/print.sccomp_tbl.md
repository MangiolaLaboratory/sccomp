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
#> Path [1] :Initial log joint density = -481519.132435 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      2.963e-03   2.543e-01    6.736e-01  6.736e-01      3566 -3.701e+03 -3.714e+03                   
#> Path [1] :Best Iter: [55] ELBO (-3700.551911) evaluations: (3566) 
#> Path [2] :Initial log joint density = -482505.926561 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.592e-02   2.819e-01    1.000e+00  1.000e+00      3411 -3.701e+03 -3.702e+03                   
#> Path [2] :Best Iter: [57] ELBO (-3701.021956) evaluations: (3411) 
#> Path [3] :Initial log joint density = -481539.444908 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      3.040e-02   3.682e-01    1.000e+00  1.000e+00      3420 -3.697e+03 -3.701e+03                   
#> Path [3] :Best Iter: [58] ELBO (-3696.555171) evaluations: (3420) 
#> Path [4] :Initial log joint density = -481955.974928 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.040e-03   3.419e-01    5.267e-01  1.000e+00      3160 -3.699e+03 -3.708e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3698.651697) evaluations: (3160) 
#> Path [5] :Initial log joint density = -481813.510514 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.285e-02   2.787e-01    1.000e+00  1.000e+00      3512 -3.700e+03 -3.704e+03                   
#> Path [5] :Best Iter: [59] ELBO (-3699.544436) evaluations: (3512) 
#> Path [6] :Initial log joint density = -481621.525693 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.287e-02   3.510e-01    1.000e+00  1.000e+00      2794 -3.706e+03 -3.714e+03                   
#> Path [6] :Best Iter: [50] ELBO (-3705.662382) evaluations: (2794) 
#> Path [7] :Initial log joint density = -481594.090403 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      2.338e-02   1.716e-01    1.000e+00  1.000e+00      3689 -3.703e+03 -3.700e+03                   
#> Path [7] :Best Iter: [62] ELBO (-3699.564861) evaluations: (3689) 
#> Path [8] :Initial log joint density = -482356.448871 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.576e-02   4.088e-01    1.000e+00  1.000e+00      3302 -3.700e+03 -3.711e+03                   
#> Path [8] :Best Iter: [55] ELBO (-3700.357055) evaluations: (3302) 
#> Path [9] :Initial log joint density = -481982.345651 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.055e-02   3.255e-01    9.356e-01  9.356e-01      2997 -3.710e+03 -3.718e+03                   
#> Path [9] :Best Iter: [52] ELBO (-3710.413108) evaluations: (2997) 
#> Path [10] :Initial log joint density = -481520.921594 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      8.646e-03   1.893e-01    1.000e+00  1.000e+00      3506 -3.701e+03 -3.703e+03                   
#> Path [10] :Best Iter: [58] ELBO (-3701.449643) evaluations: (3506) 
#> Path [11] :Initial log joint density = -481541.580387 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.021e-03   1.925e-01    7.618e-01  7.618e-01      3077 -3.707e+03 -3.712e+03                   
#> Path [11] :Best Iter: [48] ELBO (-3707.475593) evaluations: (3077) 
#> Path [12] :Initial log joint density = -481528.633352 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.919e-03   2.616e-01    9.263e-01  9.263e-01      3248 -3.700e+03 -3.713e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3700.353286) evaluations: (3248) 
#> Path [13] :Initial log joint density = -481560.234891 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.512e-03   1.832e-01    9.357e-01  9.357e-01      3039 -3.708e+03 -3.723e+03                   
#> Path [13] :Best Iter: [52] ELBO (-3708.086687) evaluations: (3039) 
#> Path [14] :Initial log joint density = -481499.826279 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.942e-02   3.370e-01    1.000e+00  1.000e+00      3300 -3.699e+03 -3.705e+03                   
#> Path [14] :Best Iter: [56] ELBO (-3698.705898) evaluations: (3300) 
#> Path [15] :Initial log joint density = -481377.436067 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.131e-02   2.113e-01    1.000e+00  1.000e+00      3794 -3.702e+03 -3.705e+03                   
#> Path [15] :Best Iter: [62] ELBO (-3702.407726) evaluations: (3794) 
#> Path [16] :Initial log joint density = -481533.485779 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.536e-02   3.338e-01    1.000e+00  1.000e+00      2950 -3.708e+03 -3.707e+03                   
#> Path [16] :Best Iter: [53] ELBO (-3706.820636) evaluations: (2950) 
#> Path [17] :Initial log joint density = -481680.716203 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      9.163e-03   1.798e-01    1.000e+00  1.000e+00      3529 -3.700e+03 -3.702e+03                   
#> Path [17] :Best Iter: [59] ELBO (-3700.228285) evaluations: (3529) 
#> Path [18] :Initial log joint density = -481522.277112 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.002e-03   2.176e-01    1.000e+00  1.000e+00      3154 -3.704e+03 -3.701e+03                   
#> Path [18] :Best Iter: [56] ELBO (-3700.806459) evaluations: (3154) 
#> Path [19] :Initial log joint density = -482960.472689 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.050e-03   2.572e-01    1.000e+00  1.000e+00      3085 -3.709e+03 -3.721e+03                   
#> Path [19] :Best Iter: [46] ELBO (-3708.699345) evaluations: (3085) 
#> Path [20] :Initial log joint density = -481884.634542 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.154e-02   2.570e-01    1.000e+00  1.000e+00      3356 -3.705e+03 -3.710e+03                   
#> Path [20] :Best Iter: [57] ELBO (-3704.778768) evaluations: (3356) 
#> Path [21] :Initial log joint density = -481762.022200 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      7.395e-03   2.114e-01    1.000e+00  1.000e+00      2972 -3.708e+03 -3.715e+03                   
#> Path [21] :Best Iter: [43] ELBO (-3708.056675) evaluations: (2972) 
#> Path [22] :Initial log joint density = -481781.691048 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.065e-02   1.943e-01    1.000e+00  1.000e+00      3163 -3.707e+03 -3.703e+03                   
#> Path [22] :Best Iter: [56] ELBO (-3702.698073) evaluations: (3163) 
#> Path [23] :Initial log joint density = -481600.717290 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      4.222e-03   2.074e-01    7.778e-01  7.778e-01      3321 -3.701e+03 -3.713e+03                   
#> Path [23] :Best Iter: [57] ELBO (-3701.281674) evaluations: (3321) 
#> Path [24] :Initial log joint density = -481729.964525 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.280e-02   3.900e-01    9.986e-01  9.986e-01      2713 -3.707e+03 -3.722e+03                   
#> Path [24] :Best Iter: [49] ELBO (-3707.272483) evaluations: (2713) 
#> Path [25] :Initial log joint density = -481738.427441 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.138e-02   1.853e-01    1.000e+00  1.000e+00      3299 -3.704e+03 -3.702e+03                   
#> Path [25] :Best Iter: [57] ELBO (-3701.761491) evaluations: (3299) 
#> Path [26] :Initial log joint density = -481759.287770 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      6.094e-03   2.148e-01    7.878e-01  7.878e-01      3555 -3.697e+03 -3.709e+03                   
#> Path [26] :Best Iter: [59] ELBO (-3696.594408) evaluations: (3555) 
#> Path [27] :Initial log joint density = -483429.718767 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.017e-02   1.968e-01    1.000e+00  1.000e+00      3439 -3.703e+03 -3.700e+03                   
#> Path [27] :Best Iter: [58] ELBO (-3699.716308) evaluations: (3439) 
#> Path [28] :Initial log joint density = -481615.704268 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.302e-02   2.694e-01    1.000e+00  1.000e+00      3305 -3.701e+03 -3.702e+03                   
#> Path [28] :Best Iter: [58] ELBO (-3700.761654) evaluations: (3305) 
#> Path [29] :Initial log joint density = -482146.180662 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      4.583e-03   1.912e-01    4.776e-01  1.000e+00      3705 -3.701e+03 -3.712e+03                   
#> Path [29] :Best Iter: [60] ELBO (-3700.565735) evaluations: (3705) 
#> Path [30] :Initial log joint density = -481654.572753 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.777e-03   2.060e-01    1.000e+00  1.000e+00      3305 -3.701e+03 -3.702e+03                   
#> Path [30] :Best Iter: [55] ELBO (-3700.895979) evaluations: (3305) 
#> Path [31] :Initial log joint density = -481933.407499 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.319e-02   2.637e-01    1.000e+00  1.000e+00      3688 -3.699e+03 -3.700e+03                   
#> Path [31] :Best Iter: [59] ELBO (-3699.036890) evaluations: (3688) 
#> Path [32] :Initial log joint density = -481467.003508 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.065e-02   2.485e-01    1.000e+00  1.000e+00      3076 -3.708e+03 -3.707e+03                   
#> Path [32] :Best Iter: [54] ELBO (-3707.004209) evaluations: (3076) 
#> Path [33] :Initial log joint density = -481621.916356 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.223e-02   1.762e-01    1.000e+00  1.000e+00      3499 -3.698e+03 -3.697e+03                   
#> Path [33] :Best Iter: [60] ELBO (-3697.185929) evaluations: (3499) 
#> Path [34] :Initial log joint density = -482212.964546 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.395e-02   3.045e-01    9.768e-01  9.768e-01      3363 -3.699e+03 -3.709e+03                   
#> Path [34] :Best Iter: [57] ELBO (-3698.549092) evaluations: (3363) 
#> Path [35] :Initial log joint density = -483129.929131 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      4.541e-03   2.334e-01    6.694e-01  6.694e-01      3644 -3.697e+03 -3.714e+03                   
#> Path [35] :Best Iter: [58] ELBO (-3696.890837) evaluations: (3644) 
#> Path [36] :Initial log joint density = -481774.389743 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.169e-03   2.509e-01    1.000e+00  1.000e+00      3262 -3.709e+03 -3.704e+03                   
#> Path [36] :Best Iter: [55] ELBO (-3704.194206) evaluations: (3262) 
#> Path [37] :Initial log joint density = -481379.886200 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.755e-02   3.969e-01    9.917e-01  9.917e-01      3109 -3.701e+03 -3.709e+03                   
#> Path [37] :Best Iter: [55] ELBO (-3700.763283) evaluations: (3109) 
#> Path [38] :Initial log joint density = -483802.891988 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.527e-03   2.239e-01    1.000e+00  1.000e+00      3381 -3.699e+03 -3.703e+03                   
#> Path [38] :Best Iter: [57] ELBO (-3698.689275) evaluations: (3381) 
#> Path [39] :Initial log joint density = -481601.015853 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.779e-02   3.222e-01    1.000e+00  1.000e+00      3188 -3.712e+03 -3.705e+03                   
#> Path [39] :Best Iter: [55] ELBO (-3705.046014) evaluations: (3188) 
#> Path [40] :Initial log joint density = -481405.640036 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      9.103e-03   2.326e-01    9.653e-01  9.653e-01      2783 -3.707e+03 -3.716e+03                   
#> Path [40] :Best Iter: [51] ELBO (-3707.484713) evaluations: (2783) 
#> Path [41] :Initial log joint density = -481629.384881 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.884e-02   2.552e-01    9.168e-01  9.168e-01      2990 -3.708e+03 -3.719e+03                   
#> Path [41] :Best Iter: [51] ELBO (-3707.782735) evaluations: (2990) 
#> Path [42] :Initial log joint density = -481495.803019 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      9.727e-03   1.689e-01    1.000e+00  1.000e+00      2673 -3.710e+03 -3.709e+03                   
#> Path [42] :Best Iter: [50] ELBO (-3709.336967) evaluations: (2673) 
#> Path [43] :Initial log joint density = -481683.905872 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.874e-03   1.578e-01    8.247e-01  8.247e-01      3255 -3.705e+03 -3.712e+03                   
#> Path [43] :Best Iter: [55] ELBO (-3705.012040) evaluations: (3255) 
#> Path [44] :Initial log joint density = -481449.217773 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.951e-03   2.225e-01    1.000e+00  1.000e+00      2877 -3.706e+03 -3.710e+03                   
#> Path [44] :Best Iter: [49] ELBO (-3705.849003) evaluations: (2877) 
#> Path [45] :Initial log joint density = -481740.578828 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      2.586e-03   3.710e-01    5.638e-01  5.638e-01      3117 -3.709e+03 -3.718e+03                   
#> Path [45] :Best Iter: [52] ELBO (-3708.593670) evaluations: (3117) 
#> Path [46] :Initial log joint density = -481787.181334 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.469e-02   3.082e-01    1.000e+00  1.000e+00      3312 -3.705e+03 -3.703e+03                   
#> Path [46] :Best Iter: [58] ELBO (-3702.776413) evaluations: (3312) 
#> Path [47] :Initial log joint density = -482721.260156 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.677e-03   2.574e-01    1.000e+00  1.000e+00      3327 -3.700e+03 -3.702e+03                   
#> Path [47] :Best Iter: [57] ELBO (-3699.840746) evaluations: (3327) 
#> Path [48] :Initial log joint density = -482042.322279 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.006e-02   3.246e-01    1.000e+00  1.000e+00      3796 -3.700e+03 -3.708e+03                   
#> Path [48] :Best Iter: [59] ELBO (-3700.494030) evaluations: (3796) 
#> Path [49] :Initial log joint density = -481977.261111 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.203e-02   3.663e-01    1.000e+00  1.000e+00      3220 -3.700e+03 -3.708e+03                   
#> Path [49] :Best Iter: [56] ELBO (-3700.226490) evaluations: (3220) 
#> Path [50] :Initial log joint density = -481379.579678 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      4.954e-03   1.247e-01    1.000e+00  1.000e+00      3016 -3.708e+03 -3.717e+03                   
#> Path [50] :Best Iter: [41] ELBO (-3708.086923) evaluations: (3016) 
#> Finished in  13.6 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
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
#>    cell_group parameter   factor c_lower c_effect c_upper c_rhat c_ess_bulk
#>    <chr>      <chr>       <chr>    <dbl>    <dbl>   <dbl>  <dbl>      <dbl>
#>  1 B1         (Intercept) NA      0.941     1.20   1.45    1.000      3545.
#>  2 B1         typecancer  type   -0.925    -0.614 -0.310   1.00       3933.
#>  3 B2         (Intercept) NA      0.506     0.768  1.02    1.000      3645.
#>  4 B2         typecancer  type   -0.975    -0.672 -0.357   1.00       3866.
#>  5 B3         (Intercept) NA     -0.580    -0.328 -0.0810  1.000      3957.
#>  6 B3         typecancer  type   -0.590    -0.276  0.0330  1.00       4077.
#>  7 BM         (Intercept) NA     -1.23     -0.964 -0.723   1.000      3867.
#>  8 BM         typecancer  type   -0.581    -0.283  0.0153  1.00       4071.
#>  9 CD4 1      (Intercept) NA      0.206     0.369  0.530   1.00       3789.
#> 10 CD4 1      typecancer  type   -0.0113    0.205  0.422   1.00       3771.
#> # ℹ 62 more rows
#> # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
#> #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>
# }
```
