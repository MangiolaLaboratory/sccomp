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
#> Path [1] :Initial log joint density = -481536.076859 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.931e-03   2.066e-01    1.000e+00  1.000e+00      3119 -3.691e+03 -3.686e+03                   
#> Path [1] :Best Iter: [55] ELBO (-3686.453735) evaluations: (3119) 
#> Path [2] :Initial log joint density = -481503.280653 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.366e-03   2.039e-01    1.000e+00  1.000e+00      3259 -3.690e+03 -3.692e+03                   
#> Path [2] :Best Iter: [43] ELBO (-3690.115423) evaluations: (3259) 
#> Path [3] :Initial log joint density = -481737.828469 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.591e-02   3.788e-01    1.000e+00  1.000e+00      3414 -3.684e+03 -3.686e+03                   
#> Path [3] :Best Iter: [58] ELBO (-3684.335552) evaluations: (3414) 
#> Path [4] :Initial log joint density = -482359.203672 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.469e-03   2.074e-01    8.067e-01  8.067e-01      3363 -3.683e+03 -3.698e+03                   
#> Path [4] :Best Iter: [57] ELBO (-3683.164345) evaluations: (3363) 
#> Path [5] :Initial log joint density = -481529.481534 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.737e-02   2.765e-01    1.000e+00  1.000e+00      3451 -3.678e+03 -3.683e+03                   
#> Path [5] :Best Iter: [59] ELBO (-3678.409514) evaluations: (3451) 
#> Path [6] :Initial log joint density = -481616.271056 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.990e-02   2.354e-01    1.000e+00  1.000e+00      2879 -3.692e+03 -3.696e+03                   
#> Path [6] :Best Iter: [44] ELBO (-3691.566449) evaluations: (2879) 
#> Path [7] :Initial log joint density = -481678.388610 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      7.334e-03   2.146e-01    9.188e-01  9.188e-01      2826 -3.692e+03 -3.704e+03                   
#> Path [7] :Best Iter: [38] ELBO (-3692.425802) evaluations: (2826) 
#> Path [8] :Initial log joint density = -481630.073491 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.049e-02   2.906e-01    1.000e+00  1.000e+00      2830 -3.687e+03 -3.701e+03                   
#> Path [8] :Best Iter: [42] ELBO (-3687.358694) evaluations: (2830) 
#> Path [9] :Initial log joint density = -483532.928508 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.381e-03   2.505e-01    1.000e+00  1.000e+00      3073 -3.689e+03 -3.694e+03                   
#> Path [9] :Best Iter: [52] ELBO (-3689.314853) evaluations: (3073) 
#> Path [10] :Initial log joint density = -481533.356727 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.583e-02   2.796e-01    4.687e-01  1.000e+00      3279 -3.680e+03 -3.690e+03                   
#> Path [10] :Best Iter: [57] ELBO (-3679.874259) evaluations: (3279) 
#> Path [11] :Initial log joint density = -481339.036616 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.481e-02   2.435e-01    1.000e+00  1.000e+00      3470 -3.680e+03 -3.682e+03                   
#> Path [11] :Best Iter: [58] ELBO (-3679.777948) evaluations: (3470) 
#> Path [12] :Initial log joint density = -481520.834698 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      5.664e-03   2.581e-01    1.000e+00  1.000e+00      2802 -3.691e+03 -3.695e+03                   
#> Path [12] :Best Iter: [43] ELBO (-3691.142234) evaluations: (2802) 
#> Path [13] :Initial log joint density = -483045.058615 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.181e-02   3.255e-01    1.000e+00  1.000e+00      3740 -3.681e+03 -3.693e+03                   
#> Path [13] :Best Iter: [60] ELBO (-3681.045358) evaluations: (3740) 
#> Path [14] :Initial log joint density = -481717.648784 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      2.785e-03   3.063e-01    5.326e-01  5.326e-01      3478 -3.681e+03 -3.697e+03                   
#> Path [14] :Best Iter: [58] ELBO (-3681.178917) evaluations: (3478) 
#> Path [15] :Initial log joint density = -481442.402326 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.402e-03   1.551e-01    7.999e-01  7.999e-01      3109 -3.691e+03 -3.696e+03                   
#> Path [15] :Best Iter: [55] ELBO (-3690.602283) evaluations: (3109) 
#> Path [16] :Initial log joint density = -481705.132991 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.095e-02   2.013e-01    1.000e+00  1.000e+00      2783 -3.687e+03 -3.692e+03                   
#> Path [16] :Best Iter: [50] ELBO (-3686.528181) evaluations: (2783) 
#> Path [17] :Initial log joint density = -481911.675929 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      4.626e-03   2.329e-01    8.283e-01  8.283e-01      2875 -3.687e+03 -3.701e+03                   
#> Path [17] :Best Iter: [49] ELBO (-3686.988614) evaluations: (2875) 
#> Path [18] :Initial log joint density = -481646.360795 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.347e-02   2.739e-01    1.000e+00  1.000e+00      3363 -3.685e+03 -3.684e+03                   
#> Path [18] :Best Iter: [59] ELBO (-3683.764982) evaluations: (3363) 
#> Path [19] :Initial log joint density = -481574.468170 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      9.751e-03   2.247e-01    1.000e+00  1.000e+00      2818 -3.693e+03 -3.703e+03                   
#> Path [19] :Best Iter: [40] ELBO (-3693.230235) evaluations: (2818) 
#> Path [20] :Initial log joint density = -488437.183249 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.946e-03   2.265e-01    1.000e+00  1.000e+00      3173 -3.689e+03 -3.705e+03                   
#> Path [20] :Best Iter: [51] ELBO (-3689.499928) evaluations: (3173) 
#> Path [21] :Initial log joint density = -481758.669848 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.012e-03   2.517e-01    1.000e+00  1.000e+00      3249 -3.687e+03 -3.686e+03                   
#> Path [21] :Best Iter: [58] ELBO (-3686.156822) evaluations: (3249) 
#> Path [22] :Initial log joint density = -481777.887932 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      4.542e-03   1.797e-01    1.000e+00  1.000e+00      3486 -3.686e+03 -3.689e+03                   
#> Path [22] :Best Iter: [57] ELBO (-3686.336666) evaluations: (3486) 
#> Path [23] :Initial log joint density = -481561.072315 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.009e-02   2.735e-01    1.000e+00  1.000e+00      3254 -3.688e+03 -3.685e+03                   
#> Path [23] :Best Iter: [56] ELBO (-3685.466762) evaluations: (3254) 
#> Path [24] :Initial log joint density = -481526.912460 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      5.948e-03   2.784e-01    4.171e-01  1.000e+00      3270 -3.686e+03 -3.696e+03                   
#> Path [24] :Best Iter: [56] ELBO (-3685.785458) evaluations: (3270) 
#> Path [25] :Initial log joint density = -481518.176608 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.207e-02   2.328e-01    8.976e-01  8.976e-01      3074 -3.691e+03 -3.702e+03                   
#> Path [25] :Best Iter: [47] ELBO (-3690.650565) evaluations: (3074) 
#> Path [26] :Initial log joint density = -481598.849530 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.936e-02   3.386e-01    1.000e+00  1.000e+00      3603 -3.682e+03 -3.691e+03                   
#> Path [26] :Best Iter: [60] ELBO (-3682.325293) evaluations: (3603) 
#> Path [27] :Initial log joint density = -481594.376726 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.442e-02   3.667e-01    9.828e-01  9.828e-01      3159 -3.690e+03 -3.692e+03                   
#> Path [27] :Best Iter: [54] ELBO (-3689.806345) evaluations: (3159) 
#> Path [28] :Initial log joint density = -485276.698314 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.102e-02   2.466e-01    1.000e+00  1.000e+00      3510 -3.681e+03 -3.679e+03                   
#> Path [28] :Best Iter: [60] ELBO (-3678.773464) evaluations: (3510) 
#> Path [29] :Initial log joint density = -481222.196437 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      4.935e-03   1.780e-01    8.784e-01  8.784e-01      2863 -3.689e+03 -3.705e+03                   
#> Path [29] :Best Iter: [50] ELBO (-3689.038992) evaluations: (2863) 
#> Path [30] :Initial log joint density = -484893.162484 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      4.659e-03   1.682e-01    1.000e+00  1.000e+00      3433 -3.686e+03 -3.694e+03                   
#> Path [30] :Best Iter: [56] ELBO (-3686.120999) evaluations: (3433) 
#> Path [31] :Initial log joint density = -481556.751159 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.998e-03   1.785e-01    8.630e-01  8.630e-01      3200 -3.689e+03 -3.693e+03                   
#> Path [31] :Best Iter: [55] ELBO (-3689.047124) evaluations: (3200) 
#> Path [32] :Initial log joint density = -481404.320052 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.768e-02   2.577e-01    1.000e+00  1.000e+00      3038 -3.691e+03 -3.695e+03                   
#> Path [32] :Best Iter: [45] ELBO (-3691.237634) evaluations: (3038) 
#> Path [33] :Initial log joint density = -481807.494190 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.603e-03   1.990e-01    1.000e+00  1.000e+00      3320 -3.685e+03 -3.694e+03                   
#> Path [33] :Best Iter: [55] ELBO (-3684.613471) evaluations: (3320) 
#> Path [34] :Initial log joint density = -481703.958627 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.665e-03   2.664e-01    1.000e+00  1.000e+00      3220 -3.682e+03 -3.690e+03                   
#> Path [34] :Best Iter: [55] ELBO (-3682.168163) evaluations: (3220) 
#> Path [35] :Initial log joint density = -481634.008623 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.954e-03   2.355e-01    1.000e+00  1.000e+00      3089 -3.693e+03 -3.701e+03                   
#> Path [35] :Best Iter: [41] ELBO (-3693.017855) evaluations: (3089) 
#> Path [36] :Initial log joint density = -481564.345620 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.135e-02   2.005e-01    1.000e+00  1.000e+00      3391 -3.685e+03 -3.688e+03                   
#> Path [36] :Best Iter: [55] ELBO (-3684.749865) evaluations: (3391) 
#> Path [37] :Initial log joint density = -481851.035872 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      4.954e-03   2.248e-01    9.200e-01  9.200e-01      2913 -3.689e+03 -3.703e+03                   
#> Path [37] :Best Iter: [50] ELBO (-3688.783510) evaluations: (2913) 
#> Path [38] :Initial log joint density = -481491.315933 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      3.404e-03   2.639e-01    7.058e-01  7.058e-01      2995 -3.689e+03 -3.712e+03                   
#> Path [38] :Best Iter: [51] ELBO (-3688.590139) evaluations: (2995) 
#> Path [39] :Initial log joint density = -481917.035347 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.893e-02   3.815e-01    1.000e+00  1.000e+00      3793 -3.680e+03 -3.691e+03                   
#> Path [39] :Best Iter: [62] ELBO (-3680.479928) evaluations: (3793) 
#> Path [40] :Initial log joint density = -482552.321197 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.303e-03   2.669e-01    1.000e+00  1.000e+00      3209 -3.692e+03 -3.692e+03                   
#> Path [40] :Best Iter: [44] ELBO (-3691.524760) evaluations: (3209) 
#> Path [41] :Initial log joint density = -482487.482142 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.852e-03   2.463e-01    7.402e-01  7.402e-01      3453 -3.683e+03 -3.701e+03                   
#> Path [41] :Best Iter: [57] ELBO (-3683.375343) evaluations: (3453) 
#> Path [42] :Initial log joint density = -481762.925986 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      5.049e-03   2.534e-01    7.926e-01  7.926e-01      2991 -3.690e+03 -3.704e+03                   
#> Path [42] :Best Iter: [38] ELBO (-3690.292988) evaluations: (2991) 
#> Path [43] :Initial log joint density = -482062.894667 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.404e-03   2.253e-01    1.000e+00  1.000e+00      3325 -3.684e+03 -3.686e+03                   
#> Path [43] :Best Iter: [55] ELBO (-3684.298786) evaluations: (3325) 
#> Path [44] :Initial log joint density = -481450.556836 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      4.988e-03   1.713e-01    7.968e-01  7.968e-01      2990 -3.690e+03 -3.703e+03                   
#> Path [44] :Best Iter: [47] ELBO (-3690.405973) evaluations: (2990) 
#> Path [45] :Initial log joint density = -481656.479170 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.560e-03   2.286e-01    1.000e+00  1.000e+00      3364 -3.685e+03 -3.685e+03                   
#> Path [45] :Best Iter: [58] ELBO (-3684.886361) evaluations: (3364) 
#> Path [46] :Initial log joint density = -482454.832977 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.051e-03   1.903e-01    9.735e-01  9.735e-01      3343 -3.685e+03 -3.690e+03                   
#> Path [46] :Best Iter: [56] ELBO (-3684.972859) evaluations: (3343) 
#> Path [47] :Initial log joint density = -481687.100378 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.108e-02   1.673e-01    1.000e+00  1.000e+00      2883 -3.690e+03 -3.691e+03                   
#> Path [47] :Best Iter: [38] ELBO (-3690.479384) evaluations: (2883) 
#> Path [48] :Initial log joint density = -481708.402264 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.786e-03   3.231e-01    9.350e-01  9.350e-01      2832 -3.691e+03 -3.704e+03                   
#> Path [48] :Best Iter: [51] ELBO (-3690.690864) evaluations: (2832) 
#> Path [49] :Initial log joint density = -481512.736458 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.767e-03   2.532e-01    1.000e+00  1.000e+00      3128 -3.690e+03 -3.688e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3688.046233) evaluations: (3128) 
#> Path [50] :Initial log joint density = -481938.378064 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.041e-03   3.112e-01    8.465e-01  8.465e-01      3263 -3.683e+03 -3.697e+03                   
#> Path [50] :Best Iter: [55] ELBO (-3682.501777) evaluations: (3263) 
#> Finished in  13.3 seconds.
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
#>    cell_group parameter   factor  c_lower c_effect c_upper c_rhat c_ess_bulk
#>    <chr>      <chr>       <chr>     <dbl>    <dbl>   <dbl>  <dbl>      <dbl>
#>  1 B1         (Intercept) NA      0.956      1.21   1.46    1.000      4193.
#>  2 B1         typecancer  type   -0.936     -0.616 -0.307   1.00       3649.
#>  3 B2         (Intercept) NA      0.507      0.775  1.03    1.00       3923.
#>  4 B2         typecancer  type   -0.961     -0.665 -0.362   1.00       3644.
#>  5 B3         (Intercept) NA     -0.585     -0.326 -0.0743  1.00       3556.
#>  6 B3         typecancer  type   -0.578     -0.277  0.0308  1.00       4158.
#>  7 BM         (Intercept) NA     -1.22      -0.968 -0.702   1.00       4056.
#>  8 BM         typecancer  type   -0.594     -0.282  0.0126  1.00       4156.
#>  9 CD4 1      (Intercept) NA      0.201      0.370  0.534   1.000      3814.
#> 10 CD4 1      typecancer  type   -0.00472    0.209  0.425   1.00       4222.
#> # ℹ 62 more rows
#> # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
#> #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>
# }
```
