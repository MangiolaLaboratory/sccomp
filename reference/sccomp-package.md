# sccomp: Differential Composition and Variability Analysis for Single-Cell Data

Comprehensive R package for differential composition and variability
analysis in single-cell RNA sequencing, CyTOF, and microbiome data.
Provides robust Bayesian modeling with outlier detection, random
effects, and advanced statistical methods for cell type proportion
analysis. Features include probabilistic outlier identification,
mixed-effect modeling, differential variability testing, and
comprehensive visualization tools. Perfect for cancer research,
immunology, developmental biology, and single-cell genomics
applications.

The sccomp package provides comprehensive tools for differential
composition and variability analysis in single-cell RNA sequencing,
CyTOF, and microbiome data. It implements robust Bayesian modeling with
outlier detection, random effects, and advanced statistical methods for
cell type proportion analysis.

## Details

The main functions are:

- [`sccomp_estimate`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_estimate.md) -
  Perform differential composition and variability analysis

- [`sccomp_test`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_test.md) -
  Test for differential composition with significance testing

- [`sccomp_remove_outliers`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_remove_outliers.md) -
  Identify and remove outlier samples

- [`sccomp_predict`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_predict.md) -
  Predict cell type proportions for new samples

- [`sccomp_remove_unwanted_variation`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_remove_unwanted_variation.md) -
  Remove unwanted variation from data

- [`sccomp_proportional_fold_change`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_proportional_fold_change.md) -
  Calculate proportional fold changes

- Plotting functions:
  [`plot.sccomp_tbl`](https://mangiolalaboratory.github.io/sccomp/reference/plot.sccomp_tbl.md),
  [`sccomp_boxplot`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_boxplot.md),
  [`plot_1D_intervals`](https://mangiolalaboratory.github.io/sccomp/reference/plot_1D_intervals.md),
  [`plot_2D_intervals`](https://mangiolalaboratory.github.io/sccomp/reference/plot_2D_intervals.md)

For detailed information on usage, see the package vignettes:

`browseVignettes("sccomp")`

All software-related questions should be posted to the GitHub Issues
page:

<https://github.com/MangiolaLaboratory/sccomp/issues>

The code can be viewed at the GitHub repository:

<https://github.com/MangiolaLaboratory/sccomp>

## References

Mangiola, S., Roth-Schulze, A.J., Trussart, M., Zozaya-Valdés, E., Ma,
M., Gao, Z., Rubin, A.F., Speed, T.P., Shim, H., & Papenfuss, A.T.
(2023). sccomp: Robust differential composition and variability analysis
for single-cell data. Proceedings of the National Academy of Sciences,
120(33), e2203828120.
[doi:10.1073/pnas.2203828120](https://doi.org/10.1073/pnas.2203828120)

## See also

Useful links:

- <https://github.com/MangiolaLaboratory/sccomp>

- Report bugs at <https://github.com/MangiolaLaboratory/sccomp/issues>

Useful links:

- <https://github.com/MangiolaLaboratory/sccomp>

- Report bugs at <https://github.com/MangiolaLaboratory/sccomp/issues>

## Author

**Maintainer**: Stefano Mangiola <stefano.mangiola@unimelb.edu.au>

Authors:

- Alexandra J. Roth-Schulze

- Marie Trussart

- Enrique Zozaya-Valdés

- Mengyao Ma

- Zijie Gao

- Alan F. Rubin

- Terence P. Speed

- Heejung Shim

- Anthony T. Papenfuss
