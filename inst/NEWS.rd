\name{NEWS}
\title{News for Package \pkg{sccomp}}

\section{News in version 2.1.21}{
\itemize{
    \item NOTE: This version was developed on a feature branch (allow-FDR-text-in-boxplot) with a commit date of July 28, 2025, but version numbering occurred after v2.1.19 and v2.1.20 were created (October 27, 2025).
    \item Enhanced sccomp_boxplot flexibility. Added tests to ensure sccomp_boxplot can accept additional ggplot layers, improving flexibility for users. Updated vignettes to reflect changes in function usage and added examples for custom ggplot layers.
    \item Commit ID: 3d4a9ef72b64ecc2db81347f8b34e16ce67f6477
}}

\section{News in version 2.1.20}{
\itemize{
    \item Major enhancements to sccomp_proportional_fold_change function. Added support for random effects in composition and variability models. Enhanced handling of complex interaction categories in from/to parameters, supporting both two-factor and three-factor interactions.
    \item Improved interaction model support. Function now dynamically handles sample column names for interaction models and ensures correct new_data structure creation across various model types (simple, two-factor, and interaction models).
    \item Added extensive unit tests covering various scenarios, including complex random effects and error handling for invalid inputs. Improved parsing of random effects formulas and new_data creation for robustness.
    \item Bug fix: Fixed contrasts_to_parameter_list function handling of variables containing digits (PR #234, #240).
    \item Between commits: 6359065f9872bd9bafe69297c418778cb329dba7..896b9e50120aab4cca5c08ac2d776552fc6bd495
}}

\section{News in version 2.1.19}{
\itemize{
    \item Added package help man page for sccomp. Improved package documentation with comprehensive help page.
    \item Function naming updates. Renamed function references from remove_unwanted_variation to remove_unwanted_effects for consistency across documentation.
    \item Improved package man documentation.
    \item Between commits: 4c5774cdb1932ae02e4c17fe0387dff1f07dbccc..6359065f9872bd9bafe69297c418778cb329dba7
}}

\section{News in version 2.1.18, Bioconductor 3.22 Release}{
\itemize{
    \item Removed deprecated .sample argument handling. Cleaned up sccomp_estimate function by removing deprecated argument handling and added tests for old argument style compatibility.
    \item Commit ID: c6f2afe (tag v2.1.18)
}}

\section{News in version 2.1.17}{
\itemize{
    \item Added cache_stan_model parameter to relevant functions with updated documentation for better model caching control.
    \item Fixed vb_iterative functionality for variational inference.
    \item Comprehensive documentation improvements. Updated README and vignettes with new figures, corrected figure paths to use relative paths, removed obsolete images, and enhanced documentation clarity.
    \item Refactored README and plotting functions to enhance clarity and consistency. Updated table formatting, improved plot caption handling, and refined boxplot aesthetics.
    \item Adjusted parameter handling in plots for better visualization and consistency.
    \item Enhanced documentation for better user experience with corrected image rendering logic in RMarkdown.
    \item Between commits: e9f64e6272e466d188e64b28929450b980394d41..49574b41ea46a607c93d35263eb9b496aa352f46
}}

\section{News in version 2.1.16}{
\itemize{
    \item Major code cleanup and refactoring. Removed unused and deprecated functions from codebase and tests.
    \item Refactored data handling in multiple functions to replace deprecated as_matrix with column_to_rownames and pivot_wider for improved data manipulation.
    \item Updated import statements and enhanced data selection logic in add_formula_columns. Improved conditional handling in get_specific_annotation_columns.
    \item Enhanced NA handling. Improved tests to validate handling of NA values in various scenarios.
    \item Cache directory improvements. Refactored cache directory handling with dedicated function for cache directory retrieval. Updated global variable initialization.
    \item Enhanced comparison logic in validation functions for clarity and consistency.
    \item Fixed documentation warnings.
    \item Between commits: c369dbe0909376e1257792c85f21fc025780f0cc..e9f64e6272e466d188e64b28929450b980394d41
}}

\section{News in version 2.1.15}{
\itemize{
    \item Enhanced visualization capabilities. Major improvements to plot_2D_intervals function including addition of corrected/adjusted intercept as a new facet. Updated data processing to include corrected variability effects and adjusted facet ordering.
    \item Added regression lines to plot_2D_intervals based on prec_coeff parameters. Added calculations for intercept and slope, and updated plotting logic to display the regression line with improved aesthetics.
    \item Improved color scales in plot_1D_intervals for significance thresholds. Enhanced traceplot aesthetics in the introduction vignette.
    \item Replaced 'corrected' terminology with 'adjusted' for intercept data throughout documentation and code. Added horizontal line for the adjusted intercept facet.
    \item Documentation updates. Updated README with logo and video link. Removed outdated intro.Rmd file. Adjusted references to the package in badges and descriptions for clarity.
    \item Updated README figures to reflect recent changes in analysis results and improved plot aesthetics. Adjusted parameter handling in plot_2D_intervals to exclude intercepts from significance coloring.
    \item Data handling refactoring. Replaced deprecated as_matrix with column_to_rownames and pivot_wider for improved data manipulation across multiple functions.
    \item Added tests to verify the inclusion of the regression line in the generated plots. Enhanced tests to validate handling of NA values in various scenarios.
    \item Between commits: ae4c9a67064c119253aaf176513d2a9dcd0bf163..493b7ebbb6d2855fdee68f9a9458937ee0b7e43c
}}

\section{News in version 2.1.14, Bioconductor 3.22 Release}{
\itemize{
    \item **Major improvement:** Transitioned to using built-in sum_to_zero_vector in Stan for improved efficiency and clarity. Dropped post-model correction calculations in favor of native Stan implementation.
    \item Removed deprecated QR-based sum-to-zero functions and replaced with new implementation. This update introduces a new way of handling certain statistical constraints in the underlying model, using a "sum-to-zero" variable type.
    \item Why was this change made? The previous method for enforcing sum-to-zero constraints (a requirement for compositional data) relied on QR decomposition, which could be numerically unstable and inefficient, especially for complex models or large datasets. The new method leverages recent advances in the Stan modeling language (requiring cmdstanr version 0.9.0 or higher).
    \item What does this mean for users? The model now converges faster and produces more reliable results, especially for analyses involving many cell types or complex experimental designs. The improvement in effective sample size (ESS) can be dramatic - across real-world datasets, the median improvement is in the range of 500--1000 additional effective samples per parameter, with many parameters seeing their ESS doubled or more.
    \item Refactored random effect handling in Stan models. Updated from matrices to arrays for improved flexibility. Adjusted array dimensions for random effects and their corresponding sigma parameters to ensure consistency. Applied proper sum-to-zero constraints.
    \item Updated cmdstanr version check to require version 0.9.0 or higher. Added error handling for cmdstanr version checks.
    \item Corrected beta_raw matrix dimensions and improved error handling during HMC inference. Adjusted verbosity for exception display and refined error reporting for model compilation issues.
    \item Updated priors and random effects in glm_multi_beta_binomial Stan model to use correct scaling for sum-to-zero constraints. Adjusted normal distributions for random effects and their corresponding sigma parameters.
    \item Added compare_branches_job.R script for parallel branch comparison to validate improvements.
    \item Enhanced NEWS.rd documentation with detailed explanations of improvements. See https://github.com/MangiolaLaboratory/sccomp/pull/211 for comparison plots and benchmarks.
    \item Between commits: 4c1b7541c18ad60806481e8d5d5fb97cac80fd9a..ae4c9a67064c119253aaf176513d2a9dcd0bf163
}}

\section{News in version 2.1.13}{
\itemize{
    \item Improved error handling in cmdstanr installation. Refactored error handling in check_and_install_cmdstanr function for better user experience.
    \item Commit ID: 4c1b7541c18ad60806481e8d5d5fb97cac80fd9a
}}

\section{News in version 2.1.12, Bioconductor 3.22 Release}{
\itemize{
    \item Completely removed the deprecated old framework (methods_OLD_framework), including all its functions and files.
    \item The function \code{sccomp_remove_unwanted_variation} is still present as it is a recent deprecation.
    \item Commit ID: f7df6ce1241bfc593b1a30a9559c8d8970835177
}}

\section{News in version 2.1.11}{
\itemize{
    \item Introduced create_multipanel_theme function for creating consistent multi-panel plot themes.
    \item Enhanced convergence metrics for random effects. Added convergence checking on random effects and comprehensive tests for convergence metrics in both fixed and random effects models.
    \item Improved plotting functionality. Enhanced plotting functions with significance_statistic parameter for flexible significance testing display. Added show_fdr_message parameter for better control over FDR messaging in visualizations.
    \item Documentation improvements. Updated introduction vignette with mixed linear models examples.
    \item Refactored tests for replicate_data to properly handle X_random_effect_unseen scenarios.
    \item Fixed tests and documentation for improved clarity and accuracy.
    \item Between commits: 4e21d30a21ef887d9a10b5d11904237951941370..a98d9ffb4c94b70149b63fe541679fb64d6ac3e9
}}

\section{News in version 2.1.10}{
\itemize{
    \item Enhanced sccomp_predict and summary_to_tibble functions to support robust statistics for more reliable predictions.
    \item Bug fixes. Fixed issues with unmodelled random effect combinations. Fixed design matrix with rownames for proper data handling.
    \item Function deprecation. Deprecated sccomp_remove_unwanted_variation function in favor of sccomp_remove_unwanted_effects for clearer terminology. Added comprehensive tests to ensure deprecation warning and result consistency between the two functions.
    \item Updated documentation with new Rd file for sccomp_remove_unwanted_effects and clarified usage examples.
    \item Removed redundant import of percent-dollar-percent from magrittr in NAMESPACE file for clarity and consistency.
    \item Solved CHECK issues for better package compliance.
    \item Between commits: f1005dde0020b0dce98dc00c13249ac0134770d3..4e21d30a21ef887d9a10b5d11904237951941370
}}

\section{News in version 2.1.0}{
\itemize{
    \item Bioconductor 3.21 Development Release. Version bumped for development following creation of RELEASE_3_21 branch.
    \item Commit ID: c8b900dc8d9498f9577514380d1d53f220b1f487
}}

\section{News in version 2.0.0}{
\itemize{
    \item Bioconductor 3.21 Release. Major version milestone for Bioconductor 3.21 release branch creation.
    \item Commit ID: a4668822705280dfc67d029e9b1faadfe504ca00
}}

\section{News in version 1.99.15}{
\itemize{
    \item Fixed parallel chains for pathfinder inference method. Improved pathfinder algorithm to properly support parallel chain execution for faster and more reliable inference.
    \item Commit ID: 505abc244e81bcc146284c3c464f15b0f273ac10
}}

\section{News in version 1.7.14}{
\itemize{
    \item Upgraded pathfinder with support for multiple chains. Enhanced the pathfinder variational inference method to allow for multiple chains, improving robustness and convergence assessment.
    \item Commit ID: db009538ee645febf33ff883c09b9a29daf5ab62
}}

\section{News in version 1.7.10}{
\itemize{
    \item **Major feature: Added pathfinder as an inference method.** Introduced pathfinder, a fast variational inference algorithm, as an alternative to HMC and variational Bayes. Pathfinder provides a good balance between speed and accuracy, offering faster inference than full HMC while maintaining better accuracy than traditional variational Bayes. Modified inference_method argument to support "hmc", "vb", and "pathfinder" options.
    \item Tweaked pathfinder algorithm parameters for improved performance and stability.
    \item Commit ID: 0018ea1759c29ef682dfc39d62158cdce57a0101 (initial), f77c9cd14fe69ff0cc8f62902904aedbe615f3c9 (tweaks)
}}

\section{News in version 1.7.4, Bioconductor 3.19 Release}{
\itemize{
    \item Single-cell transcriptomics allows the unbiased characterisation of the cellular composition of tissues. The cellular composition can be compared between biological or clinical conditions to identify potential cellular drivers. This strategy has been critical to unveil drivers of immune response in cancer and pathogen infection from single-cell data. Developing a robust statistical method for differential composition analyses from single-cell data is crucial for driving discoveries. The compositional data from single-cell experiments has four main properties. The data is in count form; counts underlie inversely correlated proportions that sum to one; larger cell groups are more variable across samples than small groups; real-world data is rich in outlier observation. A model that covers more than two of these properties is currently lacking. Here, we present a robust and outlier-aware method for testing differential tissue composition from single-cell data. This model can also transfer knowledge from a large set of integrated datasets to increase accuracy further. We present how this model can be applied to identify novel compositional and heterogeneity changes in existing studies.
    \item Commit ID: 1e0de5b229ddd9e95214a3097ff7650eb8bdfff4
}}

\section{News in version 1.0.0, Bioconductor 3.14 Release}{
\itemize{
    \item Blog post here: https://tidyomics.github.io/tidyomicsBlog/post/2023-12-07-tidy-sccomp/. Introduction of a new modular, tidy framework. This modularise estimation |> outlier exclusion |> hypothesis testing |> visualisation |> simulation(). New functionalities have also been added such as, removal of unwanted variation.
    \item Commit ID: dcc1b2624aaba990b882ed9a2c0e5142bc9aa037
}}


