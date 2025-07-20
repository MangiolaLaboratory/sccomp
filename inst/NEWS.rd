\name{NEWS}
\title{News for Package \pkg{sccomp}}

\section{News in version 2.1.18, Bioconductor 3.22 Release}{
\itemize{
    \item {Fixed issue with sccomp_proportional_fold_change function for interaction models.} The function now properly handles complex models with interactions by correctly creating the new_data structure that sccomp_predict expects. This resolves GitHub issue #193 where the function would fail with the error "your new_data might be malformed" when working with interaction models. The fix ensures that the sample column name is dynamically retrieved from the model attributes rather than hardcoded, making the function more robust and flexible.
    \item {Enhanced interaction category support in sccomp_proportional_fold_change.} The function now properly handles interaction categories in the from/to parameters, supporting both two-factor interactions (e.g., "treatment:followup") and three-factor interactions (e.g., "treatment:followup:B"). The function automatically parses interaction strings and creates the correct new_data structure with individual factor columns.
    \item {Added comprehensive unit tests for sccomp_proportional_fold_change.} A new test file (test-sccomp_proportional_fold_change.R) has been added with 20 different test scenarios covering simple models, interaction models, error handling, and edge cases. The tests progress from simple single-factor models to complex three-factor interaction models, ensuring the function works correctly across all use cases. New tests specifically cover interaction categories in from/to parameters.
}}

\section{News in version 2.1.14, Bioconductor 3.22 Release}{
\itemize{
    \item {Improved model efficiency and reliability with sum-to-zero variable constraints.} This update introduces a new way of handling certain statistical constraints in the underlying model, using a "sum-to-zero" variable type. Previously, these constraints were managed with a workaround that could be less efficient and less stable. The new approach is more mathematically direct and robust, leading to a substantial improvement in the model's ability to draw independent samples from the data (known as the "effective sample size" or ESS).
    \item {Why was this change made?} The previous method for enforcing sum-to-zero constraints (a requirement for compositional data) relied on QR decomposition, which could be numerically unstable and inefficient, especially for complex models or large datasets. The new method leverages recent advances in the Stan modeling language, allowing these constraints to be handled natively and more accurately.
    \item {What does this mean for users?} In practical terms, this change means that the model now converges faster and produces more reliable results, especially for analyses involving many cell types or complex experimental designs. The improvement in effective sample size (ESS) can be seen in the comparison plot: for the same amount of computation, the new method yields more independent samples, making statistical estimates more trustworthy and reducing the risk of misleading results (see https://github.com/MangiolaLaboratory/sccomp/pull/211).
    \item {Summary of improvement:} Across a wide range of real-world datasets, the new sum-to-zero variable approach consistently increases the effective sample size, sometimes dramatically, compared to the previous method. Visual inspection of the comparison plot shows that the median improvement in ESS is in the range of 500--1000 additional effective samples per parameter, with many parameters seeing their ESS doubled or more. This makes sccomp analyses faster, more robust, and more reproducible.
}}

\section{News in version 2.1.12, Bioconductor 3.22 Release}{
\itemize{
    \item Completely removed the deprecated old framework (methods_OLD_framework), including all its functions and files.
    \item The function \code{sccomp_remove_unwanted_variation} is still present as it is a recent deprecation.
}}

\section{News in version 1.7.4, Bioconductor 3.19 Release}{
\itemize{
    \item Single-cell transcriptomics allows the unbiased characterisation of the cellular composition of tissues. The cellular composition can be compared between biological or clinical conditions to identify potential cellular drivers. This strategy has been critical to unveil drivers of immune response in cancer and pathogen infection from single-cell data. Developing a robust statistical method for differential composition analyses from single-cell data is crucial for driving discoveries. The compositional data from single-cell experiments has four main properties. The data is in count form; counts underlie inversely correlated proportions that sum to one; larger cell groups are more variable across samples than small groups; real-world data is rich in outlier observation. A model that covers more than two of these properties is currently lacking. Here, we present a robust and outlier-aware method for testing differential tissue composition from single-cell data. This model can also transfer knowledge from a large set of integrated datasets to increase accuracy further. We present how this model can be applied to identify novel compositional and heterogeneity changes in existing studies.
}}

\section{News in version 1.0.0, Bioconductor 3.14 Release}{
\itemize{
    \item Blog post here: https://tidyomics.github.io/tidyomicsBlog/post/2023-12-07-tidy-sccomp/. Introduction of a new modular, tidy framework. This modularise estimation |> outlier exclusion |> hypothesis testing |> visualisation |> simulation(). New functionalities have also been added such as, removal of unwanted variation.
}}


