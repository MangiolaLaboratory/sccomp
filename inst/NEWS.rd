\name{NEWS}
\title{News for Package \pkg{sccomp}}

\section{News in version 2.1.12}{
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


