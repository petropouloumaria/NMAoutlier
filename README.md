NMAoutlier: Forward Search Algorithm in Network Meta-Analysis to identify outlying and influential studies
================

<img src="NMAoutlier-Logo.png" width=300 align="right" style="margin-left:20px; margin-right: 20px;"/>

Description
-----------

A package that provides forward search algorithm for detecting outlying or influential studies in network meta-analysis.

-   Provides the length of the initial oultying-free clean subset for forward search algorithm.
-   Iterations of forward search algorithm.
-   Basic set of studies in each iteration of forward search algorithm.
-   Summary estimates and their confidence intervals in each iteration of forward search algorithm.
-   Outlier and influential case diagnostics measures.
-   Ranking measures.
-   Heterogeneity and inconsistency measures.
-   Forward plots for summary estimates and their confidence intervals.
-   Forward plots for monitored measures: outlier and influential case diagnostics measures, ranking measures, heterogeneity and inconsistency measures.

Installation
------------

You can install the **NMAoutlier** package from GitHub repository as follows:

Installation using R package **[ghit](https://cran.r-project.org/package=ghit)** (without [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for Windows):

``` r
install.packages("ghit")
ghit::install_github("petropouloumaria/NMAoutlier")
```

Installation using R package **[devtools](https://cran.r-project.org/package=devtools)** (with [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for Windows):

``` r
install.packages("devtools")
devtools::install_github("petropouloumaria/NMAoutlier")
```
