NMAoutlier: Forward Search Algorithm in Network Meta-Analysis to identify outlying and influential studies
================

[![Build Status](https://travis-ci.org/petropouloumaria/NMAoutlier.svg?branch=master)](https://travis-ci.org/petropouloumaria/NMAoutlier)

<img src="man/figures/NMAoutlier_logo.png" width=300 align="right" style="margin-left:20px; margin-right: 20px;"/>

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

Installation using R package **[devtools](https://cran.r-project.org/package=devtools)** (with [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for Windows):

``` r
install.packages("devtools")
devtools::install_github("petropouloumaria/NMAoutlier")
```

Usage
-----

Example of outlying detection in network meta-analysis comparing the relative effects of four smoking cessation counseling programs, no contact (A), self-help (B), individual counseling (C) and group counseling (D). The outcome is the number of individuals with successful smoking cessation at 6 to 12 months. The dataset used as example in Dias et al.(2013). The dataset for smoking cessation is a part of data in netmeta package.

Reference: Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G and Ades AE (2013). Evidence Synthesis for Decision Making 4: Inconsistency in networks of evidence based on randomized controlled trials. Medical Decision Making 33, 641–656.

You can load the **NMAoutlier** library 

``` r
library(NMAoutlier)
```

Load the dataset smoking cessation from netmeta package. 
``` r
data(smokingcessation, package = "netmeta")
```

Transform data from arm-based format to contrast-based format using the function pairwise from netmeta package.

``` r
library(netmeta)
p1 <- pairwise(list(treat1, treat2, treat3),
              list(event1, event2, event3),
              list(n1, n2, n3),              
              data=smokingcessation,
              sm="OR")
```

You can conduct the forward search algorithm with function **NMAoutlier** as follows:

``` r
FSresult <- NMAoutlier(p1, small.values = "bad")
```

You can see the forward plots with function **fwdplot** for monitoring measures. For example, you can plot the influential diagnostic measure Cook distance as follows:

``` r
fwdplot(FSresult,"cook")
```

<img src="man/figures/fwdplot-cook-distance.png" width=450 style="margin-left: auto; margin-right: auto; display: block;"/>

Or you can plot the ratio of variances as follows:

``` r
fwdplot(FSresult,"ratio")
```

<img src="man/figures/fwdplot-ratio.png" width=450 style="margin-left: auto; margin-right: auto; display: block;"/>

You can plot the differences of direct and indirect estimates (z-values) as follows:

``` r
fwdplot(FSresult,"nsplit")
```

<img src="man/figures/fwdplot-diff-direct-indirect.png" width=700 style="margin-left: auto; margin-right: auto; display: block;"/>

You can see the forward plots for summary estimates for each treatment B, C and D with function **fwdplotest** as follows:

``` r
fwdplotest(FSresult)
```

<img src="man/figures/fwdplot-summary-estimates.png" style="margin-left: auto; margin-right: auto; display: block;"/>
