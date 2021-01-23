NMAoutlier: Detecting Outliers in Network Meta-Analysis
================

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/NMAoutlier)](https://cran.r-project.org/package=NMAoutlier)
[![CRAN\_time\_from\_release](https://www.r-pkg.org/badges/ago/NMAoutlier)](https://cran.r-project.org/package=NMAoutlier)
[![Monthly
Downloads](http://cranlogs.r-pkg.org/badges/NMAoutlier)](http://cranlogs.r-pkg.org/badges/NMAoutlier)
[![Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/NMAoutlier)](http://cranlogs.r-pkg.org/badges/grand-total/NMAoutlier)
[![Build
Status](https://travis-ci.org/petropouloumaria/NMAoutlier.svg?branch=master)](https://travis-ci.org/petropouloumaria/NMAoutlier)

<img src="man/figures/NMAoutlier_logo.png" width=300 align="right" style="margin-left:20px; margin-right: 20px;"/>

## Description

A package that provides for detecting outlying studies in network
meta-analysis.

-   **Several outlier detection measures provided:** Raw, Standardized,
    Studentized residuals; Mahalanobis distance and leverage.
-   **Plots for outlier and influence measures:** Raw, Standardized,
    Studentized residuals; Mahalanobis distance and leverage.
-   **Several outlier and influence detection measures considered
    deletion:** Raw, Standardized,Studentized deleted residuals; Cook
    distance; COVRATIO; weight “leave one out”; leverage “leave one
    out”; heterogeneity “leave one out”; R heterogeneity; R Qtotal; R
    Qheterogeneity; R Qinconsistency and DFBETAS.
-   **Plots for outlier and influence detection measures considered
    deletion:** the aboved measures.
-   **Q-Q plot for network meta-analysis.**
-   **Forward search algorithm in network meta-analysis (FS).**
-   **Forward plots (fwdplot) for the monitoring statistics** in each
    step of forward search algorithm P-scores; z-values for difference
    of direct and indirect evidence with back-calculation method;
    Standardized residuals; heterogeneity variance estimator; cook
    distance; ratio of variances; Q statistics.
-   **Forward plot for summary estimates** and their confidence
    intervals for each treatment in each step of forward search
    algorithm.  
-   **Random shift variance NMA model (RVSOM NMA).**
-   **Plots for the monitoring measures for random shift variance
    model.**
-   **Plots for the for summary estimates** of random shift variance
    model.

## Installation

You can install the **NMAoutlier** package from GitHub repository as
follows:

Installation using R package
**[devtools](https://cran.r-project.org/package=devtools)**:

``` r
install.packages("devtools")
devtools::install_github("petropouloumaria/NMAoutlier")
```

## Usage

Example of outlying detection in network meta-analysis comparing the
relative effects of four smoking cessation counseling programs, no
contact (A), self-help (B), individual counseling (C) and group
counseling (D). The outcome is the number of individuals with successful
smoking cessation at 6 to 12 months. These data are in contrast format
with effect size odds ratio (OR) and its standard error. Arm-level data
can be found in Dias et al.(2013).

Reference:

Higgins D, Jackson JK, Barrett G, Lu G, Ades AE, and White IR.
Consistency and inconsistency in network meta-analysis: concepts and
models for multi-arm studies. Research Synthesis Methods 2012, 3(2):
98–110.

Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G and Ades AE (2013).
Evidence Synthesis for Decision Making 4: Inconsistency in networks of
evidence based on randomized controlled trials. Medical Decision Making
33, 641–656.

You can load the **NMAoutlier** library

``` r
library(NMAoutlier)
```

Load the dataset smoking cessation from netmeta package.

``` r
data(smokingcessation, package = "netmeta")
```

Transform data from arm-based format to contrast-based format using the
function pairwise from netmeta package.

``` r
library(netmeta)
p1 <- pairwise(list(treat1, treat2, treat3),
              list(event1, event2, event3),
              list(n1, n2, n3),              
              data=smokingcessation,
              sm="OR")
```

**Part 1: Simply outlier detection measures**

You can calculate some simply outlier detection measures with function
**NMAoutlier.measures** as follows:

``` r
measures <- NMAoutlier.measures(p1)
```

You can see the Mahalanobis distance for each study

``` r
measures$Mahalanobis.distance
```

You can plot the Mahalanobis distance for each study with function
**measplot** as follows:

``` r
measplot(measures, "mah")
```

You can plot the Q-Q plot for network meta-analysis with function
**Qnetplot** as follows:

``` r
Qnetplot(measures)
```

**Part 2: Outlier detection measures considered deletion**

You can calculate some outlier detection measures considered deletion
with function **NMAoutlier.measures** as follows:

``` r
deletion <- NMAoutlier.measures(p1, measure = "deletion")
```

You can see the standardized deleted residuals for each study

``` r
deletion$estand.deleted
```

We can see the values of COVRATIO when considering deletion for each
study

``` r
deletion$Covratio
```

We can plot the R statistic for Qinconsistency with function
**measplot** as follows:

``` r
measplot(deletion, "rqinc", measure = "deletion")
```

**Part 3: Forward Search Algorithm - Detection Methodology**

You can conduct the forward search algorithm with function
**NMAoutlier** as follows:

``` r
FSresult <- NMAoutlier(p1, small.values = "bad")
```

You can see the forward plots with function **fwdplot** for monitoring
measures. For example, you can plot the diagnostic measure Cook distance
as follows:

``` r
fwdplot(FSresult,"cook")
```

<img src="man/figures/fwdplot-cook-distance.png" width=450 style="margin-left: auto; margin-right: auto; display: block;"/>

Or you can plot the ratio of variances as follows:

``` r
fwdplot(FSresult,"ratio")
```

<img src="man/figures/fwdplot-ratio.png" width=450 style="margin-left: auto; margin-right: auto; display: block;"/>

You can plot the differences of direct and indirect estimates (z-values)
as follows:

``` r
fwdplot(FSresult,"nsplit")
```

<img src="man/figures/fwdplot-diff-direct-indirect.png" width=700 style="margin-left: auto; margin-right: auto; display: block;"/>

You can see the forward plots for summary estimates for each treatment
B, C and D with function **fwdplotest** as follows:

``` r
fwdplotest(FSresult)
```

<img src="man/figures/fwdplot-summary-estimates.png" style="margin-left: auto; margin-right: auto; display: block;"/>

**Part 4: Shift Variance Network Meta-analysis – Detection methodology
and sensitivity analysis downweighing outlier**

You can conduct the random shift variance model for each study with
function **NMAoutlier.rsv** as follows:

``` r
RSVresult <- NMAoutlier.rsv(p1, small.values = "bad")
```

You can see the the Likelihood Ratio Test (LRT) with random shift
variance model of each study

``` r
RSVresult$LRT 
```

Or you can see the over-dispersion with random shift variance model of
each study

``` r
RSVresult$over_disp
```

You can plot the of Likelihood Ratio Test (LRT) of random shift variance
model for each study with function **rsvplot** as follows:

``` r
rsvplot(RSVresult, "LRT")
```

You can see the plots for summary estimates for each treatment B, C and
D with function **rsvplotest** as follows:

``` r
rsvplotest(RSVresult)
```
