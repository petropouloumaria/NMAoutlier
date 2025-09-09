NMAoutlier: Detecting Outliers in Network Meta-Analysis
================

[![License: GPL
(\>=2)](https://img.shields.io/badge/license-GPL-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![CRAN
Version](https://www.r-pkg.org/badges/version/NMAoutlier)](https://cran.r-project.org/package=NMAoutlier)
[![CRAN_time_from_release](https://www.r-pkg.org/badges/ago/NMAoutlier)](https://cran.r-project.org/package=NMAoutlier)
[![Monthly
Downloads](https://cranlogs.r-pkg.org/badges/NMAoutlier)](https://cranlogs.r-pkg.org/badges/NMAoutlier)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/NMAoutlier)](https://cranlogs.r-pkg.org/badges/grand-total/NMAoutlier)

<img src="man/figures/NMAoutlier_logo.png" width=300 align="right" style="margin-left:20px; margin-right: 20px;"/>

## Description

A package that provides measures and methodologies for detecting
outlying and influential studies in network meta-analysis.

- **1) Simply outlier and influential detection measures:** Raw,
  Standardized, Studentized residuals; Mahalanobis distance and
  leverage.
- **2) Outlier and influential detection measures by considering a study
  deletion (Shift the mean):** Raw, Standardized, Studentized deleted
  residuals; Cook’s distance; COVRATIO; weight “leave one out”; leverage
  “leave one out”; heterogeneity “leave one out”; R heterogeneity; R
  Qtotal; R Qheterogeneity; R Qinconsistency and DFBETAS.
- Plots for all the above outlier and influential detection measures
  (simple and deletion measures) and Q-Q plot for network meta-analysis.
- **3) Forward search algorithm in network meta-analysis (FS).**
- Forward plots (fwdplot) for the monitoring measures in each step of
  forward search algorithm. Monitoring measures: P-scores; z-values for
  difference of direct and indirect evidence with back-calculation
  method; Standardized residuals; heterogeneity variance estimator;
  Cook’s distance; ratio of variances; Q statistics.
- Forward plot for summary estimates and their confidence intervals for
  each treatment in each step of forward search algorithm.

## Installation

You can install the **NMAoutlier** package from GitHub repository as
follows:

Installation using R package
**[remotes](https://cran.r-project.org/package=remotes)**:

``` r
install.packages("remotes")
remotes::install_github("petropouloumaria/NMAoutlier")
```

## Usage

Example of network meta-analysis comparing the relative effects of four
smoking cessation counseling programs, no contact (A), self-help (B),
individual counseling (C) and group counseling (D). The outcome is the
number of individuals with successful smoking cessation at 6 to 12
months. The data are in contrast format with odds ratio (OR) and its
standard error. Arm-level data can be found in Dias et al. (2013).

References:

Higgins D, Jackson JK, Barrett G, Lu G, Ades AE, and White IR.
Consistency and inconsistency in network meta-analysis: concepts and
models for multi-arm studies. Research Synthesis Methods 2012, 3(2):
98–110.

Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, and Ades AE. Evidence
Synthesis for Decision Making 4: Inconsistency in networks of evidence
based on randomized controlled trials. Medical Decision Making 2013, 33:
641–656.

You can load the **NMAoutlier** library

``` r
library(NMAoutlier)
```

Load the dataset smoking cessation from **netmeta** package.

``` r
data(smokingcessation, package = "netmeta")
```

Transform data from arm-based to contrast-based format using the
function **pairwise** from **netmeta** package.

``` r
library(netmeta)
p1 <- pairwise(list(treat1, treat2, treat3),
               list(event1, event2, event3),
               list(n1, n2, n3),
               data = smokingcessation,
               sm = "OR")
```

**Part 1: Simply outlier detection measures**

You can calculate simply outlier and influential detection measures with
**NMAoutlier.measures** function as follows:

``` r
measures <- NMAoutlier.measures(p1)
```

You can see the Mahalanobis distance for each study

``` r
measures$Mahalanobis.distance
```

You can plot the Mahalanobis distance for each study with **measplot**
function as follows:

``` r
measplot(measures, "mah")
```

<img src="man/figures/Mahalanobis.png" width=450 style="margin-left: auto; margin-right: auto; display: block;"/>

You can figure out the Q-Q plot for network meta-analysis with
**Qnetplot** function as follows:

``` r
Qnetplot(measures)
```

<img src="man/figures/Q-Q.png" width=450 style="margin-left: auto; margin-right: auto; display: block;"/>

**Part 2: Outlier detection measures considered deletion (Shift the
mean)**

You can calculate outlier and influential detection measures considered
study deletion with **NMAoutlier,measures** function as follows:

``` r
deletion <- NMAoutlier,measures(p1, measure = "deletion")
```

You can see the standardized deleted residuals for each study

``` r
deletion$estand.deleted
```

You can see the COVRATIO for each study

``` r
deletion$Covratio
```

You can plot the R statistic for Qinconsistency with function
**measplot** as follows:

``` r
measplot(deletion, "rqinc", measure = "deletion")
```

<img src="man/figures/Qinc.png" width=450 style="margin-left: auto; margin-right: auto; display: block;"/>

**Part 3: Forward Search Algorithm - (Outlier detection Methodology)**

You can conduct the Forward Search algorithm with **NMAoutlier**
function as follows:

``` r
FSresult <- NMAoutlier(p1, small.values = "bad")
```

You can see the forward plots with **fwdplot** function for Cook’s
distance as follows:

``` r
fwdplot(FSresult,"cook")
```

<img src="man/figures/fwdplot-cook-distance.png" width=450 style="margin-left: auto; margin-right: auto; display: block;"/>

Or you can plot the Ratio of variances as follows:

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

You can see the forward plots for summary relative treatment estimates
of B, C and D versus the reference A with **fwdplotest** function as
follows:

``` r
fwdplotest(FSresult)
```

<img src="man/figures/fwdplot-summary-estimates.png" style="margin-left: auto; margin-right: auto; display: block;"/>
