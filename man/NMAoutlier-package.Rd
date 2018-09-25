\name{NMAoutlier-package}

\alias{NMAoutlier-package}

\docType{package}

\title{NMAoutlier: Brief overview of methodology for detection of influential and outlying studies.}

\description{R package \bold{NMAoutlier} provides a diagnostic methodology of influential and outlying studies in network meta-analysis}

\details{
  R package \bold{NMAoutlier} providing the methodology for
  detection of extreme study effects by Petropoulou et al. (2018).
  The methodology based on forward search algorithm for detection on
  oultying and influential studies fitted in the
  network meta-analysis model from graph-theory as provided by Rücker (2012).
  \bold{NMAoutlier} package implements the following methodological aspects:
  \itemize{
  \item forward search algorithm in network meta-analysis (function
    \code{\link{NMAoutlier}}) based on Petropoulou et al. (2018).
  \item forward plots (\code{\link{fwdplot}}) for the monitoring statistics
  in each step of forward search algorithm (a) P-scores, Rücker G & Schwarzer G (2015);
  (b) z-values for difference of direct and indirect evidence with back-calculation method, König (2013), Dias (2010); (c) Standardized residuals; (d) heterogeneity variance estimator; (e) cook distance;
  (f) ratio of variances; (g) Q statistics, Krahn et al. (2013).
  \item forward plot (\code{\link{fwdplotest}}) for summary estimates and their confidence intervals for each
  treatment in each step of forward search algorithm as provided by Petropoulou et al. (2018).
  }

  Type \code{help(package = "NMAoutlier")} for a listing of R functions
  available in \bold{NMAoutlier}.

  Type \code{citation("NMAoutlier")} on how to cite \bold{NMAoutlier} in
  publications.

  To report problems and bugs
  \itemize{
  \item please send an email to Maria Petropoulou
    \email{mpetrop@cc.uoi.gr}.
  }
  The development version of \bold{NMAoutlier} is available on GitHub
  \url{https://github.com/petropouloumaria/NMAoutlier}.
}
\references{
  Rücker G (2012),
  Network meta-analysis, electrical networks and graph theory.
  \emph{Research Synthesis Methods},
  \bold{3}, 312--24.

  Rücker G & Schwarzer G (2015),
  Ranking treatments in frequentist network meta-analysis works
  without resampling methods.
  \emph{BMC Medical Research Methodology},
  \bold{15}, 58, DOI:10.1186/s12874-015-0060-8.

  Dias S, Welton NJ, Caldwell DM, Ades AE (2010).
  Checking consistency in mixed treatment comparison meta-analysis.
  \emph{Statistics in Medicine}, \bold{29}, 932--44.

  Krahn U, Binder H, König J (2013),
  A graphical tool for locating inconsistency in network
  meta-analyses.
  \emph{BMC Medical Research Methodology}, \bold{13}, 35.

  Petropoulou M, Salanti G, Rücker G, Schwarzer G, Moustaki I, Mavridis D (2018),
  A forward search algorithm for detection of extreme study effects
  in network meta-analysis. Manuscript.

  König J, Krahn U, Binder H (2013). Visualizing the flow of evidence in network
  meta-analysis and characterizing mixed treatment comparisons.
  \emph{Statistics in Medicine}, \bold{32}(30), 5414–29.

}
\author{Petropoulou Maria \email{mpetrop@cc.uoi.gr}}
