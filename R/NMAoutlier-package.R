#' NMAoutlier: Brief overview of methodology for detection of outlying studies in network meta-analysis.
#'
#' @description
#' R package \bold{NMAoutlier} provides diagnostic methods to detect
#' outlying studies in network meta-analysis.
#'
#' @details
#' R package \bold{NMAoutlier} implements the forward search (FS)
#' algorithm for the detection of outlying studies (studies with extreme results) in
#' network meta-analysis (NMA) by Petropoulou et al. (2019). The
#' underlying model considered is the frequentist NMA approach based
#' on graph theory by Rücker (2012) which is implemented in R package
#' \bold{netmeta}.
#'
#' The \bold{NMAoutlier} package implements the following:
#' \itemize{
#' \item forward search algorithm in network meta-analysis (function
#'   \code{\link{NMAoutlier}}) based on Petropoulou et al. (2019);
#' \item forward plots (\code{\link{fwdplot}}) with monitoring
#'   statistics in each step of the FS algorithm:
#' \enumerate{
#' \item P-scores (Rücker & Schwarzer, 2015),
#' \item z-values for difference of direct and indirect evidence
#'   with back-calculation method (König et al., 2013; Dias et al.,
#'   2010),
#' \item standardized residuals,
#' \item heterogeneity variance estimator,
#' \item Cook's distance,
#' \item ratio of variances,
#' \item Q statistics (Krahn et al., 2013);
#' }
#' \item forward plot (\code{\link{fwdplotest}}) for summary estimates
#'   and their confidence intervals for each treatment in each step of
#'   the FS algorithm as provided by Petropoulou et
#'   al. (2019).
#' }
#'
#' Type \code{help(package = "NMAoutlier")} for a listing of R functions
#' available in \bold{NMAoutlier}.
#'
#' Type \code{citation("NMAoutlier")} on how to cite \bold{NMAoutlier} in
#' publications.
#'
#' To report problems and bugs, please send an email to Maria
#' Petropoulou \email{mpetrop@cc.uoi.gr}.
#'
#' The development version of \bold{NMAoutlier} is available on GitHub
#' \url{https://github.com/petropouloumaria/NMAoutlier}.
#'
#' @name NMAoutlier-package
#'
#' @docType package
#'
#' @author Petropoulou Maria \email{mpetrop@cc.uoi.gr}
#'
#' @references
#' Dias S, Welton NJ, Caldwell DM, Ades AE (2010):
#' Checking consistency in mixed treatment comparison meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{29}, 932--44
#'
#' König J, Krahn U, Binder H (2013):
#' Visualizing the flow of evidence in network meta-analysis and
#' characterizing mixed treatment comparisons.
#' \emph{Statistics in Medicine},
#' \bold{32}, 5414--29
#'
#' Krahn U, Binder H, König J (2013):
#' A graphical tool for locating inconsistency in network meta-analyses.
#' \emph{BMC Medical Research Methodology},
#' \bold{13}, 35
#'
#' Petropoulou M, Salanti G, Rücker G, Schwarzer G, Moustaki I,
#' Mavridis D (2019):
#' A forward search algorithm for detection of extreme study effects
#' in network meta-analysis.
#' \emph{Manuscript}
#'
#' Rücker G (2012):
#' Network meta-analysis, electrical networks and graph theory.
#' \emph{Research Synthesis Methods},
#' \bold{3}, 312--24
#'
#' Rücker G, Schwarzer G (2015):
#' Ranking treatments in frequentist network meta-analysis works
#' without resampling methods.
#' \emph{BMC Medical Research Methodology},
#' \bold{15}, 58
#'
#' @keywords package


NULL