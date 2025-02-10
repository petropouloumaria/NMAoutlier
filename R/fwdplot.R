#' Forward plot(s) to monitor selected statistic(s)/method(s)
#'
#' @description
#' This function generates forward plot(s) to monitor selected
#' statistic(s) and/or method(s). The function creates a plot of the
#' selected monitoring measure throughout the iterations of the Forward Search
#' algorithm. Candidate statistics to be monitored can be: P-score;
#' z-values by back-calculation method to derive indirect estimates
#' from direct pairwise comparisons and network estimates;
#' standardized residuals; heterogeneity variance estimator; Cook's
#' distance; ratio of variances; Q statistics (Overall heterogeneity /
#' inconsistency Q statistic (\code{Q}), overall heterogeneity Q
#' statistic (\code{Q}), between-designs Q statistic (\code{Q}), based
#' on a random effects design-by-treatment interaction model).
#'
#' @param x an object of class NMAoutlier (mandatory).
#' @param stat statistical measure to be monitored in forward plot(s)
#'   (mandatory), available choices are: "pscore", "nsplit", "estand",
#'   "heterog", "cook", "ratio", or "Q" (can be abbreviated).
#' @param select.st selected statistic (pscore/nsplit/estand) for
#'   selected treatment(s)/comparison(s)/study
#'
#' @details
#' Plot of statistical measures for each iteration of search.
#' Vertical axis provides the FS iterations. Horizontal axis
#' provides the values of the monitoring statistical measure.
#'
#' @keywords hplot
#'
#' @examples
#' \dontrun{
#' library("netmeta")
#' data(smokingcessation)
#' smokingcessation$id <- 1:nrow(smokingcessation)
#'
#' study912 <- subset(smokingcessation, id %in% 9:12)
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   list(event1, event2, event3), list(n1, n2, n3),
#'   data = study912, sm = "OR")
#'
#' # Forward search algorithm
#' #
#' FSresult <- NMAoutlier(p1, P = 1, small.values = "bad", n_cores = 2)
#'
#' # forward plot for Cook's distance
#' fwdplot(FSresult, "cook")
#'
#' data(smokingcessation)
#'
#' # Transform data from arm-based to contrast-based format
#' # Use 'sm' argument for odds ratios.
#' # Use function pairwise from netmeta package
#'
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   list(event1, event2, event3), list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#'
#' # Forward Search algorithm
#' FSresult <- NMAoutlier(p1, small.values = "bad")
#' FSresult
#'
#' # forward plot for Cook's distance
#' fwdplot(FSresult, "cook")
#'
#' # forward plot for ratio of variances
#' fwdplot(FSresult, "ratio")
#'
#' # forward plot for heterogeneity estimator
#' fwdplot(FSresult, "heterog")
#'
#' # forward plot for Q statistics
#' fwdplot(FSresult, "Q")
#'
#' # forward plot for P-scores
#' fwdplot(FSresult, "pscore")
#'
#' # forward plot monitoring P-scores for treatment A
#' fwdplot(FSresult,"pscore", "A")
#'
#' # forward plot for z-values of disagreement of direct and indirect evidence
#' fwdplot(FSresult, "nsplit")
#'
#' # forward plot for z-values of disagreement of direct and indirect evidence
#' # monitoring treatment comparison A versus B
#' fwdplot(FSresult, "nsplit", "A:B")
#'
#' # forward plot for standardized residual for study 4
#' fwdplot(FSresult, "estand", 4)
#' }
#'
#' @export
#'
#' @author Maria Petropoulou <maria.petropoulou@uniklinik-freiburg.de>



fwdplot <- function(x, stat, select.st = NULL) {

  ## Check class
  ##
  chkclass(x, "NMAoutlier")

  stat <- setchar(stat, c("pscore", "nsplit", "estand", "heterog",
                          "cook", "ratio", "Q"))

  plothelper(x, stat, select.st)

}

