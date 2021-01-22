#' Forward plots of summary estimates
#'
#' @description
#' Forward plots of summary estimate with 95 percent confidence
#' interval for each treatment.
#'
#' @param x object of class NMAoutlier (mandatory).
#'
#' @details
#' Plot of summary estimates and their confidence intervals for each
#' iteration of search.  Vertical axis provides iterations of
#' search. Horizontal axis provides summary estimate of a treatment.
#'
#' @keywords hplot
#'
#' @examples
#' data(smokingcessation, package = "netmeta")
#' smokingcessation$id <- 1:nrow(smokingcessation)
#'
#' study912 <- subset(smokingcessation, id %in% 9:12)
#' p1 <- netmeta::pairwise(list(treat1, treat2, treat3),
#'                         list(event1, event2, event3),
#'                         list(n1, n2, n3),
#'                         data = study912,
#'                         sm = "OR")
#'
#' # Forward search algorithm
#' #
#' FSresult <- NMAoutlier(p1, P = 1, small.values = "bad", n_cores = 2)
#'
#' # Forward plot for summary estimates for each treatment
#' # and their confidence intervals
#' fwdplotest(FSresult)
#' \dontrun{
#' data(smokingcessation, package = "netmeta")
#'
#' # Transform data from arm-based format to contrast-based format
#' # We use 'sm' argument for odds ratios.
#' # We use function pairwise from netmeta package
#' #
#' p1 <- netmeta::pairwise(list(treat1, treat2, treat3),
#'                         list(event1, event2, event3),
#'                         list(n1, n2, n3),
#'                         data=smokingcessation,
#'                         sm="OR")
#'
#' # forward search algorithm
#' FSresult <- NMAoutlier(p1, small.values = "bad")
#'
#' # Forward plot for summary estimates for each treatment
#' # and their confidence intervals
#' fwdplotest(FSresult)
#' }
#'
#' @export
#'
#' @author Maria Petropoulou <petropoulou@imbi.uni-freiburg.de>


fwdplotest <- function(x) {

  ## Check class
  ##
  chkclass(x, "NMAoutlier")

  plotesthelper(x, lower = x$lb, upper = x$ub, estimate = x$estb, xdata = x$dat, xtitle = "Iterations", method = "fs")

}
