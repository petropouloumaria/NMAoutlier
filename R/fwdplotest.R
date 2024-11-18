#' Forward plots of summary treatment estimates
#'
#' @description
#' Forward plots of summary treatment estimates with their 95 percent
#' confidence intervals.
#'
#' @param x object of class NMAoutlier (mandatory).
#'
#' @details
#' Plot of summary treatment estimates and their confidence intervals
#' for each FS iteration. Vertical axis provides the FS
#' iterations. Horizontal axis provides summary treatment estimates.
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
#' # Forward plot for summary treatment estimates and their confidence
#' # intervals
#' #
#' fwdplotest(FSresult)
#'
#' data(smokingcessation)
#'
#' # Transform data from arm-based format to contrast-based format
#' # We use 'sm' argument for odds ratios.
#' # We use function pairwise from netmeta package
#' #
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   list(event1, event2, event3), list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#'
#' # forward search algorithm
#' FSresult <- NMAoutlier(p1, small.values = "bad")
#'
#' # Forward plot for summary treatment estimates
#' # and their confidence intervals
#' fwdplotest(FSresult)
#' }
#'
#' @export
#'
#' @author Maria Petropoulou <maria.petropoulou@uniklinik-freiburg.de>


fwdplotest <- function(x) {

  ## Check class
  ##
  chkclass(x, "NMAoutlier")

  plotesthelper(x, lower = x$lb, upper = x$ub, estimate = x$estb, xdata = x$dat, xtitle = "Iterations")

}
