#' Plots for summary treatment estimates in RSV NMA model.
#'
#' @description
#' Plot for summary treatment estimates for RSV NMA model.

#'
#' @param x object of class NMAoutlier.rsv (mandatory).
#'
#' @details
#' Plot of summary treatment estimates and their confidence intervals
#'  RSV NMA model fitted for each study.
#' Vertical axis provides study. Horizontal axis
#' provides summary treatment estimates.
#'
#' @keywords hplot
#'
#' @examples
#' data(smokingcessation, package = "netmeta")
#'
#' # Transform data from arm-based to contrast-based format
#' # We use 'sm' argument for odds ratios.
#' # We use function pairwise from netmeta package
#' #
#' p1 <- netmeta::pairwise(list(treat1, treat2, treat3),
#'                         list(event1, event2, event3),
#'                         list(n1, n2, n3),
#'                         data=smokingcessation,
#'                         sm="OR")
#'
#' # RSV NMA model for study 1 of smoking cessation data
#' RSVresult <- NMAoutlier.rsv(p1, small.values = "bad", study = c(1), n_cores = 2)
#'
#' # Plot for summary treatment estimates
#' # and their confidence intervals with RSV NMA model
#' # fitted for study 1 (shift the variance of study 1)
#' rsvplotest(RSVresult)
#'
#' \dontrun{
#' data(smokingcessation, package = "netmeta")
#'
#' # Transform data from arm-based to contrast-based format
#' # We use 'sm' argument for odds ratios.
#' # We use function pairwise from netmeta package
#' #
#' p1 <- netmeta::pairwise(list(treat1, treat2, treat3),
#'                         list(event1, event2, event3),
#'                         list(n1, n2, n3),
#'                         data=smokingcessation,
#'                         sm="OR")
#'
#' # RSV NMA model for each study of smoking cessation data
#' RSVresult <- NMAoutlier.rsv(p1, small.values = "bad")
#'
#' # Plot for summary treatment estimates and their confidence intervals
#' # for RSV NMA model fitted for each study
#' rsvplotest(RSVresult)
#' }
#'
#' @export
#'
#' @author Maria Petropoulou <petropoulou@imbi.uni-freiburg.de>


rsvplotest <- function(x) {

  ## Check class
  ##
  chkclass(x, "NMAoutlier.rsv")

  plotesthelper(x, lower = x$l, upper = x$u, estimate = x$b, xdata = x$dat, xtitle = "Study", method = "rsv")

}
