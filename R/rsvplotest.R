#' Plot for summary estimates for each treatment
#' and their confidence intervals
#' for Random Shift Variace Model
#' fitted for selected/all studies of a dataset.
#'
#' @description
#' Plot for summary estimates for each treatment
#' and their 95 percent confidence intervals for Random Shift Variace Model
#' fitted for selected or all studies of a dataset
#'
#' @param x object of class NMAoutlier.rsv (mandatory).
#'
#' @details
#' Plot of summary estimates and their confidence intervals for each treatment
#' for Random Shift Variance Model fitted for each study.
#' Vertical axis provides study. Horizontal axis
#' provides summary estimate of a treatment.
#'
#' @keywords hplot
#'
#' @examples
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
#' # Random Shift Variace Model for study 1 of smoking cessation data
#' RSVresult <- NMAoutlier.rsv(p1, small.values = "bad", study = c(1), n_cores = 2)
#'
#' # Plot for summary estimates for each treatment
#' # and their confidence intervals for Random Shift Variace Model
#' # fitted for study 1
#' rsvplotest(RSVresult)
#'
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
#' # Random Shift Variace Model for each study of smoking cessation data
#' RSVresult <- NMAoutlier.rsv(p1, small.values = "bad")
#'
#' # Plot for summary estimates for each treatment
#' # and their confidence intervals for Random Shift Variace Model
#' # fitted for each study
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
