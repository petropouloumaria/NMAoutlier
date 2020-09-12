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
#' @param x object of class NMAsvr (mandatory).
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
#' SVRresult <- NMAsvr(p1, small.values = "bad", study = c(1), n_cores = 2)
#'
#' # Plot for summary estimates for each treatment
#' # and their confidence intervals for Random Shift Variace Model
#' # fitted for study 1
#' svrplotest(SVRresult)
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
#' SVRresult <- NMAsvr(p1, small.values = "bad")
#'
#' # Plot for summary estimates for each treatment
#' # and their confidence intervals for Random Shift Variace Model
#' # fitted for each study
#' svrplotest(SVRresult)
#' }
#'
#' @export
#'
#' @author Maria Petropoulou <mpetrop@cc.uoi.gr>


svrplotest <- function(x) {

  ## Check class
  ##
  chkclass(x, "NMAsvr")

  plotesthelper(x, lower = x$l, upper = x$u, estimate = x$b, xdata = x$dat, xtitle = "Study", method = "svr")

}
