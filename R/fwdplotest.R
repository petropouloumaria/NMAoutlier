#' Forward plots of summary estimates.
#'
#' Forward plots of summary estimate with 95 percent confidence interval for each treatment.
#' @param x object of class NMAoutlier (results of FS algorithm).
#' @return forward plots of summary estimates.
#' @export
#'
#' @author
#' Maria Petropoulou <mpetrop@cc.uoi.gr>
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics par lines abline
#' @importFrom plotrix plotCI


fwdplotest <- function(x) {

  ## Check class
  ##
  chkclass(x, "NMAoutlier")

  dev.new()
  ##
  ## Forward plots of summary estimate with 95% confidence interval for each treatment
  ## the reference treatment (j=1)
  ##
  lim.u <- round(max(x$ub), 2)
  lim.l <- round(min(x$lb), 2)

  ##
  nt <- length(unique(c(x$dat[,4], x$dat[,5])))
  ##
  par(mfrow = c(2, round(nt/2)))
  ##
  for (j in 2:nt) {

     plotrix::plotCI(x$estb[j,], y = NULL, ui = x$ub[j,],
                   li = x$lb[j,], ylim = c(lim.l,lim.u),
                   xlab = "iteration", ylab = rownames(x$estb)[j])
    ##
    lines(x$estb[j,], y = NULL)
    abline(h = 0, lty = 3)
  }

}




