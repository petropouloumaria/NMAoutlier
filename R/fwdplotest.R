#'
#' @title Forward plots of summary estimates.
#'
#' @description Forward plots of summary estimate with 95 percent confidence interval for each treatment.
#'
#' @usage
#' fwdplotest(x)
#'
#' @param x object of class NMAoutlier (mandatory).
#' @details Plot of summary estimates and their confidence intervals for each iteration of search.
#' Vertical axis provides iterations of search. Horizontal axis provides summary estimate of a treatment.
#' @return forward plots of summary estimates.
#'
#' @examples
#' \dontrun{
#' data(smokingdata)
#'
#' # forward search algorithm
#' FSresult <- NMAoutlier(TE, seTE, treat1, treat2,
#'                        studlab, data = smokingdata,
#'                        small.values = "bad")
#'
#' # forward plot for summary estimates for each treatment
#' # and their confidence intervals
#' fwdplotest(FSresult)
#'
#' }
#'
#' @export
#'
#' @author
#' Maria Petropoulou <mpetrop@cc.uoi.gr>
#'
#' @importFrom ggplot2 ggplot aes theme element_rect element_line geom_line geom_point geom_errorbar geom_line coord_cartesian geom_hline labs
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange

fwdplotest <- function(x) {

  ## Check class
  ##
  chkclass(x, "NMAoutlier")

  ## setting the lower and upper y limit for the graphs, so that we have the same limits for all treatments' graphs
  lim.u <- round(max(x$ub), 2)
  lim.l <- round(min(x$lb), 2)

  nt <- length(unique(c(x$dat[,4], x$dat[,5])))

  graphs <- vector("list", nt-1)

  ## generating the plots for every treatment. We ignore treatment A (j counts from 2)
  for (j in 2:nt) {
    local({ ## localizing variables - workaround for ggplot problem with handling variables in multiple plots (variable environment scope problem)
      melt_data <- melt(x$estb[j,])
      j <- j
      g <- eval(substitute(ggplot(data=melt_data, aes(x=1:length(x$estb[j,]), y=melt_data$value)) + ## eval(substitute) is another workaround for the aforementioned problem
         theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
         geom_point(color='#016FB9', size=2) +
         geom_errorbar(aes(ymin=x$lb[j,], ymax=x$ub[j,]), width=.5, colour="#8486F4") +
         geom_line(color='#016FB9') +
         coord_cartesian(ylim = c(lim.l,lim.u)) +
         labs(y=rownames(x$estb)[j], x="Iterations"), list(j=j))) +
         geom_hline(yintercept=0, linetype="dashed", color = "#EF2917", size=0.5)
      graphs[[j-1]] <<- g ## j-1 as we count from 2 (treatment B) in the loop
    })
  }

  ## creating a multiple graph consisting of a graph for each treatment
  grid.arrange(grobs = graphs, ncol = 2)

}




