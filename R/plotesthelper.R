#' Helper function to create the plots of summary estimates and their intervals for FS and RVSOM.
#'
#' @param x object of class NMAoutlier or class NMAoutlier.rsv (mandatory).
#' @param lower lower bundary of confidence interval of summary estimate.
#' @param upper upper bundary of confidence interval of summary estimate.
#' @param estimate summary estimate.
#' @param xdata data.
#' @param xtitle title for plots in x-axis.
#' @param method The method that was used.
#' Select "fs" to monitor statistics during the forward search algorithm.
#' Select "rsv" to monitor statistics by fitting the data with the Shift Variance Model.
#'
#' @details
#' Plot of summary estimate and its confidence interval for each treatment for FS and RVSOM methodologies.
#' Vertical axis provides the iterations of FS methodology or the study for RVSOM methodology.
#' Horizontal axis provides a summary estimate of a treatment.
#'
#'
#' @keywords internal
#'
#' @author Maria Petropoulou <petropoulou@imbi.uni-freiburg.de>
#'
#' @importFrom ggplot2 ggplot aes theme element_rect element_line
#'   geom_line geom_point geom_errorbar geom_line coord_cartesian
#'   geom_hline labs
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange



plotesthelper <- function(x, lower, upper, estimate, xdata, xtitle, method) {

  ## Setting the lower and upper y limit for the graphs, so that we
  ## have the same limits for all treatments' graphs
  lim.u <- round(max(upper), 2)
  lim.l <- round(min(lower), 2)

  nt <- length(unique(c(xdata[, 4], xdata[, 5])))

  graphs <- vector("list", nt - 1)


  if (tolower(method) == "rsv") {
    if (length(x$z[-1]) == 1) {
      estimate <- as.data.frame(estimate)
      lower <- as.data.frame(lower)
      upper <- as.data.frame(upper)
    }
    xlabels <- factor(as.character(x$z[-1]), levels = as.character(x$z[-1])) # as factor to prevent ggplot from reordering the x labels in alphabetical order
  } else {
    xlabels <- factor(as.character(1:x$index), levels = as.character(1:x$index))
  }


  ## Generating the plots for every treatment. We ignore treatment A (j counts from 2)
  for (j in 2:nt) {

    ## Localizing variables - workaround for ggplot problem with
    ## handling variables in multiple plots (variable environment
    ## scope problem)
    local({

      melt_data <- melt(estimate[j, ], id.vars = 0)
      j <- j
      g <- eval(substitute(ggplot(data = melt_data, aes(x = xlabels, y = melt_data$value)) + # eval(substitute) is another workaround for the aforementioned problem
                             theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
                             geom_point(color = '#016FB9', size = 2) +
                             geom_errorbar(aes(ymin = lower[j,], ymax = upper[j,]), width = .5, colour = "#8486F4") +
                             geom_line(color = '#016FB9') +
                             coord_cartesian(ylim = c(lim.l,lim.u)) +
                             labs(y = rownames(estimate)[j], x = xtitle), list(j = j))) +
                             geom_hline(yintercept = 0, linetype = "dashed", color = "#EF2917", size = 0.5)
      graphs[[j - 1]] <<- g # j - 1 as we count from 2 (treatment B) in the loop
    })
  }

  ## creating a multiple graph consisting of a graph for each treatment
  grid.arrange(grobs = graphs, ncol = 2)
}
