#' Q-Q plot for network meta-analysis (Q-Q netplot).
#'
#' @description
#' This function generates the Q-Q plot for network meta-analyisis model.
#'
#' @param data object of class NMAoutlier.measures (mandatory).
#'
#' @details
#' Plot of Q-squared Mahalanobis distance for each study included in the network meta-analysis.
#' Vertical axis provides the Q-squared Mahalanobis distance for each i study included in the network meta-analysis.
#' Horizontal axis provides Q estimated quantiles (theoretical quantiles from the normal distribution).
#' A reference line is fitted from the cartesian points of the two measures.
#' The Q-Q plot can visualize studies that are away from the reference line (potiential outliers).
#'
#' Q-Q plot for network meta-analysis has been introduced by Petropoulou (2020).
#'
#' @references
#' Petropoulou M (2020):
#' Exploring methodological challenges in network meta-analysis models and
#' developing methodology for outlier detection.
#' \emph{PhD dissertation}
#'
#'
#' @examples
#' data(smokingcessation, package = "netmeta")
#'
#' p1 <- netmeta::pairwise(list(treat1, treat2, treat3),
#'                         list(event1, event2, event3),
#'                         list(n1, n2, n3),
#'                         data = smokingcessation,
#'                         sm = "OR")
#'
#' # Outlier and influential detection measures
#' measures <- NMAoutlier.measures(p1)
#'
#' # Mahalanobis distance values for each study in the network
#' measures$Mah
#'
#' # Q-Q netplot for the network of smoking cessation dataset
#' Qnetplot(measures)
#'
#' @keywords hplot
#'
#' @export
#'
#' @author Maria Petropoulou <petropoulou@imbi.uni-freiburg.de>
#'
#' @importFrom ggplot2 ggplot aes theme element_rect element_line
#'   geom_point geom_abline labs
#' @importFrom reshape2 melt
#' @importFrom stats qchisq quantile

Qnetplot <- function(data){

  chkclass(data, "NMAoutlier.measures")

  data <- data$Mah
  data <- data[order(data)]
  data <- t(data)
  melt_data <- melt(data)


  ## ith estimated quantile
  qch_chi <- c()
  for (i in 1:length(data)) {

    pr <- (i - 0.5)/length(data)
    qch_chi[i] <- qchisq(pr, 1)

  }

  yline <- quantile(data, c(.25, .75))
  pline <- quantile(qch_chi , c(.25, .75))

  slope_line <- (yline[2] - yline[1])/(pline[2] - pline[1])
  slope_line

  y_values <- melt_data$value
  ggplot(data = melt_data, aes(x = qch_chi, y = y_values)) +
    theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
    geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
    labs(title = "Qi-Qi netplot", y = "Qi-Squared Mahalanobis distance", x = "Qi-estimated quantile") +
    geom_abline(slope = 1, intercept = slope_line)

}


