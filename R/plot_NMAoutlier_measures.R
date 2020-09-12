
#' Plot(s) to monitor selected outlier and influential  statistical measure(s) for each study.
#'
#' @description
#' This function generates plot(s) of the selected outlier detection measure(s)
#' for each study included in the network. The function creates a plot of the
#' selected detection statistic for each study in the network.
#' Candidate statistics to be monitored can be Standardized residual;
#' Studentized residual; Mahalanobis distance and
#' leverage for each study.
#'
#' @param object an object of class NMAoutlier_measures (mandatory).
#' @param stat selected statistical outlier and influential detection measure
#' (mandatory), available choices are: ("estand"/ "estud"/ "mah"/ "leverage").
#'
#' @details
#' Plot of outlier and influential detection measures for each study included in the network.
#' Vertical axis provides the studies in the network. Horizontal axis
#' provides a monitoring outlier and influential detection measure.
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
#' # outlier and influential detection measures for each study in the network
#' measures <- NMAoutlier_measures(p1)
#'
#' # plot of Standardized residuals for each study
#' plot_NMAoutlier_measures(measures, "estand")
#'
#' # plot of Mahalanobis distance values for each study
#' plot_NMAoutlier_measures(measures, "mah")
#'
#' # plot of leverage values for each study
#' plot_NMAoutlier_measures(measures, "leverage")
#'
#'
#' @export
#'
#' @author Maria Petropoulou <mpetrop@cc.uoi.gr>
#'
#' @importFrom ggplot2 ggplot aes theme element_rect element_line
#'   geom_point labs
#' @importFrom reshape2 melt



plot_NMAoutlier_measures <- function(object, stat){

  chkclass(object, "NMAoutlier_measures")

  stat <- setchar(tolower(stat), c("estand", "estud",
                          "mah", "leverage"))

  stlab <- unique(object$dat[ ,3])
  xlabel <- "study"

  if (tolower(stat) == "estand") {

    title = "Standardized residuals plot"
    ylabel = "Standardized residuals"
    data <- object$estand
    help_plot(data, stlab, title, xlabel, ylabel, "residual")

  } else if (tolower(stat) == "estud") {

    title = "Studentized residuals plot"
    ylabel = "Studentized residuals"
    data <- object$estud
    help_plot(data, stlab, title, xlabel, ylabel, "residual")

  } else if (tolower(stat) == "mah") {

    title = "Mahalanobis plot"
    ylabel = "Qi-Mahalanobis distance"
    data <- object$Mahalanobis.distance
    help_plot(data, stlab, title, xlabel, ylabel, "mah")

  } else if (tolower(stat) == "leverage") {

    title = "Leverage plot"
    ylabel = "Leverage"
    data <- object$leverage
    help_plot(data, stlab, title, xlabel, ylabel, "leverage")
  }


}



