#' Plot(s) to monitor selected outlier and influential measure(s).
#'
#' @description
#' This function generates plot(s) of the selected outlier detection
#' measure(s) for each study included in the network. Candidate
#' statistics to be monitored are: Standardized residual; Studentized
#' residual; Mahalanobis distance and leverage.
#'
#' The function also generates plot(s) of the selected outlier
#' detection measure(s) considering a deletion of each study included
#' in the network (Shift the mean measures).  Candidate statistics to
#' be monitored are: Standardized deleted residual; Studentized
#' deleted residual; Cook distance between the treatment estimates for
#' study j and treatment estimates when study j is removed; Ratio of
#' determinants of variance-covariance matrix of treatment estimates
#' for study j to treatment estimates when study j is removed; weight
#' leave one out;leverage leave one out; heterogeneity estimator leave
#' one out; R statistic for heterogeneity; R statistic for Q
#' (\code{Qtotal}), R statistic for heterogeneity Q (\code{Qhet}), R
#' statistic for Qinconsistency (\code{Qinc}), DFbetas.
#'
#'
#' @param object an object of class NMAoutlier.measures (mandatory).
#' @param stat selected statistical outlier and influential detection
#'   measure (mandatory), For simply outlier and influential measures
#'   available choices are: ("estand"/ "estud"/ "mah"/ "leverage").
#'   For outlier and influential deletion measures available choices
#'   are: ("estand.deleted", "estud.deleted", "leverage.leaveoneout",
#'   "weight.leaveoneout", "heterog.leaveoneout", "covratio", "cook",
#'   "rheterogeneity", "restimates", "rqhet", "rqinc", "rqtotal",
#'   "dfbetas")
#' @param measure Outlier and influential detection measures. Simple
#'   measures (default: "simple") and measures considered study
#'   deletion (measure = "deletion").
#'
#' @details
#' Plot of outlier and influential (simple or/and deletion) detection
#' measures for each study included in the network.  Vertical axis
#' provides each study included in the network (or the study deleted
#' for outlier deletion measures). Horizontal axis provides a
#' monitoring outlier and influential detection measure.
#'
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
#' # Outlier and influential detection measures for each study in the
#' # network
#' measures <- NMAoutlier.measures(p1)
#'
#' # plot of standardized residuals for each study
#' measplot(measures, "estand")
#' 
#' # plot of Mahalanobis distance values for each study
#' measplot(measures, "mah")
#'
#' # plot of leverage values for each study
#' measplot(measures, "leverage")
#'
#' \dontrun{
#' # Outlier detection measures considered deletion each time of an
#' # included study
#' deletion <- NMAoutlier.measures(p1, measure = "deletion")
#'
#' # plot for R statistic for heterogeneity estimator
#' measplot(deletion, "rheterogeneity", measure = "deletion")
#'
#' # plot for R statistic for Qinconsistency
#' measplot(deletion, "rqinc", measure = "deletion")
#'
#' # plot of COVRATIO values when considering deletion for each study
#' measplot(deletion, "covratio", measure = "deletion")
#' }
#'
#' @export
#'
#' @author Maria Petropoulou <petropoulou@imbi.uni-freiburg.de>
#'
#' @importFrom ggplot2 ggplot aes theme element_rect element_line
#'   geom_point labs
#' @importFrom reshape2 melt



measplot <- function(object, stat, measure = "simple"){

  chkclass(object, "NMAoutlier.measures")

  if (measure == "simple") {

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
      if (!is.null(object$leverage))
        data <- object$leverage
      else
        data <- object$leverage.fixed
      help_plot(data, stlab, title, xlabel, ylabel, "leverage")

    }

  } else if (measure == "deletion") {

    stat <- setchar(stat, c("estand.deleted", "estud.deleted",
                            "leverage.leaveoneout", "weight.leaveoneout", "heterog.leaveoneout", "covratio", "cook",
                            "rheterogeneity", "restimates", "rqhet", "rqinc", "rqtotal", "dfbetas"))

    stlab <- unique(object$dat[ ,3])
    xlabel <- "study deleted"

   if (tolower(stat) == "estand.deleted") {

      title = "Standardized study deleted residual"
      ylabel = "Standardized study deleted residual"
      data <- object$estand.deleted
      help_plot(data, stlab, title, xlabel, ylabel, "residual")

    } else if (tolower(stat) == "estud.deleted") {

      title = "Studentized study deleted residual"
      ylabel = "Studentized study deleted residual"
      data <- object$estud.deleted
      help_plot(data, stlab, title, xlabel, ylabel, "residual")

    }  else if (tolower(stat) == "leverage.leaveoneout") {

      title = "leverage leave-one-out"
      ylabel = "leverage"
      data <- object$H.leaveoneout
      help_plot(data, stlab, title, xlabel, ylabel, "H.leaveoneout")

    } else if (tolower(stat) == "heterog.leaveoneout") {

      title = "Heterogeneity leave-one-out"
      ylabel = "heterogeneity estimator"
      data <- object$heterog.leaveoneout
      help_plot(data, stlab, title, xlabel, ylabel, "heterog.leaveoneout")

    } else if (tolower(stat) == "weight.leaveoneout") {

      title = "Weight leave-one-out"
      ylabel = "weight"
      data <- object$w.leaveoneout
      help_plot(data, stlab, title, xlabel, ylabel, "weight.leaveoneout")

    } else if (tolower(stat) == "covratio") {

      title = "Covratio leave-one-out"
      ylabel = "Covtatio"
      data <- object$Covratio
      help_plot(data, stlab, title, xlabel, ylabel, "covratio")

    } else if (tolower(stat) == "cook") {

      title = "Cook's distance leave-one-out"
      ylabel = "Cook's distance"
      data <- object$Cooks.distance
      help_plot(data, stlab, title, xlabel, ylabel, "cook")


    } else if (tolower(stat) == "rheterogeneity") {

      title = "R statistic for heterogeneity leave-one-out"
      ylabel = "R statistic for heterogeneity"
      data <- object$Rheterogeneity
      help_plot(data, stlab, title, xlabel, ylabel, "rheterogeneity")

    } else if (tolower(stat) == "rqhet") {

      title = "R statistic for Qheterogeneity leave-one-out"
      ylabel = "R statistic for Qheterogeneity"
      data <- object$RQhet
      help_plot(data, stlab, title, xlabel, ylabel, "rqhet")

    } else if (tolower(stat) == "rqinc") {

      title = "R statistic for Qinconsistency leave-one-out"
      ylabel = "R statistic for inconsistency"
      data <- object$RQinc
      help_plot(data, stlab, title, xlabel, ylabel, "rqinc")

    } else if (tolower(stat) == "rqtotal") {

      title = "R statistic for Qtotal leave-one-out"
      ylabel = "R statistic for Qtotal"
      data <- object$RQtotal
      help_plot(data, stlab, title, xlabel, ylabel, "rqtotal")

    } else if (tolower(stat) == "dfbetas" || tolower(stat) == "restimates") {

      if (tolower(stat) == "dfbetas") {

        title = "DFbetas"
        ylabel = "BFbetas"
        obj <- object$DFbetas
        lab <- "DFbetas for treatment"

      } else if (tolower(stat) == "restimates") {

        title = "R statistic for Qestimates leave-one-out"
        ylabel = "R statistic for Qestimates"
        obj <- object$Restimates
        lab <- "R statistic for Qestimates of treatment"

      }

      xlabels <- factor(as.character(stlab), levels = as.character(stlab)) # as factor to prevent ggplot from reordering the x labels in alphabetical order
      nt <- length(unique(c(object$dat[, 4], object$dat[, 5])))

      graphs <- vector("list", nt - 1)


      upper <- max(obj)
      lower <- min(obj)

      limu <- round(max(upper), 2)
      liml <- round(min(lower), 2)


      for (j in 1:(nt - 1)) {

        ## Localizing variables - workaround for ggplot problem with
        ## handling variables in multiple plots (variable environment
        ## scope problem)
        local({

          melt_data <- melt(obj[j, ], id.vars = 0)
          y_values <- melt_data$value
          j <- j
          g <- eval(substitute(ggplot(data = melt_data, aes(x = xlabels, y = y_values)) + # eval(substitute) is another workaround for the aforementioned problem
                                 theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
                                 coord_cartesian(ylim = c(liml,limu)) +
                                 geom_point(color = '#016FB9', size = 2) +
                                 labs(y = paste(lab, rownames(obj)[j]), x = "Study deleted"), list(j = j)))

          graphs[[j]] <<- g
        })
      }

      ## creating a multiple graph consisting of a graph for each treatment
      grid.arrange(grobs = graphs, ncol = 3)
    }

  }

}
