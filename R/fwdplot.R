#' Forward plot(s) to monitor selected statistic(s)/method(s)
#'
#' @description
#' This function generates forward plot(s) to monitor selected
#' statistic(s) and/or method(s).  The function creates a plot of the
#' selected statistic throughout the iterations of the forward search
#' algorithm.  Candidate statistics to be monitored can be P-score;
#' z-values by back-calculation method to derive indirect estimates
#' from direct pairwise comparisons and network estimates;
#' standardized residuals; heterogeneity variance estimator; Cook's
#' distance; ratio of variances; Q statistics (Overall heterogeneity /
#' inconsistency Q statistic (\code{Q}), overall heterogeneity Q
#' statistic (\code{Q}), between-designs Q statistic (\code{Q}), based
#' on a random effects design-by-treatment interaction model).
#'
#' @param x an object of class NMAoutlier (mandatory).
#' @param stat statistical measure to be monitored in forward plot(s)
#'   (mandatory), available choice are: "pscore", "nsplit", "estand",
#'   "heterog", "cook", "ratio", or "Q" (can be abbreviated).
#' @param select.st selected statistic (pscore/nsplit/estand) for
#'   selected treatment(s)/comparison(s)/study
#'
#' @details
#' Plot of statistical measures for each iteration of search.
#' Vertical axis provides iterations of search. Horizontal axis
#' provides a monitoring statistical measure.
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
#' # Forward search algorithm
#' #
#' FSresult <- NMAoutlier(p1, P = 1, small.values = "bad", n_cores = 2)
#'
#' # forward plot for Cook's distance
#' fwdplot(FSresult, "cook")
#'
#' \dontrun{
#' data(smokingcessation, package = "netmeta")
#'
#' # Transform data from arm-based format to contrast-based format
#' # We use 'sm' argument for odds ratios.
#' # We use function pairwise from netmeta package
#'
#' p1 <- netmeta::pairwise(list(treat1, treat2, treat3),
#'                         list(event1, event2, event3),
#'                         list(n1, n2, n3),
#'                         data=smokingcessation,
#'                         sm="OR")
#'
#' # forward search algorithm
#' FSresult <- NMAoutlier(p1, small.values = "bad")
#'
#' FSresult
#'
#' # forward plot for Cook's distance
#' fwdplot(FSresult, "cook")
#'
#' # forward plot for ratio of variances
#' fwdplot(FSresult, "ratio")
#'
#' # forward plot for heterogeneity variance estimator
#' fwdplot(FSresult, "heterog")
#'
#' # forward plot for Q statistics
#' fwdplot(FSresult, "Q")
#'
#' # forward plot for P-scores
#' fwdplot(FSresult, "pscore")
#'
#' # forward plot monitoring P-scores for treatment A
#' fwdplot(FSresult,"pscore", "A")
#'
#' # forward plot for z-values of disagreement of direct and indirect evidence
#' fwdplot(FSresult, "nsplit")
#'
#' # forward plot for z-values of disagreement of direct and indirect evidence
#' # monitoring treatment comparison A versus B
#' fwdplot(FSresult, "nsplit", "A:B")
#'
#' # forward plot for standardized residuals for study 4
#' fwdplot(FSresult, "estand", 4)
#' }
#'
#' @export
#'
#' @author Maria Petropoulou <mpetrop@cc.uoi.gr>
#'
#' @importFrom ggplot2 ggplot aes theme element_rect element_line
#'   geom_line geom_point labs guides guide_legend scale_x_discrete
#'   scale_linetype_manual scale_shape_manual scale_y_continuous
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange


fwdplot <- function (x, stat, select.st = NULL) {

  ## Check class
  ##
  chkclass(x, "NMAoutlier")

  stat <- setchar(stat, c("pscore", "nsplit", "estand", "heterog",
                          "cook", "ratio", "Q"))

  if (tolower(stat) == "pscore") {
    data<-getSelected(x$p.score, select.st)
    melt_data <- melt(data) # melt formats our data in a tall format which is proper for the ggplot function.
    var1_factors <- as.factor(melt_data$Var1)
    ggplot(data=melt_data, aes(x=melt_data$Var2, y=melt_data$value, colour=var1_factors)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_line(aes(group=var1_factors, color=var1_factors, linetype = var1_factors), size=1, na.rm=TRUE) +
      geom_point(aes(shape=var1_factors, color=var1_factors), size=3, na.rm=TRUE) +
      labs(title="Forward plot for P-score", y="P-Score", x="Iterations") +
      guides(colour = guide_legend("Treatments"), shape = guide_legend("Treatments"), linetype = guide_legend("Treatments")) +
      scale_x_discrete(labels=1:length(factor(melt_data$Var2))) +
      scale_shape_manual(values=seq(1,length(var1_factors))) +
      scale_linetype_manual(values=seq(1,length(var1_factors)))
  }
  else if (tolower(stat) == "nsplit") {
    data<-getSelected(x$dif, select.st)
    melt_data <- melt(data)
    var1_factors <- as.factor(melt_data$Var1)
    ggplot(data=melt_data, aes(x=melt_data$Var2, y=melt_data$value, colour=var1_factors)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_line(aes(group=var1_factors, color=var1_factors, linetype = var1_factors), size=1, na.rm=TRUE) +
      geom_point(aes(shape=var1_factors, color=var1_factors), size=3, na.rm=TRUE) +
      labs(title="Forward plot for difference of direct and indirect estimate (z-values)", y="Difference of direct and indirect estimate ", x="Iterations") +
      guides(colour = guide_legend("Comparisons"), shape = guide_legend("Comparisons"), linetype = guide_legend("Comparisons")) +
      scale_x_discrete(labels=1:length(factor(melt_data$Var2))) +
      scale_shape_manual(values=seq(1,length(var1_factors))) +
      scale_linetype_manual(values=seq(1,length(var1_factors)))
  }
  else if (tolower(stat) == "estand") {
    data<-getSelected(x$estand, select.st)
    melt_data <- melt(data)
    var1_factors <- as.factor(melt_data$Var1)
    ggplot(data=melt_data, aes(x=melt_data$Var2, y=melt_data$value, colour=var1_factors)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_line(aes(group=var1_factors, color=var1_factors, linetype = var1_factors), size=1, na.rm=TRUE) +
      geom_point(aes(shape=var1_factors, color=var1_factors), size=3, na.rm=TRUE) +
      labs(title="Forward plot for standardized residuals", y="Standardized residuals", x="Iterations") +
      guides(colour = guide_legend("Studies"), shape = guide_legend("Studies"), linetype = guide_legend("Studies")) +
      scale_x_discrete(labels=1:length(factor(melt_data$Var2))) +
      scale_shape_manual(values=seq(1,length(var1_factors))) +
      scale_linetype_manual(values=seq(1,length(var1_factors)))
  }
  else if (tolower(stat) == "heterog") {
    data<-getSelected(x$tau, select.st)
    melt_data <- melt(data)
    upper_limit <- NA
    if (any(melt_data$value == 0)) {
      upper_limit <- 1
    }
    ggplot(data=melt_data, aes(x=1:length(data), y=melt_data$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color='#016FB9', size=3, na.rm=TRUE) +
      labs(title="Forward plot for heterogeneity", y="Heterogeneity", x="Iterations") +
      scale_x_discrete(limits = c(1:length(data))) +
      scale_y_continuous(limits = c(0,upper_limit))
  }
  else if (tolower(stat) == "cook") {
    data<-getSelected(x$cook_d, select.st)
    melt_data <- melt(data)
    ggplot(data=melt_data, aes(x=2:(length(data)+1), y=melt_data$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color='#016FB9', size=3, na.rm=TRUE) +
      labs(title="Forward plot for Cook's distance", y="Cook's distance", x="Iterations") +
      scale_x_discrete(limits = c(2:(length(data)+1)))
  }
  else if (tolower(stat) == "ratio") {
    data<-getSelected(x$Ratio, select.st)
    melt_data <- melt(data)
    ggplot(data=melt_data, aes(x=2:(length(data)+1), y=melt_data$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color='#016FB9', size=3, na.rm=TRUE) +
      labs(title="Forward plot for ratio of variances", y="Ratio of variances", x="Iterations") +
      scale_x_discrete(limits = c(2:(length(data)+1)))
  }
  else if (tolower(stat) == "q") {
    data<-getSelected(x$Qb, select.st)
    melt_data1 <- melt(data)
    p1 <- ggplot(data=melt_data1, aes(x=1:length(x$Qb), y=melt_data1$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color='#016FB9', size=3, na.rm=TRUE) +
      labs(title="Forward plot for Qtotal", y="Qtotal", x="Iterations")

    data2<-getSelected(x$Qhb, select.st)
    melt_data2 <- melt(data2)
    p2 <- ggplot(data=melt_data2, aes(x=1:length(x$Qhb), y=melt_data2$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color='#016FB9', size=3, na.rm=TRUE) +
      labs(title="Forward plot for Qheterogeneity", y="Qheterogeneity", x="Iterations")

    data3<-getSelected(x$Qib, select.st)
    melt_data3 <- melt(data3)
    p3 <- ggplot(data=melt_data3, aes(x=1:length(x$Qib), y=melt_data3$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color='#016FB9', size=3, na.rm=TRUE) +
      labs(title="Forward plot for Qinconsistency", y="Qinconsistency", x="Iterations")

    grid.arrange(p1, p2, p3, ncol=3)
  }
}


getSelected <- function (dataSet, select.st) {

  if (!is.null(select.st)) {
    tr.names <- rownames(dataSet)
    select <- list()
    select <- match(select.st, tr.names)

    newDataSet <- as.matrix(t(dataSet[select[1],]))
    if(length(select)>1) {
      for (i in 2:length(select)){
        newDataSet <- rbind(newDataSet, dataSet[select[i],])
      }
    }
    toreturn <- as.matrix(newDataSet)
    rownames(toreturn) <- select.st
    return(toreturn)
  } else {
    return(dataSet)
  }
}
