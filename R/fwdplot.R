#' @title Forward plot(s) to monitor selected statistic(s)/method(s).
#'
#' @description This function generates forward plot(s) to monitor selected statistic(s) and/or method(s).
#' The function creates the plot of the selected statistic(s) for each iteration of forward search algorithm.
#' The selected statistic(s) to be monitored can be P-score; z-values by back-calculation method to derive indirect estimates
#' from direct pairwise comparisons and network estimates approach;
#' standardized residuals; heterogeneity variance estimator; cook distance; ratio of variances;
#' Q statistics (Overall heterogeneity / inconsistency Q statistic (\code{Q}), overall heterogeneity Q statistic (\code{Q}),
#' between-designs Q statistic (\code{Q}), based on a random effects design-by-treatment interaction model).
#'
#' @usage
#' fwdplot(x, stat, select.st = "NULL")
#'
#' @param x an object of class NMAoutlier (mandatory).
#' @param stat statistical measure(s) to be monitored in forward plot(s) (mandatory), available choices pscore; nsplit; estand; heterog; cook; ratio; Q.
#' @param select.st choice to monitor seleced statistic(s) (P-scores/z-values of disagreement of direct and indirect evidence/standardized resiluals) for selected item(s) (treatment(s)/comparison(s)/study(studies))
#' @return forward plot(s) for selected statistic(s)/method(s).
#'
#' @examples
#' \dontrun{
#' data(Dias2013)
#'
#' # forward search algorithm
#' FSresult <- NMAoutlier(TE, seTE, treat1, treat2,
#'                        studlab, data = Dias2013)
#'
#' FSresult
#'
#' # forward plot for cook distance
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
#' # forward plot for p-scores
#' fwdplot(FSresult, "pscore")
#'
#' # forward plot monitoring p-scores for treatment A
#' fwdplot(FSresult,"pscore", "A")
#'
#' # forward plot for z-values of disagreement of direct and indirect evidence
#' fwdplot(FSresult, "nsplit")
#'
#' # forward plot for z-values of disagreement of direct and indirect evidence
#' # monitoring treatment comparison A versus B
#' fwdplot(FSresult, "nsplit", "A:B")
#'
#' # forward plot for standardized residuals for study 3
#' fwdplot(FSresult, "estand", 3)
#'
#' }
#'
#' @export
#'
#' @author
#' Maria Petropoulou <mpetrop@cc.uoi.gr>
#'
#' @importFrom ggplot2 ggplot aes theme element_rect element_line geom_line geom_point labs guides guide_legend scale_x_discrete scale_linetype_manual scale_shape_manual
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange

fwdplot <- function (x, stat, select.st = NULL) {

  ## Check class
  ##
  chkclass(x, "NMAoutlier")

  if(stat == "pscore") {
    data<-getSelected(x$p.score, select.st)
    melt_data <- melt(data) ## melt formats our data in a tall format which is proper for the ggplot function.
    var1_factors <- as.factor(melt_data$Var1)
    ggplot(data=melt_data, aes(x=melt_data$Var2, y=melt_data$value, colour=var1_factors)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_line(aes(group=var1_factors, color=var1_factors, linetype = var1_factors), size=1, na.rm=TRUE) +
      geom_point(aes(shape=var1_factors, color=var1_factors), size=3, na.rm=TRUE) +
      labs(title="forward plot for P-score", y="P-Score", x="Iterations") +
      guides(colour = guide_legend("Treatments"), shape = guide_legend("Treatments"), linetype = guide_legend("Treatments")) +
      scale_x_discrete(labels=1:length(factor(melt_data$Var2))) +
      scale_shape_manual(values=seq(1,length(var1_factors))) +
      scale_linetype_manual(values=seq(1,length(var1_factors)))
  }
  else if (stat == "nsplit") {
    data<-getSelected(x$dif, select.st)
    melt_data <- melt(data)
    var1_factors <- as.factor(melt_data$Var1)
    ggplot(data=melt_data, aes(x=melt_data$Var2, y=melt_data$value, colour=var1_factors)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_line(aes(group=var1_factors, color=var1_factors, linetype = var1_factors), size=1, na.rm=TRUE) +
      geom_point(aes(shape=var1_factors, color=var1_factors), size=3, na.rm=TRUE) +
      labs(title="forward plot for node-splitting (z-values)", y="Node-splitting", x="Iterations") +
      guides(colour = guide_legend("Comparisons"), shape = guide_legend("Comparisons"), linetype = guide_legend("Comparisons")) +
      scale_x_discrete(labels=1:length(factor(melt_data$Var2))) +
      scale_shape_manual(values=seq(1,length(var1_factors))) +
      scale_linetype_manual(values=seq(1,length(var1_factors)))
  }
  else if (stat == "estand") {
    data<-getSelected(x$estand, select.st)
    melt_data <- melt(data)
    var1_factors <- as.factor(melt_data$Var1)
    ggplot(data=melt_data, aes(x=melt_data$Var2, y=melt_data$value, colour=var1_factors)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_line(aes(group=var1_factors, color=var1_factors, linetype = var1_factors), size=1, na.rm=TRUE) +
      geom_point(aes(shape=var1_factors, color=var1_factors), size=3, na.rm=TRUE) +
      labs(title="forward plot for standardized residuals", y="Standardized residuals", x="Iterations") +
      guides(colour = guide_legend("Studies"), shape = guide_legend("Studies"), linetype = guide_legend("Studies")) +
      scale_x_discrete(labels=1:length(factor(melt_data$Var2))) +
      scale_shape_manual(values=seq(1,length(var1_factors))) +
      scale_linetype_manual(values=seq(1,length(var1_factors)))
  }
  else if (stat == "heterog") {
    data<-getSelected(x$tau, select.st)
    melt_data <- melt(data)
    ggplot(data=melt_data, aes(x=1:length(x$tau), y=melt_data$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color='#016FB9', size=3, na.rm=TRUE) +
      labs(title="forward plot for heterogeneity", y="Heterogeneity", x="Iterations")
  }
  else if (stat == "cook") {
    data<-getSelected(x$cook_d, select.st)
    melt_data <- melt(data)
    ggplot(data=melt_data, aes(x=1:length(x$cook_d), y=melt_data$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color='#016FB9', size=3, na.rm=TRUE) +
      labs(title="forward plot for cook distance", y="Cook distance", x="Iterations")
  }
  else if (stat == "ratio") {
    data<-getSelected(x$Ratio, select.st)
    melt_data <- melt(data)
    ggplot(data=melt_data, aes(x=1:length(x$Ratio), y=melt_data$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color='#016FB9', size=3, na.rm=TRUE) +
      labs(title="forward plot for ratio of variances", y="Ratio of variances", x="Iterations")
  }
  else if (stat == "Q") {
    data<-getSelected(x$Qb, select.st)
    melt_data1 <- melt(data)
    p1 <- ggplot(data=melt_data1, aes(x=1:length(x$Qb), y=melt_data1$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color='#016FB9', size=3, na.rm=TRUE) +
      labs(title="forward plot for Qtotal", y="Qtotal", x="Iterations")

    data2<-getSelected(x$Qhb, select.st)
    melt_data2 <- melt(data2)
    p2 <- ggplot(data=melt_data2, aes(x=1:length(x$Qhb), y=melt_data2$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color='#016FB9', size=3, na.rm=TRUE) +
      labs(title="forward plot for Qheterogeneity", y="Qheterogeneity", x="Iterations")

    data3<-getSelected(x$Qib, select.st)
    melt_data3 <- melt(data3)
    p3 <- ggplot(data=melt_data3, aes(x=1:length(x$Qib), y=melt_data3$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color='#016FB9', size=3, na.rm=TRUE) +
      labs(title="forward plot for Qinconsistency", y="Qinconsistency", x="Iterations")

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


