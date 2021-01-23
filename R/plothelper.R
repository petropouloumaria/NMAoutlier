#' Helper function to create the plots for FS and RVSOM methodologies.
#'
#' @param x object of class NMAoutlier or class NMAoutlier.rsv (mandatory).
#' @param method The method that was used.
#' Select "fs" to monitor statistics during the forward search algorithm.
#' Select "rsv" to monitor statistics by fitting the data with the Shift Variance Model.
#' @param stat The monitored statistic to generate the plot.
#' @param select.st selected statistic (pscore/nsplit/estand) for
#'   selected treatment(s)/comparison(s)/study
#'
#' @details
#' Plot of several monitoring measures for FS and RVSOM methodologies.
#' Vertical axis provides the iterations of FS methodology or the study for RVSOM methodology.
#' Horizontal axis provides a monitoring statistical measure in the methodology.
#'
#'
#' @keywords internal
#'
#' @author Maria Petropoulou <petropoulou@imbi.uni-freiburg.de>
#'
#' @importFrom ggplot2 ggplot aes theme element_rect element_line
#'   geom_line geom_point labs guides guide_legend scale_x_discrete
#'   scale_linetype_manual scale_shape_manual
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange



plothelper <- function(x, method, stat, select.st){

  if (tolower(method) == "fs") {
    xlabel <- "Iterations"
  } else {
    # RVSOM
    xlabel <- "Study"
  }

  if (tolower(stat) == "pscore") {

    if (tolower(method) == "fs") {
      title <- "Forward plot for P-score"
    } else {
      # RVSOM
      title <- "P-score for Random Shift Variance Model"
    }

    data <- getSelected(x$p.score, select.st)
    melt_data <- melt(data) # melt formats our data in a tall format which is proper for the ggplot function.
    var1_factors <- as.factor(melt_data$Var1)
    var2_factors <- as.factor(melt_data$Var2)

    labels_factors <- fixLabelsAndFactors(x, method, var2_factors)
    xlabels <- unlist(labels_factors[1])
    var2_factors <- unlist(labels_factors[2])

    if (length(var1_factors) == 0) {
      var1_factors <- as.factor(rownames(melt_data))
    }

    ggplot(data = melt_data, aes(x = var2_factors, y = value, colour = var1_factors)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_line(aes(group = var1_factors, color = var1_factors, linetype = var1_factors), size = 1, na.rm = TRUE) +
      geom_point(aes(shape = var1_factors, color = var1_factors), size = 3, na.rm = TRUE) +
      labs(title = title, y = "P-Score", x = xlabel) +
      guides(colour = guide_legend("Treatments"), shape = guide_legend("Treatments"), linetype = guide_legend("Treatments")) +
      scale_x_discrete(labels = xlabels) +
      scale_shape_manual(values = seq(1,length(var1_factors))) +
      scale_linetype_manual(values = seq(1,length(var1_factors)))

  } else if (tolower(stat) == "nsplit") {

    if (tolower(method) == "fs") {
      title <- "Forward plot for difference of direct and indirect estimate (z-values)"
    } else {
      # RVSOM
      title <- "Difference of direct and indirect estimate (z-values) for Random Shift Variance Model"
    }

    data <- getSelected(x$dif, select.st)
    melt_data <- melt(data)
    var1_factors <- as.factor(melt_data$Var1)
    var2_factors <- as.factor(melt_data$Var2)

    labels_factors <- fixLabelsAndFactors(x, method, var2_factors)
    xlabels <- unlist(labels_factors[1])
    var2_factors <- unlist(labels_factors[2])

    if (length(var1_factors) == 0) {
      var1_factors <- as.factor(rownames(melt_data))
    }

    ggplot(data = melt_data, aes(x = var2_factors, y = value, colour = var1_factors)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_line(aes(group = var1_factors, color = var1_factors, linetype = var1_factors), size = 1, na.rm = TRUE) +
      geom_point(aes(shape = var1_factors, color = var1_factors), size = 3, na.rm = TRUE) +
      labs(title = title, y = "Difference of direct and indirect estimate", x = xlabel) +
      guides(colour = guide_legend("Comparisons"), shape = guide_legend("Comparisons"), linetype = guide_legend("Comparisons")) +
      scale_x_discrete(labels = xlabels) +
      scale_shape_manual(values = seq(1,length(var1_factors))) +
      scale_linetype_manual(values = seq(1,length(var1_factors)))

  } else if (tolower(stat) == "estand") {

    if (tolower(method) == "fs") {
      title <- "Forward plot for standardized residuals"
    } else {
      # RVSOM
      title <- "Standardized residuals for Random Shift Variance Model"
    }

    data <- getSelected(x$estand, select.st)
    melt_data <- melt(data)
    var1_factors <- as.factor(melt_data$Var1)
    var2_factors <- as.factor(melt_data$Var2)

    labels_factors <- fixLabelsAndFactors(x, method, var2_factors)
    xlabels <- unlist(labels_factors[1])
    var2_factors <- unlist(labels_factors[2])

    if (length(var1_factors) == 0) {
      var1_factors <- as.factor(rownames(melt_data))
    }

    ggplot(data = melt_data, aes(x = var2_factors, y = value, colour = var1_factors)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_line(aes(group = var1_factors, color = var1_factors, linetype = var1_factors), size = 1, na.rm = TRUE) +
      geom_point(aes(shape = var1_factors, color = var1_factors), size = 3, na.rm = TRUE) +
      labs(title = title, y = "Standardized residuals", x = xlabel) +
      guides(colour = guide_legend("Studies"), shape = guide_legend("Studies"), linetype = guide_legend("Studies")) +
      scale_x_discrete(labels = xlabels) +
      scale_shape_manual(values = seq(1,length(var1_factors))) +
      scale_linetype_manual(values = seq(1,length(var1_factors)))

  } else if (tolower(stat) == "heterog") {

    if (tolower(method) == "fs") {
      title <- "Forward plot for heterogeneity"
    } else {
      # RVSOM
      title <- "Heterogeneity for Random Shift Variance Model"
    }

    data <- getSelected(x$tau, select.st)
    melt_data <- melt(data)
    ggplot(data = melt_data, aes(x = 1:length(x$tau), y = value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      labs(title = title, y = "Heterogeneity", x = xlabel)

  } else if (tolower(stat) == "cook" & tolower(method) == "fs") {

    data <- getSelected(x$cook_d, select.st)
    melt_data <- melt(data)
    ggplot(data = melt_data, aes(x = 2:(length(data)+1), y = value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(title = "Forward plot for Cook's distance", y = "Cook's distance", x = xlabel) +
      scale_x_discrete(limits = as.factor(c(2:(length(data)+1))))

  } else if (tolower(stat) == "ratio" & tolower(method) == "fs") {

    data <- getSelected(x$Ratio, select.st)
    melt_data <- melt(data)
    ggplot(data = melt_data, aes(x = 2:(length(data)+1), y = value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(title = "Forward plot for ratio of variances", y = "Ratio of variances", x = xlabel) +
      scale_x_discrete(limits = as.factor(c(2:(length(data)+1))))

  } else if (tolower(stat) == "over_disp" & tolower(method) == "rsv") {

    data <- getSelected(x$over_disp, select.st)
    melt_data <- melt(data)
    ggplot(data = melt_data, aes(x = 1:length(x$over_disp), y = value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      labs(title = "Shift variance estimator for Random Shift Variance Model", y = "Shift variance estimator", x = xlabel)

  } else if (tolower(stat) == "lrt" & tolower(method) == "rsv") {

    data <- getSelected(x$LRT, select.st)
    melt_data <- melt(data)
    ggplot(data = melt_data, aes(x = 1:length(x$LRT), y = value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      labs(title = "Likelihood Ratio Test (LRT)", y = "Likelihood Ratio Test (LRT)", x = xlabel)

  } else if (tolower(stat) == "leverage" & tolower(method) == "rsv") {

    data <- getSelected(x$leverage, select.st)
    melt_data <- melt(data)
    var1_factors <- as.factor(melt_data$Var1)
    var2_factors <- as.factor(melt_data$Var2)

    labels_factors <- fixLabelsAndFactors(x, method, var2_factors)
    xlabels <- unlist(labels_factors[1])
    var2_factors <- unlist(labels_factors[2])

    if (length(var1_factors) == 0) {
      var1_factors <- as.factor(rownames(melt_data))
    }

    ggplot(data = melt_data, aes(x = var2_factors, y = value, colour = var1_factors)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      #geom_line(aes(group=var1_factors, color=var1_factors, linetype = var1_factors), size=1, na.rm=TRUE) +
      geom_point(aes(shape = var1_factors, color = var1_factors), size = 3, na.rm = TRUE) +
      labs(title = "Leverage for each pairwise comparison in Shift Variance Model", y = "Leverage", x = "Study") +
      guides(colour = guide_legend("Studies"), shape = guide_legend("Studies"), linetype = guide_legend("Studies")) +
      scale_x_discrete(labels = xlabels) +
      scale_shape_manual(values = seq(1,length(var1_factors))) +
      scale_linetype_manual(values = seq(1,length(var1_factors)))

  } else if (tolower(stat) == "q") {

    if (tolower(method) == "fs") {
      # FS algorithm
      title1 <- "Forward plot for Qtotal"
      title2 <- "Forward plot for Qheterogeneity"
      title3 <- "Forward plot for Qinconsistency"


    } else {
      # RVSOM
      title1 <- "Qtotal for Random Shift Variance Model"
      title2 <- "Qheterogeneity for Random Shift Variance Model"
      title3 <- "Qinconsistency for Random Shift Variance Model"

      x$Qb <- x$Q
      x$Qhb <- x$Qhet
      x$Qib <- x$Qinc
    }


    data <- getSelected(x$Qb, select.st)

    melt_data1 <- melt(data)
    p1 <- ggplot(data = melt_data1, aes(x = 1:length(x$Qb), y = melt_data1$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      labs(title = title1, y = "Qtotal", x = xlabel)

    data2 <- getSelected(x$Qhb, select.st)
    melt_data2 <- melt(data2)
    p2 <- ggplot(data = melt_data2, aes(x = 1:length(x$Qhb), y = melt_data2$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      labs(title = title2, y = "Qheterogeneity", x = xlabel)

    data3 <- getSelected(x$Qib, select.st)
    melt_data3 <- melt(data3)
    p3 <- ggplot(data = melt_data3, aes(x = 1:length(x$Qib), y = melt_data3$value)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      labs(title = title3, y = "Qinconsistency", x = xlabel)

    grid.arrange(p1, p2, p3, ncol = 3)
  } else {
    stop("The stat ", stat, " for the method ", method, " is not available")
  }
}

getSelected <- function(dataSet, select.st) {

  if (!is.null(select.st)) {
    tr.names <- rownames(dataSet)
    select <- list()
    select <- match(select.st, tr.names)

    newDataSet <- as.matrix(t(dataSet[select[1],]))
    if (length(select) > 1) {
      for (i in 2:length(select)) {
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

fixLabelsAndFactors <- function(x, method, var2_factors) {
  if (tolower(method) == "rsv") {
    xlabels <- factor(as.character(x$z[-1]), levels = as.character(x$z[-1])) # as factor to prevent ggplot from reordering the x labels in alphabetical order
    if (length(var2_factors) == 0) {
      var2_factors <- as.factor(x$z[-1])
    }
  } else {
    xlabels <- factor(as.character(1:x$index), levels = as.character(1:x$index))
    if (length(var2_factors) == 0) {
      var2_factors <- as.factor(1:x$index)
    }
  }

  return(list(xlabels, var2_factors))
}
