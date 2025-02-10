#' Helper function to create the plots for FS algorithm.
#'
#' @param x object of class NMAoutlier (mandatory).
#' @param stat The monitored statistic to generate the plot.
#' @param select.st selected statistic (pscore/nsplit/estand) for
#'   selected treatment(s)/comparison(s)/study
#'
#' @details
#' Plot of several monitoring measures for FS algorithm.
#' Vertical axis provides the iterations of FS algorithm.
#' Horizontal axis provides a monitoring statistical measure in the methodology.
#'
#'
#' @keywords internal
#'
#' @author Maria Petropoulou <maria.petropoulou@uniklinik-freiburg.de>



plothelper <- function(x, stat, select.st){

    xlabel <- "Iterations"


  if (tolower(stat) == "pscore") {

    title <- "Forward plot for P-score"


    data <- getSelected(x$p.score, select.st)
    melt_data <- melt(data) # melt formats our data in a tall format which is proper for the ggplot function.
    var1_factors <- as.factor(melt_data$Var1)
    var2_factors <- as.factor(melt_data$Var2)

    labels_factors <- fixLabelsAndFactors(x, var2_factors)
    xlabels <- unlist(labels_factors[1])
    var2_factors <- unlist(labels_factors[2])

    if (length(var1_factors) == 0) {
      var1_factors <- as.factor(rownames(melt_data))
    }

    y_values <- melt_data$value
    ggplot(data = melt_data, aes(x = var2_factors, y = y_values, colour = var1_factors)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_line(aes(group = var1_factors, color = var1_factors, linetype = var1_factors), size = 1, na.rm = TRUE) +
      geom_point(aes(shape = var1_factors, color = var1_factors), size = 3, na.rm = TRUE) +
      labs(title = title, y = "P-Score", x = xlabel) +
      guides(colour = guide_legend("Treatments"), shape = guide_legend("Treatments"), linetype = guide_legend("Treatments")) +
      scale_x_discrete(labels = xlabels) +
      scale_shape_manual(values = seq(1,length(var1_factors))) +
      scale_linetype_manual(values = seq(1,length(var1_factors)))

  } else if (tolower(stat) == "nsplit") {


    title <- "Forward plot for difference of direct and indirect estimate (z-values)"

    data <- getSelected(x$dif, select.st)
    melt_data <- melt(data)
    var1_factors <- as.factor(melt_data$Var1)
    var2_factors <- as.factor(melt_data$Var2)

    labels_factors <- fixLabelsAndFactors(x, var2_factors)
    xlabels <- unlist(labels_factors[1])
    var2_factors <- unlist(labels_factors[2])

    if (length(var1_factors) == 0) {
      var1_factors <- as.factor(rownames(melt_data))
    }

    y_values <- melt_data$value
    ggplot(data = melt_data, aes(x = var2_factors, y = y_values, colour = var1_factors)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_line(aes(group = var1_factors, color = var1_factors, linetype = var1_factors), size = 1, na.rm = TRUE) +
      geom_point(aes(shape = var1_factors, color = var1_factors), size = 3, na.rm = TRUE) +
      labs(title = title, y = "Difference of direct and indirect estimate", x = xlabel) +
      guides(colour = guide_legend("Comparisons"), shape = guide_legend("Comparisons"), linetype = guide_legend("Comparisons")) +
      scale_x_discrete(labels = xlabels) +
      scale_shape_manual(values = seq(1,length(var1_factors))) +
      scale_linetype_manual(values = seq(1,length(var1_factors)))

  } else if (tolower(stat) == "estand") {


    title <- "Forward plot for standardized residuals"

    data <- getSelected(x$estand, select.st)
    melt_data <- melt(data)
    var1_factors <- as.factor(melt_data$Var1)
    var2_factors <- as.factor(melt_data$Var2)

    labels_factors <- fixLabelsAndFactors(x, var2_factors)
    xlabels <- unlist(labels_factors[1])
    var2_factors <- unlist(labels_factors[2])

    if (length(var1_factors) == 0) {
      var1_factors <- as.factor(rownames(melt_data))
    }
    y_values <- melt_data$value
    ggplot(data = melt_data, aes(x = var2_factors, y = y_values, colour = var1_factors)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_line(aes(group = var1_factors, color = var1_factors, linetype = var1_factors), size = 1, na.rm = TRUE) +
      geom_point(aes(shape = var1_factors, color = var1_factors), size = 3, na.rm = TRUE) +
      labs(title = title, y = "Standardized residuals", x = xlabel) +
      guides(colour = guide_legend("Studies"), shape = guide_legend("Studies"), linetype = guide_legend("Studies")) +
      scale_x_discrete(labels = xlabels) +
      scale_shape_manual(values = seq(1,length(var1_factors))) +
      scale_linetype_manual(values = seq(1,length(var1_factors)))

  } else if (tolower(stat) == "heterog") {


    title <- "Forward plot for heterogeneity"

    data <- getSelected(x$tau, select.st)
    melt_data <- melt(data)

    y_values <- melt_data$value
    ggplot(data = melt_data, aes(x = 1:length(x$tau), y = y_values)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      labs(title = title, y = "Heterogeneity", x = xlabel)

  } else if (tolower(stat) == "cook") {

    data <- getSelected(x$cook_d, select.st)
    melt_data <- melt(data)

    y_values <- melt_data$value
    ggplot(data = melt_data, aes(x = 1:length(data), y = y_values)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(title = "Forward plot for Cook's distance", y = "Cook's distance", x = xlabel) +
      scale_x_discrete(limits = as.factor(c(2:(length(data)+1))))

  } else if (tolower(stat) == "ratio") {

    data <- getSelected(x$Ratio, select.st)
    melt_data <- melt(data)

    y_values <- melt_data$value
    ggplot(data = melt_data, aes(x = 1:length(data), y = y_values)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(title = "Forward plot for ratio of variances", y = "Ratio of variances", x = xlabel) +
      scale_x_discrete(limits = as.factor(c(2:(length(data)+1))))

  } else if (tolower(stat) == "q") {

    # FS algorithm
    title1 <- "Forward plot for Qtotal"
    title2 <- "Forward plot for Qheterogeneity"
    title3 <- "Forward plot for Qinconsistency"





    data <- getSelected(x$Qb, select.st)

    melt_data1 <- melt(data)
    y_values1 <- melt_data1$value
    p1 <- ggplot(data = melt_data1, aes(x = 1:length(x$Qb), y = y_values1)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      labs(title = title1, y = "Qtotal", x = xlabel)

    data2 <- getSelected(x$Qhb, select.st)
    melt_data2 <- melt(data2)
    y_values2 <- melt_data2$value
    p2 <- ggplot(data = melt_data2, aes(x = 1:length(x$Qhb), y = y_values2)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      labs(title = title2, y = "Qheterogeneity", x = xlabel)

    data3 <- getSelected(x$Qib, select.st)
    melt_data3 <- melt(data3)
    y_values3 <- melt_data3$value
    p3 <- ggplot(data = melt_data3, aes(x = 1:length(x$Qib), y = y_values3)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      labs(title = title3, y = "Qinconsistency", x = xlabel)

    grid.arrange(p1, p2, p3, ncol = 3)
  } else {
    stop("The stat ", stat, " for the FS algorithm is not available")
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

fixLabelsAndFactors <- function(x, var2_factors) {

    xlabels <- factor(as.character(1:x$index), levels = as.character(1:x$index))
    if (length(var2_factors) == 0) {
      var2_factors <- as.factor(1:x$index)
    }

  return(list(xlabels, var2_factors))
}
