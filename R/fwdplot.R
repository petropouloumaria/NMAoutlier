#' Forward plot(s) to monitor selected statistic(s)/method(s).
#'
#' Forward plot(s) to monitor selected statistic(s) and/or method(s) (P-score, Node-splitting, Standardized residuals, heterogeneity, cook distance, ratio of variances, Q statistics).
#' @param x an object of class NMAoutlier (results from forward search algorithm in network meta-analysis).
#' @param stat measure(s) to be monitored in forward plot(s), available choices pscore, nsplit, estand, heterog, cook, ratio, Q.
#' @param select.st choice to monitor seleced statistic(s) (P-scores/node-splitting z-value/standardized resiluals) for selected item(s) (treatment(s)/comparison(s)/study(studies))
#' @return forward plot(s) for selected statistic(s)/method(s).
#' @export
#'
#' @author
#' Maria Petropoulou <mpetrop@cc.uoi.gr>
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics par matplot layout legend plot


fwdplot <- function (x, stat, select.st = "NULL") {

  ## Check class
  ##
  chkclass(x, "NMAoutlier")

  ##
  if (stat == "pscore") {

    dataSet <- x$p.score
    plotName <- "P-score"
    nameTitle <- "treatment"
    ##
    dev.new()
    ##
    ## Add extra space to right of plot area; change clipping to figure
    ##
    par(mar = c(5.1, 4.1, 4.1, 9), xpd = TRUE)
    ##
    tr.names <- rownames(dataSet)
    tr.n <- length(tr.names)
    ##
    if (is.character(select.st) && select.st != "NULL") {
      select <- list()
      select <- match(select.st, tr.names)
      ##
      dataSet <- dataSet[select,]
      ##
      tr.names <- select.st
      tr.n <- length(tr.names)
         if (length(select.st) == 1)
           dataSet <- t(dataSet)
    ##
    }
    ##
    matplot(t(dataSet), xlab = "iteration", ylab = plotName,
                      lwd = 2, lty = c(1:tr.n), col = c(1:tr.n), pch = 0:25 , type = "o",
                      main = "forward plot for p-score")
    ##
    legend("topright", inset = c(-0.3, 0), legend = tr.names,
                     ncol = 1, lty = c(1:tr.n), col = c(1:tr.n), pch = 0:25,
                     bg = ("white"), lwd = 2, horiz = F, cex = 0.7, title = nameTitle)

  }
  ##
  if (stat == "nsplit") {

    dataSet <- x$dif
    plotName <- "Node-splitting"
    nameTitle <- "comparison"
    ##
    dev.new()
    ##
    ## Add extra space to right of plot area; change clipping to figure
    ##
    par(mar = c(5.1, 4.1, 4.1, 11), xpd = TRUE)
    ##
    tr.names <- rownames(dataSet)
    tr.n <- length(tr.names)
    ##
    if (is.character(select.st) && select.st != "NULL") {
      select <- list()
      select <- match(select.st, tr.names)
      ##
      dataSet <- dataSet[select,]
      ##
      tr.names <- select.st
      tr.n <- length(tr.names)
      ##
      if (length(select.st) == 1)
        dataSet <- t(dataSet)
    ##
    }
    ##
    matplot(t(dataSet), xlab = "iteration", ylab = plotName,
                      lwd = 2, lty = c(1:tr.n), col = c(1:tr.n), pch = 0:25 , type = "o",
                      main ="forward plot for node-splitting (z-values)")
    ##
    legend("topright", inset = c(-0.5,0), legend = tr.names, ncol = 2,
                     lty = c(1:tr.n), col = c(1:tr.n), pch = 0:25, bg = ("white"),
                     lwd = 2, horiz = F, cex = 0.7, title = nameTitle)

  }
  ##
  if (stat == "estand") {

    dataSet <- x$estand
    plotName <- "Standardized residuals"
    nameTitle <- "study"
    ##
    dev.new()
    ##
    ## Add extra space to right of plot area; change clipping to figure
    ##
    par(mar = c(5.1, 4.1, 4.1, 9), xpd = TRUE)
    ##
    tr.names <- rownames(dataSet)
    tr.n <- length(tr.names)
    ##
    ##
    if (is.numeric(select.st)) {

      select <- list()
      select <- match(select.st, tr.names)
      ##
      dataSet <- dataSet[select,]
      ##
      tr.names <- select.st
      tr.n <- length(tr.names)
      ##
      if (length(select.st) == 1)
        dataSet <- t(dataSet)
    ##
    }
    ##
    matplot(t(dataSet), xlab = "iteration", ylab = plotName,
                      lwd = 2, lty = c(1:tr.n), col = c(1:tr.n), pch = 0:25,
                      type = "o", main = "forward plot for standardized residuals")
    ##
    legend("topright", inset = c(-0.3, 0), legend = tr.names, ncol = 2,
                     lty = c(1:tr.n), col = c(1:tr.n), pch = 0:25, bg = ("white"),
                     lwd = 2, horiz = F, cex = 0.7, title = nameTitle)

  }
  ##
  if (stat == "heterog") {
    ## forward plot for heterogeneity
    dev.new()
    ##
    plot(x$tau, xlab = "iteration", ylab = "heterogeneity",
                   main ="forward plot for heterogeneity")
  }
  ##
  if (stat == "cook") {
    ## forward plot for cook distance
    dev.new()
    ##
    plot(x$cook_d, xlab = "iteration", ylab = "cook distance",
                   main = "forward plot for cook distance")
  }
  ##
  if (stat == "ratio") {
    ## forward plot for ratio of variances
    dev.new()
    ##
    plot(x$Ratio, xlab = "iteration", ylab = "ratio of variances",
                   main = "forward plot for ratio of variances")
  }
  ##
  if (stat == "Q") {
    ## forward plot for Q statistics
    dev.new()
    ##
    layout(matrix(c(1:3), 1, byrow = TRUE))
    ##
    for (j in 1:3) {
      plot(list(x$Qb, x$Qhb, x$Qib)[[j]], xlab = "iteration",
                     ylab = c("Qtotal", "Qheterogeneity", "Qinconsistency")[j],
                     main = paste("forward plot for",
                     c("Qtotal", "Qheterogeneity", "Qinconsistency")[j]))
    }
  }
}


