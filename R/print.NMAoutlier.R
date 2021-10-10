#' @method print NMAoutlier
#' @export


print.NMAoutlier <- function(x, digits = 4, ...) {

  ## Check class
  ##
  chkclass(x, "NMAoutlier")


  cat("Original data:\n")
  ##
  Mydata <- x$dat
  datamatrix <- cbind(formatN(as.numeric(Mydata[, 1]), digits),
                      formatN(as.numeric(Mydata[, 2]), digits),
                      as.character(Mydata[, 3]),
                      as.character(Mydata[, 4]),
                      as.character(Mydata[, 5]))
  ##
  prmatrix(datamatrix,
           rowlab = paste(c(1:length(Mydata[, 1])), ""),
           collab = c("TE", "seTE", "studylab", "treat1", "treat2"),
           quote = FALSE, right = TRUE)
  ##
  cat("\n")
  cat("Length of the initial clean dataset:", x$length.initial, "studies\n")
  cat("Number of forward search iterations:", x$index ,"\n\n")
  cat(paste("Study entered into the basic set in each step of",
            "forward search algorithm:\n"))
  prmatrix(x$basic)
  ##
  cat("\n")
  cat(paste("Monitored measures of basic set in each step of",
            "forward search algorithm:\n"))
  cat("Heterogeneity and inconsistency measures:\n")
  ##
  cat("\n")
  cat("Q statistics and heterogeneity:\n")
  prmatrix(cbind(formatN(x$Qb),
                 formatN(x$Qhb),
                 formatN(x$Qib),
                 formatN(x$taub, digits)),
           rowlab = paste("it=", iteration = c(1:x$index)),
           collab = c("Qtotal", "Qheterogeneity", "Qinconsistency",
                      "heterogeneity"),
           quote = FALSE, right = TRUE)
  ##
  cat("\n")
  cat(paste("Z-values from difference of direct and",
            "indirect evidence (Node-splitting)\n"))
  prmatrix(formatN(t(x$dif), digits), quote = FALSE, right = TRUE)
  ##
  cat("\n")
  cat("Outlying measures:\n")
  prmatrix(cbind(formatN(x$Ratio,digits),
                 formatN(x$cook_d, digits)),
           rowlab = paste("it=", iteration = c(1:x$index)),
           collab = c("Ratio of variances", "Cook's distance"),
           quote = FALSE, right = TRUE)
  ##
  cat("\n")
  cat("Ranking measures, P-score for each treatment:\n")
  prmatrix((formatN(x$p.score, digits)), quote = FALSE, right = TRUE)
  
  invisible(NULL)

}
