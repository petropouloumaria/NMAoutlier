


print.NMAoutlier_measures <- function(x, digits = 4) {

  ## Check class
  ##
  chkclass(x,"NMAoutlier_measures")


  cat("Original data:\n")
  ##
  Mydata <- x$dat
  datamatrix <- cbind(formatN(as.numeric(Mydata[, 1]), digits), formatN(as.numeric(Mydata[, 2]), digits),
                      as.character(Mydata[, 3]), as.character(Mydata[, 4]), as.character(Mydata[, 5]))
  ##
  prmatrix(datamatrix, rowlab = paste(c(1:length(Mydata[, 1])), ""), collab = c("TE","seTE","studylab","treat1","treat2"),
           quote = FALSE, right = TRUE)

  cat("\n")
  ##
  ##
  cat("Simply outlier detection measures:\n")
  ##
  ##
  cat("\n")
  cat("Raw, standardized, Studentized residuals, Mahalanobis diastance and leverage for each included study in the network:\n")
  cat("\n")
  prmatrix(cbind(formatN(x$eraw, digits), formatN(x$estand, digits),
                 formatN(x$estud, digits), formatN(x$Mahalanobis.distance, digits), formatN(x$leverage, digits)),
           collab = c("Raw residual", "Standardized residual", "Studentized residual", "Mahalanobis distance", "leverage"),
           rowlab = c(unique(as.character(Mydata[, 3]))),
           quote = FALSE, right = TRUE)


  invisible(NULL)

}

