

print.NMAsvr <- function(x, digits = 4) {

  ## Check class
  ##
  chkclass(x,"NMAsvr")


  cat("Original data:\n")
  ##
  Mydata <- x$dat
  datamatrix<-cbind(formatN(as.numeric(Mydata[, 1]), digits), formatN(as.numeric(Mydata[, 2]), digits),
                    as.character(Mydata[, 3]), as.character(Mydata[, 4]), as.character(Mydata[, 5]))
  ##
  prmatrix(datamatrix, rowlab = paste(c(1:length(Mydata[, 1])), ""), collab = c("TE","seTE","studylab","treat1","treat2"),
           quote = FALSE, right = TRUE)

  cat("\n")
  ##
  ##
  cat("Results with Standard Random Variance Model:\n")
  ##
  ##
  cat("\n")
  cat("Summary estimates with their 95% confidence intervals (ci) and heterogeneity from Standard Random Variance Model:\n")
  cat("\n")
  cat("Summary estimates from Standard Random Variance Model:\n")
  prmatrix(cbind(formatN(x$b_standard, digits), formatN(x$l_standard, digits),
                 formatN(x$u_standard, digits)),
           collab = c("Summary estimates", "Lower ci", "Upper ci"),
           quote = FALSE, right = TRUE)
  cat("\n")
  cat("Heterogeneity maximum estimator from Standard Random Variance Model:\n", formatN(x$tau_standard, digits))
  cat("\n")
  ##
  ##
  cat("\n")
  ##
  cat("Results with Random Shift Variance Model:\n")
  ##
  ##
  cat("\n")
  cat("Summary estimates and their 95% confidence intervals from Random Shift Variance Model:\n")
  cat("\n")
  cat("Summary estimates from Random Shift Variance Model:\n")
  prmatrix(formatN(x$b, digits), quote = FALSE, right = TRUE)
  ##
  cat("Lower confidence interval from Random Shift Variance Model:\n")
  prmatrix(formatN(x$l, digits), quote = FALSE, right = TRUE)
  ##
  ##
  cat("Upper confidence interval from Random Shift Variance Model:\n")
  prmatrix(formatN(x$u, digits), quote = FALSE, right = TRUE)
  ##
  ##
  cat("\n")
  cat("Variance estimator of Shift Variance Model:\n")
  cat("\n")
  prmatrix(cbind(formatN(x$tau, digits), formatN(x$over_disp, digits)),
  rowlab = paste("study=", iteration = x$z[!is.na(x$z)]),
  collab = c("heterogeneity", "over-dispresion"),
  quote = FALSE, right = TRUE)

  ##
  ##
  cat("\n")
  cat("Likelihood statistics:\n")
  prmatrix(cbind(formatN(c(x$twiceloglik),digits), formatN(x$converge, digits), formatN(c(x$LRT),digits)),
           collab = c("2*(maximum log-likelihood)", "Convergence diagnostic", "Likelihood Ratio Test (LRT)"),
           rowlab = paste("study=", iteration = x$z[!is.na(x$z)]),
           quote = FALSE, right = TRUE)

  cat("\n")
  cat("Monitored measures for Random Shift Variance Model for selected study:\n")
  cat("Heterogeneity and inconsistency measures:\n")
  ##
  ##
  cat("\n")
  cat("Q statistics and heterogeneity:\n")
  prmatrix(cbind(formatN(x$Q), formatN(x$Qhet), formatN(x$Qinc), formatN(x$tau, digits)),
           rowlab = paste("study=", iteration = x$z[!is.na(x$z)]),
           collab = c("Qtotal", "Qheterogeneity", "Qinconsistency", "heterogeneity"),
           quote = FALSE, right = TRUE)
  ##
  ##
  cat("\n")
  cat("Z-values from difference of direct and indirect evidence (Node-splitting)\n")
  prmatrix(formatN(t(x$dif), digits), quote = FALSE, right = TRUE)

  ##
  ##
  cat("\n")
  cat("Ranking measures, P-score for each treatment:\n")
  prmatrix((formatN(x$p.score, digits)), quote = FALSE, right = TRUE)


  invisible(NULL)

}

