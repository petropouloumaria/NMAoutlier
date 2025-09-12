#' @method print NMAoutlier.measures
#' @export


print.NMAoutlier.measures <- function(x, digits = 4, details = TRUE, ...) {

  ## Check class and arguments
  ##
  chkclass(x, "NMAoutlier.measures")
  chknumeric(digits, min = 0, single = TRUE)
  chklogical(details)
  ##
  Mydata <- x$dat


  if (details) {
    cat("Original data:\n")
    ##
    datamatrix <- cbind(formatN(as.numeric(Mydata[, 1]), digits),
                        formatN(as.numeric(Mydata[, 2]), digits),
                        as.character(Mydata[, 3]),
                        as.character(Mydata[, 4]),
                        as.character(Mydata[, 5]),
                        as.character(Mydata[, 6]))
    ##
    prmatrix(datamatrix,
             rowlab = paste(c(1:length(Mydata[, 1])), ""),
             collab = c("TE", "seTE", "studylab.orig", "studylab", "treat1", "treat2"),
             quote = FALSE, right = TRUE)
  }


  ids  <- as.numeric(x$dat[, 4])      # studylab (numeric id)
  labs <- as.character(x$dat[, 3])    # studylab.orig (label)

  dfu <- data.frame(id = ids, lab = labs, stringsAsFactors = FALSE)
  dfu <- dfu[!duplicated(dfu$id), ]   # unique per study id

  labs_in_order <- dfu$lab[ match(x$stud_id, dfu$id) ]


  if (x$measure == "simple") {
    ##
    cat("\n")
    ##
    cat("Simple outlier detection measures (for each included study):\n")
    ##
    cat("\n")
    prmatrix(cbind(x$stud_id, formatN(x$eraw, digits),
                   formatN(x$estand, digits),
                   formatN(x$estud, digits),
                   formatN(x$Mahalanobis.distance, digits),
                   formatN(x$leverage, digits)),
             collab = c("studylab", "raw", "standardized",
                        "studentized", "mahalanobis",
                        "leverage"),
             rowlab = c(as.character(labs_in_order)),
             quote = FALSE, right = TRUE)

  }
  else if (x$measure == "deletion") {
    ##
    cat("\n")
    ##
    cat("Outlier detection measures considering deletion:\n")
    ##
    cat("\n")
    cat(paste("Raw, standardized, Studentized deleted residuals for",
              "each included study in the network:\n"))
    cat("\n")
    prmatrix(cbind(c(x$stud_id),
                   formatN(x$eraw.deleted, digits),
                   formatN(x$estand.deleted, digits),
                   formatN(x$estud.deleted, digits)),
             collab = c("studylab", "Raw deleted residual",
                        "Standardized deleted residual",
                        "Studentized deleted residual"),
             rowlab = c(as.character(labs_in_order)),
             quote = FALSE, right = TRUE)
    ##
    cat("\n")
    cat(paste("'Leave one out' measures for each included study",
              "in the network:\n"))
    cat("\n")
    prmatrix(cbind(c(x$stud_id),
                   formatN(x$w.leaveoneout, digits),
                   formatN(x$H.leaveoneout, digits),
                   formatN(x$heterog.leaveoneout, digits)),
             collab = c("studylab", "Weight 'leave one out'",
                        "Leverage 'leave one out'",
                        "heterogeneity 'leave one out'"),
             rowlab = c(as.character(labs_in_order)),
             quote = FALSE, right = TRUE)
    ##
    cat("\n")
    cat(paste("Cook's distance and COVRATIO considered deletion of",
              "study in the network:\n"))
    cat("\n")
    prmatrix(cbind(c(x$stud_id),
                   formatN(x$Cooks.distance, digits),
                   formatN(x$Covratio, digits)),
             collab = c("studylab", "Cook's distance", "COVRATIO"),
             rowlab = c(unique(labs_in_order)),
             quote = FALSE, right = TRUE)
    ##
    cat("\n")
    cat("R statistics considered deletion of study in the network:\n")
    cat("\n")
    prmatrix(cbind(c(x$stud_id),
                   formatN(x$Rheterogeneity, digits),
                   formatN(x$RQtotal, digits),
                   formatN(x$RQhet, digits),
                   formatN(x$RQinc, digits)),
             collab = c("studylab",
                        "R statistic for heterogeneity",
                        "R statistic for Qtotal",
                        "R statistic for Qheterogeneity",
                        "R statistic for Qinconsistency"),
             rowlab = c(unique(labs_in_order)),
             quote = FALSE, right = TRUE)
    ##
    cat("\n")
    cat(paste("DFbetas for each treatment versus the reference",
              "considered deletion of study in the network:\n"))
    cat("\n")
    prmatrix(formatN(t(x$DFbetas), digits),
             rowlab = c(unique(labs_in_order)),
             collab = rownames(x$DFbetas),
             quote = FALSE, right = TRUE)
  }

  invisible(NULL)

}
