#' Outlier and influential detection measures in network meta-analysis.
#'
#' @description
#' Employs the computation of several measures for detection of outlying
#' studies (studies with extreme results) and influential studies fitted in network
#' meta-analysis model from graph-theory. Calculation of these measures for each study
#' included in the network can detect evidence of outliers. It can also be used to detect
#' studies that are potential sources for heterogeneity and
#' inconsistency.
#'
#' Statistical outlier and influential measures are:
#' \itemize{
#' \item outlying and influential measures for each study (raw residuals, standardized residuals,
#'  studentized residuals, Mahalanobis distance, leverage).
#' }
#'
#' A description of the several outlier detection measures in the context of network meta-analysis
#' can be found in Petropoulou et al. (2019).
#'
#' @param TE Estimate of treatment effect, i.e. difference between
#'   first and second treatment (e.g. log odds ratio, mean difference,
#'   or log hazard ratio).
#' @param seTE Standard error of treatment estimate.
#' @param treat1 Label/Number for first treatment.
#' @param treat2 Label/Number for second treatment.
#' @param studlab Study labels (important when multi arm studies are
#'   included).
#' @param data A data frame containing the study information.
#' @param sm A character string indicating underlying summary measure,
#'   e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
#'   \code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param reference Reference treatment group.
#'
#' @details
#' Description of several outlier and influential measures are calculated for
#' the network meta-analysis model from graph theory (Rücker, 2012)
#' fitted with (\code{netmeta} function) of R package \bold{netmeta} (Rücker et al., 2015).
#' The researcher can choose the reference treatment \code{reference} fitted in NMA model.
#'
#' An overview of the several outlier detection measures is described in Petropoulou et al. 2019.
#'
#' Let \emph{n} be the number of treatments in a network and let
#' \emph{m} be the number of pairwise treatment comparisons.  If there
#' are only two-arm studies, \emph{m} is the number of studies. Let
#' \code{TE} and \code{seTE} be the vectors of observed effects and their standard
#' errors.  Comparisons belonging to multi-arm studies are identified
#' by identical study labels (argument \code{studlab}). It is
#' therefore important to use identical study labels for all
#' comparisons belonging to the same multi-arm study.
#'
#' The function calculates several outlier detection measures for each study.
#' The statistical measures calculated are:
#' Raw residuals, Standardized residuals, Studentized residuals, Mahalanobis distance
#' and leverage for each study.
#'
#' @return
#' An object of class \code{NMAoutlier_measures}; a list containing the
#' following components:
#'    \item{dat}{Matrix containing the data \code{"TE"}, \code{"seTE"}, \code{"studlab"}, \code{"treat1"}, \code{"treat2"} as defined above.}
#'    \item{eraw}{Raw residual for each study included in the network.}
#'    \item{estand}{Standardized residual for each study included in the network.}
#'    \item{estud}{Studentized residual for each study included in the network.}
#'    \item{Mah}{Mahalanobis distance for each study included in the network.}
#'    \item{leverage}{Leverage for each study included in the network.}
#'    \item{call}{Function call}
#'
#' @references
#' Rücker G (2012):
#' Network meta-analysis, electrical networks and graph theory.
#' \emph{Research Synthesis Methods},
#' \bold{3}, 312--24
#'
#' Rücker G, Schwarzer G (2015):
#' Ranking treatments in frequentist network meta-analysis works
#' without resampling methods.
#' \emph{BMC Medical Research Methodology},
#' \bold{15}, 58
#'
#' Petropoulou M (2019):
#' Outlier detection measures in network meta-analysis.
#' \emph{Manuscript}
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
#'
#' #' # outlier and influential detection measures for studies 9, 10, 11, 12
#' measures <- NMAoutlier_measures(p1)
#'
#' # Standardized residual for each study included in the network
#' measures$estand
#'
#'
#' \dontrun{
#' data(smokingcessation, package = "netmeta")
#'
#' # Transform data from arm-based format to contrast-based format
#' # We use 'sm' argument for odds ratios.
#' # We use function pairwise from netmeta package
#' #
#' p1 <- netmeta::pairwise(list(treat1, treat2, treat3),
#'                         list(event1, event2, event3),
#'                         list(n1, n2, n3),
#'                         data = smokingcessation,
#'                         sm = "OR")
#'
#' # outlier and influential detection measures for each study in the network
#' measures <- NMAoutlier_measures(p1)
#'
#' # Mahalanobis distance for each study included in the network
#' measures$Mah
#' }
#'
#' @export
#'
#' @author Maria Petropoulou <mpetrop@cc.uoi.gr>
#'
#' @importFrom netmeta netmeta



NMAoutlier_measures <- function(TE, seTE, treat1, treat2, studlab,
                                data = NULL,
                                sm,
                                reference = ""){

  ## Check arguments
  ##
  chkchar(reference)


  ## Read data
  ##
  ##
  nulldata <- is.null(data)
  ##
  if (nulldata)
    data <- sys.frame(sys.parent())
  ##
  mf <- match.call()
  ##
  ## Catch TE, treat1, treat2, seTE, studlab from data:
  ##
  TE <- eval(mf[[match("TE", names(mf))]],
             data, enclos = sys.frame(sys.parent()))
  ##
  if (inherits(TE, "pairwise")) {
    sm <- attr(TE, "sm")
    ##
    seTE <- TE$seTE
    treat1 <- TE$treat1
    treat2 <- TE$treat2
    studlab <- TE$studlab
    ##
    if (!is.null(TE$n1))
      n1 <- TE$n1
    if (!is.null(TE$n2))
      n2 <- TE$n2
    if (!is.null(TE$event1))
      event1 <- TE$event1
    if (!is.null(TE$event2))
      event2 <- TE$event2
    ##
    is.pairwise <- TRUE
    pairdata <- TE
    data <- TE
    ##
    TE <- TE$TE
  }
  else {
    is.pairwise <- FALSE
    if (missing(sm))
      if (!is.null(data) && !is.null(attr(data, "sm")))
        sm <- attr(data, "sm")
      else
        sm <- ""
      ##
      seTE <- eval(mf[[match("seTE", names(mf))]],
                   data, enclos = sys.frame(sys.parent()))
      ##
      treat1 <- eval(mf[[match("treat1", names(mf))]],
                     data, enclos = sys.frame(sys.parent()))
      ##
      treat2 <- eval(mf[[match("treat2", names(mf))]],
                     data, enclos = sys.frame(sys.parent()))
      ##
      studlab <- eval(mf[[match("studlab", names(mf))]],
                      data, enclos = sys.frame(sys.parent()))
      ##
      n1 <- eval(mf[[match("n1", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
      ##
      n2 <- eval(mf[[match("n2", names(mf))]],
                 data, enclos = sys.frame(sys.parent()))
      ##
      event1 <- eval(mf[[match("event1", names(mf))]],
                     data, enclos = sys.frame(sys.parent()))
      ##
      event2 <- eval(mf[[match("event2", names(mf))]],
                     data, enclos = sys.frame(sys.parent()))
  }
  ##
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  ##
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  ##
  if (!is.numeric(studlab))
    studlab <- as.numeric(as.factor(studlab))



  ## Additional checks
  ##
  ## Check NAs and zero standard errors
  ##
  excl <- is.na(TE) | is.na(seTE) | seTE <= 0
  ##
  if (any(excl)) {
    dat.NAs <- data.frame(studlab = studlab[excl],
                          treat1 = treat1[excl],
                          treat2 = treat2[excl],
                          TE = format(round(TE[excl], 4)),
                          seTE = format(round(seTE[excl], 4))
    )
    warning("Comparison",
            if (sum(excl) > 1) "seTE",
            " with missing TE / seTE or zero seTE not considered in network meta-analysis.",
            call. = FALSE)
    cat(paste("Comparison",
              if (sum(excl) > 1) "s",
              " not considered in network meta-analysis:\n", sep = ""))
    prmatrix(dat.NAs, quote = FALSE, right = TRUE,
             rowlab = rep("", sum(excl)))
    ##
    studlab <- studlab[!(excl)]
    treat1  <- treat1[!(excl)]
    treat2  <- treat2[!(excl)]
    TE      <- TE[!(excl)]
    seTE    <- seTE[!(excl)]
  }
  ## Check for correct number of comparisons (after removing
  ## comparisons with missing data)
  ##
  is.wholenumber <-
    function(x, tol = .Machine$double.eps ^ 0.5)
      abs(x - round(x)) < tol
  ##
  tabnarms <- table(studlab)
  sel.narms <- !is.wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  ##
  if (sum(sel.narms) == 1)
    stop(paste("After removing comparisons with missing treatment effects",
               " or standard errors,\n  study '",
               names(tabnarms)[sel.narms],
               "' has a wrong number of comparisons.",
               " Please check data and\n  consider to remove study",
               " from network meta-analysis.",
               sep = ""))
  if (sum(sel.narms) > 1)
    stop(paste("After removing comparisons with missing treatment effects",
               " or standard errors,\n  the following studies have",
               " a wrong number of comparisons: ",
               paste(paste("'", names(tabnarms)[sel.narms], "'", sep = ""),
                     collapse = ", "),
               "\n  Please check data and consider to remove studies",
               " from network meta-analysis.",
               sep = ""))
  ##
  ## Check number of subgraphs
  ##
  n.subnets <- netconnection(treat1, treat2, studlab)$n.subnets
  ##
  if (n.subnets > 1)
    stop(paste("After removing comparisons with missing treatment effects",
               " or standard errors,\n  network consists of ",
               n.subnets, " separate sub-networks.\n  ",
               "Please check data and consider to remove studies",
               " from network meta-analysis.",
               sep = ""))
  ##
  ## Check for correct treatment order within comparison
  ##
  wo <- treat1 > treat2
  ##
  if (any(wo)) {
    warning("Note, treatments within a comparison have been re-sorted in increasing order.",
            call. = FALSE)
    TE[wo] <- -TE[wo]
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
  }


  # names of treatments
  names.treat <- sort(unique(c(treat1, treat2)))

  ## if no option exist, set reference treatment as the first in
  ## alphabetic / numeric order
  ##
  if (reference == "")
    reference <- names.treat[1]


  ## Conduct network meta-analysis (NMA) with random effects model,
  ## Rücker model
  model <- netmeta(TE, seTE, treat1, treat2, studlab,
                   comb.random = TRUE, reference.group = reference)


  ## Model objects
  ##
  ## number of treatments
  ##
  nt <- model$n


  ## treatment positions
  tr1 <- model$treat1.pos
  tr2 <- model$treat2.pos

  ## effect
  y.m <- model$TE


  ## B is the design matrix, the edge-vertex incidence matrix (mxn)
  ##
  #B <- createB(tr1, tr2, nt)

  #b <- model$TE.random[,reference]
  ##

  ## predicted estimate
  y.m.est <- model$TE.nma.random

  ## Outlier and influence diagnostics measures
  ##
  ## Raw residuals for each pairwise comparison
  ##
  rawres <- y.m - y.m.est

  ## Raw residuals for each study
  ##
  eraw <- res_multi(studlab, rawres)$res


  ##
  ## Standardized residuals for each pairwise comparison
  ##
  standres <- sqrt(model$w.random) * rawres


  ## Standardized residuals for each study
  ##
  estand <- res_multi(studlab, standres)$res


  ##
  ## Studentized residuals for each pairwise comparison
  ##
  studres <- 1/sqrt(1 - diag(model$H.matrix)) * sqrt(model$w.random) * rawres


  ## Studentized residuals for each study
  ##
  estud <- res_multi(studlab, studres)$res


  #Qi <- t(y.m[1:3] - model$TE.nma.fixed[1:3]) %*% diag(model$w.fixed[1:3]) %*% (y.m[1:3] - model$TE.nma.fixed[1:3])
  #t(y.m[7] - B[7,] %*% b) %*% model$w.fixed[7] %*% (y.m[7] - B[7,] %*% b)



  #t(y.m[j] - model$TE.nma.fixed[j]) %*% model$w.fixed[j] %*% (y.m[j] -  model$TE.nma.fixed[j])

  #j=11
 # model$TE.nma.fixed
 # model$w.fixed[j] * (y.m[j] - B[j,] %*% b)^2
  #Qi/3
  #model$TE.nma.fixed[j]
  #B[j,] %*% b

  ## right formula for Qi contribution
  #model$w.fixed[j] * (y.m[j] - model$TE.nma.fixed[j])^2
  #t(y.m[j] - model$TE.nma.fixed[j]) %*% model$w.fixed[j] %*% (y.m[j] -  model$TE.nma.fixed[j])

  #model$Q.heterogeneity

  ## Mahalanobis distance for each pairwise comparison
  ##
  Mah <- model$Q.fixed
  # Mahalanobis random
  #Mah <- model$w.random * (model$TE - model$TE.nma.random)^2

  ## right formula for Qi contribution
  Q.pooled <- model$w.fixed * (model$TE - model$TE.nma.fixed)^2
  Q.random <- model$w.random * (model$TE - model$TE.nma.random)^2

  ## Mahalanobis distance for each study
  ##
  Mahalanobis.distance <- res_multi(studlab, Mah)$res


  # leverage for each pairwise comparison
  lev <- as.numeric(diag(model$H.matrix))

  # leverage for each study
  leverage <- res_multi(studlab, lev)$res


  ##
  dat <- noquote(cbind(TE, seTE, studlab, treat1, treat2))
  ##
  rownames(dat) <- c(1:length(TE))
  ##


   res <- list(dat = dat,
               eraw = eraw,
               estand = estand,
               estud = estud,
               Mah = Mah,
               Mahalanobis.distance = Mahalanobis.distance,
               lev = lev,
               leverage = leverage)

  class(res) <- "NMAoutlier_measures"

  res


}
