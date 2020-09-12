#' Outlier and influential detection measures in network meta-analysis considered deletion of study.
#'
#' @description
#' Employs the computation of several measures for detection of outlying
#' studies (studies with extreme results) and influential studies considering study deletion fitted in network
#' meta-analysis model from graph-theory. Calculation of these measures for each study
#' included in the network can detect evidence of outliers. It can also be used to detect
#' studies that are potential sources for heterogeneity and
#' inconsistency.
#'
#' Statistical outlier and influential deletion measures are:
#' \itemize{
#' \item outlying and influential measures for each study are:
#' Standardized deleted residual;Studentized deleted residual; Cook
#' distance between the treatment estimates for study j and treatment estimates when study j is removed;
#' ratio of determinants of variances of variance-covariance matrix
#' of treatment estimates for study j to treatment estimates when study j is removed;
#' weight leave one out; leverage leave one out; heterogeneity estimator leave one out;
#' R statistic for heterogeneity;  R statistic for Q (\code{Qtotal}),  R statistic for  heterogeneity Q
#' (\code{Qhet}), R statistic for Qinconsistency (\code{Qinc}), DFbetas.
#' }
#'
#' A description of the several outlier detection measures considering study deletion in the context of network meta-analysis
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
#' The function calculates several outlier detection measures for each study.
#' The statistical measures calculated such as are raw, standardized and studentized deleted residuals,
#' Cook distance, COVRATIO, weight leave one out etc.
#'
#' @return
#' An object of class \code{NMAoutlier_measures}; a list containing the
#' following components:
#'    \item{dat}{Matrix containing the data \code{"TE"}, \code{"seTE"}, \code{"studlab"}, \code{"treat1"}, \code{"treat2"} as defined above.}
#'    \item{eraw.deleted}{Raw deleted residual for each study.}
#'    \item{estand.deleted}{Standardized deleted residual for each study.}
#'    \item{estud.deleted}{Studentized deleted residual for each study.}
#'    \item{Cooks.distance}{Cook's distance.}
#'    \item{Covratio}{ratio of determinants of variances of variance-covariance matrix of treatment estimates.}
#'    \item{w.leaveoneout}{weight leave one out.}
#'    \item{H.leaveoneout}{leverage leave one out.}
#'    \item{heterog.leaveoneout}{heterogeneity estimator leave one out.}
#'    \item{Rheterogeneity}{R statistic for heterogeneity.}
#'    \item{RQtotal}{R statistic for Qtotal.}
#'    \item{RQhet}{R statistic for Qheterogeneity.}
#'    \item{RQinc}{R statistic for Qinconsistency.}
#'    \item{DFbetas}{DFbetas.}
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
#' #outlier and influential detection measures for studies 9, 10, 11, 12.
#' deletion <- NMAoutlier_deletion_measures(p1)
#'
#' # Standardized deleted residual for studies 9, 10, 11, 12.
#' deletion$estand.deleted
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
#' deletion <- NMAoutlier_deletion_measures(p1)
#'
#' # Cook's distance between the treatment estimates for study j
#' # and treatment estimates when study j is removed.
#' deletion$Cooks.distance
#'
#'
#' }
#'
#' @export
#'
#' @author Maria Petropoulou <mpetrop@cc.uoi.gr>
#'
#' @importFrom netmeta netmeta
#' @importFrom MASS ginv



NMAoutlier_deletion_measures <- function(TE, seTE, treat1, treat2, studlab,
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
            " with missing TE / seTE or zero seTE not considered in computations of deletion measures.",
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
               " from network meta-analysis dataset.",
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
  nt <- model$n                     # number of treatments

  tr1 <- model$treat1.pos           # treatment positions
  tr2 <- model$treat2.pos
  y.m <- model$TE
  s.m <-  model$seTE

  ## B is the design matrix, the edge-vertex incidence matrix (mxn)
  ##
  B <- createB(tr1, tr2, nt)
  ##
  t <- (model$tau)^2                   # heterogeneity
  ##
  b <- model$TE.random[, reference]    # summary estimate of treatment effects



  ## Computations for variance-covariance matrix (random effects model)
  ##
  ind <- which(names.treat == reference) # index of reference treatment
  t.pos1 <- rep(ind, (nt - 1))
  t.pos2 <- setdiff(1:nt, ind)
  ##
  B.r <- createB(t.pos1, t.pos2, ncol = nt) # Restricted matrix B

  # Laplacian matrix

  L.random <- t(B) %*% diag(model$w.random) %*% B
  Lplus.random <- solve(L.random - 1 / nt) + 1 / nt


  # Variance-covariance matrix (random effects model)
  Cov <- B.r %*% Lplus.random %*% t(B.r)


  ## Q statistics
  ##
  Q.standard <- model$Q                                       # overall
  Qh.standard <- model$Q.heterogeneity                        # within-designs
  Qi.standard <- decomp.design(model)$Q.inc.random$Q          # between-designs


  ## Outlier and influence diagnostics measures considered deletion


  studies <- unique(studlab)


  heterog.leaveoneout <- w.leaveoneout <- H.leaveoneout <- eraw.deleted <- estand.deleted <- estud.deleted <- Cooks.distance <- Covratio <- Rstat.heterogeneity <- RQtotal <- RQhet <- RQinc <- list()
  DFbetas <- Rstat.estimates <- NULL

  for (i in 1:length(studies)) {

    # deleted study
    deleted <- studies[i]
    remaining <- studies[studies != deleted]
    remaining.ind <- which(studlab %in% remaining)
    netmeta.res <- netmeta(TE, seTE, treat1, treat2, studlab, comb.random = TRUE, reference.group = reference, subset = remaining.ind)


    estimate <- netmeta.res$TE.random[,reference]                # summary estimates

    heterog <- (netmeta.res$tau)^2                               # heterogeneity
    heterog.leaveoneout[[i]] <- heterog


    ## Q statistics
    Qt <- netmeta.res$Q                                          # overall
    Qhe <- netmeta.res$Q.heterogeneity                           # within-designs
    Qin <- decomp.design(netmeta.res)$Q.inc.random$Q             # between-designs

    # index of study deleted
    ind.deleted <- which(studlab == deleted)


    # weight
    w.leave <- 1/(s.m[ind.deleted]^2 + heterog)

    ## Standardized study deleted residuals
    w.leaveoneout[[i]] <- res_multi(studlab[ind.deleted], w.leave)$res

    # hat values
    #
    Bi <- B[ind.deleted,]
    Bi.matrix <- matrix(Bi, ncol = nt)


    # leverage "leave-one-out"
    wi.matrix <- diag(w.leave,  nrow = length(w.leave), ncol = length(w.leave))
    hii <- diag(Bi.matrix %*% ginv(t(Bi.matrix) %*% wi.matrix %*% Bi.matrix) %*% t(Bi.matrix) %*% wi.matrix)

    H.leaveoneout[[i]] <- res_multi(studlab[ind.deleted], hii)$res

    n <- netmeta.res$n                                           # number of treatments

    t1 <- netmeta.res$treat1.pos                                 # treatment positions
    t2 <- netmeta.res$treat2.pos                                 # treatment positions

    ## B is the design matrix, the edge-vertex incidence matrix (mxn)
    Brem <- createB(t1, t2, n)

    ## Laplacian matrix
    ##
    L.r <- t(Brem) %*% diag(netmeta.res$w.random) %*% Brem
    Lplus <- solve(L.r - 1 / n) + 1 / n



    ## Computations for variance-covariance matrix (random effects model)
    ##
    ind <- which(names.treat == reference)    # index of reference treatment
    t.pos1 <- rep(ind, (n - 1))
    t.pos2 <- setdiff(1:n, ind)
    ##
    Br.remove <- createB(t.pos1, t.pos2, ncol = n)  # Restricted matrix B
    Cov.remove <- Br.remove %*% Lplus %*% t(Br.remove)    # Variance-covariance matrix (random effects model)


    ## Raw pairwise deleted residuals
    rawres <- c(y.m[ind.deleted] - B[ind.deleted,] %*% estimate)

    ## Raw study deleted residuals
    eraw.deleted[[i]] <- res_multi(studlab[ind.deleted], rawres)$res


    ## Standardized pairwise deleted residuals
    standres <- sqrt(w.leave) * rawres

    ## Standardized study deleted residuals
    estand.deleted[[i]] <- res_multi(studlab[ind.deleted], standres)$res


    ## Studentized pairwise deleted residuals
    studres <- 1/sqrt(s.m[ind.deleted]^2 + t + hii *  w.leave) * rawres

    ## Studentized study deleted residuals
    estud.deleted[[i]] <-  res_multi(studlab[ind.deleted], studres)$res


    ## Cook's statistic considered deletion
    Cooks.distance[[i]] <- c(t(b[2:length(b)] - estimate[2:length(estimate)]) %*% ginv(Cov) %*% (b[2:length(b)] - estimate[2:length(estimate)]))

    ## Covratio considered deletion

    ## Ratio of the determinants of the variance-covariance matrix
    Covratio[[i]] <- det(Cov.remove) / det(Cov)

    ## R statistic for heterogeneity
    Rstat.heterogeneity[[i]] <- ((t - heterog) / t) * 100

    ## R statistic for summary estimates
    Rstat.estimate <- ((b[2:length(b)] - estimate[2:length(estimate)]) / b[2:length(b)]) * 100
    Rstat.estimates  <- cbind(Rstat.estimates, Rstat.estimate)

    ## R statistic total
    RQtotal[[i]] <- ((Q.standard - Qt) / Q.standard) * 100

    ## R statistic hererogeneity
    RQhet[[i]] <- ((Qh.standard - Qhe) / Qh.standard) * 100

    ## R statistic inconsistency
    RQinc[[i]] <- ((Qi.standard - Qin) / Qi.standard) * 100

    ## DFbetas
    DFbeta <- (b[2:length(b)] - estimate[2:length(estimate)]) * sqrt(sum(w.leave)/length(w.leave))
    DFbetas <- cbind(DFbetas, DFbeta)
  }

  ##
  dat <- noquote(cbind(TE, seTE, studlab, treat1, treat2))
  ##
  rownames(dat) <- c(1:length(TE))
  ##


  res <- list(dat = dat,
              eraw.deleted = unlist(eraw.deleted),
              estand.deleted = unlist(estand.deleted),
              estud.deleted = unlist(estud.deleted),
              Cooks.distance = unlist(Cooks.distance),
              Covratio = unlist(Covratio),
              w.leaveoneout = unlist(w.leaveoneout),
              H.leaveoneout = unlist(H.leaveoneout),
              heterog.leaveoneout = unlist(heterog.leaveoneout),
              Rheterogeneity = unlist(Rstat.heterogeneity),
              Restimates = Rstat.estimates,
              RQtotal = unlist(RQtotal),
              RQhet = unlist(RQhet),
              RQinc = unlist(RQinc),
              DFbetas = DFbetas)

  class(res) <- "NMAoutlier_deletion_measures"

  res

}




