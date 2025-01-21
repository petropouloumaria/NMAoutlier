#' Forward Search algorithm in network meta-analysis
#'
#' @description
#' This function employs the Forward Search algorithm to detect
#' outliers and influential studies fitted in network meta-analysis
#' model from graph-theory.  This is an outlying diagnostic tool to
#' detect outliers and studies that are potential sources for
#' heterogeneity and inconsistency in network meta-analysis.
#'
#' Monitoring measures during the search are:
#' \itemize{
#' \item outlier detection measures (standardized residuals, Cook's
#'   distance, ratio of variance);
#' \item ranking measures (P-scores);
#' \item heterogeneity and inconsistency measures (Q statistics for
#'   overall heterogeneity / inconsistency, inconsistency by
#'   design-by-treatment interaction model, z-values for comparison
#'   between direct and indirect evidence by back-calculation method).
#' }
#'
#' A description of the outlier detection methodology can be found in
#' Petropoulou et al. (2021).
#'
#' @param TE Estimate of treatment effect, i.e. difference between
#'   first and second treatment (e.g. log odds ratio, mean difference,
#'   or log hazard ratio). This can also be a pairwise object
#'   (i.e. the result of pairwise function of netmeta package).  In
#'   this case, the pairwise object should include the following: TE,
#'   seTE, treat1, treat2, studlab
#' @param seTE Standard error of treatment estimate.
#' @param treat1 Label/Number for first treatment.
#' @param treat2 Label/Number for second treatment.
#' @param studlab Study labels (important when multi arm studies are
#'   included).
#' @param data A data frame containing the study information.
#' @param crit1 A character string indicating the criterion to be used
#'   for selecting the initial subset, this criterion may be the
#'   minimum of median absolute residuals ("R") or the maximum of
#'   median absolute likelihood contributions ("L"). Default value is
#'   "R".
#' @param crit2 A character string indicating the criterion to be used
#'   for selecting the study entered from non-basic set to basic set,
#'   this criterion may be the minimum of absolute residuals ("R") or
#'   the maximum of absolute likelihood contributions ("L"). Default
#'   value is "R".
#' @param studies An optional vector specifying the number of the
#'   initial subset of studies. The default value is the maximum of
#'   the number of treatments and the 20 percent of the total number
#'   of studies.
#' @param P An optional vector specifying the number of candidate
#'   sample of studies (with length equal to \code{studies}) for the
#'   choice of the initial subset. Default value is 100.
#' @param sm A character string indicating underlying summary measure,
#'   e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"ASD"},
#'   \code{"HR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param Isub A vector for the studies to be included in the initial
#'   subset (default: NULL, the initial subset not specified by the
#'   user).
#' @param reference Reference treatment group.
#' @param small.values A character string indicating if small values
#'   are considered beneficial (option:"good") or harmfull
#'   (option:"bad") on outcome.  This is requirement for p-scores
#'   computation. The default value is considered benefial outcome
#'   ("good").
#' @param n_cores The number of cores that the process is running
#'   using the parallel (default: NULL, the process is running using
#'   all the available cores)
#' @param \dots Additional arguments passed on to
#'   \code{\link[netmeta]{netmeta}}.
#'
#' @details
#' FS algorithm for network meta-analysis model from graph theory is
#' described in Petropoulou et al. (2021).
#'
#' Let \emph{n} be the number of treatments and let \emph{m} be the
#' number of pairwise treatment comparisons.  If there are only
#' two-arm studies, \emph{m} is equal to the number of studies. Let TE
#' and seTE be the vectors of observed effects and their standard
#' errors.  Comparisons belonging to multi-arm studies are identified
#' by identical study labels (argument \code{studlab}).
#'
#' The FS algorithm is an outlier diagnostic iterative procedure. FS
#' algorithm apart from three steps. It starts with a subset of
#' studies and it gradually adds studies until all studies entered.
#' After the search, statistical measures are monitored for sharp
#' changes.
#'
#' In more detail, the FS algorithm starts with an initial subset of
#' the dataset with size \emph{l}. Let (argument \code{P})
#' (eg. \emph{P} = 100) a large number of candidate subset of studies
#' with size \emph{l}. The candidate subset that optimize the
#' criterion (argument \code{crit1}) is taken as the initial subset
#' (considered ideally to be outlying-free).  Criterion (\code{crit1})
#' to be used for selecting the initial subset, can be the minimum of
#' median absolute residuals \code{"R"} or the maximum of median
#' absolute likelihood contributions \code{"L"}. It is conventionally
#' refer this subset as basic set, whereas the remaining studies
#' constitute the non-basic set.
#'
#' The FS algorithm gradually adds studies from the non-basic to the
#' basic subset based on how close the former studies are to the
#' hypothesized model fit in the basic set.  A study from non-basic
#' set entered into the basic set if optimize the criterion (argument
#' \code{crit2}).  Criterion (\code{crit2}) for selecting the study
#' from non-basic to basic set may be the minimum of median absolute
#' residuals \code{"R"} or the maximum of median absolute likelihood
#' contributions \code{"L"}.  The algorithm order the studies
#' according to their closeness to the basic set by adding the study
#' that optimize the criterion from non-basic set to basic set.
#'
#' The process is repeated until all studies are entered into the
#' basic set.  The number of iterations of algorithm \emph{index} is
#' equal to the total number of studies minus the number of studies
#' entered into the initial subset. Through the FS procedure,
#' parameter estimates (summary estmates, heterogeneity estimator) and
#' other statistics of interest (outlying measures, heterogeneity and
#' inconsistency measures, ranking measures) are monitored. In each
#' iteration, network meta-analysis model from graph theory (Rücker,
#' 2012) is fitted (\code{netmeta} function) with R package
#' \bold{netmeta}.
#'
#' Monitoring statistical measures for each FS iteration can be:
#'
#' \bold{Outlying detection measures:}
#' Standardized residuals (arithmetic mean in case of multi-arm
#' studies); Cook's statistic; Ratio of determinants of
#' variance-covariance matrix
#'
#' \bold{Ranking measures:}
#' P-scores for ranking of treatments (Rücker G & Schwarzer G, 2015)
#' for each basic set with implementation of (\code{netrank} function)
#' from R package \bold{netmeta}.
#'
#' \bold{Heterogeneity and inconsistency measures:}
#' Overall heterogeneity / inconsistency Q statistic (\code{Q}) This
#' is the design-based decomposition of Cochran Q as provided by Krahn
#' et al. (2013); Overall heterogeneity Q statistic (\code{Q});
#' Between-designs Q statistic (\code{Q}), based on a random effects
#' model with square-root of between-study variance estimated embedded
#' in a full design-by-treatment interaction model.  Implementation
#' with (\code{decomp.design} function) from R package \bold{netmeta};
#' Z-values (Dias et al., 2010; König et al., 2013) for comparison
#' between direct and indirect evidence in
#' each iteration of forward search algorithm.  By monitoring
#' difference of direct and indirect evidence, potential sources of
#' consistency can be detected with the implementation of
#' (\code{netsplit} function) from R package \bold{netmeta} for each
#' iteration of the search.
#'
#' @return
#' An object of class \code{NMAoutlier}; a list containing the
#' following components:
#'    \item{dat}{Matrix containing the data \code{"TE"}, \code{"seTE"}, \code{"studlab"}, \code{"treat1"}, \code{"treat2"} as defined above.}
#'    \item{length.initial}{The number of studies that constitute the initial (outlying-clean) subset of studies.}
#'    \item{index}{The number of iterations of forward search algorithm.}
#'    \item{basic}{Studies entered into the basic set in each iteration of the search.
#'    At the first iteration, basic set constitute the studies that are included in the basic-initial subset.
#'    The number of studies in the first iteration is equal to length.initial.}
#'    \item{taub}{Heterogeneity estimator variance for basic set in each iteration of forward search algorithm.}
#'    \item{Qb}{Overall heterogeneity - inconsistency Q statistic (\code{Q}) for the basic set in each iteration of forward search algorithm.}
#'    \item{Qhb}{Overall heterogeneity Q statistic (\code{Q}) for the basic set in each iteration of forward search algorithm.}
#'    \item{Qib}{Overall inconsistency Q statistic (\code{Q}) from design-by-treatment interaction model for the basic set in each iteration of forward search algorithm.}
#'    \item{estb}{Summary estimates for each treatment for the basic set in each iteration of forward search algorithm.}
#'    \item{lb}{Lower 95\% confidence interval of summary estimates for the basic set in each iteration of forward search algorithm.}
#'    \item{ub}{Upper 95\% confidence interval of summary estimates for the basic set in each iteration of forward search algorithm.}
#'    \item{Ratio}{Ratio of determinants (\code{COVRATIOj}) of variance-covariance matrix of treatment estimates at iteration j to that iteration at (j-1).}
#'    \item{cook_d}{Cook's statistic (\code{Cj}) at iteration j of forward search algorithm.}
#'    \item{p.score}{P-score for ranking each treatment for the basic set in each iteration of forward search algorithm.}
#'    \item{dif}{Z-values for comparison between direct and indirect evidence for each iteration of forward search algorithm.
#'     Based on back-calculation method to derive indirect estimates from direct pairwise comparisons and network estimates.}
#'    \item{estand}{Standardized residuals for each study for the basic set in each iteration of forward search algorithm.}
#'    \item{call}{Function call}
#'
#' @references
#' Dias S, Welton NJ, Caldwell DM, Ades AE (2010):
#' Checking consistency in mixed treatment comparison meta-analysis.
#' \emph{Statistics in Medicine},
#' \bold{29}, 932--44
#'
#' König J, Krahn U, Binder H (2013):
#' Visualizing the flow of evidence in network meta-analysis and
#' characterizing mixed treatment comparisons.
#' \emph{Statistics in Medicine},
#' \bold{32}, 5414--29
#'
#' Krahn U, Binder H, König J (2013):
#' A graphical tool for locating inconsistency in network meta-analyses.
#' \emph{BMC Medical Research Methodology},
#' \bold{13}, 35
#'
#' Petropoulou M, Salanti G, Rücker G, Schwarzer G, Moustaki I,
#' Mavridis D (2021):
#' A forward search algorithm for detecting extreme study effects in
#' network meta-analysis.
#' \emph{Statistics in Medicine}
#'
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
#' @examples
#' \dontrun{
#' library("netmeta")
#' data(smokingcessation)
#' smokingcessation$id <- 1:nrow(smokingcessation)
#'
#' study912 <- subset(smokingcessation, id %in% 9:12)
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'                         list(event1, event2, event3),
#'                         list(n1, n2, n3),
#'                         data = study912,
#'                         sm = "OR")
#'
#' # Forward search algorithm
#' #
#' FSresult <- NMAoutlier(p1, P = 1, small.values = "bad", n_cores = 2)
#' FSresult
#'
#' data(smokingcessation)
#'
#' # Transform data from arm-based to contrast-based format
#' # We use 'sm' argument for odds ratios.
#' # We use function pairwise from netmeta package
#' #
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'                         list(event1, event2, event3),
#'                         list(n1, n2, n3),
#'                         data = smokingcessation,
#'                         sm = "OR")
#'
#' # Forward search algorithm
#' #
#' FSresult1 <- NMAoutlier(p1, small.values = "bad")
#'
#' # Basic set for each iteration of forward search algorithm
#' #
#' FSresult1$basic
#'
#' # Forward search algorithm using the criteria (crit1, crit2)
#' # with the maximum of absolute likelihood contributions ("L")
#' #
#' FSresult2 <- NMAoutlier(p1, crit1 = "L", crit2 = "L",
#'                         small.values = "bad")
#' FSresult2
#' }
#'
#' @export
#'
#' @author Maria Petropoulou <maria.petropoulou@uniklinik-freiburg.de>


NMAoutlier <- function(TE, seTE, treat1, treat2, studlab,
                       data = NULL,
                       crit1 = "R", crit2 = "R",
                       studies = NULL,
                       P = 100,
                       sm,
                       Isub = NULL,
                       reference = "", small.values = "good", n_cores = NULL,
                       ...) {


  ## Check arguments
  ##
  crit1 <- setchar(crit1, c("R", "L"))
  crit2 <- setchar(crit2, c("R", "L"))
  chkchar(reference)
  chkchar(small.values)
  ##
  ##
  chknumeric(P, min = 0, single = TRUE)
  ##
  ##
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
  if (inherits(TE, "pairwise") ||
      is.data.frame(TE) & !is.null(attr(TE, "pairwise"))) {
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
            " with missing TE / seTE or zero seTE not considered in Forward Search algorithm.",
            call. = FALSE)
    cat(paste("Comparison",
              if (sum(excl) > 1) "s",
              " not considered in Forward Search algorithm:\n", sep = ""))
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
               " from Forward Search algorithm.",
               sep = ""))
  if (sum(sel.narms) > 1)
    stop(paste("After removing comparisons with missing treatment effects",
               " or standard errors,\n  the following studies have",
               " a wrong number of comparisons: ",
               paste(paste("'", names(tabnarms)[sel.narms], "'", sep = ""),
                     collapse = ", "),
               "\n  Please check data and consider to remove studies",
               " from Forward Search algorithm.",
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


  ##
  ##
  ##
  m <- length(TE)                                # number of pairwise comparisons (edges)
  n <- length(unique(c(treat1, treat2)))         # number of treatments (vertices)
  nk <- length(unique(studlab))                  # number of studies
  names.treat <- sort(unique(c(treat1, treat2))) # names of treatments
  ##
  ##

  ## Number of studies for the initial subset (default choice)
  ##
  nullstudies <- is.null(studies)
  ##
  if (nullstudies)
    studies <- max(round(nk * 0.20), n)
  ##
  ## Check argument
  ##
  chknumeric(studies, min = max(round(nk * 0.20), n),
             single = TRUE)
  ##

  ## Adapt numbers to treatment IDs
  ##
  t1.label <- match(treat1, names.treat)
  t2.label <- match(treat2, names.treat)


  ## Names of comparisons
  ##
  ##
  names.comp <- c() # rep("", )
  k <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      k <- k + 1
      names.comp[k] <- paste(names.treat[i], names.treat[j], sep = ":")
    }
  }
  ##
  ##
  ##
  ## At least number of studies equal to the number of treatments
  ##
  if (nk < n)
    stop("Network meta-analysis has fewer studies (",
         nk,
         ") than number of treatments (", n,
         "). Note, Forward Search algorithm requires a network-meta-analysis ",
         "with at least as many studies as number of treatments.",
         call. = FALSE)
  else {
    ##
    ## if no option exist, set reference treatment as the first in
    ## alphabetic / numeric order
    ##
    if (reference == "")
      reference <- names.treat[1]


    ## NMA for the whole dataset
    ##
    res <- nma(TE, seTE, treat1, treat2, studlab, reference, names.treat, ...)

    ## initialize lists for variance-covariance matrix of basic set
    ## for each step of FS algorithm
    ##
    liv <- list()

    ## number of iterations of FS algorithm
    index <- 1

    ## study that inputed from non-basic to basic set
    subset <- NULL

    if (is.null(Isub)){
      ## Take the initial subset
      Isub <- c(InitialSubset(TE, seTE, treat1, treat2, studlab,
                              crit1, studies, P, reference,
                              t1.label, t2.label, n_cores, ...)$set)
    }

    ## Define the initial basic set
    ##
    bs <- Isub
    ##
    ## indices of basic set
    ##
    ind.bs <- which(studlab %in% bs)
    ##
    S1 <- netmet(TE, seTE, treat1, treat2, studlab,
                 ind.bs, reference, small.values, names.treat, ...)
    ##
    ## Initial basic set
    ##
    taub <- S1$t           # heterogeneity
    estb <- S1$b           # summary estimate
    ##
    lb <- S1$l             # confidence interval of summary estimate
    ub <- S1$u
    ##
    Qb <- S1$Q             # Q statistic
    Qhb <- S1$Qhet
    Qib <- S1$Qinc
    ##
    dif <- S1$diff          # difference of direct and indirect evidence (node-splitting)
    ##
    p.score <- S1$p.rank    # P-score
    ##
    estand <- S1$st.res    # standardized residuals
    length(estand) <- nk
    ##
    liv[[index]] <- S1$Cov # Variance-covariance matrix
    ##
    cook_d <- NULL
    Ratio <- NULL
    ##
    ##
    ## Do Forward Search algorithm for the number of studies minus the
    ## number of studies of basic set (nk - length(bs)) steps
    ##
    while (nk - length(bs) > 0) {
      ##
      condition1 <- nk - length(bs) > 1
      condition2 <- nk - length(bs) == 1
      ##
      ##
      ## Add to the basic set, the study with the optimal function
      ##
      if (condition1) { # big loop if
        ##
        mi <- ma <- c()
        ##
        ## Conduct network meta-analysis (NMA) with random effects model, Rucker model
        ##
        model <- netmeta(TE, seTE, treat1, treat2, studlab,
                         comb.random = TRUE, reference.group = reference,
                         subset = ind.bs, ...)
        ##
        t.basic <- (model$tau) ^ 2              # heterogeneity
        e.basic <- model$TE.random[, reference] # summary estimate
        ##
        ##
        ## Non-basic set
        ##
        c <- setdiff(1:m, ind.bs)      # indices of non-basic set
        ##
        ## treatment labels of non-basic set
        non_b.lab <- sort(unique(c(as.character(treat1[c]), as.character(treat2[c]))))
        t_ind <- match(non_b.lab, names.treat)
        ##
        ## summary estimate of basic set for the corresponding treatments in non-basic set
        ##
        m.basic <- e.basic[t_ind]
        ##
        ## design matrix of non-basic set
        ##
        B.nb <- res$X[c, t_ind]
        ##
        ## weights of random effects model
        ##
        wr <- prepare(TE[c], seTE[c], treat1[c], treat2[c], studlab[c],
                      t.basic)$weights
        ##
        ## residuals
        ##
        rs <- c(TE[c] - B.nb %*% m.basic)
        ##
        ## standardized residuals
        ##
        st.res <- sqrt(wr) * rs
        ##
        rn <- Multi(studlab[c], st.res, wr)
        ##
        if (crit2 == "R") {
          mi <- abs(rn$res)
          pick2 <- which.min(mi)
          subset <- rn$study[pick2]
        }
        else if (crit2 == "L") {
          ma <- abs(rn$logl)
          pick2 <- which.max(ma)
          subset <- rn$study[pick2]
        }
      }
      else if (condition2)
        subset <- setdiff(studlab, bs)
      ##
      ## end big loop if (end of adding)
      ##
      index <- index + 1
      ##

      ## Re-define the "basic set: D(m + index)"
      ##
      bs <- c(bs, subset) # basic set
      ##
      ind.bs <- which(studlab %in% bs)
      ##

      ##
      Si <- netmet(TE, seTE, treat1, treat2, studlab,
                   ind.bs, reference, small.values, names.treat, ...)
      ##
      ## Basic set
      ##
      taub <- c(taub, Si$t)                # heterogeneity
      estb <- cbind(estb, Si$b)            # summary estimate
      ##
      lb <- cbind(lb, Si$l)                # confidence interval of summary estimate
      ub <- cbind(ub, Si$u)
      ##
      Qb <- c(Qb, Si$Q)                    # Q statistic
      Qhb <- c(Qhb, Si$Qhet)
      Qib <- c(Qib, Si$Qinc)
      ##
      dif <- cbind(dif, Si$diff)           # difference of direct and indirect evidence (node-splitting)
      ##
      p.score <- cbind(p.score, Si$p.rank) # P-score
      ##
      length(Si$st.res) <- nk
      estand <- cbind(estand, Si$st.res)   # standardized residuals
      ##
      liv[[index]] <- Si$Cov                # Variance-covariance matrix
      ##
      ## Cook's statistic at index step
      ##
      ind <- 2:length(estb[, index])
      ##
      cookj <- t(estb[ind, index] - estb[ind, index - 1]) %*%
        ginv(liv[[index]]) %*% (estb[ind, index] - estb[ind, index - 1])
      cook_d <- c(cook_d, cookj)
      ##
      ## Ratio of the determinants of the variance-covariance matrix
      ## at step (index) to that at step (index - 1)
      ##
      Rj <- (det(liv[[index]])) / (det(liv[[index - 1]]))
      Ratio <- c(Ratio, Rj)

    } # End while


    ##
    ## length of the initial basic set
    ##
    length.initial <- length(Isub)
    int <- rep(1, length.initial)
    iteration <- c(int, 2:index)
    ##
    basic <- matrix(bs, nrow = 1)
    colnames(basic) <- paste("it=", iteration)
    rownames(basic) <- paste("study entered")
    ##
    dat <- noquote(cbind(TE, seTE, studlab, treat1, treat2))
    ##
    rownames(dat) <- c(1:length(TE))
    ##
    colnames(estb) <- colnames(lb) <- colnames(ub) <- colnames(dif) <-
    colnames(p.score) <- colnames(estand) <- paste("it=", 1:index)
    ##
    rownames(dif) <- names.comp
    rownames(estand) <- paste(bs)
    ##
    ## remove z-values with NA treatment comparisons
    ##
    for (j in 1:ncol(dif))
      excl <- is.na(dif[, j])
    ##
    ##

    dif <- dif[which(as.numeric(!excl) == 1), ]


    res2 <- list(dat = dat, length.initial = length.initial,
                 index = index, basic = basic,
                 taub = taub, Qb = Qb, Qhb = Qhb, Qib = Qib,
                 estb = estb, lb = lb, ub = ub, Ratio = Ratio,
                 cook_d = cook_d, p.score = p.score,
                 dif = dif, estand = estand,
                 call = match.call())


  } # end ("if" requirement to be equal the number of studies with the number of treatments)


  class(res2) <- "NMAoutlier"

  res2

} # End function
