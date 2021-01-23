#' Random shift variance network meta-analysis (RVSOM NMA) model.
#'
#' Fits the RVSOM NMA model for detecting outliers.
#'
#' @description
#' Fits the Random Shift Variance model for detection of outlying
#' studies (studies with extreme results) fitted in network
#' meta-analysis model from graph-theory. This is a diagnostic tool
#' for detection evidence of outliers. It can also be used to detect
#' studies that are potential sources for heterogeneity and
#' inconsistency.
#'
#' Monitoring measures during the search are:
#' \itemize{
#' \item outlying measures (standardized residuals, leverage);
#' \item ranking measures (P-scores);
#' \item heterogeneity and inconsistency measures (Q statistics for
#'   overall heterogeneity / inconsistency, inconsistency by
#'   design-by-treatment interaction model, z-values for comparison
#'   between direct and indirect evidence by back-calculation method).
#'  \item over dispresion parameter
#'  \item Likelihood Ratio Test (LRT)
#' }
#'
#' A description of the methodology can be found in Petropoulou et
#' al. (2020).
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
#' @param small.values A character string indicating if small values
#'   are considered beneficial (option:"good") or harmfull
#'   (option:"bad") on outcome.  This is requirement for p-scores
#'   computation. The default value is considered benefial outcome
#'   ("good").
#' @param study The study (or studies) to fit by shift its (their) variance.
#'   By default, the standard random effects model is fitted, value NULL.
#' @param n_cores The number of cores that the process is running
#'   using the parallel (default: NULL, the process is running
#'   using all the available cores)
#'
#' @details
#' Description of methodology by fitting random shift variance model (RVSOM) in
#' network meta-analysis. Methodology of RVSOM fitted in NMA
#' model from graph theory is described in Petropoulou (2020).
#'
#'
#' NMAoutlier.rsv function employs the random shift variance model (RVSOM)
#' for detection of outlying and influential studies fitted in network meta-analysis model
#' from graph-theory.The random effects variance shift outlier model (RVSOM)
#' allows inflated variances for the study j to detect if study j is an outlier.
#' The RVSOM model is fitted by shifting the variance for each study
#' and sharp changes in statistical measures between RVSOM and standard NMA model
#' can indicate which study is an outlier. RVSOM NMA model is fitted in network meta-analysis model
#' from graph theory (Rücker, 2012) (\code{netmeta} function) with
#' R package \bold{netmeta} (Rücker et al., 2015) using the restricted maximum likelihood
#' estimation method for heterogeneity variance. The variance estimator of RVSOM NMA model
#' decomposed by the heterogeneity variance and the over-dispersion.
#' The researcher can choose the reference treatment (\code{reference}) fitted in NMA model.
#' Monitoring is helpful to identify outlying and/or influential studies.
#'
#' Monitoring statistical measures (outlying measures,
#' heterogeneity and inconsistency measures, ranking measures)
#' and likelihood statistics can be:
#'
#'
#' Monitoring is helpful to identify outlying studies:
#'
#' \bold{Likelihood statistics:}
#' Heterogeneity estimation method is conducted under
#' the Restricted Maximum likelihood estimator.
#' The variance estimator of shift variance model.
#' Likelihood statistics offered from calculation:
#' the twice of maximum log-likelihood;
#' the convergence diagnostic and the Likelihood Ratio test (LRT) test.
#'
#' \bold{Outlier and influential case diagnostics measures:}
#' Standardized residuals (arithmetic mean in case of multi-arm
#' studies); leverage of hat matrix
#'
#' \bold{Ranking measures:}
#' P-scores for ranking of treatments (Rücker G & Schwarzer G (2015))
#' for each basic set with implementation of (\code{netrank} function)
#' from R package \bold{netmeta}.
#'
#' \bold{Heterogeneity and inconsistency measures:}
#' Overall heterogeneity / inconsistency Q statistic (\code{Q}) This
#' is the design-based decomposition of Cochran Q as provided by Krahn
#' et al.(2013); Overall heterogeneity Q statistic (\code{Q});
#' Between-designs Q statistic (\code{Q}), based on a random effects
#' model with square-root of between-study variance estimated embedded
#' in a full design-by-treatment interaction model. Implementation
#' with (\code{decomp.design} function) from R package \bold{netmeta};
#' Z-values for comparison between direct and indirect evidence for
#' each iteration of forward search algorithm.  By monitoring
#' difference of direct and indirect evidence, potential sources of
#' consistency can be detected with the implementation of
#' (\code{netsplit} function) from R package \bold{netmeta} for each
#' basic set of the search.  Based on the methodology with
#' back-calculation method to derive indirect estimates from direct
#' pairwise comparisons and network estimates (Dias et al., 2010;
#' König et al., 2013).
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
#' Petropoulou M (2020):
#' Exploring methodological challenges in network meta-analysis models and
#' developing methodology for outlier detection.
#' \emph{PhD dissertation}.
#'
#' Krahn U, Binder H, König J (2013):
#' A graphical tool for locating inconsistency in network meta-analyses.
#' \emph{BMC Medical Research Methodology},
#' \bold{13}, 35
#'
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
#' @examples
#' data(smokingcessation, package = "netmeta")
#' smokingcessation$id <- 1:nrow(smokingcessation)
#' p1 <- netmeta::pairwise(list(treat1, treat2, treat3),
#'                         list(event1, event2, event3),
#'                         list(n1, n2, n3),
#'                         data = smokingcessation,
#'                         sm="OR")
#'
#' # Random Shift Variance NMA model implementation for study 1
#' #
#' RVSOMres <- NMAoutlier.rsv(p1,  small.values = "bad", study = c(1) , n_cores = 2)
#' #
#'
#' \dontrun{
#' # Random Shift Variance Model
#' #
#' RVSOMresult <- NMAoutlier.rsv(p1, small.values = "bad")
#'
#' # Random shift variance estimator
#' RVSOMresult$over-disp
#' #
#' # Plot of Likelihood Ratio Test
#' RVSOMresult$LRT
#' #
#' # Summary estimators and their 95% confidence intervals for each treatment
#' # by downweighting each study
#' RVSOMresult$b
#' RVSOMresult$l
#' RVSOMresult$u
#' }
#'
#' @export
#'
#' @author Maria Petropoulou <petropoulou@imbi.uni-freiburg.de>
#'
#' @importFrom parallel clusterExport detectCores makeCluster parLapply stopCluster



NMAoutlier.rsv <- function(TE, seTE, treat1, treat2, studlab,
                   data = NULL,
                   sm,
                   reference = "", small.values = "good",
                   study = NULL, n_cores =  NULL){


   ## Check arguments
   ##
   chkchar(reference)
   chkchar(small.values)

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
             " with missing TE / seTE or zero seTE not considered in Random Shift Variance model.",
             call. = FALSE)
     cat(paste("Comparison",
               if (sum(excl) > 1) "s",
               " not considered in Random Shift Variance model:\n", sep = ""))
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
                " from Random Shift Variance model.",
                sep = ""))
   if (sum(sel.narms) > 1)
     stop(paste("After removing comparisons with missing treatment effects",
                " or standard errors,\n  the following studies have",
                " a wrong number of comparisons: ",
                paste(paste("'", names(tabnarms)[sel.narms], "'", sep = ""),
                      collapse = ", "),
                "\n  Please check data and consider to remove studies",
                " from data.",
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
                " from Random Shift Variance model.",
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



   if (is.null(study) == "TRUE") {

     ## take all the studies from dataset
     study <- c(NA, 1:length(unique(studlab)))

   }else{

     study <- c(NA, study)

   }

   ## fit RVSOM model for NMA parallel for each study of NMA
   if (is.null(n_cores)) {
     ## Use all available cores
     n_cores <- max(1, detectCores())
   }



    cl <- makeCluster(n_cores)


    clusterExport(cl = cl,
                  varlist = c("TE", "seTE", "treat1", "treat2",
                        "studlab","reference", "small.values","res_multi",
                        "NMA_RVSOM", "netmeta", "createB", "NMA_REML",
                        "prepare", "multiarm",
                        "decomp.design", "netrank", "netsplit"),
                   envir = environment())


      ## fit RVSOM model for NMA parallel for each study of NMA
      paral <- parLapply(cl, study, function(z) {



       RVmodel <- NMA_RVSOM(TE = TE, seTE = seTE,
                           treat1 = treat1, treat2 = treat2, studlab = studlab,
                           reference = reference, small.values = small.values,
                           observations = z)


       tau <- RVmodel$tau2est
       Q <- RVmodel$Q
       Qhet <- RVmodel$Qhet
       Qinc <- RVmodel$Qinc

       b <- RVmodel$b
       l <- RVmodel$l
       u <- RVmodel$u

       leverage <- RVmodel$leverage
       estand <- RVmodel$estand

       p.score <- RVmodel$p.score

       dif <- RVmodel$dif

       over_disp <- RVmodel$over_disp
       converge <- RVmodel$converge
       twiceloglik <- -1/2*RVmodel$twiceloglik

       result <- list(z = z,
                     tau = tau,
                     Q = Q, Qhet = Qhet, Qinc = Qinc,
                     b = b, l = l, u = u,
                     leverage = leverage, estand = estand,
                     p.score = p.score, dif = dif,
                     over_disp = over_disp,
                     converge = converge, twiceloglik = twiceloglik)



   })

   ## close the cluster after done
   stopCluster(cl)



   b <- z <- tau <- Q <- Qhet <- Qinc <- b <- l <- c()
   u <- leverage <- estand <- p.score <- dif <- c()
   over_disp <- converge <-  twiceloglik <- c()

   result <- as.data.frame(do.call(rbind, paral))


   for (i in 1:length(paral)){

     z <- cbind(z, result$z[[i]])

     tau <- cbind(tau, result$tau[[i]])

     Q <- cbind(Q, result$Q[[i]])

     Qhet <- cbind(Qhet, result$Qhet[[i]])

     Qinc <- cbind(Qinc, result$Qinc[[i]])

     b <- cbind(b, result$b[[i]])

     l <- cbind(l, result$l[[i]])

     u <- cbind(u, result$u[[i]])

     leverage <- cbind(leverage, result$leverage[[i]])

     estand <- cbind(estand, result$estand[[i]])

     p.score <- cbind(p.score, result$p.score[[i]])

     dif <- cbind(dif, result$dif[[i]])

     over_disp <- cbind(over_disp, result$over_disp[[i]])

     converge <- cbind(converge, result$converge[[i]])

     twiceloglik <- cbind(twiceloglik, result$twiceloglik[[i]])

   }

   ## Loglikelihood test
   twiceloglikNULL <- twiceloglik[, 1]
   twiceloglikREST <- twiceloglik[, 2:length(twiceloglik)]

   LRT <- 2 * (rep(twiceloglikNULL,length(twiceloglikREST)) - twiceloglikREST)


   names.treat <- sort(unique(c(treat1, treat2))) # names of treatments

   ## Names of comparisons

   ## number of treatments
   nt <- length(names.treat)
   names.comp <- c()
   k <- 0
   for (i in 1:(nt - 1)) {
     for (j in (i + 1):nt) {
       k <- k + 1
       names.comp[k] <- paste(names.treat[i], names.treat[j], sep = ":")
     }
   }

   ##
   dat <- noquote(cbind(TE, seTE, studlab, treat1, treat2))
   ##
   rownames(dat) <- c(1:length(TE))
   ##
   colnames(b) <- colnames(l) <- colnames(u) <- colnames(dif) <-
     colnames(p.score) <- colnames(estand) <- colnames(leverage) <- paste("RVSOM:",study)
   ##


   rownames(dif) <- names.comp
   rownames(leverage) <- paste(studlab)
   rownames(estand) <- paste(unique(studlab))
   ##
   ## remove z-values with NA treatment comparisons
   ##
   for (j in 1:ncol(dif))
     excl <- is.na(dif[, j])
   ##
   ##

   dif <- dif[which(as.numeric(!excl) == 1), ]

   ## objects for standard random variance model

   tau_standard <- tau[1]
   Q_standard <- Q[1]
   Qhet_standard <- Qhet[1]
   Qinc_standard <- Qinc[1]
   b_standard <- b[, 1]
   l_standard <- l[ , 1]
   u_standard <- u[ , 1]
   leverage_standard <- leverage[ , 1]
   estand_standard <- estand[ , 1]
   p.score_standard <- p.score[ , 1]
   dif_standard <- dif[ ,1]
   over_disp_standard <- over_disp[1]
   converge_standard <- converge[1]
   twiceloglik_standard <- twiceloglik[1]


   tau <- tau[!tau %in% tau[ , 1]]
   Q <- Q[!Q %in% Q[ , 1]]
   Qhet <- Qhet[!Qhet %in% Qhet[ , 1]]
   Qinc <- Qinc[!Qinc %in% Qinc[ , 1]]
   b <- b[ , 2:length(study)]
   l <- l[ , 2:length(study)]
   u <- u[ , 2:length(study)]
   leverage <- leverage[ , 2:length(study)]
   estand <- estand[ , 2:length(study)]
   p.score <- p.score[, 2:length(study)]
   dif <- dif[ , 2:length(study)]
   over_disp <- over_disp[ , 2:length(study)]
   converge <- converge[ , 2:length(study)]
   twiceloglik <- twiceloglik[ , 2:length(study)]


   res <- list(dat = dat, z = z,
               tau_standard = tau_standard,
               Q_standard = Q_standard,
               Qhet_standard = Qhet_standard,
               Qinc_standard = Qinc_standard,
               b_standard = b_standard,
               l_standard = l_standard,
               u_standard = u_standard,
               leverage_standard = leverage_standard,
               estand_standard = estand_standard,
               p.score_standard = p.score_standard,
               dif_standard = dif_standard,
               over_disp_standard = over_disp_standard,
               converge_standard = converge_standard,
               twiceloglik_standard = twiceloglik_standard,

               tau = tau, Q = Q, Qhet = Qhet,
               Qinc = Qinc, b = b, l = l, u = u,
               leverage = leverage, estand = estand, p.score = p.score,
               dif = dif, over_disp = over_disp, converge = converge,
               twiceloglik = twiceloglik, LRT = LRT)



   class(res) <- "NMAoutlier.rsv"

   res

}




