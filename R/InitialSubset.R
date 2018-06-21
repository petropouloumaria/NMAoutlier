#' InitialSubset
#'
#' Choose an initial clean (i.e. likely outlier-free) subset of studies.
#' @param TE   Estimate of treatment effect, i.e. difference between first and second treatment (e.g. log odds ratio, mean difference, or log hazard ratio).
#' @param seTE Standard error of treatment estimate.
#' @param treat1 Label/Number for first treatment.
#' @param treat2 Label/Number for second treatment.
#' @param studlab Study labels (important when multi arm studies are included).
#' @param crit1 A character string indicating the criterion to be used for selecting the initial subset,
#'              this criterion could be the minimum of median absolute residuals ("R") or
#'              the maximum of median absolute likelihood contributions ("L"). Default value is "R".
#' @param studies An optional vector specifying a subset of studies to be used.
#'                The default value is the number of studies.
#' @param P       An optional vector specifying the number of samples for the choice of the initial subset.
#' @param reference Reference treatment group.
#' @param m number of pairwise comparisons (edges)
#' @param n  number of treatments (vertices)
#' @param t1.label numbers to treatment 1 IDs.
#' @param t2.label numbers to treatment 2 IDs.
#' @return An initial clean subset of studies.
#'


InitialSubset <- function(TE, seTE, treat1, treat2, studlab,
                          crit1, studies, P, reference,
                          m, n, t1.label, t2.label) {


  ##
  ## Generate a list of length P
  ##
  ran <- vector("list", P)


  ## Set the number of studies to be equal with the number of treatments
  ## studies = n
  ##
  ## Examine a large number (P) of candidate initial subsets
  ##
  ## Create random subsets with size equal to the number of treatments
  ## that include all treatments

  ## Set up a cluster

  ## Calculate the number of cores
  no_cores <- max(1, parallel::detectCores())
  cl <- parallel::makeCluster(no_cores)

  parallel::clusterExport(cl=cl, varlist=c("Subset", "m", "n", "t1.label", "t2.label", "studlab", "studies",
                                           "Indices", "ran",
                                           "netconnection", "treat1", "treat2", "sub",
                                           "netmeta", "TE", "seTE", "reference",
                                           "createB", "prepare",
                                           "Multi", "multiarm",
                                           "crit1"
  ))
  ##
  ## Create P initial subset with parallel computations
  ##
  ##

  to_parallel <- function(z) {
    ##
    ##
    minmax <- NA


    ## studies of the subset
    ran <- Subset(m, n, t1.label, t2.label, studlab, studies)$random

    ## indices of the subset
    sub <- Indices(studlab, ran)

    ## Check if the subset is a connected network
    codit <- netconnection(treat1, treat2, studlab, subset = sub)$n.subnets == 1

    ## Take the subset if it is a connected network
    ##
    if (codit) {

      ## Conduct network meta-analysis (NMA) with random effects
      ## model, Rucker model
      ##
      netm <- netmeta(TE, seTE, treat1, treat2, studlab,
                      reference.group = reference, subset = sub)

      B <- createB(netm$treat1.pos, netm$treat2.pos, netm$n)
      b <- netm$TE.random[, reference]
      st.m <- sqrt(netm$w.random) * (netm$TE - B %*% b)

      ## multi-arm studies
      studlb <- studlab[sub]                      # studylab of set
      multi <- unique(studlb[duplicated(studlb)]) # studylab of multi-arm studies

      ## weights of random effects model
      wr.m <- prepare(TE[sub], seTE[sub], treat1[sub], treat2[sub],
                      studlab[sub], netm$tau^2)$weights

      ## Compute standardized residuals and log-likelihood contributions
      ##
      r <- Multi(studlab[sub], st.m, wr.m)
      ##
      ## Compute the median of the absolute standardized residuals for
      ## each subset
      ##,
      if (crit1 == "R")
        minmax <- stats::median(abs(r$res))
      else if (crit1 == "L")
        minmax <- stats::median(abs(r$logl))

      ##
      resul <- list(ran=ran, minmax=minmax)
      return(resul)
      ##

    } # End if

  }


  paral <- parallel::parLapply(cl, 1:P, to_parallel )

  ## close the cluster after done
  parallel::stopCluster(cl)

  ## Generate a vector of length P
  mnax <- rep(NA, P)


  if(any(sapply(paral, is.null))){
    par <- paral[-which(sapply(paral, is.null))]
    valid_length <- length(par)
    to_merge <- list()
    for(i in 1:(P-valid_length)){
      to_merge[[i]] <- to_parallel()
    }
    paral <- c(par, to_merge)
  }


  for(j in 1:P){

    ran[[j]] <- paral[[j]]$ran

    mnax[j] <- paral[[j]]$minmax

  }


  ##
  ##
  if (crit1 == "R")
    pick <- which.min(mnax)
  else if (crit1 == "L")
    pick <- which.max(mnax)


  res <- list(set = ran[[pick[1]]])


  ##
  res


}
