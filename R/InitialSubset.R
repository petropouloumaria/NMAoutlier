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
#' @param t1.label numbers to treatment 1 IDs.
#' @param t2.label numbers to treatment 2 IDs.
#' @return An initial clean subset of studies.
#'


InitialSubset <- function(TE, seTE, treat1, treat2, studlab,
                          crit1, studies, P, reference,
                          t1.label, t2.label) {


  ## Set the number of studies (n) to be equal to the number of treatments
  ## Examine a large number (P) of candidate initial subsets
  ## Create random subsets which include all treatments with size equal to the number of treatments


  ## Calculate the number of cores
  no_cores <- max(1, parallel::detectCores())
  cl <- parallel::makeCluster(no_cores)

  parallel::clusterExport(cl=cl, varlist=c("t1.label", "t2.label", "studlab", "studies",
                                           "Indices",
                                           "netconnection", "treat1", "treat2", "sub",
                                           "netmeta", "TE", "seTE", "reference",
                                           "createB", "prepare",
                                           "Multi", "multiarm",
                                           "crit1"), envir=environment())


  ## create P initial subsets with parallel computations
  paral <- parallel::parLapply(cl, 1:P, function(z) {

    ## ensuring that the network is connected
    repeat{
      ## getting random subset from studies

      ## ensuring that random subset includes studies that compare all treaments
      repeat {
        ## get a random subset of studies with length equal to variable studies
        random <- sample(unique(studlab), studies, replace=FALSE)

        ## indices of the random subset
        chosen_ind <- Indices(studlab, random)

        ## check if random subset includes studies that compare all treaments
        if(all(unique(c(t1.label, t2.label)) %in% c(t1.label[chosen_ind], t2.label[chosen_ind])))
          break
      }

      ran <- random

      ## indices of the subset (studlab)
      sub <- Indices(studlab, ran)

      ## Check if the subset is a connected network
      is_connected <- netconnection(treat1, treat2, studlab, subset = sub)$n.subnets == 1

      if(is_connected) break
    }

    ## Conduct network meta-analysis (NMA) with random effects model, Ruecker model
    netm <- netmeta(TE, seTE, treat1, treat2, studlab, reference.group = reference, subset = sub)

    ## create design matrix
    B <- createB(netm$treat1.pos, netm$treat2.pos, netm$n)

    ## summary estimate
    b <- netm$TE.random[, reference]

    ## standardized residual
    st.m <- sqrt(netm$w.random) * (netm$TE - B %*% b)

    ## multi-arm studies
    studlb <- studlab[sub]                      # studylab of set
    multi <- unique(studlb[duplicated(studlb)]) # studylab of multi-arm studies

    ## weights of random effects model
    wr.m <- prepare(TE[sub], seTE[sub], treat1[sub], treat2[sub], studlab[sub], netm$tau^2)$weights

    ## Compute standardized residuals and log-likelihood contributions
    r <- Multi(studlab[sub], st.m, wr.m)

    ## Compute the median of the absolute standardized residuals for each subset
    if (crit1 == "R") {
      minmax <- stats::median(abs(r$res))
    } else if (crit1 == "L") {
      minmax <- stats::median(abs(r$logl))
    }

    result <- list(ran=ran, minmax=minmax)
  })

  ## close the cluster after done
  parallel::stopCluster(cl)

  ## Generate a list of length P
  rand <- vector("list", P)
  ## Generate a vector of length P
  mimax <- rep(NA, P)

  result <- as.data.frame(do.call(rbind, paral))
  rand <- result$ran
  mimax <- result$minmax

  if (crit1 == "R")
    pick <- which.min(mimax)
  else if (crit1 == "L")
    pick <- which.max(mimax)

  res <- list(set = rand[[pick[1]]])

  res

}
