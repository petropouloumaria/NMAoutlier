#' Network meta-analysis for the whole dataset.
#'
#' Conduct network meta-analysis (Rücker model) for the whole dataset.
#'
#' @param TE Estimate of treatment effect, i.e. difference between
#'   first and second treatment (e.g. log odds ratio, mean difference,
#'   or log hazard ratio).
#' @param seTE Standard error of treatment estimate.
#' @param treat1 Label/Number for first treatment.
#' @param treat2 Label/Number for second treatment.
#' @param studlab Study labels (important when multi arm studies are
#'   included).
#' @param reference Reference treatment group
#' @param names.treat names of treatments
#'
#' @return
#' Within-study standard error, design matrix, summary estimate,
#' heterogeneity from network meta-analysis.
#'
#' @keywords internal


nma <- function(TE, seTE, treat1, treat2, studlab,
                reference, names.treat, ...) {

  ## NMA for the whole dataset
  ##
  met <- netmeta(TE, seTE, treat1, treat2, studlab,
                 reference.group = reference, ...)

  ## Use adjusted standard errors from random effects model for
  ## multi-arm studies
  ##
  if (!is.null(met$seTE.adj.random))
    within.se <- met$seTE.adj.random
  else
    within.se <- met$seTE.adj

  ## X is the design matrix, the edge-vertex incidence matrix (mxn)
  ##
  X <- createB(met$treat1.pos, met$treat2.pos, met$n)
  colnames(X) <- names.treat
  rownames(X) <- studlab

  ## Summary estimates of treatment effects
  ##
  est <- met$TE.random[, reference]
  het <- met$tau^2

  list(within.se = within.se, X = X, est = est, het = het)
}
