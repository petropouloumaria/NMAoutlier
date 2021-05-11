#' Network meta-analysis for a set of studies
#'
#' Conduct network meta-analysis (Rücker model) for a set of studies.
#'
#' @param TE Estimate of treatment effect, i.e. difference between
#'   first and second treatment (e.g. log odds ratio, mean difference,
#'   or log hazard ratio).
#' @param seTE Standard error of treatment estimate.
#' @param treat1 Label/Number for first treatment.
#' @param treat2 Label/Number for second treatment.
#' @param studlab Study labels (important when multi arm studies are
#'   included).
#' @param ind.bs indices of studies in basic set
#' @param reference Reference treatment group
#' @param small.values small.values for P-score. A character string
#'   indicating if small values are considered bad ("bad") or good
#'   ("good").
#' @param names.treat names of treatments
#' @param \dots Additional arguments passed on to
#'   \code{\link{netmeta}}.
#'
#' @return results and statistics from network meta-analysis.
#'
#' @keywords internal
#'
#' @importFrom netmeta decomp.design netmeta netrank netsplit



netmet <- function(TE, seTE, treat1, treat2, studlab,
                   ind.bs, reference, small.values, names.treat,
                   ...) {

  ## Conduct network meta-analysis (NMA) with random effects model,
  ## Rücker model
  model <- netmeta(TE, seTE, treat1, treat2, studlab,
                   comb.random = TRUE, reference.group = reference,
                   subset = ind.bs, ...)

  nt <- model$n # number of treatments

  tr1 <- model$treat1.pos # treatment positions
  tr2 <- model$treat2.pos
  y.m <- model$TE

  ## B is the design matrix, the edge-vertex incidence matrix (mxn)
  ##
  B <- createB(tr1, tr2, nt)
  ##
  t <- (model$tau)^2                   # heterogeneity
  ##
  b <- model$TE.random[, reference]    # summary estimate of treatment effects
  l <- model$lower.random[, reference] # lower confidence interval of treatment effects
  u <- model$upper.random[, reference] # upper confidence interval of treatment effects

  ## Laplacian matrix
  ##
  L.random <- t(B) %*% diag(model$w.random) %*% B
  Lplus.random <- solve(L.random - 1 / nt) + 1 / nt
  G.random <- B %*% Lplus.random %*% t(B)

  ## Hat matrix (random-effects model)
  ##
  H.random <- G.random %*% diag(model$w.random)


  ## Computations for variance-covariance matrix (random effects model)
  ##
  ind <- which(names.treat == reference) # index of reference treatment
  t.pos1 <- rep(ind, (nt - 1))
  t.pos2 <- setdiff(1:nt, ind)
  ##
  B.r <- createB(t.pos1, t.pos2, ncol = nt) # Restricted matrix B
  Cov <- B.r %*% Lplus.random %*% t(B.r)    # Variance-covariance matrix (random effects model)


  ## Standardized residuals
  ##
  standres <- sqrt(model$w.random) * (y.m - B %*% b)

  ## P-scores (random effects model)
  ##
  p.rank <- netrank(model, small.values)$Pscore.random

  ## Node splitting
  ## Split direct from indirect evidence
  ## Take the difference between direct and indirect evidence
  ## z-values of test for disagreement (direct versus indirect)
  ##
  diff <- netsplit(model)$compare.random$statistic


  ## Q statistics
  ##
  Q <- model$Q                                         # overall
  Qhet <- model$Q.heterogeneity                        # within-designs
  Qinc <- decomp.design(model)$Q.inc.random$Q          # between-designs
  ##
  st.res <- mult(studlab, ind.bs, standres, model)


  res <- list(y.m = y.m, t = t, Q = Q, Qhet = Qhet, Qinc = Qinc,
              b = b, l = l, u = u, Cov = Cov, H.random = H.random,
              st.res = st.res, B = B, p.rank = p.rank, diff = diff)

  res
}
