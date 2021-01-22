#' Restricted Maximum log-likelihood (REML) function for NMA model.
#'
#' Conduct REML log-likelihood function for NMA model.
#'
#' @param TE Estimate of treatment effect, i.e. difference between
#'   first and second treatment (e.g. log odds ratio, mean difference,
#'   or log hazard ratio).
#' @param seTE Standard error of treatment estimate.
#' @param treat1 Label/Number for first treatment.
#' @param treat2 Label/Number for second treatment.
#' @param studlab Study labels (important when multi arm studies are
#'   included).
#' @param Xdesign design matrix.
#' @param tr1 positions for the first treatment.
#' @param tr2 positions for the second treatment.
#' @param nt number of treatments.
#' @param observations (optional) (default: NA) The study (or studies) to fit by shift its (their) variance.
#'   the default value is NA and shift variance model is fitted for each study.
#' @return REML log-likelihood function.
#'
#' @keywords internal



NMA_REML <- function(x, TE, seTE, treat1, treat2, studlab,
                     Xdesign, tr1, tr2, nt, observations)
{


  if (length(observations) > 0)
  {

    id_study <- which(studlab == observations)

    # shift variance for study
    for (i in 1:length(id_study))

      seTE[id_study[i]] <- seTE[id_study[i]] + exp(x[2])


  }

  # heterogeneity variance estimator
  tau2 <- exp(x[1])


  p2 <- prepare(TE, seTE, treat1, treat2, studlab, tau = sqrt(tau2))


  w.pooled <- p2$weights

  W <- diag(w.pooled,                    # Weighted degree diagonal matrix
            nrow = length(w.pooled))     #


  ## L is the weighted Laplacian (Kirchhoff) matrix (n x n)
  L.matrix <- t(Xdesign) %*% W %*% Xdesign

  ## Lplus is its Moore-Penrose pseudoinverse
  ##
  Lplus <- solve(L.matrix - 1 / nt) + 1 / nt

  ##
  G <- Xdesign %*% Lplus %*% t(Xdesign)
  H <- G %*% W
  ##
  ##
  ## Resulting effects at numbered edges
  ##
  v <- as.vector(H %*% TE)

  ##
  ## Resulting effects, all edges, as a nt x nt matrix:
  ##
  all <- matrix(NA, nrow = nt, ncol = nt)
  for (i in 1:length(TE)) {
       all[tr1[i], tr2[i]] <- v[i]
  }
  for (i in 1:nt) {
      for (j in 1:nt) {
           for (k in 1:nt) {
              if (!is.na(all[i, k]) & !is.na(all[j, k])) {
                 all[i, j] <- all[i, k] - all[j, k]
                 all[j, i] <- all[j, k] - all[i, k]
               }
              if (!is.na(all[i, j]) & !is.na(all[k, j])) {
                 all[i, k] <- all[i, j] - all[k, j]
                 all[k, i] <- all[k, j] - all[i, j]
              }
              if (!is.na(all[i, k]) & !is.na(all[i, j])) {
                 all[j, k] <- all[i, k] - all[i, j]
                 all[k, j] <- all[i, j] - all[i, k]
              }
           }
       }
   }



   ## summary estimate of RE model
   theta_hat <- all[, 1]


   ## REML log-likelihood function
   - log(det(abs(W))) + t(TE -  Xdesign %*% c(theta_hat)) %*% W %*% (TE -  Xdesign %*% c(theta_hat)) + log(det(abs(L.matrix)))


}

