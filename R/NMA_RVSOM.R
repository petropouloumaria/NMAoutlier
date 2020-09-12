#' Random shift variance network meta-analysis (RVSOM NMA) model.
#'
#' Conduct RVSOM NMA model.
#'
#' @param TE Estimate of treatment effect, i.e. difference between
#'   first and second treatment (e.g. log odds ratio, mean difference,
#'   or log hazard ratio).
#' @param seTE Standard error of treatment estimate.
#' @param treat1 Label/Number for first treatment.
#' @param treat2 Label/Number for second treatment.
#' @param studlab Study labels (important when multi arm studies are
#'   included).
#' @param reference Reference treatment group.
#' @param small.values A character string indicating if small values
#'   are considered beneficial (option:"good") or harmfull
#'   (option:"bad") on outcome.  This is requirement for p-scores
#'   computation. The default value is considered benefial outcome
#'   ("good").
#' @param observations (optional) (default: NA) The study (or studies) to fit by shift its (their) variance.
#'   the default value is NA and shift variance model is fitted for each study.
#' @return A list containing the following components:
#'
#' @keywords internal
#'
#' @importFrom netmeta netmeta decomp.design netrank netsplit
#' @importFrom stats optim optimize



NMA_RVSOM <- function(TE, seTE, treat1, treat2, studlab,
                       reference, small.values, observations)
{

## network meta-analysis
  model <- netmeta(TE, seTE, treat1, treat2, studlab,
                   reference.group = reference)


  tr1 <- model$treat1.pos # treatment positions
  tr2 <- model$treat2.pos

  nt <-  model$n # number of treatments


  ## Xdesign is the design matrix, the edge-vertex incidence matrix (mxn)
  ##
  Xdesign <- createB(tr1, tr2, nt)

  ## Generalized Dersimonial and Laird heterogeneity estimator
  tau2hat <- model$tau^2

  ## starting s^2? not understand this
  ## s^2=0.1

  ## variance
  variance <- tau2hat + 0.1


  if (is.na(observations)) {

    ## transform variance to log value for starting value
    start <- rep(log(variance), 1)

    # why is missing x, what is x as argument in function REML_MA?
    inference <- optimize(NMA_REML, TE = TE, seTE = seTE, treat1 = treat1,
                          treat2 = treat2, studlab = studlab, Xdesign = Xdesign,
                          tr1 = tr1, tr2 = tr2, nt = nt,
                          reference = reference,
                          observations = observations, interval = c(-10,3))

   # inferenceML <- optimize(NMA_ML, TE = TE, seTE = seTE, treat1 = treat1,
   #                       treat2 = treat2, studlab = studlab, Xdesign = Xdesign,
   #                       tr1 = tr1, tr2 = tr2, nt = nt,
    #                      reference = reference,
    #                      observations = observations, interval = c(-10,3))
    #exp(inferenceML$minimum)

    ## from inference we have the REML estimator in log scale
    ## the estimate of the between-study variance
    tau2est <- exp(inference$minimum)


    ## network meta-analysis with specified heterogeneity, tau2est
    model <- netmeta(TE, seTE, treat1, treat2, studlab,
                     comb.random = TRUE, reference.group = reference,
                     tau.preset = sqrt(tau2est))


    y.m <- model$TE

    ##
    ##
    b <- model$TE.random[, reference]    # summary estimate of treatment effects
    l <- model$lower.random[, reference] # lower confidence interval of treatment effects
    u <- model$upper.random[, reference] # upper confidence interval of treatment effects

    ## Laplacian matrix
    ##
    L.random <- t(Xdesign) %*% diag(model$w.random) %*% Xdesign
    Lplus.random <- solve(L.random - 1 / nt) + 1 / nt
    G.random <- Xdesign %*% Lplus.random %*% t(Xdesign)


    ## Hat matrix (random-effects model)
    ##
    H.random <- G.random %*% diag(model$w.random)
    ##
    ## leverage
    leverage <- diag(H.random)


    ## Q statistics
    ##
    Q <- model$Q                                         # overall
    Qhet <- model$Q.heterogeneity                        # within-designs
    Qinc <- decomp.design(model)$Q.inc.random$Q          # between-designs


    ## Standardized residuals of pairwise comparisons
    ##
    stand <- sqrt(model$w.random) * (y.m - Xdesign %*% b)

    ## standardized residulas for each study
    estand <- res_multi(studlab, stand)$res



    ## P-scores (random effects model)
    ##
    p.score <- netrank(model, small.values)$Pscore.random

    ## Node splitting
    ## Split direct from indirect evidence
    ## Take the difference between direct and indirect evidence
    ## z-values of test for disagreement (direct versus indirect)
    ##
    dif <- netsplit(model)$compare.random$z

    ## over-dispersion parameter
    over_disp <- NA

    ## convergence diagnostic (should be zero)
    converge <- 0

    ## twice the maximum log-likelihood
    twiceloglik <- -inference$objective


    res <- list(tau2est = tau2est,
                Q = Q, Qhet = Qhet, Qinc = Qinc,
                b = b, l = l, u = u,
                leverage = leverage,
                estand = estand,
                p.score = p.score,
                dif = dif,
                over_disp =  over_disp,
                converge = converge, twiceloglik = twiceloglik,
                call = match.call())

  } else if (length(observations) > 0) {

    ## transform variance to log value for starting value
    start <- rep(log(variance), length(observations) + 1)



      inference <- optim(start, NMA_REML, TE = TE, seTE = seTE, treat1 = treat1,
                        treat2 = treat2, studlab = studlab,  Xdesign = Xdesign,
                        tr1 = tr1, tr2 = tr2, nt = nt, reference = reference,
                        observations = observations
                        # lower = -10, upper = 3
                        )


          ## find studylab of the study observation
          id_study <- which(studlab == observations)

          ## shift variance for the selected observation-study
          for (i in 1:length(id_study))
            seTE[id_study[i]] <- seTE[id_study[i]] + exp(inference$par[2])


      ## from inference we have the REML estimator in log scale
      ## the estimate of the between-study variance
      tau2est <- exp(inference$par[1])

      ## over-dispersion parameter
      over_disp <- exp(inference$par)[-1]

      ## converge
      converge <- inference$convergence

      ## twice the restricted maximum log-likelihood
      twiceloglik <- -inference$value

      ## network meta-analysis with specified heterogeneity, tau2est
      model <- netmeta(TE, seTE, treat1, treat2, studlab,
                       comb.random = TRUE, reference.group = reference,
                        tau.preset = sqrt(tau2est))


      ##
      ##
      b <- model$TE.random[, reference]    # summary estimate of treatment effects
      l <- model$lower.random[, reference] # lower confidence interval of treatment effects
      u <- model$upper.random[, reference] # upper confidence interval of treatment effects

      ## Laplacian matrix
      ##
      L.random <- t(Xdesign) %*% diag(model$w.random) %*% Xdesign
      Lplus.random <- solve(L.random - 1 / nt) + 1 / nt
      G.random <- Xdesign %*% Lplus.random %*% t(Xdesign)

      ## Hat matrix (random-effects model)
      ##
      H.random <- G.random %*% diag(model$w.random)

      ## leverage
      leverage <- diag(H.random)


      ## Q statistics
      ##
      Q <- model$Q                                         # overall Q statistic
      Qhet <- model$Q.heterogeneity                        # within-designs Q statistic
      Qinc <- decomp.design(model)$Q.inc.random$Q          # between-designs Q statistic


      ## Standardized residuals
      ##
      stand <- sqrt(model$w.random) * (model$TE - Xdesign %*% b)

      estand <- res_multi(studlab, stand)$res


      ## P-scores (random effects model)
      ##
      p.score <- netrank(model, small.values)$Pscore.random

      ## Node splitting
      ## Split direct from indirect evidence
      ## Take the difference between direct and indirect evidence
      ## z-values of test for disagreement (direct versus indirect)
      ##
      dif <- netsplit(model)$compare.random$z



      res <- list(tau2est = tau2est,
                  Q = Q, Qhet = Qhet, Qinc = Qinc,
                  b = b, l = l, u = u,
                  leverage = leverage,
                  estand = estand,
                  p.score = p.score,
                  dif = dif,
                  over_disp = over_disp,
                  converge = converge, twiceloglik = twiceloglik,
                  call = match.call())

  }

  res
}


