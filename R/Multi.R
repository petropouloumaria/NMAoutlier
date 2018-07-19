#' Residuals and log-likelihood contributions in multi-arm trials
#'
#' Computes the arithmetic mean of standardized and log-likelihood contributions in the case of multi-arm trials.
#' @param st.lab studylab of set of studies
#' @param ei standardized residuals for each comparison from a set of studies
#' @param wr weigths of random effects model
#' @return the arithmetic mean of standardized and log-likelihood contributions.
#'
#' @keywords internal

Multi <- function(st.lab, ei, wr) {


  ## initialise the position and standardized residuals,log-likelihood contributions
  ## in case of multi-arm studies
  ##
  em <- lm <- m.arm <- mult <- list()

  ## studylab of multi-arm studies
  ##
  multi <- unique(st.lab[duplicated(st.lab)])

  ## number of multi-arm studies
  ##
  n.multi <- length(multi)

  ## In case of multi-arm studies
  if (n.multi >= 1) {


     ## Compute the arithmetic mean of the standardized residuals
     ## which correspond to multi-arm study

     for (d in 1:n.multi) {

        ##
        m.arm[[d]] <- which(st.lab == multi[d])
        em[[d]] <- mean(ei[m.arm[[d]]])
        mult[[d]] <- m.arm[[d]][1]

        ##
        ##
        loglikcon <- log(wr) - (ei)^2

        ##
        lm[[d]] <- mean(loglikcon[m.arm[[d]]])

      } #end for


      ## find the position of 2-arm studies
      ##
      pos.m <- setdiff(1:length(st.lab), unlist(m.arm))

      ##
      res <- unlist(c(em, ei[pos.m]))
      ##
      logl <- unlist(c(lm, loglikcon[pos.m]))
      ##

      ##
      ## find the studylab corresponding in residuals and log-likelihood contributions
      ##
      study <- c(st.lab[unlist(mult)], st.lab[pos.m])

  } else {

     ## In case that there are not exist multi-arm studies
     ##
     res <- ei
     logl <- log(wr) - (ei)^2
     ##
     ##
     study <- st.lab

  } #end if

  result <- list(res = res, logl = logl, study = study)

}


