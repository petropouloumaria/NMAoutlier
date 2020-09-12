#' Statistical measure in multi-arm studies
#'
#' Computes the arithmetic mean of a statistical measure for the pairwise comparisons to study level data
#' (this can be done in the case of multi-arm trials).
#'
#' @param st.lab studylab of set of studies.
#' @param ei statistical measure for each pairwise comparison of studies.
#' @return the arithmetic mean of statistical measure.
#'
#' @keywords internal


res_multi <- function(st.lab, ei) {


  ## initialise the position of statistical measure calculated for each pairwise comparison
  ## in case of multi-arm studies
  ##
  em <- m.arm <- mult <- list()

  ## studylab of multi-arm studies
  ##
  multi <- unique(st.lab[duplicated(st.lab)])

  ## number of multi-arm studies
  ##
  n.multi <- length(multi)

  ## In case of multi-arm studies
  if (n.multi >= 1) {


    ## Compute the arithmetic mean of statistical measure
    ## which correspond to multi-arm study

    for (d in 1:n.multi) {

      ##
      m.arm[[d]] <- which(st.lab == multi[d])
      em[[d]] <- mean(ei[m.arm[[d]]])
      mult[[d]] <- m.arm[[d]][1]


    } #end for


    ## find the position of 2-arm studies
    ##
    pos.m <- setdiff(1:length(st.lab), unlist(m.arm))

    ##
    res <- unlist(c(em, ei[pos.m]))


    ##
    ## find the studylab corresponding in the statistical measure
    ##
    study <- c(st.lab[unlist(mult)], st.lab[pos.m])

  } else {

    ## In case that there are not exist multi-arm studies
    ##
    res <- ei

    ##
    ##
    study <- st.lab

  } #end if

  result <- list(res = res)

  result
}
