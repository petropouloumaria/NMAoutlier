#' mult
#'
#' Computes the arithmetic mean of standardized residuals in case of multi-arm trials.
#' @param studlab Study labels (important when multi arm studies are included).
#' @param ind.bs indices of basic set.
#' @param standres standardized residuals of each pairwise comparison.
#' @param model model assigned from netmeta
#' @return standardized residuals and arithmetic mean of standardized residuals in case of multi-arm trials.
#'
#' @keywords internal

mult <- function(studlab, ind.bs, standres, model) {


  ## studylab of set
  studlb <- studlab[ind.bs]
  ##
  ## studylab of multi-arm studies
  ##
  multi <- unique(studlb[duplicated(studlb)])
  ## number of multi-arm studies
  ##
  n.multi <- length(multi)


  ## In case of multi-arm studies
  ##

  if (n.multi >= 1) {


     ## initialize indicies for the position of multi-arm studies
     ## initialize standardized residuals of multi-arm comparisons
     ##
     em <- m.arm <- ind <- list()
     st.res <- matrix(NA, length(unique(studlb)))

     ## Compute the arithmetic mean of standardized residuals which corresponding to multi-arm study
     for (d in 1:n.multi) {

         ## position of multi arm study
         m.arm[[d]] <- which(studlb == multi[d])

         ## studylab of multi-arm study
         study <- studlb[m.arm[[d]]][1]

         ## re-define the position of multi-arm study
         ##
         ind[[d]] <- match(study, unique(studlb))

         ##
         em <- mean(standres[m.arm[[d]]])

         ##
         st.res[ind[[d]],] <- em

     } #end for

     ## find the position of 2-arm studies
     pos.m <- setdiff(1:length(studlb), unlist(m.arm))

     ## re-define the position of 2-arm studies
     stud <- setdiff(1:length(unique(studlb)), unlist(ind))

     ##
     st.res[stud,] <- standres[pos.m]


 } else {


    ## standardized residuals
    ##
    st.res <- standres
    ##


 }

 st.res


}

