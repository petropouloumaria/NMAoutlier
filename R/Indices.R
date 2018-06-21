#' Indices
#'
#' Finds the indices of a set of studies (this is helpful to find the indices for the case of a multi-arm trial).
#' @param studlab Study labels (important when multi arm studies are included).
#' @param set set of studies
#' @return indices for a set of studies
#'

Indices <- function(studlab, set) {

  ##
  A <- st <- list()
  ##
  ## Studylab in case of multi-arm study
  ##
  for (j in 1:length(set)) {

      ##
      A <- which(studlab == set[j], arr.ind = TRUE)
      ##
      for (z in 1:length(A))
           st <- cbind(st, A[z])

      }

  ## Studylab for the set
  ##
  indices <- unlist(st)
  ##

  c(indices)

}


