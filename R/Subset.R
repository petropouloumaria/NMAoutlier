#' Subset
#'
#' Finds a subset of studies with size equal to the initial set that includes all treatments.
#' @param m number of pairwise comparisons (edges)
#' @param n  number of treatments (vertices)
#' @param t1.label numbers to treatment 1 IDs.
#' @param t2.label numbers to treatment 2 IDs.
#' @param studlab Study labels (important when multi arm studies are included).
#' @param studies An optional vector specifying a subset of studies to be used.
#'                The default value is the number of studies.
#' @return subset of studies that includes all treatments.
#'


Subset <- function(m, n, t1.label, t2.label,
                   studlab, studies) {

  t <- NULL
  t1 <- NULL
  t2 <- NULL
  multiple_pos <- NULL
  chosen_studylab <- NULL
  chosen_pos <- NULL

  add1 <- FALSE
  add2 <- FALSE
  must_add_random <- FALSE
  multiple_already_added <- FALSE

  k <- sample.int(m, m)

  ## random sample with size equal to the number of studies
  ## we get <n because in the first run length of unique is 0

  while (length(unique(as.vector(chosen_studylab))) < studies) {

      counter <- 1
      k <- setdiff(k, chosen_pos)

      while(counter <= length(k) && (length(t) < n ||
          ( must_add_random && length(unique(as.vector(chosen_pos)))
            < studies ))) {

          ## if t does not contain t1
          add1 =! (t1.label[k[counter]] %in% t)
          add2 =! (t2.label[k[counter]] %in% t)


          if((add1 || add2) || must_add_random){

             ## checking for multiple studies with same studylab again
             ## is a multiple
             if (length(which(studlab == studlab[k[counter]])) > 1){


               ## checking if we already added this to multiple_pos before
               multiple_already_added <- k[counter] %in% multiple_pos

               ##
               if (!multiple_already_added) {


                ## keeping the other positions
                multiple_pos <- c(multiple_pos,
                                  setdiff(which(studlab == studlab[k[counter]]),
                                  k[counter]))
              }
            }

            if (!multiple_already_added) {

              ## adding the current one to the chosen
              t1 <- cbind(t1, t1.label[k[counter]])
              t2 <- cbind(t2, t2.label[k[counter]])

              ##
              if (add1) {
                t <- cbind(t, t1.label[k[counter]])
              }
              ##
              if (add2) {
                t <- cbind(t, t2.label[k[counter]])
              }

            ##
            chosen_studylab <- cbind(chosen_studylab, studlab[k[counter]])
            chosen_pos <- cbind(chosen_pos, k[counter])
            }

          } # end if add1||add2

          ##
          counter <- counter + 1
          multiple_already_added <- FALSE

       } # end while t is not full
       ##
       must_add_random <- length(as.vector(chosen_studylab)) < studies

  } # end while

  ## Add remaining comparisons of multi-arm studies
  for (i in 1:length(multiple_pos)) {

     t1 <- cbind(t1, t1.label[multiple_pos[i]])
     t2 <- cbind(t2, t2.label[multiple_pos[i]])
     ##
     ##

     chosen_studylab <- cbind(chosen_studylab, studlab[multiple_pos[i]])
     chosen_pos <- cbind(chosen_pos, multiple_pos[i])
  }


  ## random sample of studies that includes all treatments
  random <- unique(as.vector(chosen_studylab))

  list(random = random)

}

