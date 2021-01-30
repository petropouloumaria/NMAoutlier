#' Network meta-analysis comparing interventions for actinic keratosis
#'
#' @description
#' Network meta-analysis dataset for comparing interventions for actinic keratosis.
#'
#' @details
#' The dataset compares the relative effects of nine interventions:
#' \itemize{
#' \item placebo / vehicle (including placebo-PDT) (treatment 1),
#' \item diclofenac 3 percent in 2.5 percent hyaluronic acid (DCF/HA)
#'   (treatment 2),
#' \item 5-fluorouracil (5-FU) 0.5 percent (treatment 3),
#' \item imiquimod (IMI) 5 percent (treatment 4),
#' \item methyl aminolaevulinate (MAL)-PDT (treatment 5),
#' \item 5-aminolaevulinic acid (ALA)-photodynamic therapy (PDT)
#'   (treatment 6),
#' \item 5-fluorouracil (5-FU) 5.0 percent (treatment 7),
#' \item cryotherapy (treatment 8),
#' \item ingenol mebutate (IMB) 0.015-0.05 percent (treatment 9).
#' }
#'
#' The outcome is the number of individuals with participant complete
#' clearance or equivalent efficacy. These data are in contrast
#' format with effect size the odds ratio (OR).
#' The arm-level data were used in Gupta and Paquet (2013).
#'
#' @name Gupta2013
#' @aliases Gupta2013
#'
#' @docType data
#'
#' @format
#' A data frame in contrast format with the following columns:
#' \tabular{rl}{
#' \bold{\emph{logOR}}\tab log odds ratio \cr
#' \bold{\emph{selogOR}}\tab standard error of log odds ratio \cr
#' \bold{\emph{id}}\tab study ID \cr
#' \bold{\emph{t1}}\tab first treatment \cr
#' \bold{\emph{t2}}\tab second treatment
#' }
#'
#' @source
#' Gupta AK, Paquet M (2013):
#' Network meta-analysis of the outcome participant complete
#' clearance in nonimmunosuppressed participants of eight
#' interventions for actinic keratosis: a follow-up on a Cochrane
#' review.
#' \emph{British Journal of Dermatology},
#' \bold{169}, 250--9
#'
#' @keywords datasets
#'
#' @examples
#' \donttest{
#' data(Gupta2013)
#' # Conduct forward search algorithm for the network of actinic keratosis
#' #
#' FSresult <- NMAoutlier(logOR, selogOR, t1, t2, id, data = Gupta2013, n_cores = 2)
#'
#' # Plovide the forward plot for z-values from difference of direct and
#' # indirect evidence
#' #
#' fwdplot(FSresult, "nsplit")
#'
#' # Provide forward plot for Q statistic
#' #
#' fwdplot(FSresult, "Q")
#' }


NULL
