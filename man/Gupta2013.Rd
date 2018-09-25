\name{Gupta2013}
\docType{data}
\alias{Gupta2013}
\title{
Network meta-analysis comparing interventions for actinic keratosis.
}
\description{
Example of network meta-analysis with inconsistency. The dataset comparing the relative effects of  interventions, placebo/vehicle (including placebo-PDT)(treatment 1), diclofenac 3 percent in 2.5 percent hyaluronic acid (DCF/HA) (treatment 2), 5-fluorouracil (5-FU) 0.5 percent (treatment 3), imiquimod (IMI) 5 percent (treatment 4), methyl aminolaevulinate (MAL)-PDT (treatment 5), 5-aminolaevulinic acid (ALA)-photodynamic therapy (PDT) (treatment 6), 5-fluorouracil (5-FU) 5.0 percent (treatment 7), cryotherapy (treatment 8), ingenol mebutate (IMB) 0.015–0.05 percent (treatment 9).
The outcome is the number of individuals with participant complete clearance or an equivalent efficacy.
These data are in contrast format with effect size odds ratio (OR) and its standard error.
Arm-level data can be found in Gupta and Paquet (2013).
}
\usage{data(Gupta2013)}

\format{A data frame in contrast format with the following columns:
  \describe{
  \item{TE}{Observed effect size, OR.}
  \item{seTE}{Standard error of the observed effects OR.}
  \item{studlab}{Study label of studies.}
  \item{treat1}{Treatment 1.}
  \item{treat2}{Treatment 2.}
  }
}

\source{
Gupta AK, Paquet M (2013). Network meta-analysis of the outcome “participant complete clearance” in nonimmunosuppressed participants of eight interventions for actinic keratosis: a follow-up on a Cochrane review.
Evidence Synthesis for Decision Making 4: Inconsistency in networks of evidence based on randomized controlled trials. British Journal of Dermatology 169(2):250-9.
}

\examples{
data(Gupta2013)

# Conduct forward search algorithm
# for the network of smoking cessation
\donttest{FSresult <- NMAoutlier(TE, seTE, treat1, treat2,
                                 studlab, data = Gupta2013)}

# Draw forward plot for z-values
# from difference of direct and indirect evidence
\donttest{fwdplot(FSresult,"nsplit")}

# Draw forward plot for Q statistic
\donttest{fwdplot(FSresult,"Q")}
}

\keyword{datasets}

