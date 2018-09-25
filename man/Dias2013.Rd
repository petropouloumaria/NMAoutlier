\name{Dias2013}
\docType{data}
\alias{Dias2013}
\title{
Network meta-analysis comparing smoking cessation counseling programs.
}
\description{
Example of outlying detection in network meta-analysis comparing the relative effects of four smoking cessation counseling programs, no contact (A), self-help (B), individual counseling (C) and group counseling (D).
The outcome is the number of individuals with successful smoking cessation at 6 to 12 months.
These data are in contrast format with effect size odds ratio (OR) and its standard error.
Arm-level data can be found in Dias et al.(2013).
}
\usage{data(Dias2013)}

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
Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G and Ades AE (2013). Evidence Synthesis for Decision Making 4: Inconsistency in networks of evidence based on randomized controlled trials. Medical Decision Making \bold{33}, 641–656.
}

\examples{
data(Dias2013)

# Conduct forward search algorithm
# for the network of smoking cessation
\donttest{FSresult <- NMAoutlier(TE, seTE, treat1, treat2,
                                 studlab, data = Dias2013,
                                 small.values = "bad")}

# Draw forward plot for Cook statistic measure
\donttest{fwdplot(FSresult,"cook")}
}

\keyword{datasets}
