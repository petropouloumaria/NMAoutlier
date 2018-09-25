\name{Laparoscopic}
\docType{data}
\alias{Laparoscopic}
\title{
Network meta-analysis comparing the effects after Laparoscopic Heller myotomy (LHM).
}
\description{
Example of network meta-analysis with inconsistency. The dataset comparing the effects after Laparoscopic Heller myotomy (LHM).
The outcome is the number of individuals with successful rates at 12 months.
These data are in contrast format with effect size odds ratio (OR) and its standard error.
Arm-level data can be found in Schoenberg et al. (2013).
}
\usage{data(Laparoscopic)}

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
Schoenberg MB, Marx S, Kersten JF, Rösch T, Belle S, Kähler G, Vassiliou MC, Lüth S, von Renteln D (2013). Laparoscopic Heller myotomy versus endoscopic balloon dilatation for the treatment of achalasia: a network meta-analysis. Annals of Surgery \bold{258}(6), 943-52.
}

\examples{
data(Laparoscopic)

# Conduct forward search algorithm
# for the network of Laparoscopic Heller myotomy
\donttest{FSresult <- NMAoutlier(TE, seTE, treat1, treat2,
                                 studlab, data = Laparoscopic)}

# Draw forward plot for z-values
# from difference of direct and indirect evidence
\donttest{fwdplot(FSresult,"nsplit")}
}

\keyword{datasets}
