## NMAoutlier, version 0.1.18 (2021-10-11)

* Use hat matrix from random effects model (if available) to calculate
  leverages

* Extended check whether first argument of NMAoutlier() or
  NMAoutlier.measures() is R object created with pairwise() from R
  package **netmeta**

* Add argument '...' to provide additional arguments to internal calls
  of netmeta()

* Use Markdown for NEWS


## NMAoutlier, version 0.1.17

* code refactoring

* added n_cores to examples to comply with CRAN limitation to 2 cores


## NMAoutlier, version 0.1.16

* added initial subset argument to NMAoutlier function


## NMAoutlier, version 0.1.15

* new functions added: NMAoutlier.measures, measplot, Qnetplot. The
  functions provide the calculation of several outlier and influential
  measures


## NMAoutlier, version 0.1.14

* new functions added: Qnetplot. The functions provide the calculation
  of several outlier and influential measures


## NMAoutlier, version 0.1.13

* added n_cores argument to limit cores to 2 for small toy examples
  that run in < 5 sec for CRAN submission


## NMAoutlier, version 0.1.12

* added maintainer to the description file


## NMAoutlier, version 0.1.11

* DESCRIPTION: ORCIDs added

* help pages updated (fix error for keratosis)

* datasets renamed


## NMAoutlier, version 0.1.10

* added small examples for CRAN submission


## NMAoutlier, version 0.1.7

* new datasets Laparoscopic and Gupta2013 added


## NMAoutlier, version 0.1.6

* changed matplot to ggplot2 for the graphs of fwdplot and fwdplotest


## NMAoutlier, version 0.1.5

* added logo

* added README file

*added data (Dias 2013)

* deleted Indices function and implemented a solution inline where
  needed


## NMAoutlier, version 0.1.4

* performance improvements for Subset

* removed the Subset function and integrated the code to the
  InitialSubset file


## NMAoutlier, version 0.1.3

* initial release on GitHub
