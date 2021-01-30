##
## (1) Make R packages available
##
library(devtools)
library(roxygen2)


##
## (2) Create documentation file(s)
##
document("../NMAoutlier") # Also considers datasets in subdirectory NMAoutlier/data


##
## (3) Build R package and PDF file with help pages
##
build("../NMAoutlier")
build_manual("../NMAoutlier")


##
## (4) Install R package
##
install("../NMAoutlier")


##
## (5) Check R package
##
check("../NMAoutlier")


##
## (6) Check R package (with dontrun examples)
##
check("../NMAoutlier", run_dont_test = TRUE)
