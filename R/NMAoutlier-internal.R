
setchar <- function(x, val, text, list = FALSE, name = NULL) {

  ## Bug fixes and code cleaning in 'default' functions
  ##
  ## Based on R function setchar from R package meta
  ## Author: G. Schwarzer <sc@imbi.uni-freiburg.de>
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  nval <- length(val)
  ##
  idx <- charmatch(tolower(x), tolower(val), nomatch = NA)
  ##
  if (anyNA(idx) || any(idx == 0)) {
    if (list)
      first <- "List element '"
    else
      first <- "Argument '"
    ##
    if (missing(text)) {
      if (nval == 1)
        vlist <- paste('"', val, '"', sep = "")
      else if (nval == 2)
        vlist <- paste('"', val, '"', collapse = " or ", sep = "")
      else
        vlist <- paste(paste('"', val[-nval],
                             '"', collapse = ", ", sep = ""),
                       ', or ', '"', val[nval], '"', sep = "")
      ##
      stop(first, name, "' should be ", vlist, ".",
           call. = FALSE)
    }
    else
      stop(first, name, "' ", text, ".", call. = FALSE)
  }
  ##
  val[idx]
}


chkchar <- function(x, single = TRUE, name = NULL, nchar = NULL) {

  ## Check character variable
  ##
  ## Based on R function chkchar from R package meta
  ## Author: G. Schwarzer <sc@imbi.uni-freiburg.de>
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (single) {
    if (length(x) != 1 || !is.character(x))
      stop(paste("Argument '", name,
                 "' must be a character string.", sep = ""))
    ##
    if (!is.null(nchar) && nchar(x) != nchar)
      if (nchar == 1)
        stop(paste("Argument '", name,
                   "' must be a single character.", sep = ""))
      else
        stop(paste("Argument '", name,
                   "' must be a character string of length ",
                   nchar, ".", sep = ""))
  }
  else if (!is.character(x))
    stop(paste("Argument '", name,
               "' must be a character vector.", sep = ""))
  else {
    if (!is.null(nchar) & any(nchar(x) != nchar))
      if (nchar == 1)
        stop(paste("Argument '", name,
                   "' must be a vector of single characters.", sep = ""))
      else
        stop(paste("Argument '", name,
                   "' must be a character vector where each element has ", nchar,
                   " characters.", sep = ""))
  }
}


chknumeric <- function(x, min, max, zero = FALSE, single = FALSE,
                       name = NULL) {

  ## Check numeric variable
  ##
  ## Based on R function chknumeric from R package meta
  ## Author: G. Schwarzer <sc@imbi.uni-freiburg.de>
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  x <- x[!is.na(x)]
  if (length(x) == 0)
    return(invisible(NULL))
  ##
  if (!is.numeric(x))
    stop("Non-numeric value for argument '", name, "'.",
         call. = FALSE)
  ##
  if (single & length(x) != 1)
    stop("Argument '", name, "' must be a numeric of length 1.",
         call. = FALSE)
  ##
  if (!missing(min) & missing(max)) {
    if (zero & min == 0 & any(x <= min, na.rm = TRUE))
      stop("Argument '", name, "' must be positive.",
           call. = FALSE)
    else if (any(x < min, na.rm = TRUE))
      stop("Argument '", name, "' must be larger equal ",
           min, ".", call. = FALSE)
  }
  ##
  if (missing(min) & !missing(max)) {
    if (zero & max == 0 & any(x >= max, na.rm = TRUE))
      stop("Argument '", name, "' must be negative.",
           call. = FALSE)
    else if (any(x > max, na.rm = TRUE))
      stop("Argument '", name, "' must be smaller equal ",
           min, ".", call. = FALSE)
  }
  ##
  if ((!missing(min) & !missing(max)) &&
      (any(x < min, na.rm = TRUE) | any(x > max, na.rm = TRUE)))
    stop("Argument '", name, "' must be between ",
         min, " and ", max, ".", call. = FALSE)
  ##
  invisible(NULL)
}

chkclass <- function(x, class, name = NULL) {

  ##
  ## Check class of R object
  ##
  ## Based on R function chkclass from R package meta
  ## Author: G. Schwarzer <sc@imbi.uni-freiburg.de>
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (!inherits(x, class))
    stop("Argument '", name,
         "' must be an object of class \"",
         class, "\".", call. = FALSE)
  ##
  invisible(NULL)
}

formatN <- function(x, digits = 2, text.NA = "--", big.mark = "") {

  ##
  ## format of a vector
  ##
  ## Based on R function fortatN from R package meta
  ## Author: G. Schwarzer <sc@imbi.uni-freiburg.de>
  ##
  outdec <- options()$OutDec

  res <- format(ifelse(is.na(x),
                       text.NA,
                       formatC(x, decimal.mark = outdec,
                               format = "f", digits = digits,
                               big.mark = big.mark)
                       )
                )
  ##
  res <-  rmSpace(res, end = TRUE)
  ##
  res
}


rmSpace <- function(x, end = FALSE, pat = " ") {

  ##
  ## Based on R function rmSpace from R package meta
  ## Author: G. Schwarzer <sc@imbi.uni-freiburg.de>
  ##
  if (!end) {
    while (any(substring(x, 1, 1) == pat, na.rm = TRUE)) {
      sel <- substring(x, 1, 1) == pat
      x[sel] <- substring(x[sel], 2)
    }
  }
  else {
    last <- nchar(x)

    while (any(substring(x, last, last) == pat, na.rm = TRUE)) {
      sel <- substring(x, last, last) == pat
      x[sel] <- substring(x[sel], 1, last[sel] - 1)
      last <- nchar(x)
    }
  }

  x
}


createB <- function (pos1, pos2, ncol) {

  ##
  ## Create full edge-vertex incidence matrix
  ##
  ## Based on R function createB from R package netmeta
  ## Author: G. Schwarzer <sc@imbi.uni-freiburg.de>
  ##
  if (missing(pos1) | missing(pos2)) {

    nrow <- choose(ncol, 2)
    B <- matrix(0, nrow = nrow, ncol = ncol)
    ##
    i <- 0
    ##
    for (pos1.i in 1:(ncol - 1)) {
      for (pos2.i in (pos1.i + 1):ncol) {
        i <- i + 1
        B[i, pos1.i] <-  1
        B[i, pos2.i] <- -1
      }
    }
  }
  else {
    ##
    ## Create edge-vertex incidence matrix
    ##
    nrow <- length(pos1)
    ncol <- length(unique(c(pos1, pos2)))
    ##
    B <- matrix(0, nrow = nrow, ncol = ncol)
    ##
    for (i in 1:nrow) {
      B[i, pos1[i]] <-  1
      B[i, pos2[i]] <- -1
    }
  }


  B
}

prepare <- function(TE, seTE, treat1, treat2, studlab, tau ) {

  ##
  ## Ordering data set
  ##
  ## Based on R function prepare from R package netmeta
  ## Author: G. Schwarzer <sc@imbi.uni-freiburg.de>
  ##
  weights <- 1 / (seTE^2 + tau^2)

  data <- data.frame(studlab,
                     treat1, treat2,
                     treat1.pos = NA, treat2.pos = NA,
                     TE, seTE, weights,
                     narms = NA, stringsAsFactors = FALSE)
  ##
  ## Ordering data set
  ##
  o <- order(data$studlab, data$treat1, data$treat2)
  ##
  data <- data[o, ]
  ##
  ## Adapt numbers to treatment IDs
  ##
  names.treat <- sort(unique(c(data$treat1, data$treat2)))
  data$treat1.pos <- match(data$treat1, names.treat)
  data$treat2.pos <- match(data$treat2, names.treat)

  newdata <- data[1, ][-1, ]
  ##
  sl <- unique(data$studlab)
  ##
  ## Determining number of arms and adjusting weights of
  ## multi-arm studies
  ##
  for (s in sl) {
    subgraph <- data[data$studlab == s, ]
    subgraph$narms <- (1 + sqrt(8 * dim(subgraph)[1] + 1)) / 2
    if (dim(subgraph)[1] > 1)
      subgraph$weights <- 1 / multiarm(1 / subgraph$weights)$v # Reciprocal new weights
    ##
    newdata <- rbind(newdata, subgraph)
  }
  res <- newdata
  ##
  res$order <- o
  ##
  res
}

multiarm <- function(r) {

  ##
  ## Measures in case of multi arm studies
  ##
  ## Based on R function multiarm from R package netmeta
  ## Author: G. Schwarzer <sc@imbi.uni-freiburg.de>
  ##
  m <- length(r)                 # Number of edges
  n <- (1 + sqrt(8 * m + 1)) / 2 # Number of vertices
  ##
  ## Construct adjacency matrix and edge vertex incidence matrix of
  ## complete graph of dimension n
  ##
  A <- 1 - diag(rep(1, n))
  B <- createB(ncol = n)
  ##
  ## Distribute the edge variances on a symmetrical n x n matrix, R
  ##
  R <- diag(diag(t(B) %*% diag(r) %*% B)) - t(B) %*% diag(r) %*% B
  ##
  ## Construct pseudoinverse Lt from given variance (resistance) matrix R
  ## using a theorem equivalent to Theorem 7 by Gutman & Xiao
  ## Lt <- -0.5 * (R - (R %*% J + J %*% R) / n + J %*%R %*% J / n^2)
  ##
  Lt <- -0.5 * t(B) %*% B %*% R %*% t(B) %*% B / n^2
  ##
  ## Compute Laplacian matrix L from Lt
  ##
  L <- solve(Lt - 1 / n) + 1 / n
  ##
  ## Compute weight matrix W and variance matrix V from Laplacian L
  ##
  W <- diag(diag(L)) - L
  V <- 1 / W
  ##
  ## Compute original variance vector v from V
  ##
  v <- rep(0, n)
  edge <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      edge <- edge + 1
      v[edge] <- V[i, j]
    }
  }
  ##
  ## Result
  ##
  res <- list(n = n, r = r, R = R, Lt = Lt, L = L, W = W, V = V, v = v)
  res
}
