## Auxiliary functions to conduct network meta-analysis
##
## Package: netmeta
## Authors: Gerta Rücker <ruecker@imbi.uni-freiburg.de>,
##          Guido Schwarzer <sc@imbi.uni-freiburg.de>
## License: GPL (>= 2)
##
createB <- function(pos1, pos2, ncol, aggr = FALSE) {
  
  
  if (!aggr) {
    if (missing(pos1) | missing(pos2)) {
      ##
      ## Create full edge-vertex incidence matrix
      ##
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
  }
  else {
    nrow <- 0
    ##
    ## Determine number of edges (no. of rows of B)
    ##
    for (i in 1:(ncol - 1)) {
      for (j in (i + 1):ncol) {
        ij.count <- 0
        ## Cycle through every possible edge ij
        ## Search pos1 and pos2 to see if at least one of these
        ## combinations is ij
        for (k in seq_along(pos1)) {
          if (pos1[k] == i & pos2[k] == j) {
            ij.count <- ij.count + 1
          }
          else {
            ij.count <- ij.count
          }
        }
        if (ij.count > 0)
          nrow <- nrow + 1
        else
          nrow <- nrow
      }
    }
    ##
    ## Create aggregate B matrix with dimensions e x n
    ##
    B <- matrix(0, nrow = nrow, ncol = ncol)
    ##
    r <- 0
    ## Cycle through each possible pairwise comparison ij
    for (i in 1:(ncol - 1)) {
      for (j in (i + 1):ncol) {
        ij.count <- 0
        for (k in 1:length(pos1)) {
          ## If there is an edge for that pairwise comparison ...
          if (pos1[k] == i & pos2[k] == j)
            ij.count <- ij.count + 1 # ...then ij.count is no longer = 0 ...
          else
            ij.count <- ij.count
        }
        if (ij.count > 0) {
          ## ...and we add this row to B
          r <- r + 1
          B[r, i] <-  1
          B[r, j] <- -1
        }
      }
    }
  }
  
  B
}
multiarm <- function(r) {
  ##
  ## Dimension of r and R
  ##
  m <- length(r) # Number of edges
  ##
  k <- (1 + sqrt(8 * m + 1)) / 2 # Number of vertices
  if (!(abs(k - round(k)) < .Machine$double.eps^0.5))
    stop("Wrong number of comparisons in multi-arm study.", call. = FALSE)
  ##
  ## Construct edge-vertex incidence matrix of complete graph of
  ## dimension k
  ##
  B <- createB(ncol = k)
  ##
  ## Distribute the edge variances on a symmetrical matrix R of
  ## dimension k x k
  ##
  R <- diag(diag(t(B) %*% diag(r, nrow = m) %*% B)) -
    t(B) %*% diag(r, nrow = m) %*% B
  ##
  ## Construct pseudoinverse Lt from given variance (resistance) matrix R
  ## using a theorem equivalent to Theorem 7 by Gutman & Xiao
  ## Lt <- -0.5 * (R - (R %*% J + J %*% R) / k + J %*%R %*% J / k^2)
  ##
  Lt <- -0.5 * t(B) %*% B %*% R %*% t(B) %*% B / k^2
  ##
  ## Compute Laplacian matrix L from Lt
  ##
  L <- invmat(Lt)
  ##
  ## Compute weight matrix W and variance matrix V from Laplacian L
  ## 
  W <- diag(diag(L)) - L
  ##
  ## Replace small negative weights with zeros
  ## (i.e., if an absolute negative weight contributes less than 0.1%)
  ##
  W[W < 0 & (abs(W) / sum(abs(W)[lower.tri(W)])) < 0.001] <- 0
  ##
  V <- 1 / W
  ##
  ## Compute original variance vector v from V
  ##
  v <- rep(0, m)
  edge <- 0
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      edge <- edge + 1
      v[edge] <- V[i, j]
    }
  }
  ##
  ## Result
  ##
  res <- list(k = k, r = r, R = R, Lt = Lt, L = L, W = W, V = V, v = v)
  res
}
invmat <- function(X) {
  n <- nrow(X)
  J <- matrix(1, nrow = n, ncol = n)
  ##
  res <- solve(X - J / n) + J / n
  ##
  res
}
prepare <- function(TE, seTE, treat1, treat2, studlab, tau = 0) {
  
  if (is.na(tau))
    tau <- 0
  
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
