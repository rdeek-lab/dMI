#' Gower's centering
#'
#' Function to convert a distance matrix into a double centered similarity matrix
#' using Gower's centering: \eqn{G = -\frac{1}{2} (I - \frac{11^T}{n}) D^2 (I - \frac{11^T}{n})}.
#'
#' @param D An \eqn{n \times n} symmetric distance matrix of class `dist` or `matrix`.
#' \eqn{n} refers to the total number of observations, or samples, in the data set.
#'
#' @returns An \eqn{n \times n} symmetric double centered similarity matrix.
#' @export
#' @examplesIf requireNamespace("vegan", quietly = TRUE)
#' mb.bray = vegan::vegdist(microbiomeData, method = "bray")
#' D.x = gower.centering(mb.bray)
gower.centering = function(D){
  if(!("dist" %in% class(D) | "matrix" %in% class(D)) ) {
    stop("ERROR: D must be a distance matrix of class 'dist' or 'matrix'.")
  } else if(nrow(as.matrix(D)) != ncol(as.matrix(D)) ) {
    stop("ERROR: D must be a square (n by n) distance matrix.")
  } else if(isSymmetric(as.matrix(D)) == FALSE) {
    stop("ERROR: D must be a symmetric distance matrix.")
  }

  # if D is of class dist, make matrix
  if("dist" %in% class(D)){
    D = as.matrix(D)
  }

  # convert D to A
  A = (-1/2)*D^2

  # centering
  n = nrow(as.matrix(D))                                   # sample size (nrow/ncol of D)
  I.mat = diag(n)                                          # n x n identify matrix
  one.vec = rep(1, n)                                      # n-length vector of 1s
  centering.mat = I.mat - (1/n) * (one.vec %*% t(one.vec)) # centering matrix

  # Gower's centered matrix
  G = centering.mat %*% A %*% centering.mat

  return(G)
}
