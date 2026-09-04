#' Distance-based mutual information
#'
#' Function to perform estimation of distance-based mutual information using kNN
#' density estimation and a permutation test for the global null hypothesis of
#' independence.
#'
#' @param D.x An \eqn{n \times n} symmetric distance matrix of class `dist` or
#' `matrix` for modality \eqn{X}. \eqn{n} refers to the total number of observations
#' (samples) in the data set. Individuals must have observations in both \eqn{X}
#' and \eqn{Y}.
#' @param D.y An \eqn{n \times n} symmetric distance matrix of class `dist` or
#' `matrix` for modality \eqn{Y}. \eqn{n} refers to the total number of observations
#' (samples) in the data set.
#' @param kbar A vector of nearest neighbors to be used in density estimation.
#' It is recommend to use \eqn{1:k}, where \eqn{k} greater than or equal to
#' \eqn{\sqrt{n(n-1)/2}}. \eqn{k} cannot be larger than \eqn{n(n-1)/2}.
#' @param n.perm Number of permutations to be used in hypothesis testing. Default
#' is set as 999.
#' @param n.cores Number of cores to be used in parallel processing. Default is
#' set to 1 (no parallelization).
#' @param MI Logical (`TRUE`/`FALSE`) indicating if the mutual information
#' estimate should be returned. Default is `FALSE` and suggested if independence
#' testing is of primary interest.
#'
#' @returns A list containing the permutation test p-value for the global null
#' hypothesis of independence, by default. If `MI = TRUE`, the list also contains
#' the KL estimate of distance-based mutual information.
#' @export
#' @examplesIf requireNamespace("vegan", quietly = TRUE)
#' mb.bray = vegan::vegdist(microbiomeData, method = "bray")
#' host.euclidean = vegan::vegdist(hostData, method = "euclidean")
#' dMI(D.x = mb.bray, D.y = host.euclidean, kbar = 1:71)
dMI = function(D.x, D.y, kbar, n.perm = 999, n.cores = 1, MI = FALSE){
  if(!("dist" %in% class(D.x) | "matrix" %in% class(D.x)) ) {
    stop("ERROR: D.x must be a distance matrix of class 'dist' or 'matrix'.")
  } else if(nrow(as.matrix(D.x)) != ncol(as.matrix(D.x)) ) {
    stop("ERROR: D.x must be a square (n by n) distance matrix.")
  } else if(isSymmetric(as.matrix(D.x)) == FALSE) {
    stop("ERROR: D.x must be a symmetric distance matrix.")
  } else if(!("dist" %in% class(D.y) | "matrix" %in% class(D.y)) ) {
    stop("ERROR: D.y must be a distance matrix of class 'dist' or 'matrix'.")
  } else if(nrow(as.matrix(D.y)) != ncol(as.matrix(D.y)) ) {
    stop("ERROR: D.y must be a square (n by n) distance matrix.")
  } else if(isSymmetric(as.matrix(D.y)) == FALSE) {
    stop("ERROR: D.y must be a symmetric distance matrix.")
  } else if(length(kbar) == 1) {
    stop("ERROR: kbar must be a vector denoting a sequence of nearest neighbors to use for density estimation, typically 1 to k.")
  } else if(nrow(as.matrix(D.x)) != nrow(as.matrix(D.y))) {
    stop("ERROR: D.x and D.y must be distance matrices of the same dimension (i.e., are paired observations across modalities).")
  } else if(max(kbar) > nrow(as.matrix(D.x)) * (nrow(as.matrix(D.x)) - 1) / 2) {
    stop("ERROR: kbar cannot contain values larger than n(n-1)/2, where n is the number of paired observations in the data set.")
  }

  # 1. distance matrix to Gower's centered matrix
  G.x = gower.centering(D.x)
  G.y = gower.centering(D.y)

  # 2. matrix to vector
  g.x.lt.vect = as.vector(G.x[lower.tri(G.x)])
  g.y.lt.vect = as.vector(G.y[lower.tri(G.y)])

  # 3. rank and probit transformation
  g.x.lt.vect.rank = rank(g.x.lt.vect)/(length(g.x.lt.vect)+1)
  g.y.lt.vect.rank = rank(g.y.lt.vect)/(length(g.y.lt.vect)+1)

  g.x.lt.vect.probit = stats::qnorm(g.x.lt.vect.rank)
  g.y.lt.vect.probit = stats::qnorm(g.y.lt.vect.rank)

  # 4. dMI test statistic & p-value
  dist.data = cbind(g.x.lt.vect.probit, g.y.lt.vect.probit)
  n = dim(dist.data)[1]

  ## joint entropy -- marginal entropy is the unchanged by permutation
  H = IndepTest::KLentropy(dist.data, k = max(kbar), weights = F)[[1]][kbar]

  if(MI == TRUE) {
    Hgx = IndepTest::KLentropy(g.x.lt.vect.probit, k = max(kbar), weights = F)[[1]][kbar]
    Hgy = IndepTest::KLentropy(g.y.lt.vect.probit, k = max(kbar), weights = F)[[1]][kbar]

    mutinf.per.k = Hgx + Hgy - H
    mutinf.per.k = ifelse(mutinf.per.k < 0, 0, mutinf.per.k)
    mutinf.kbar = sum(mutinf.per.k)/max(kbar)
  }

  ## permuted entropy
  future::plan(future::multisession, workers = n.cores)

  H.perm.ls = future.apply::future_lapply(
    X = 1:n.perm,
    FUN = function(b) {
      y.perm = as.matrix(g.y.lt.vect.probit)[sample(n), ]
      data.perm = cbind(g.x.lt.vect.probit, y.perm)
      H.perm.b = IndepTest::KLentropy(data.perm, k = max(kbar), weights = F)[[1]][kbar]
      return(H.perm.b)
    },
    future.seed = TRUE
  )
  H.perm = do.call(rbind, H.perm.ls)

  ## k_bar (averaged) test statistic & p-value
  teststat = sum(H)
  nullstats = apply(H.perm, 1, sum)
  p = (1 + sum(nullstats <= teststat))/(n.perm + 1)

  # output
  out = list(
    "p.value.kbar" = p
  )

  if(MI == TRUE) {
    out = list(
      "dMI.kbar" = mutinf.kbar,
      "p.value.kbar" = p
    )
  }

  return(out)
}
