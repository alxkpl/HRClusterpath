#' Functions to navigate between \eqn{R} and \eqn{\Theta} for given clusters
#'
#' The block matrix models is defined from two elements : the cluster's partition
#' \eqn{\{C_1, \dots, C_K\}}, included in \eqn{V} and the \eqn{R} matrix which belongs
#' to \eqn{\mathcal S_K(\mathbb R)}, the set of symmetric \eqn{K \times K} matrix.
#'The expression of the precision \eqn{\Theta} which statisfies the block matrix model is 
#' \deqn{
#'    \Theta = U R U^t + A
#' }
#' where \eqn{U} is the cluster matrix and \eqn{A} a diagonal matrix such that the rows of
#' theta sum to zero :
#' \eqn{
#' a_{ii} = - \sum_l p_l r_{kl}
#' }
#' for \eqn{i} in cluster \eqn{C_k}.
#'
#' @name theta-r
#'
#' @param Theta For `extract_R_matrix()`, a \eqn{d \times d} block matrix which can be factorizable.
#' @param R For `build_theta()`, the \eqn{K \times K} matrix of the clusters.
#' @param clusters A list of indices associated to a partition of \eqn{V}.
#'
#' @return For `extract_R_matrix()`, a matrix on size the number of clusters (i.e. \eqn{K})
#' corresponding to the reduced matrix \eqn{R} of the original matrix \eqn{\Theta}.
#'
#' For `build_theta()`, a matrix on size the number of variables (i.e. \eqn{d})
#' corresponding to the precision matrix \eqn{\Theta} induced by \eqn{R}.
#'
#' @examples
#' ##############################################################
#' #                       FROM THETA TO R
#' ##############################################################
#' Theta <- matrix(
#'  c(4.5, .5, .5, .5, -2, -2, -2,
#'    .5, 4.5, .5, .5, -2, -2, -2,
#'    .5, .5, 4.5, .5, -2, -2, -2,
#'    .5, .5, .5, 4.5, -2, -2, -2,
#'    -2, -2, -2, -2, 6, 1, 1,
#'    -2, -2, -2, -2, 1, 6, 1,
#'    -2, -2, -2, -2, 1, 1, 6),
#'  nc = 7
#' )
#'
#' clusters <- list(c(1,2,3,4), c(5,6,7))
#'
#' extract_R_matrix(Theta, clusters)
#'
#' ##############################################################
#' #                       FROM R TO THETA
#' ##############################################################
#'
#' R <- matrix(c(1, -3, 0,
#'               -3, 2, -2,
#'               0, -2, 1), nc = 3)
#'
#' clusters <- list(1:4, 5:8, 9:12)
#'
#' build_theta(R, clusters)
#'
#'
NULL

#' @rdname theta-r
#'
#' @export
build_theta <- function(R, clusters) {
  U <- U_matrix(clusters)
  URUt <- U %*% R %*% t(U)
  a <-  - rowSums(URUt)

  URUt + diag(a)

}

#' @rdname theta-r
#'
#' @export
extract_R_matrix <- function(Theta, clusters) {

  K <- length(clusters)
  indx <- 1:K
  R <- matrix(rep(NA, K * K), nc = K)

  for (i in indx) {
    k <- clusters[[i]][1]
    if (length(clusters[[i]]) > 1) {
      l <- clusters[[i]][2]
      R[i, i] <- Theta[k, l]
    }
    for (j in indx[-i]){
      l <- clusters[[j]][1]
      R[i, j] <- Theta[k, l]
    }
  }

  R
}

#' ChiToGamma
#'
#' Transform a \eqn{\chi} matrix to the corresponding variogram \eqn{\Gamma}
#'
#' @param Chi_matrix The matrix with \eqn{\chi_{ij}} entries.
#'
#' @return Gives the variogram \eqn{\Gamma} according to the \eqn{\chi} matrix for
#' Hüsler-Reiss MGDP. In a such case, there exists a closed equation which
#' link the variogram and the extremal coefficients, given by :
#' \deqn{
#'    \chi_{ij} = 2 - 2 \phi(\sqrt{\Gamma_{ij}/2})
#' }
#' where \eqn{\phi} is the standard normal distribution function.
#'
#' @examples
#' ChiToGamma(matrix(c(1, 0.7,
#'                     0.7, 1),
#'            nrow = 2))
#'
#' @export
ChiToGamma <- function(Chi_matrix) {

  (2 * qnorm((2 - Chi_matrix) / 2))**2

}

# internal --------------------------------------------------------------------------

#' Semi-definite checker
#'
#' @param M A /eqn{n/times n} symmetric matrix.
#'
#' @returns `TRUE` if the matrix is semi-definite positive, `FALSE` otherwise. The
#' verification is done by checking the sign of the eigen-values of the matrix.
#'
#' @keywords internal
#'
semi_def <- function(M) {
  !sum(Re(eigen(M)$values) < -1e-10)        # for numerical errors
}

#' Computation of the matrix of clusters
#'
#' @param clusters a list of vector : each vector gives the element of a
#' cluster.
#'
#' @returns The d x K matrix U such that u_jk = 1 if the variable j belongs to
#' cluster C_k and 0 otherwise.
#'
#' @keywords internal
#' @examples
#' clusters <- list(c(1,2,3),c(4,5))
#' U_matrix(clusters)
#'
U_matrix <- function(clusters) {
  K <- length(clusters)
  d <- length(unlist(clusters))
  U <- matrix(rep(0, d * K), nc = K)
  for (k in 1:K){
    for (j in 1:d){
      if (j %in% clusters[[k]]) {
        U[j, k] <- 1
      }
    }
  }
  U
}



#' Variogram transformation application gamma
#'
#' @param sigma A d x d numeric matrix.
#'
#' @returns For a symmetric positive matrix sigma (covariance matrix),
#' return the corresponding variogram matrix. Can be used for other
#' but with no interpretation.
#'
#' @keywords internal
#' @examples
#' s_sigma <- matrix(rnorm(16, 2), nc = 4)
#' gamma_function(s_sigma %*% t(s_sigma))
gamma_function <- function(sigma) {
  indic <- rep(1, nrow(sigma))

  tcrossprod(diag(sigma), indic) + tcrossprod(indic, diag(sigma)) - 2 * sigma
}

#' Moore-Penrose pseudo inverse
#'
#' @param A a d x d symmetric positive semi-definite matrix.
#'
#' @returns Computes the Moore-Penrose inverse of a matrix. The calculation is
#' done thanks to an article and if  :
#'                                A = L L^t       (e.g. Crout decomposition)
#' (with L having no zero-columns) then we have :
#'                        A^+ = L (L^t L)^-1 (L^t L)^-1 L^t
#'
#' @keywords internal
#' @examples
#' A <- matrix(c(1, 2, 3,
#'               2, 5, 6,
#'               3, 6, 9), nc = 3)
#' psolve(A)
psolve <- function(A, tol = 1e-12) {

  A <- unname(as.matrix(A))

  n <- nrow(A)
  if (n != ncol(A)) {
    stop("no square matrix.")
  }

  S <- crout_decomposition_rcpp(A, tol = tol)

  L <- S[, which(diag(S) != 0)]         # to get no null columns

  psolve_rcpp(L)

}

#' Computation of the first weight matrix
#' @param data the data
#'
#' @returns The initial wieght matrix without fused variables.
#'
#' @keywords internal
#'
compute_W <- function(data) {
  Gamma_est <- graphicalExtremes::emp_vario(data)
  d <- ncol(data)
  R.init <- graphicalExtremes::Gamma2Theta(Gamma_est)
  # Exponential weights construction
  D <- distance_matrix(R.init, as.list(1:d))
  W <- matrix(rep(0, d * d), nc = d)

  for (k in 1:(d - 1)) {
    for (l in (k + 1):d){
      W[k, l] <- exp(-1 * D[k, l])
    }
  }
  W + t(W)
}

#' Computation of the clustered weight matrix
#'
#' @param weights a d x d symmetric matrix with a zero diagonal.
#'
#' @returns A function of clusters : the weight matrix for each pair
#' of clusters :
#'                          W_kl = sum_(i in C_k) sum_(j in C_l) w_ij
#'
#' @keywords internal
#' @examples
#' clusters <- list(c(1,2,3), c(4,5))
#' W <- matrix(1:25, nc = 5)
#' W_c <- weight_clustered(W)
#' W_c(clusters)
weight_clustered <- function(weights) {
  function(clusters) {
    U <- U_matrix(clusters)

    t(U) %*% weights %*% U
  }
}


#' Generalised determinant
#'
#' @param A a d x d real valued matrix.
#' @param tol a positive value.
#'
#' @returns Compute the generalised determinant of a matrix A. We recall
#' that the generalised determinant is an extension of the determinant for
#' singular matrix. It corresponds to the product of all the non zero eigen
#' values.
#'
#' @keywords internal
#' @examples
#' A <- matrix(c(1,2,3,
#'               2,5,6,
#'               3,6,9), nc = 3)
#' gen_det(A)
gen_det <- function(A, tol = 1e-10) {
  res <- Re(eigen(A, only.values = TRUE)$values)      # avoid complex number (impossible for symmetric matrices)
  return(
    prod(res[res > tol])
  )
}


sub_theta <- function(R, clusters) {
  p <- sapply(clusters, length)           # Vector of cluster's size
  tilde_R <- - as.numeric(R %*% p)

  # Subtheta matrix :
  R %*% diag(p) + diag(tilde_R)

}