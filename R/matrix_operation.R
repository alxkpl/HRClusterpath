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

#' Cluster distance squared function
#'
#' @param R K x K symmetric matrix.
#' @param clusters a list of vector : each vector gives the element of
#' a cluster.
#'
#' @returns A function of cluster number : compute the square distance between
#' two clusters for the distance defined in section 4.2 in cluster document.
#'
#' @keywords internal
#' @examples
#'
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4))
#' D2 <- D_tilde2_r(R, clusters)
#' D2(1, 2)
D_tilde2_r <- function(R, clusters) {
  function(k, l) {
    if (k == l) {
      # 0 distance for the same cluster even with the symmetry of the matrix
      0
    }else {
      # Parameters of clusters
      K <- length(clusters)               # Number of clusters
      p <- sapply(clusters, length)       # Vector of cluster's size

      # for fixed k, l the square difference is multiplied by 1-p_k and 1-p_l
      sum((p - ((1:K) %in% c(k, l))) * (R[k, ] - R[l, ])**2)
    }
  }
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

#' From R matrix compute theta matrix
#'
#' @param R a K x K matrix : the matrix of the clusters coefficients
#' @param clusters a list of vector : each vector gives the element of a
#' cluster.
#'
#' @returns Return the theta matrix from the value of the R matrix of clusters
#' coefficient using :
#'
#'                            theta = U R U^t + A
#'
#' where U is the cluster matrix and a diagonal matrix such that the rows of
#' theta sum to zero :
#'                            a_ii = - sum_l p_l r_kl
#' for i in cluster C_k.
#'
#' @examples
#'
#' R <- matrix(c(1, 2,
#'               2, 5), nr = 2)
#' clusters <- list(c(1,2,3),c(4,5))
#' build_theta(R, clusters)
#'
#' @export
build_theta <- function(R, clusters) {
  U <- U_matrix(clusters)
  URUt <- U %*% R %*% t(U)
  a <-  - rowSums(URUt)

  URUt + diag(a)

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


#' Crout factorization algorithm.
#'
#' @param A a d x d symmetric positive matrix.
#' @param tol a positive value : tolerance for the zero diagonal
#'
#' @returns The LU Crout decomposition of a matrix A. For A a symmetric
#' positive matrix, there exists a LU decomposition such that :
#'                                   A = L' L
#'
#' @keywords internal
#' @examples
#'A <- matrix(c(1,2,3,
#'              2,5,6,
#'              3,6,9), nc = 3)
#' L <- crout_factorisation(A)
#' L %*% t(L)
crout_factorisation <- function(A, tol = 1e-12) {
  A <- unname(as.matrix(A))
  n <- nrow(A)
  if (n != ncol(A)) {
    stop("no square matrix.")
  }
  if (!semi_def(A)) {
    stop("no positive semi definite matrix.")
  }
  # d <- rep(0, n)
  # L <- diag(rep(1, n))
  # d[1] <- A[1, 1]
  # for (i in 2:n) {
  #   for (j in 1:(i - 1)) {
  #     L[i, j] <- (1 / d[j]) * (A[i, j] - sum(L[i, ] * L[j, ] * d))
  #   }
  #   d[i] <- A[i, i] - sum(d * (L[i, ])**2)
  # }

  # return(L %*% diag(sqrt(d * (d > tol))))
  crout_decomposition_rcpp(A, tol)
}

#' Moore-Penrose pseudo inverse
#'
#' @param A a d x d symmetric positive semi-definite matrix.
#'
#' @returns Computes the Moore-Penrose inverse of a matrix. The calculation is
#' done thanks to an article and if  :
#'                                A = L L^t
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

  S <- crout_factorisation(A, tol = tol)

  L <- S[, which(diag(S) != 0)]         # to get no null columns

  return(
    L %*% solve(t(L) %*% L) %*% solve(t(L) %*% L) %*% t(L)
  )
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
  D <- D_tilde2_r(R.init, as.list(1:d))
  W <- matrix(rep(0, d * d), nc = d)

  for (k in 1:(d - 1)) {
    for (l in (k + 1):d){
      W[k, l] <- exp(-1 * D(k, l))
    }
  }
  W + t(W)
}

#' Transform a chi matrix to the corresponding variogram
#'
#' @param Chi_matrix the matrix with the chi coefficient.
#'
#' @return
#' The corresponding variogram matrix \eqn{\Gamma}.
#'
#' @examples
#' ChiToGamma(matrix(c(0.5, 0.7, 0.7, 0.5), nrow = 2))
#'
#' @export
ChiToGamma <- function(Chi_matrix) {

  (2 * qnorm((2 - Chi_matrix) / 2))**2

}



#' Function to extract the coefficients of the reduced matrix R.
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
#' @param Theta A block matrix which can be factorizable.
#' @param clusters a list of vectors : the variable index per clusters.
#'
#' @return A matrix on size the number of clusters (i.e. \eqn{K}). It
#' corresponds to the reduced matrix \eqn{R} of the original matrix.
#'
#' @examples
#' # We can build a Theta which is a block matrix as described upper :
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