#' Functions to navigate between \eqn{R} and \eqn{\Theta} for given clusters
#'
#' The block matrix models are characterized by two objects : a clustering partition
#' \eqn{\{C_1, \dots, C_K\}} of \eqn{V = \{1, \dots, d\}} and a matrix \eqn{R}
#' which belongs to \eqn{\mathcal S_K(\mathbb R)}, the set of symmetric \eqn{K \times K}
#' matrix. The expression of the precision matrix \eqn{\Theta} statisfying the block matrix
#' assumption is
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
#' corresponding to the precision matrix \eqn{\Theta} induced by \eqn{R} and a list of clusters.
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
#' build_Theta(R, clusters)
#'
#'
NULL

#' @rdname theta-r
#'
#' @export
build_Theta <- function(r_matrix, clusters) {
  # Theta matrix built with the Rcpp function (see src/model.cpp)
  # ---- INITIALIZATION ----
  D_VARIABLE <- sum(sapply(clusters, length))

  # ---- OUTPUT ----
  return(
    .build_theta(
      D_VARIABLE = D_VARIABLE,
      R_matrix   = r_matrix,
      clusters   = clusters
    )
  )
}

#' @rdname theta-r
#'
#' @export
extract_R_matrix <- function(theta, clusters) {
  # ---- INITIALIZATION ----
  NB_CLUSTERS <- length(clusters)      # Number of clusters
  r_matrix <- matrix(
    data = rep(NA, NB_CLUSTERS ** 2),  # Reduced matrix R
    nc   = NB_CLUSTERS
  )

  indx <- 1:NB_CLUSTERS     # Vector of indices for the loop

  # ---- COMPUTATION ----
  for (i in indx) {
    # Select a cluster
    k <- clusters[[i]][1]

    # Take the value for the diagonal by taking another variable
    # in the selected cluster.
    if (length(clusters[[i]]) > 1) {
      l <- clusters[[i]][2]
      r_matrix[i, i] <- theta[k, l]
    }

    # Take the value for the other coefficients by selecting the
    # first variable of the other clusters.
    for (j in indx[-i]){
      l <- clusters[[j]][1]
      r_matrix[i, j] <- theta[k, l]
    }
  }

  # ---- OUTPUT ----
  return(
    r_matrix
  )
}
