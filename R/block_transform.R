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
#' build_Theta(R, clusters)
#'
#'
NULL

#' @rdname theta-r
#'
#' @export
build_Theta <- function(R, clusters) {
  return(
    .build_theta(R, clusters)
  )
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
