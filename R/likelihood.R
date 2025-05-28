#' Negative log-likelihood computation.
#'
#' @param gamma A d x d matrix: the empirical variogram matrix.
#'
#' @returns For a fixed variogram gamma, compute for a set of clusters and
#' corresponding R matrix, the value of the associated negative likelihood
#' defined by :
#'                  nllh = - log(|theta|_+) - 1/2 trace(gamma * theta)
#' where |.|_+ is the generalised determinant.
#'
#' @examples
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4))
#' gamma <- matrix(c(0, 2, 1, 0,
#'                   2, 0, 4, 1,
#'                   1, 4, 0, 7,
#'                   0, 1, 7, 0), nc = 4)
#' nllh <- neg_likelihood(gamma)
#' nllh(R, clusters)
#'
#' @export
neg_likelihood <- function(gamma) {
  function(R, clusters) {
    # Building of the theta matrix from R
    theta <- build_theta(R, clusters)

    # Computation of the log-determinant part
    log_det <- log(gen_det(theta))

    # Computation of the trace part
    tr <- sum(diag(gamma %*% theta))


    - log_det - .5 * tr
  }
}

#' Penalty function.
#'
#' @param weights a d x d symmetric matrix with a zero diagonal.
#'
#' @returns For fixed weights, returns a function which compute the value of
#' the penalty for chosen clusters and corresponding R matrix. We recall the
#' penalty is given by :
#'                    P(R) = sum_{k<l} W_kl D^2(r_.k, r_.l)
#' with W_kl the clustered weight for clusters k and l.
#'
#' @examples
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4))
#' W <- matrix(c(0, 1, 1, 1,
#'               1, 0, 1, 1,
#'               1, 1, 0, 1,
#'               1, 1, 1, 0), nc = 4)
#' P <- penalty(W)
#' P(R, clusters)
penalty <- function(weights) {
  # Fixing the weights for computing weights clustered
  get_W <- weight_clustered(weights)
  function(R, clusters) {
    # Initialization
    K <- length(clusters)              # Number of clusters
    D2 <- D_tilde2_r(R, clusters)      # Function for distance between clusters
    W <- get_W(clusters)               # Weights clustered
    D <- matrix(rep(0, K * K), nc = K) # Distance matrix for clusters

    # Computation of the distance matrix
    for (l in 2:K) {
      for (k in 1:(l - 1)) {        # we keep only the lower triangular part
        D[k, l] <- D2(k, l)
      }
    }

    sum(D * W)
  }
}

#' Computation of the penalised negative log-likelihood
#'
#' @param gamma a d x d matrix : the variogram matrix.
#' @param weights a d x d symmetric matrix with a zero diagonal.
#' @param lambda a positive number : the weight of the penalty.
#'
#' @returns A function of clusters and the R matrix which compute the penalised
#' negative log-likelihood of the model.
#'
#' @examples
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4))
#' W <- matrix(c(0, 1, 1, 1,
#'               1, 0, 1, 1,
#'               1, 1, 0, 1,
#'               1, 1, 1, 0), nc = 4)
#' gamma <- matrix(c(0,2,1,0,
#'                   2,0,4,1,
#'                   1,4,0,7,
#'                   0,1,7,0), nc = 4)
#' f <- neg_likelihood_pen(gamma, W, 0.5)
#' f(R, clusters)
#'
#' @export
neg_likelihood_pen <- function(gamma, weights, lambda) {
  nllh <- neg_likelihood(gamma)
  pen <- penalty(weights)

  function(R, clusters) {

    nllh(R, clusters) + lambda * pen(R, clusters)

  }
}