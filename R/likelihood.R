#' Negative log-likelihood computation.
#'
#' @param gamma A \eqn{d \times d} matrix: the empirical variogram matrix.
#'
#' @returns For a fixed variogram gamma, compute for a set of clusters and
#' corresponding \eqn{R} matrix, the value of the associated negative likelihood
#' defined by :
#'
#' \deqn{
#'      L(\Theta) = - \log(|\Theta|_+) - \frac 1 2 \text{tr}(\hat \Gamma \Theta)
#' }
#'
#' where \eqn{|\cdot|_+} is the generalised determinant.
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
#' @param weights a \eqn{d \times d} symmetric matrix with a zero diagonal.
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
neg_likelihood_pen <- function(gamma, weights) {
  nllh <- neg_likelihood(gamma)
  pen <- penalty(weights)

  function(R, clusters, lambda) {

    nllh(R, clusters) + lambda * pen(R, clusters)

  }
}

#' Gradient of the negative log likelihood without penalty
#'
#' @param gamma a d x d variogram matrix.
#'
#' @returns A function of the R matrix and clusters and compute the gradient
#' matrix of the negative log likelihood for a fixed variogram gamma. The
#' gradient matrix can be computed by :
#'
#'                        dnllh = dlog + dtrace
#' where :
#'      - dlog = t(U) g(Theta_p) U - 0.5 diag(t(U) g(Theta_p) U)
#'      - dtrace = - t(U) gamma U  + 0.5 diag(t(U) gamma U)
#'
#' @examples
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4))
#' gamma <- matrix(c(0,2,1,0,
#'                   2,0,4,1,
#'                   1,4,0,7,
#'                   0,1,7,0), nc = 4)
#' gradient <- nloglike_grad_np(gamma)
#' gradient(R, clusters)
#'
#' @keywords internal
nloglike_grad_np <- function(gamma) {
  function(R, clusters) {
    # Get tu U matrix of clusters indicators
    U <- U_matrix(clusters)
    p <- colSums(U)

    # Computation of gamma(Theta_p) with Theta_p the Penrose inverse of Theta
    G_theta_p <- build_theta(R, clusters) |>
      psolve() |>
      gamma_function()

    # Gradient (factorisation of the formula)
    dlog <- t(U) %*% (G_theta_p - gamma) %*% U
    diag <- - .5 * diag(diag(dlog * (p == 1)))

    dlog + diag
  }
}

#' Gradient matrix of distance between two columns
#'
#' @param R K x K symmetric matrix.
#' @param clusters a list of vector : each vector gives the element of
#' a cluster.
#'
#' @returns Return a functon of indices (k', l') computing the gradient
#' matrix of tilde D^2(r_k', r_l'). See section 4.3.3 in cluster document
#' for details.
#'
#' @examples
#' R <- matrix(c(-1,0,-2,
#'               0,-3,-1,
#'               -2,-1,-1), 3)
#' clusters <- list(c(1,3), c(2), 4)
#' grad <- gradient_D2(R, clusters)
#' grad(1, 2)
#' @keywords internal
gradient_D2 <- function(R, clusters) {
  # Initialization
  K <- length(clusters)                   # Number of clusters
  p <- sapply(clusters, length)           # Vector of cluster's size
  function(k, l) {
    A <- matrix(rep(0, K * K), nc = K)

    # Computation of the non zero row and column
    A[k, ] <- 2 * p * (R[k, ] - R[l, ])
    A[, l] <- 2 * p * (R[l, ] - R[k, ])
    A[k, l] <- 2 * ((p[k] - 1) * (R[k, l] - R[k, k]) +
                      (p[l] - 1) * (R[k, l] - R[l, l]))

    # Build the symmetry of the gradient (except for the diagonal)
    grad <- A + t(A)

    # Computation of the diagonal
    grad[k, k] <- 2 * (p[k] - 1) * (R[k, k] - R[k, l])
    grad[l, l] <- 2 * (p[l] - 1) * (R[l, l] - R[k, l])

    grad

  }
}

#' Computation of the penalty's gradient
#'
#' @param weights a d x d symmetric matrix with a zero diagonal.
#'
#' @returns A function of clusters and corresponding R matrix. Compute
#' the gradient with fixed weight. The expression of the gradient is
#' just the weighted sum of the gradient of each tilde D^2 where the
#' weights are the clustered weights.
#' See equations in section 4.3.3 for details.
#'
#' @examples
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4))
#' W <- matrix(c(0, 1, 1, 1,
#'               1, 0, 1, 1,
#'               1, 1, 0, 1,
#'               1, 1, 1, 0), nc = 4)
#' dpen <- penalty_grad(W)
#' dpen(R, clusters)
#' @keywords internal
penalty_grad <- function(weights) {
  get_W <- weight_clustered(weights)
  function(R, clusters) {
    # Initialization
    W <- get_W(clusters)                # Weight clustered
    K <- length(clusters)               # Number of clusters

    # Function of gradient of indices
    grad_D2 <- gradient_D2(R, clusters)

    res <- matrix(rep(0, K * K), nc = K)

    for (k in 1:(K - 1)) {
      for (l in (k + 1):K) {
        res <- res + W[k, l] * grad_D2(k, l)        # Weighted sum
      }
    }

    res
  }
}
