#' Step for the gradient descent
#'
#' @param gamma a d x d matrix : the variogram matrix.
#' @param weights a d x d symmetric matrix with a zero diagonal.
#' @param lambda a positive number : the weight of the penalty.
#' @param size_grid integer : size of the search grid for the optimal step.
#'
#' @returns A function of clusters and R matrix which returns the next step of
#' the optimisation for the gradient descent algorithm.
#'
#' @examples
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1, 3), c(2,4))
#' W <- matrix(c(0, 1, 1, 1,
#'               1, 0, 1, 1,
#'               1, 1, 0, 1,
#'               1, 1, 1, 0), nc = 4)
#' gamma <- matrix(c(0, 2, 1, 0,
#'                   2, 0, 4, 1,
#'                   1, 4, 0, 7,
#'                   0, 1, 7, 0), nc = 4)
#' f <- step_gradient(gamma, W, 0.5)
#' f(R, clusters, 1)
#'
#' @keywords internal
step_gradient <- function(gamma, weights, size_grid = 100) {
  # Initialization of functions
  dlog <- nloglike_grad_np(gamma)                     # Neg-lklh gradient part
  dpen <- penalty_grad(weights)                       # Penalty gradient part

  # Penalised negative log-likelihood
  nllh <- neg_likelihood_pen(gamma, weights)

  function(R, clusters, lambda) {
    # Initialization
    p <- sapply(clusters, length)           # Vector of cluster's size

    # Gradient matrix computation
    grad <- dlog(R, clusters) + lambda * dpen(R, clusters)

    # Grid line search for optimal gradient step
    # Grid line construction
    check_pos <- \(.) !semi_def(sub_theta(R - . * grad, clusters))

    if (max(p) == 1) {
      s_opt <- optim(par = 1, fn = \(.) nllh(R - . * grad, clusters, lambda),
                     method = "Brent", lower = 0, upper = 1)$par

      s_opt <- s_optimal(s_opt, check_pos)

      return(list(step = s_opt, gradient = grad))
    }
    s_max <- min(
      # Maximum step size to get positive matrix
      abs(((R %*% p) / (grad %*% p))[p > 1])
    )

    s_opt <- optim(par = 1, fn = \(.) nllh(R - . * grad, clusters, lambda),
                   method = "Brent", lower = 0, upper = min(s_max, 1))$par

    s_opt <- s_optimal(s_opt, check_pos)
    # Returning results : size step and gradient matrix
    list(step = s_opt, gradient = grad)
  }
}

#' Gradient descent algorithm for Husler-Reiss graphical models clustering
#'
#' @param gamma a d x d matrix : the variogram matrix.
#' @param weights a d x d symmetric matrix with a zero diagonal.
#' @param lambda a positive number : the weight of the penalty.
#' @param ...
#'
#' @returns A function which compute the maximum likelihood estimator using
#' cluster-path gradient descent and returns the estimation of the clusters and
#' the corresponding R matrix.
#'
#' @examples
#' W <- matrix(c(0, 1, 1, 1,
#'               1, 0, 1, 1,
#'               1, 1, 0, 1,
#'               1, 1, 1, 0), nc = 4)
#' gamma <- generate_random_Gamma(d = 4)
#' R <- matrix(c(1,0,0,-1,
#'               0,1,1,-2,
#'               0,1,1,-1,
#'               -1,-2,-1,1), nc = 4)
#' Cluster_HR <- get_cluster(gamma, W, 100)
#' Cluster_HR(R)
#'
#' @export
#'
get_cluster <- function(gamma, weights, eps_f, ...) {
  L <- neg_likelihood_pen(gamma, weights)
  step <- step_gradient(gamma, weights, ...)

  function(R.init, lambda, it_max = 1000, eps_g = 1e-3) {
    # Initialization
    d <- nrow(gamma)

    if (lambda == 0) {
      return(
        list(
          R = R.init,
          clusters = as.list(1:d),
          nllh = L(R, clusters, lambda),
          lambda = lambda
        )
      )
    }

    res <- gradient_descent_rcpp(
      R = R.init,
      clusters = as.list(1:d),
      step = step,
      lambda = lambda,
      it_max = it_max,
      eps_g = eps_g,
      eps_f = eps_f
    )

    R <- res$R
    clusters <- res$clusters

    if (length(R) == 1) {
      return(
        list(
          R = R,
          clusters = clusters,
          nllh = -(d - 1) * (d - 2) * R,
          lambda = lambda
        )
      )
    }
    return(
      list(
        R = R,
        clusters = clusters,
        nllh = L(R, clusters, lambda),
        lambda = lambda
      )
    )
  }
}

#' Optimization results from data
#'
#' @param data \eqn{n \times d} matrix : the data.
#' @param zeta a positive number : tuned parameter for the exponential weights.
#' @param lamda a numerix vector of positive number : the grid line for \eqn{\lambda}.
#' @param eps_f d
#' @param eps_g d
#' @param it_max d
#'
#' @returns Returns the optimization path results from `get_clusters()` for each
#' \eqn{\lambda} from the data.
#'
#' @importFrom graphicalExtremes emp_vario Gamma2Theta
#'
#' @export
HR_Clusterpath <- function(data, zeta, lambda, eps_g = 1e-3, eps_f = 1e-2, it_max = 1000) {
  # Initialization
  Gamma_est <- emp_vario(data)
  d <- ncol(data)
  R.init <- Gamma2Theta(Gamma_est)

  # Exponential weights construction
  D <- distance_matrix(R.init, as.list(1:d))
  W <- matrix(rep(0, d * d), nc = d)

  for (k in 1:(d - 1)){
    for (l in (k + 1):d){
      w_val <- exp(-zeta * D[k, l])
      W[k, l] <- w_val
      W[l, k] <- w_val
    }
  }

  Cluster_HR <- get_cluster(
    gamma = Gamma_est,
    weights = W,
    eps_f = eps_f
  )

  future::plan(future::multisession, workers = parallel::detectCores() - 1)

  furrr::future_map(
    lambda,
    \(.) Cluster_HR(R.init = R.init, lambda = ., it_max = it_max, eps_g = eps_g),
    .options = furrr::furrr_options(seed = TRUE)
  )

}
