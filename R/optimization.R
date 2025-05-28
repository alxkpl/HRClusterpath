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
#' f(R, clusters)
#'
#' @keywords internal
step_gradient <- function(gamma, weights, lambda, size_grid = 100) {
  # Initialization of functions
  dlog <- nloglike_grad_np(gamma)                     # Neg-lklh gradient part
  dpen <- penalty_grad(weights)                       # Penalty gradient part

  # Penalised negative log-likelihood
  nllh <- neg_likelihood_pen(gamma, weights, lambda)

  function(R, clusters) {
    # Initialization
    p <- sapply(clusters, length)           # Vector of cluster's size

    # Gradient matrix computation
    grad <- dlog(R, clusters) + lambda * dpen(R, clusters)

    # Grid line search for optimal gradient step
    # Grid line construction
    if (max(p) == 1) {
      s_opt <- optim(par = 1, fn = \(.) nllh(R - . * grad, clusters),
                     method = "Brent", lower = 0, upper = 1)$par
      while (!semi_def(sub_theta(R - s_opt * grad, clusters))) {
        s_opt <- 0.95 * s_opt
      }
      return(list(step = s_opt, gradient = grad))
    }
    s_max <- min(
      # Maximum step size to get positive matrix
      abs(((R %*% p) / (grad %*% p))[p > 1])
    )

    s_opt <- optim(par = 1, fn = \(.) nllh(R - . * grad, clusters),
                   method = "Brent", lower = 0, upper = min(s_max, 1))$par
    while (!semi_def(sub_theta(R - s_opt * grad, clusters))) {
      s_opt <- 0.95 * s_opt
    }

    # Returning results : size step and gradient matrix
    list(step = s_opt, gradient = grad)
  }
}

#' Function which merges clusters
#'
#' @param R K x K symmetric matrix.
#' @param clusters a list of vector : each vector gives the element of
#' a cluster.
#' @param eps positive value : minimal tolerance for merging clusters
#' @param cost a function : Cost function of the optimisation
#'
#' @returns Returns, if merging, a list of the new clusters and the
#' corresponding R matrix, where the coefficient of the new clustered
#' is computing by averaging the coefficient of the two previous clusters.
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
#' cost <- neg_likelihood_pen(gamma, weights, 100000)
#' merge_clusters(R, clusters, cost = cost)
#'
#' @keywords internal
merge_clusters <- function(R, clusters, eps = 1e-1, cost) {
  # Initialization
  D <- D_tilde2_r(R, clusters)               # Function of clusters distance
  K <- length(clusters)                      # Actual number of clusters

  # Computation of the distance matrix
  distance <- matrix(rep(Inf, K * K), nc = K)

  for (k in 1:(K - 1)) {
    for (l in (k + 1):K) {
      distance[k, l] <- D(k, l)
    }
  }

  # Search of the two potential clusters to merge
  index <- as.numeric(which(distance == min(distance), arr.ind = TRUE))
  k <- index[1]
  l <- index[2]

  # Checking uselessness of merging
  if (distance[k, l] > eps) {
    return(
      list(
        R = R,
        clusters = clusters
      )
    )
  }

  p <- sapply(clusters, length)           # Vector of cluster's size

  # Case when merging give only one cluster
  if (nrow(R) == 2) {
    return(
      list(
        R = (p[1] * R[1, 1] + p[2] * R[1, 2]) / (p[1] + p[2]),
        clusters = c(clusters[[1]], clusters[[2]])
      )
    )
  }

  # Merging and computation of merged coefficient in the case where K>2
  new_clusters <- clusters[-l]

  new_clusters[[k]] <- c(clusters[[k]], clusters[[l]])      # New cluster

  # Coefficient calculation
  R_new <- R[-l, -l]

  R_new[k, -k] <- ((p[k] * R[k, ] + p[l] * R[l, ]) / (p[k] + p[l]))[-c(k, l)]
  R_new[-k, k] <- ((p[k] * R[k, ] + p[l] * R[l, ]) / (p[k] + p[l]))[-c(k, l)]
  R_new[k, k] <- R[k, l]

  return(
    list(
      R = R_new,
      clusters = new_clusters
    )
  )
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
#' @keywords internal
get_cluster <- function(gamma, weights, lambda, ...) {
  L <- neg_likelihood_pen(gamma, weights, lambda)
  step <- step_gradient(gamma, weights, lambda, ...)
  function(R.init, it_max = 1000, eps_g = 1e-3) {
    # Initialization
    d <- nrow(gamma)
    R <- R.init
    clusters <- as.list(1:d)
    gradstep <- list(gradient = eps_g + 1)
    cpt <- 1
    while ((cpt < it_max) && (length(R) != 1) && (sum(gradstep$gradient**2) > eps_g)) {
      # Gradient step
      gradstep <- step(R, clusters)

      R <- R - gradstep$step * gradstep$gradient
      # Try for merging
      res.merge <- merge_clusters(R, clusters, cost = L, ...)

      if (length(res.merge$R) != length(R)) {
        R <- res.merge$R
        clusters <- res.merge$clusters
      }
      cpt <- cpt + 1
    }

    if (length(R) == 1) {
      return(
        list(
          R = R,
          clusters = clusters,
          nllh = -(d - 1) * (d - 2) * R
        )
      )
    }
    return(
      list(
        R = R,
        clusters = clusters,
        nllh = L(R, clusters)
      )
    )
  }
}

#' Search of optimal lambda for the penalty in optimization
#'
#' @param data n x d matrix : the data.
#' @param chi a positive number : tuned parameter for the exponential weights
#' @param l_grid a numerix vector of psoitive number : the grid line for lambda.
#' @param include_zero Boolean : if FALSE (default) avoid computation for
#' non-penalization setting (i.e. lambda = 0).
#'
#' @returns Returns the optimal optimization from get_clusters() with the best
#'  lambda in the grid line.
#'
#' @importFrom graphicalExtremes emp_vario Gamma2Theta
#'
#' @examples
#' @export
HR_Clusterpath <- function(data, chi, lambda, it_max = 1000) {
  # Initialization
  Gamma_est <- emp_vario(data)
  d <- ncol(data)
  R.init <- Gamma2Theta(Gamma_est)

  # Exponential weights construction
  D <- D_tilde2_r(R.init, as.list(1:d))
  W <- matrix(rep(0, d * d), nc = d)

  for (k in 1:(d - 1)){
    for (l in (k + 1):d){
      W[k, l] <- exp(-chi * D(k, l))
    }
  }
  W <- W + t(W)

  if (lambda == 0) {
    L <- neg_likelihood(Gamma_est)
    return(
      list(
        R = R.init,
        clusters = as.list(1:d),
        nllh = L(R.init, as.list(1:d)),
        lambda = 0
      )
    )
  }

  Cluster_HR <- get_cluster(gamma = Gamma_est, weights = W, lambda = lambda)

  res_base <- Cluster_HR(R.init, it_max = it_max)

  return(
    list(
      R = res_base$R,
      clusters = res_base$clusters,
      nllh = res_base$nllh,
      lambda = lambda
    )
  )

}