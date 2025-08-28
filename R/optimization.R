#' Clusterpath algorithm for Hüsler-Reiss models
#'
#' Gradient descent based on the Clusterpath algorithm adapted to the likelihood
#' of a graphical Hüsler-Reiss model. Usefull when the precision matrix \eqn{\Theta}
#' (or the variogram \eqn{\Gamma}) has a block matrix structure.
#'
#' The block matrix models is defined from two elements : the cluster's partition
#' \eqn{\{C_1, \dots, C_K\}}, included in \eqn{V} and the \eqn{R} matrix which belongs
#' to \eqn{\mathcal S_K(\mathbb R)}, the set of symmetric \eqn{K \times K} matrix (See also
#' \code{\link{build_theta}()} and \code{\link{extract_R_matrix}()}).
#'
#' The Clusterpath aims to find optimum of some penalised negative loglikelihood defined in
#' \code{\link{neg_likelihood_pen}()}.
#'
#' @name hr-clusterpath
#'
#' @param lambda For `get_clusters()`, a positive number : the weight of the penalty.
#'
#' For `HR_Clusterpath()`, a numerix vector of positive number : the grid line for \eqn{\lambda}.
#'
#' @param eps_f A positive number : tolerance threshold for merging clusters.
#' @param eps_g A positive number : tolerance threshold for the convergence of the gradient descent.
#' @param it_max An integer : maximal number of iteration of the gradient descent algorithm.
#'
#' @returns `HR_Clusterpath()` is a gradient descent from the data with default values and method for
#' the optimization : the variogram matrix \eqn{\Gamma} is the empirical variogram and the weights are
#' set as the exponential weights, which depends of only one tuning parameter \eqn{\zeta} and where 
#' we have some theoretical results.
#'
#' `get_clusters()` is a freer version of the previous function. You can use your own estimation for the
#' variogram \eqn{\Gamma} and customizable weights with the matrix \eqn{W}.
#'
#' Both of them produce a list of results :
#'  - $R : the \eqn{R} matrix of the clusters.
#'  - $clusters : a list of the variable indices, clustered.
#'  - $nllh : the value of the negative penalised negative loglikelihood.
#'  - $lambda : the value of lambda.
#'
#' In the case of `HR_Clusterpath()` it is a list of the previous list.
#'
#' @section Some results for `HR_Clusterpath()`:
#'
#' When we get replications \eqn{X_1, \dots, X_n} of a random vector \eqn{X} which
#' belongs to the attraction domain of a \eqn{K}-block Hüsler-Reiss graphical model, the
#' minimum of the penalised negative loglikelihood conveges almost surely to the
#' true precision matrix of the model \eqn{\Theta^*}, for all \eqn{\lambda}, provided that
#' we choose a sequence \eqn{(\zeta_n)_{n\in \mathbb N^*}} which grows slower than the
#' \eqn{log(n)} sequence. In general, the \eqn{\zeta \in [1,2]} is a good choice.
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
NULL

#' @rdname hr-clusterpath
#'
#' @param gamma For `get_clusters()`, a \eqn{d \times d} matrix : the variogram matrix \eqn{\Gamma}.
#' @param weights For `get_clusters()`, the \eqn{d \times d} symmetric weightsmatrix with a zero diagonal.
#' @export
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

#' @rdname hr-clusterpath
#'
#' @param data For `HR_Clusterapth()`, a \eqn{n \times d} matrix : the matrix of data.
#' @param zeta For `HR_Clusterapth()`, a positive number : tuned parameter for the exponential weights.
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

# internal ------------------------------------------------------------------------------------

#' Step for the gradient descent
#'
#' @noRd
#' @param gamma a d x d matrix : the variogram matrix.
#' @param weights a d x d symmetric matrix with a zero diagonal.
#' @param lambda a positive number : the weight of the penalty.
#' @param size_grid integer : size of the search grid for the optimal step.
#'
#' @returns A function of clusters and R matrix which returns the next step of
#' the optimisation for the gradient descent algorithm.
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
