#' Clusterpath algorithm for Hüsler-Reiss models hierarchical version
#'
#' Block gradient descent based on the Clusterpath algorithm for Hüsler-Reiss
#' graphical model. Usefull when the precision matrix \eqn{\Theta}
#' (or the variogram \eqn{\Gamma}) has a block matrix structure.
#'
#' The block matrix models are characterized by two objects : a clustering partition
#' \eqn{\{C_1, \dots, C_K\}} of \eqn{V = \{1, \dots, d\}} and a matrix \eqn{R}
#' which belongs to \eqn{\mathcal S_K(\mathbb R)}, the set of symmetric \eqn{K \times K}
#' matrix (See also \code{\link{build_theta}()} and \code{\link{extract_R_matrix}()}).
#'
#' The Clusterpath aims to find the optimum of some penalised negative loglikelihood defined in
#' \code{\link{neg_likelihood_pen}()}. The path of solution
#' is computed using the Clusterpath procedure on a mesh for \eqn{\lambda} where each
#' initialization is provided using the results of the previous \eqn{\lambda} within the grid.
#'
#' @name hrc-hierarchy
#'
#' @param Gamma A \eqn{d \times d} matrix : a variogram matrix.
#'
#' @param lambda A numeric vector of positive number : the grid line for \eqn{\lambda}.
#'
#' @param W A \eqn{d \times d} matrix of positive number : the weights for the penalty.
#' If `NULL` (default), it is set as the exponential weights, which depends of only one
#' tuning parameter \eqn{\zeta} and where we have some theoretical results.
#'
#' @param zeta A positive number : tuned parameter for the exponential weights.
#' If `W` is not `NULL`, this parameter is not used.
#'
#' @param eps_f A positive number : tolerance threshold for merging clusters. If `NULL` (default),
#' it is set as \eqn{\kappa} times the median of the non-diagonal elements of the distance matrix
#' between the column of the initial \eqn{\Theta}.
#'
#' @param kappa A positive number : tuned parameter for the fusion threshold. If `eps_f` is not `NULL`,
#' this parameter is not used.
#'
#' @param eps_conv A positive number : tolerance threshold for the convergence of the algorithm.
#'
#' @param tol_opt A positive number : tolerance for the optimal step in the gradient descent.
#'
#' @param iter_max An integer : maximal number of iterations of the procedure.
#'
#' @returns `HR_Clusterpath()` is a block gradient descent from the data with default values and method for
#' the optimization : the variogram matrix \eqn{\Gamma} is the empirical variogram and the weights are
#' set as the exponential weights, which depends of only one tuning parameter \eqn{\zeta} and where
#' some theoretical guarantees are provided.
#'
#' The output is an object of class `HRC_hierarchy` which contains a list of several datas :
#'  - $results : a list of results for each value of lambda, where each result is a
#'               list with the following elements : $R the \eqn{R} matrix of the clusters, $clusters
#'               a list of the variable indices separated per cluster, $likelihood the value of the
#'               penalised negative loglikelihood and $lambda the corresponding regularized parameter.
#'  - $Gamma : the variogram matrix used for the procedure.
#'  - $inputs : a list of the input parameters used in the procedure.
#'
#' @examples
#' # Initialization of the parameters for the simulation
#' R <- matrix(c(1, -3, 0,
#'               -3, 2, -2,
#'               0, -2, 1), nc = 3)
#' clusters <- list(1:5, 6:10, 11:15)
#'
#' # Building of the variogram matrix for the simulation
#' Theta <- build_Theta(R, clusters)
#' Gamma <- graphicalExtremes::Theta2Gamma(Theta)
#'
#' # Simulation of data
#' set.seed(123)
#' n <- 1e4                               # Sample size
#' d <- ncol(Gamma)                       # Number of variables
#' data <- graphicalExtremes::rmstable(
#'            n = n,
#'            model = "HR",
#'            d = d,
#'            par = Gamma)
#'
#' # Computation of the multivariate Pareto data
#' norm_inf <- apply(abs(data), 1, max)
#' quantile <- sort(norm_inf, decreasing = TRUE) |> as.numeric()
#' k <- floor(0.1 * nrow(data))
#' u <- quantile[k]
#'
#' data_par <- data[norm_inf > u, ] / u         # Multivariate Pareto data
#'
#' # Computation of the Hüsler-Reiss Clusterpath
#' zeta <- log(n) ** 2              # Zeta parameter for the exponential weights
#' Gamma_est <- graphicalExtremes::emp_vario(data_par) # Estimated variogram
#' lambda <- seq(0, 1, 0.1)                            # Meshgrid for lambda
#'
#' HRC <- HR_Clusterpath_hierarchy(Gamma_est, lambda = lambda, zeta = zeta)
#' HRC$results$l10
NULL

#' @rdname hrc-hierarchy
#'
#' @importFrom graphicalExtremes emp_vario Gamma2Theta
#'
#' @export
HR_Clusterpath_hierarchy <- function(
  Gamma, lambda, mu = 0, W = NULL, Z = NULL, zeta, eps_lasso = 5e-3, eps_f = NULL, kappa = 1e-2,
  eps_conv = 1e-7, tol_opt = 1e-3, iter_max = 1000
) {
  # Number of variables
  d <- ncol(Gamma)

  # Computation of the initial Theta
  res <- list()
  res$R <- Gamma2Theta(Gamma)

  # Initial clusters
  res$clusters <- as.list(1:d)

  # If no custom weights are given, we use the exponential weights
  if (is.null(W)) {
    W <- exp(-zeta * sqrt(distance_matrix(Gamma, as.list(1:d))))    # Choosen weights
  }

  # Adaptative threshold for the fusion step if no custom one is given
  if (is.null(eps_f)) {
    eps_f <-  kappa * median(sqrt(distance_matrix(res$R, as.list(1:d))) + diag(rep(NA, d)), na.rm = TRUE)
  }

  if (is.null(Z)) {
    Z <- abs(1 / res$R)
    diag(Z) <- 0
    Z <- Z / sum(Z) * d * (d - 1)
  }

  # Non singular matrix projection P for the likelihood computation
  P <- .non_singular_P(d)
  hierarchy <- list()
  hierarchy$results <- list()
  for (i in seq_along(lambda)) {
    res <- .HRClusterpath(
      res$R,
      res$clusters,
      Gamma,
      W,
      Z,
      lambda[i],
      mu,
      eps_lasso,
      eps_f,
      eps_conv,
      tol_opt,
      iter_max
    )

    hierarchy$results[[i]] <- res

    names(hierarchy$results[[i]]$clusters) <- paste0("C", seq_along(hierarchy$results[[i]]$clusters))

    # Likelihood value
    hierarchy$results[[i]]$likelihood <- .Likelihood_penalised(res$R, res$clusters, Gamma, P, W, Z, lambda[i], mu, eps_lasso)
    hierarchy$results[[i]]$lambda <- lambda[i]
    hierarchy$results[[i]]$mu <- mu
  }

  hierarchy$Gamma <- Gamma

  # Input parameters
  hierarchy$inputs <- list()
  hierarchy$inputs$eps_conv <- eps_conv
  hierarchy$inputs$eps_f <- eps_f
  hierarchy$inputs$tol_opt <- tol_opt
  hierarchy$inputs$iter_max <- iter_max

  # Naming the results with the corresponding lambda values
  names(hierarchy$results) <- paste0("l", seq_along(lambda))

  class(hierarchy) <- "HRC_hierarchy"

  return(hierarchy)

}