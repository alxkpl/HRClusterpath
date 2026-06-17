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
#' @param mu A positive number, the regularized parameter for the lasso penalty.
#'
#' @param W_cluster A \eqn{d \times d} matrix of positive number : the weights for the penalty.
#' If `NULL` (default), it is set as the exponential weights, which depends of only one
#' tuning parameter \eqn{\zeta} and where we have some theoretical results.
#'
#' @param W_lasso A \eqn{d \times d} matrix of positive number : the sparsity weights for the
#' lasso penalty. If `NULL` (default), it is set as the inverse of the absolute coefficients 
#' of the initial guess for the precision matrix.
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
#' @param EPS_CONV A positive number : tolerance threshold for the convergence of the algorithm.
#'
#' @param TOL_OPT A positive number : tolerance for the optimal step in the gradient descent.
#'
#' @param MAX_ITER An integer : maximal number of iterations of the procedure.
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
  Gamma, lambda, mu = 0, W_cluster = NULL, W_lasso = NULL,
  zeta, eps_lasso = 5e-3, eps_f = NULL, kappa = 1e-2,
  EPS_CONV = 1e-7, TOL_OPT = 1e-3, MAX_ITER = 1000
) {
  # ---- INITIALIZATION ----
  D_VARIABLE <- ncol(Gamma)     # Number of variables

  # Intermediate result
  inter_result <- list()                                  # Initialisation
  inter_result$R <- graphicalExtremes::Gamma2Theta(Gamma) # First guess for the R matrix
  inter_result$clusters <- as.list(1:D_VARIABLE)          # Initial clusters

  # Default clusterpath weights : exponential weights with parameter zeta
  if (is.null(W_cluster)) {
    W_cluster <- exp(
      - zeta * sqrt(.distance_matrix(Gamma, as.list(1:D_VARIABLE)))
    )
  }

  # Default merge threshold : data-driven threshold
  if (is.null(eps_f)) {
    # Base on data : median computed from the first guess for the precision matrix
    distance_median <- median(
      x     = sqrt(.distance_matrix(inter_result$R, as.list(1:D_VARIABLE))) + diag(rep(NA, D_VARIABLE)),
      na.rm = TRUE
    )
    eps_f <-  kappa * distance_median
  }

  # Default sparsity weights : inverse absolute coefficient of the initial guess
  if (is.null(W_lasso)) {
    W_lasso <- abs(1 / inter_result$R)
    diag(W_lasso) <- 0      # Weight matrix = null diagonal
    W_lasso <- W_lasso / sum(W_lasso) * D_VARIABLE * (D_VARIABLE - 1)   # Standardized weights
  }

  # Matrix projection P for the likelihood computation
  P <- .non_singular_P(D_VARIABLE)
  hierarchy <- list()
  hierarchy$results <- list()

  # ---- COMPUTATION ----
  for (i in seq_along(lambda)) {
    # Use the previous results as initialisation for the next lambda
    inter_result <- .HRClusterpath(
      R_init        = inter_result$R,
      clusters_init = inter_result$clusters,
      Gamma         = Gamma,
      W_cluster     = W_cluster,
      W_lasso       = W_lasso,
      lambda        = lambda[i],
      mu            = mu,
      eps_lasso     = eps_lasso,
      eps_f         = eps_f,
      EPS_CONV      = EPS_CONV,
      TOL_OPT       = TOL_OPT,
      MAX_ITER      = MAX_ITER
    )

    # Keep the results in the global variable
    hierarchy$results[[i]] <- inter_result
    names(hierarchy$results[[i]]$clusters) <- paste0("C", seq_along(hierarchy$results[[i]]$clusters))

    # Take usefull values
    hierarchy$results[[i]]$likelihood <- .Likelihood_penalised(
      R_init    = inter_result$R,
      clusters  = inter_result$clusters,
      Gamma     = Gamma,
      P         = P,
      W_cluster = W_cluster,
      W_lasso   = W_lasso,
      lambda    = lambda[i],
      mu        = mu,
      eps_lasso = eps_lasso
    )
    hierarchy$results[[i]]$lambda <- lambda[i]    # Regularised parameter
    hierarchy$results[[i]]$mu <- mu               # Lasso parameter
  }

  # ---- OUTPUT ----
  # Keep interesting computation and parameters
  # First variogram estimation in the output
  hierarchy$Gamma <- Gamma

  # Input parameters in the output
  hierarchy$inputs <- list()                 # Initialisation
  hierarchy$inputs$eps_conv <- EPS_CONV      # Convergence tolerance of the gradient descent
  hierarchy$inputs$eps_f <- eps_f            # Fusing threshold parameter
  hierarchy$inputs$tol_opt <- TOL_OPT        # Tolerance for the optimal step
  hierarchy$inputs$max_iter <- MAX_ITER      # Limit of iteration for the gradient descent

  # Naming the results with the corresponding lambda values
  names(hierarchy$results) <- paste0("l", seq_along(lambda))

  # Class of the results
  class(hierarchy) <- "HRC_hierarchy"

  return(hierarchy)
}