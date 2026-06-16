#' Clusterpath algorithm for Hüsler-Reiss models
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
#' \code{\link{neg_likelihood_pen}()}.
#'
#' @name hr-clusterpath
#'
#' @param data A \eqn{n \times d} matrix : the matrix of data.
#'
#' @param lambda If a positive number, the regularized parameter for the clsuterpath penalty. 
#' If a numerix vector of positive number : the grid line for \eqn{\lambda}.
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
#' If `W_c` is not `NULL`, this parameter is not used.
#'
#' @param p A numeric between 0 and 1, or `NULL`.  If `NULL` (default), it is
#' assumed that the data are already on multivariate Pareto scale. Else, `p` is used as the probability
#' in the function `data2mpareto()` to standardize the data (see `graphicalExtremes` documentation).
#'
#' @param eps_lasso A positive number : smoothness threshold for the absolute value.
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
#' @param tol_opt A positive number : tolerance for the optimal step in the gradient descent.
#'
#' @param MAX_ITER An integer : maximal number of iterations of the procedure.
#'
#' @returns `HR_Clusterpath()` is a block gradient descent from the data with default values and method for
#' the optimization : the variogram matrix \eqn{\Gamma} is the empirical variogram and the weights are
#' set as the exponential weights, which depends of only one tuning parameter \eqn{\zeta} and where
#' some theoretical guarantees are provided.
#'
#' The output is an object of class `HR_Clusterpath` which contains a list of several datas :
#'  - $results : a list of results for each value of lambda, where each result is a
#'               list with the following elements : $R the \eqn{R} matrix of the clusters, $clusters
#'               a list of the variable indices separated per cluster, $likelihood the value of the
#'               penalised negative loglikelihood and $lambda the corresponding regularized parameter.
#'  - $Gamma : the variogram matrix used for the procedure.
#'  - $inputs : a list of the input parameters used in the procedure.
#'
#' @section Some results for `HR_Clusterpath()`:
#'
#' When we get replications \eqn{X_1, \dots, X_n} of a random vector \eqn{X} which
#' belongs to the attraction domain of a \eqn{K}-block Hüsler-Reiss graphical model, the
#' minimum of the penalised negative loglikelihood conveges almost surely to the
#' true precision matrix of the model \eqn{\Theta^*}, for all \eqn{\lambda}, provided that
#' we choose a sequence \eqn{(\zeta_n)_{n\in \mathbb N^*}} which grows slower than the
#' \eqn{log(n)^s} sequence with \eqn{s < 8}.
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
#' HRC <- HR_Clusterpath(data_par, lambda = 10, zeta = zeta)
#' HRC
NULL


#' @rdname hr-clusterpath
#'
#' @importFrom graphicalExtremes emp_vario Gamma2Theta
#'
#' @export
HR_Clusterpath <- function(
  data, lambda, mu = 0, zeta, W_cluster = NULL,
  W_lasso = NULL, p = NULL, eps_lasso = 5e-3,
  eps_f = NULL, kappa = 1e-2, EPS_CONV = 1e-7,
  TOL_OPT = 1e-3, MAX_ITER = 1000
) {
  # ---- INITIALIZATION ----
  Gamma <- graphicalExtremes::emp_vario(
    data = data,
    p    = p
  )  # Computation of the empirical variogram

  HRC_output <- list()  # Initialisation fot the list of results

  # ---- COMPUTATION ----
  # Parallelization of the procedure for each value of lambda
  HRC_output$results <- parallel::mclapply(
    lambda,
    function(lambda_par) {
      .HRC_wrapper(
        Gamma     = Gamma,
        zeta      = zeta,
        lambda    = lambda_par,
        mu        = mu,
        W         = W_cluster,
        Z         = W_lasso,
        kappa     = kappa,
        eps_lasso = eps_lasso,
        eps_f     = eps_f,
        EPS_CONV  = EPS_CONV,
        TOL_OPT   = TOL_OPT,
        MAX_ITER  = MAX_ITER
      )
    },
    mc.cores = min(length(lambda), parallel::detectCores() - 1)
  )

  # ---- OUTPUT ----
  # Keep interesting computation and parameters
  # First variogram estimation in the output
  HRC_output$Gamma <- Gamma

  # Input parameters in the output
  HRC_output$inputs <- list()                 # Initialisation
  HRC_output$inputs$eps_conv <- EPS_CONV      # Convergence tolerance of the gradient descent
  HRC_output$inputs$eps_f <- eps_f            # Fusing threshold parameter
  HRC_output$inputs$tol_opt <- TOL_OPT        # Tolerance for the optimal step
  HRC_output$inputs$max_iter <- MAX_ITER      # Limit of iteration for the gradient descent

  # Naming the results with the corresponding lambda values
  names(HRC_output$results) <- paste0("l", seq_along(lambda))

  # Class of the results
  class(HRC_output) <- "HR_Clusterpath"

  return(HRC_output)
}
