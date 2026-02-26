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
#' @param data A \eqn{n \times d} matrix : the matrix of data.
#'
#' @param lambda If a positive number, the regularized parameter for the penalty. If a numerix
#' vector of positive number : the grid line for \eqn{\lambda}.
#'
#' @param W A \eqn{d \times d} matrix of positive number : the weights for the penalty.
#' If `NULL` (default), it is set as the exponential weights, which depends of only one
#' tuning parameter \eqn{\zeta} and where we have some theoretical results.
#'
#' @param zeta A positive number : tuned parameter for the exponential weights.
#' If `W` is not `NULL`, this parameter is not used.
#'
#' @param p A numeric between 0 and 1, or `NULL`.  If `NULL` (default), it is
#' assumed that the data are already on multivariate Pareto scale. Else, `p` is used as the probability
#' in the function `data2mpareto()` to standardize the data (see `graphicalExtremes` documentation).
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
#' @param iter_max An integer : maximal number of iteration of the procedure.
#'
#' @returns `HR_Clusterpath()` is a block gradient descent from the data with default values and method for
#' the optimization : the variogram matrix \eqn{\Gamma} is the empirical variogram and the weights are
#' set as the exponential weights, which depends of only one tuning parameter \eqn{\zeta} and where 
#' we have some theoretical results.
#'
#' The produce a list of results :
#'  - $R : the \eqn{R} matrix of the clusters.
#'  - $clusters : a list of the variable indices, clustered.
#'  - $likelihood : the value of the negative penalised negative loglikelihood.
#'  - $lambda : the value of lambda.
#'  - $inputs : a list of the input parameters used in the procedure.
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
#' ############################################################################
#' #                            With get_clusters
#' ############################################################################
#' # Customizable weights
#' W <- matrix(c(0, 1, 1, 1,
#'               1, 0, 1, 1,
#'               1, 1, 0, 1,
#'               1, 1, 1, 0), nc = 4)
#'
#' # Free choice of variogram
#' gamma <- graphicalExtremes::generate_random_Gamma(d = 4)
#'
#' # Choice of initial condition for the optimization
#' R <- matrix(c(1, 0, 0, -1,
#'               0, 1, 1, -2,
#'               0, 1, 1, -1,
#'               -1, -2, -1, 1), nc = 4)
#' lambda <- 2.2
#'
#'
#' ############################################################################
#' #                            With HR_Clusterpath
#' ############################################################################
#' # Construction of clusters and R matrix for simulation
#' R <- matrix(c(1, -3, 0,
#'               -3, 2, -2,
#'               0, -2, 1), nc = 3)
#' clusters <- list(1:5, 6:10, 11:15)
#'
#' # Construction of induced theta and corresponding variogram gamma
#' Theta <- build_Theta(R, clusters)
#' Gamma <- graphicalExtremes::Theta2Gamma(Theta)
#'
#' gr3_bal_sim_param_cluster <-
#'   list(
#'     R = R,
#'     clusters = clusters,
#'     Theta = Theta,
#'     Gamma = Gamma,
#'     chi = 1,
#'     n = 1e3,
#'     d = 15
#'   )
#'
#'
NULL


#' @rdname hr-clusterpath
#'
#'
#' @importFrom graphicalExtremes emp_vario Gamma2Theta
#'
#' @export
HR_Clusterpath <- function(
  data, lambda, zeta, W = NULL, p = NULL, eps_f = NULL, kappa = 1e-2,
  eps_conv = 1e-7, tol_opt = 1e-3, iter_max = 1000
) {
  # Parallelization of the procedure for each value of lambda
  results <- parallel::mclapply(
    lambda,
    function(l) {
      .HRC_wrapper(
        data = data,
        zeta = zeta,
        lambda = l,
        p = p,
        W = W,
        kappa = kappa,
        eps_conv = eps_conv,
        eps_f = eps_f,
        tol_opt = tol_opt,
        iter_max = iter_max
      )
    },
    mc.cores = min(length(lambda), parallel::detectCores() - 1)
  )

  # Naming the results with the corresponding lambda values
  names(results) <- paste0("l", seq_along(lambda))

  # Class of the results
  class(results) <- "HR_Clusterpath"

  return(results)
}
