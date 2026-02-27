#' Refit a Hüsler-Reiss model with a given clustering and variogram matrix
#'
#' @keywords internal
#'
.HRC_refit_wrapper <- function(clusters, Gamma) {
  # Initialization
  d <- ncol(Gamma)    # Number of variables
  P <- .non_singular_P(d) # Projection matrix
  U <- .create_U(clusters) # Cluster matrix
  K <- ncol(U) # Number of clusters

  # Problem formulation
  Theta <- CVXR::Variable(d, d, PSD = TRUE)   # The precision matrix to estimate
  R <- CVXR::Variable(K, K, symmetric = TRUE)       # The reduced matrix for Theta
  A <- CVXR::Variable(d)                      # The diagonal matrix for null sum constraint
  llh <- -CVXR::log_det(t(P) %*% Theta %*% P) - 1 / 2 * sum(CVXR::diag(Gamma %*% Theta)) # Likelihood on Theta

  # Objective and constraints
  objective <- CVXR::Minimize(llh)                    # Objective on the likelihood
  constraints <- list(                                # Constraints :
    P %*% t(P) %*% Theta %*% P %*% t(P) == Theta,     # Theta in SP_d^1
    Theta == CVXR::diag(A) + U %*% R %*% t(U)         # Theta with fixed block matrix structure
  )

  prob <- CVXR::Problem(objective, constraints)   # Problem definition
  solution <- CVXR::solve(prob)

  # Results
  results <- list()
  results$R <- solution$getValue(R)
  results$clusters <- clusters
  results$likelihood <- solution$value

  names(results$clusters) <- paste0("C", seq_along(results$clusters))

  return(results)
}


#' Hüsler-Reiss Clusterpath refitted
#'
#' The refit procedure is based on a convex optimization problem which is solved with the `CVXR` package.
#' Based on the results produced by `HR_Clusterpath`, a new problem is defined as the minimization of the
#' negative loglikelihood of a Hüsler-Reiss model with precision matrix \eqn{\Theta} under the constraints
#' of block matrix structure induced by the estimated clusters. These solutions are computed for each value of
#' \eqn{\lambda} in the output of `HR_Clusterpath`.
#'
#' @name hr-clusterpath-refit
#'
#' @param hrc_output An object of the `HR_Clusterpath` function.
#'
#' @return An object of class `HR_Clusterpath_refit` with the same structure as the `HR_Clusterpath`
#' output but with refitted R matrices for each result.
#'
#' @section Motivation :
#' The problem solution provided by `HR_clusterpath` may be biased due to the penalty in the
#' negative loglikelihood. To avoid this potential issue, the matrix \eqn{R} is refitted for
#' each solution over \eqn{\lambda} with the non-penalized negative loglikelihood constraining by the
#' estimated clusters.
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
#' n <- 1e4
#' d <- ncol(Gamma)
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
#'
#' # First estimation of R
#' HRC$results$l1$R
#'
#' # Refit of the Hüsler-Reiss Clusterpath
#' HRC_refit <- HR_Clusterpath_refit(HRC)
#'
#' # New estimation of R
#' HRC_refit$results$l1$R
NULL

#' @rdname hr-clusterpath-refit
#'
#' @export
HR_Clusterpath_refit <- function(hrc_output) {

  if (inherits(hrc_output, "HR_Clusterpath_refit")) {
    warning("Already refitted results")
    return(hrc_output)
  }

  Gamma <- hrc_output$Gamma
  results <- hrc_output$results

  for (i in seq_along(results)) {
    results[[i]] <- .HRC_refit_wrapper(results[[i]]$clusters, Gamma)
  }

  # Results of the refit procedure
  refit <- list()
  refit$results <- results
  refit$Gamma <- Gamma

  # Input parameters of the raw procedure
  refit$initial_inputs <- hrc_output$inputs


  class(refit) <- "HR_Clusterpath_refit"

  return(refit)

}
