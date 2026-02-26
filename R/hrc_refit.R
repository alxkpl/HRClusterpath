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


#' Refit the HR_clusterpath output for each results
#'
#' @param hrc_output An output of the `HR_Clusterpath` function.
#' @return An output of class `HR_Clusterpath_refit` with the same structure as the `HR_Clusterpath`
#' output but with refitted R matrices for each result.
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
