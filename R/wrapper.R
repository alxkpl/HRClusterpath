#' Wrapper function for the HRClusterpath procedure implemented in .HRClusterpath (in C++).
#'
#' @keywords internal
.HRC_wrapper <- function(
  Gamma, zeta, lambda, W, kappa,
  eps_conv, eps_f, tol_opt, iter_max
) {
  # Number of variables
  d <- ncol(Gamma)

  # Computation of the initial Theta
  R.init <- Gamma2Theta(Gamma)

  # Initial clusters
  clusters.init <- as.list(1:d)

  # If no custom weights are given, we use the exponential weights
  if (is.null(W)) {
    W <- exp(-zeta * sqrt(distance_matrix(Gamma, as.list(1:d))))    # Choosen weights
  }

  # Adaptative threshold for the fusion step if no custom one is given
  if (is.null(eps_f)) {
    eps_f <-  kappa * median(sqrt(distance_matrix(R.init, as.list(1:d))) + diag(rep(NA, d)), na.rm = TRUE)
  }

  # Non singular matrix projection P for the likelihood computation
  P <- .non_singular_P(d)

  # Results of the Clusterpath procedure
  results <- .HRClusterpath(
    R.init,
    clusters.init,
    Gamma,
    W,
    lambda,
    eps_f,
    eps_conv,
    tol_opt,
    iter_max
  )

  # Likelihood value
  results$likelihood <- .Likelihood_penalised(results$R, results$clusters, Gamma, P, W, lambda)
  results$lambda <- lambda

  names(results$clusters) <- paste0("C", seq_along(results$clusters))

  return(results)

}
