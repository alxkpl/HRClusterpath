#' Wrapper function for the HRClusterpath procedure implemented in .HRClusterpath (in C++).
#'
#' @keywords internal
.HRC_wrapper <- function(
  data, zeta, lambda, p, W, kappa,
  eps_conv, eps_f, tol_opt, iter_max
) {
  # Number of variables
  d <- ncol(data)

  # Computation of the variogram and initial Theta
  Gamma_est <- emp_vario(data, p = p)
  R.init <- Gamma2Theta(Gamma_est)

  # Initial clusters
  clusters.init <- as.list(1:d)

  # If no custom weights are given, we use the exponential weights
  if (is.null(W)) {
    W <- exp(-zeta * sqrt(distance_matrix(Gamma_est, as.list(1:d))))    # Choosen weights
  }

  # Adaptative threshold for the fusion step if no custom one is given
  if (is.null(eps_f)) {
    eps_f <-  kappa * median(sqrt(distance_matrix(R.init, as.list(1:d))) + diag(rep(NA, d)), na.rm = TRUE)
  }

  # Non singular matrix projection P for the likelihood computation
  P <- non_singular_P(d)

  # Results of the Clusterpath procedure
  results <- .HRClusterpath(
    R.init,
    clusters.init,
    Gamma_est,
    W,
    lambda,
    eps_f,
    eps_conv,
    tol_opt,
    iter_max
  )

  # Likelihood value
  results$likelihood <- Likelihood_penalised(results$R, results$clusters, Gamma_est, P, W, lambda)

  # Input parameters
  results$lambda <- lambda
  results$inputs$eps_conv <- eps_conv
  results$inputs$eps_f <- eps_f
  results$inputs$tol_opt <- tol_opt
  results$inputs$iter_max <- iter_max


  names(results$clusters) <- paste0("C", seq_along(results$clusters))

  return(results)

}
