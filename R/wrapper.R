#' Wrapper function for the HRClusterpath procedure implemented in .HRClusterpath (in C++).
#'
#' @keywords internal
.HRC_wrapper <- function(
  Gamma, zeta, lambda, mu, W_cluster, W_lasso, kappa,
  eps_lasso, eps_f, EPS_CONV, TOL_OPT, MAX_ITER
) {
  # ---- INITIALIZATION ----
  D_VARIABLE <- ncol(Gamma)               # Number of variables
  R.init <- Gamma2Theta(Gamma)            # First guess for the R matrix
  clusters.init <- as.list(1:D_VARIABLE)  # Initial clusters

  # Default clusterpath weights : exponential weights with parameter zeta
  if (is.null(W)) {
    W <- exp(
      - zeta * sqrt(.distance_matrix(Gamma, as.list(1:D_VARIABLE)))
    )
  }

  # Default merge threshold : data-driven threshold
  if (is.null(eps_f)) {
    # Base on data : median computed from the first guess for the precision matrix
    distance_median <- median(
      x     = sqrt(.distance_matrix(R.init, as.list(1:D_VARIABLE))) + diag(rep(NA, D_VARIABLE)),
      na.rm = TRUE
    )
    eps_f <-  kappa * distance_median
  }

  # Default sparsity weights : inverse absolute coefficient of the initial guess
  if (is.null(Z)) {
    Z <- 1 / abs(R.init) - diag(1 / abs(R.init))      # Null diagonal on weights matrix
    Z <- Z / sum(Z) * D_VARIABLE * (D_VARIABLE - 1)   # Standardized weights
  }

  # Matrix projection P for the likelihood computation
  P <- .non_singular_P(D_VARIABLE)

  # ---- COMPUTATION ----
  HRC_results <- .HRClusterpath(
    R_init        = R.init,
    clusters_init = clusters.init,
    Gamma         = Gamma,
    W_cluster     = W_cluster,
    W_lasso       = W_lasso,
    lambda        = lambda,
    mu            = mu,
    eps_lasso     = eps_lasso,
    eps_f         = eps_f,
    EPS_CONV      = EPS_CONV,
    TOL_OPT       = TOL_OPT,
    MAX_ITER      = MAX_ITER
  )

  # ---- OUTPUT ----
  # Keep the results in the output results
  HRC_results$likelihood <- .Likelihood_penalised(
    R_matrix  = HRC_results$R,
    clusters  = HRC_results$clusters,
    Gamma     = Gamma,
    P         = P,
    W_cluster = W_cluster,
    W_lasso   = W_lasso,
    lambda    = lambda,
    mu        = mu,
    eps_lasso = eps_lasso
  )
  HRC_results$lambda <- lambda    # Regularised parameter
  HRC_results$mu <- mu            # Lasso parameter

  names(HRC_results$clusters) <- paste0("C", seq_along(HRC_results$clusters))

  return(HRC_results)

}
