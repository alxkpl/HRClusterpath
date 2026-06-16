#' Negative log-likelihood for Hüsler-Reiss graphical models
#'
#' The penalised negative log-likelihood of the Hüsler-Reiss graphical models can be
#' written as follow:
#' \deqn{
#'  L_{\mathcal P}^{(n)}(\Theta, \lambda) = \underbrace{-\log(|\Theta|_+) -
#'  \frac 12 \text{tr}(\hat \Gamma^{(n)}\Theta)}_{L^{(n)}(\Theta)} +
#'  \underbrace{\lambda\sum_{i<j} w_{ij} d^2_{ij}(\Theta) +
#'  \mu \sum_{i \neq j} z_{ij} |\theta_{ij}|}_{\mathcal P(\Theta)},
#' }
#' where \eqn{|\cdot|_+} is the generalised determinant, \eqn{n} is the sample size,
#' \eqn{\hat \Gamma^{(n)}} an estimation of the variogram matrix \eqn{\Gamma}, \eqn{w_{ij}>0} the clusterpath
#' weights and \eqn{z_{ij}>0} the lasso weights.
#'
#' An equivalent version of the negative likelihood is used by the function which is more convenient
#' for the optimization procedure. It is defined by the following:
#'
#' \deqn{
#'  L^{(n)}(\Theta, \lambda) = -\log(|P^t \Theta P|) - \frac 12 \text{tr}(\Gamma \Theta),
#' }
#' where \eqn{P} is a matrix \eqn{d\times (d-1)} matrix such that \eqn{P P^t = \Pi}, with \eqn{\Pi}
#' the projection matrix on the orthogonal space of \eqn{<\pmb{1}_d>}.
#'
#' The block matrix structure is used to compute the negative log-likelihood (see \code{\link{build_theta}()} and
#' \code{\link{extract_R_matrix}()}) not directly the value of the precision matrix \eqn{\Theta}.
#'
#' @name likelihood
#'
#' @param r_matrix A \eqn{K \times K} matrix : the matrix of the clusters.
#'
#' @param clusters A list of indices associated to a partition of \eqn{V}.
#'
#' @returns For `NegLikelihood_HR()`, it computes the value of the negative log-likelihood from two parameters:
#' - **R** the \eqn{K \times K} clusters matrix
#' - **clusters** a list of indices associated to a partition of \eqn{V}
#' and compute the value of the non penalised negative log-likelihood \eqn{L^{(n)}(\Theta)}.
#'
#' For `penalty()`, it computes the value of the penalty with:
#' - **R** the \eqn{K \times K} clusters matrix
#' - **clusters** a list of indices associated to a partition of \eqn{V}
#' and compute the value of the penalty \eqn{\mathcal P(\Theta)}.
#'
#' For `NegLikelihood_penalised()`, it computes the value of the penalised negative log-likelihood with:
#' - **R** : the \eqn{K \times K} clusters matrix
#' - **clusters** : a list of indices associated to a partition of \eqn{V}
#' - **lambda** : a positive number, the regularized parameter.
#' and compute the value of the penalised negative log-likelihood \eqn{L_{\mathcal P}^{(n)}(\Theta, \lambda)}.
#'
#' @examples
#' ############################################################
#' #                      INITIALIZATION
#' ############################################################
#' # Block matrix structure
#' R <- matrix(c(0.5, -1,
#'               -1, -1), nr = 2)
#' clusters <- list(c(1,3), c(2,4))
#'
#' # Weights matrices
#' W_cluster <- matrix(c(0, 1, 1, 1,
#'                       1, 0, 1, 1,
#'                       1, 1, 0, 1,
#'                       1, 1, 1, 0), nc = 4)
#'
#' W_lasso <- matrix(c(0, 1, 1, 1,
#'                     1, 0, 1, 1,
#'                     1, 1, 0, 1,
#'                     1, 1, 1, 0), nc = 4)
#'
#' # Random variogram
#' gamma <- matrix(c(0, 2, 1, 0,
#'                   2, 0, 4, 1,
#'                   1, 4, 0, 7,
#'                   0, 1, 7, 0), nc = 4)
#'
#' # Regularization parameter
#' lambda <- 2.5
#' mu <- 0.1
#'
#' ############################################################
#' #                   NEGATIVE LOG-LIKEHOOD
#' ############################################################
#' NegLikelihood_HR(R, clusters, gamma)
#'
#' ############################################################
#' #                          PENALTY
#' ############################################################
#' penalty(R, clusters, W_cluster)
#'
#' ############################################################
#' #             PENALISED NEGATIVE LOG-LIKELIHOOD
#' ############################################################
#' NegLikelihood_penalised(R, clusters, gamma, W_cluster, W_lasso, lambda, mu)
#'
NULL

#' @rdname likelihood
#'
#' @param Gamma A \eqn{d \times d} matrix: the variogram matrix \eqn{\Gamma}.
#'
#' @export
NegLikelihood_HR <- function(r_matrix, clusters, Gamma) {
  # ---- INITIALIZATION ----
  D_VARIABLE <- ncol(Gamma)           # Number of variable
  P_matrix <- .non_singular_P(D_VARIABLE)    # Matrix projection P

  # --- OUTPUT ---
  # Computed with Rcpp function (see src/likelihood.cpp)
  return(
    .Likelihood_HR(
      R        = r_matrix,
      clusters = clusters,
      Gamma    = Gamma,
      P        = P_matrix
    )
  )
}

#' @rdname likelihood
#'
#' @param W_cluster The \eqn{d \times d} symmetric matrix \eqn{W} with a zero diagonal.
#'
#' @export
penalty <- function(r_matrix, clusters, W_cluster) {
  # ---- OUTPUT ----
  # Computed with Rcpp function (see src/likelihood.cpp)
  return(
    .Penalty(
      R        = r_matrix,
      clusters = clusters,
      W        = W_cluster
    )
  )
}

#' @rdname likelihood
#'
#' @param Gamma A \eqn{d \times d} matrix: the variogram matrix \eqn{\Gamma}.
#'
#' @param W_cluster The \eqn{d \times d} symmetric matrix \eqn{W} with a zero diagonal. The weights for
#' the clusterpath penalty.
#'
#' @param W_lasso The \eqn{d \times d} symmetric matrix \eqn{W_{lasso}} with a zero diagonal. The
#' weights for the lasso penalty.
#'
#' @param lambda A positive number, the regularized parameter for clusterpath penalty.
#'
#' @param mu A positive number, the regularized parameter for lasso penalty.
#'
#' @param eps_lasso A small positive number, the parameter for the smoothed Lasso penalty.
#'
#' @export
NegLikelihood_penalised <- function(
  r_matrix, clusters, Gamma, W_cluster, W_lasso = NULL,
  lambda, mu = 0, eps_lasso = 5e-3
) {
  # ---- INITIALIZATION ----
  D_VARIABLE <- ncol(Gamma)                   # Number of variable
  P_matrix <- .non_singular_P(D_VARIABLE)     # Matrix projection P

  # Default weights : uniform weights for all coefficient
  if (is.null(W_lasso)) {
    W_lasso <- matrix(1, nc = D_VARIABLE, nr = D_VARIABLE) - diag(rep(1, D_VARIABLE))
  }

  # ---- OUTPUT ----
  # Computed with Rcpp function (see src/likelihood.cpp)
  return(
    .Likelihood_penalised(
      R         = r_matrix,
      clusters  = clusters,
      Gamma     = Gamma,
      P         = P_matrix,
      W         = W_cluster,
      Z         = W_lasso,
      lambda    = lambda,
      mu        = mu,
      eps_lasso = eps_lasso
    )
  )
}
