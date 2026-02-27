#' Negative log-likelihood for Hüsler-Reiss graphical models
#'
#' The penalised negative log-likelihood of the Hüsler-Reiss graphical models can be
#' written as follow:
#' \deqn{
#'  L_{\mathcal P}^{(n)}(\Theta, \lambda) = \underbrace{-\log(|\Theta|_+) -
#'  \frac 12 \text{tr}(\hat \Gamma^{(n)}\Theta)}_{L^{(n)}(\Theta)} + \lambda
#'  \underbrace{\sum_{i<j} w_{ij} d^2_{ij}(\Theta)}_{\mathcal P(\Theta)},
#' }
#' where \eqn{|\cdot|_+} is the generalised determinant, \eqn{n} is the sample size,
#' \eqn{\hat \Gamma^{(n)}} an estimation of the variogram matrix \eqn{\Gamma} and \eqn{w_{ij}>0} the weights.
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
#' @param R A \eqn{K \times K} matrix : the matrix of the clusters.
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
#' # Weight matrix
#' W <- matrix(c(0, 1, 1, 1,
#'               1, 0, 1, 1,
#'               1, 1, 0, 1,
#'               1, 1, 1, 0), nc = 4)
#'
#' # Random variogram
#' gamma <- matrix(c(0, 2, 1, 0,
#'                   2, 0, 4, 1,
#'                   1, 4, 0, 7,
#'                   0, 1, 7, 0), nc = 4)
#'
#' # Regularization parameter
#' lambda <- 2.5
#'
#' ############################################################
#' #                   NEGATIVE LOG-LIKEHOOD
#' ############################################################
#' NegLikelihood_HR(R, clusters, gamma)
#'
#' ############################################################
#' #                          PENALTY
#' ############################################################
#' penalty(R, clusters, W)
#'
#' ############################################################
#' #             PENALISED NEGATIVE LOG-LIKELIHOOD
#' ############################################################
#' NegLikelihood_penalised(R, clusters, gamma, W, lambda)
#'
NULL

#' @rdname likelihood
#'
#' @param Gamma A \eqn{d \times d} matrix: the variogram matrix \eqn{\Gamma}.
#'
#' @export
NegLikelihood_HR <- function(R, clusters, Gamma) {
  d <- ncol(Gamma)
  P <- .non_singular_P(d)

  return(
    .Likelihood_HR(R, clusters, Gamma, P)
  )
}

#' @rdname likelihood
#'
#' @param weights The \eqn{d \times d} symmetric matrix \eqn{W} with a zero diagonal.
#'
#' @export
penalty <- function(R, clusters, weights) {
  return(
    .Penalty(R, clusters, weights)
  )
}

#' @rdname likelihood
#'
#' @param Gamma A \eqn{d \times d} matrix: the variogram matrix \eqn{\Gamma}.
#' @param weights The \eqn{d \times d} symmetric matrix \eqn{W} with a zero diagonal.
#' @param lambda A positive number, the regularized parameter.
#'
#' @export
NegLikelihood_penalised <- function(R, clusters, Gamma, weights, lambda) {
  d <- ncol(Gamma)
  P <- .non_singular_P(d)

  return(
    .Likelihood_penalised(R, clusters, Gamma, P, weights, lambda)
  )
}
