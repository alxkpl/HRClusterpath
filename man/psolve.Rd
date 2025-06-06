% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_operation.R
\name{psolve}
\alias{psolve}
\title{Moore-Penrose pseudo inverse}
\usage{
psolve(A, tol = 1e-12)
}
\arguments{
\item{A}{a d x d symmetric positive semi-definite matrix.}
}
\value{
Computes the Moore-Penrose inverse of a matrix. The calculation is
done thanks to an article and if  :
A = L L^t
(with L having no zero-columns) then we have :
A^+ = L (L^t L)^-1 (L^t L)^-1 L^t
}
\description{
Moore-Penrose pseudo inverse
}
\examples{
A <- matrix(c(1, 2, 3,
              2, 5, 6,
              3, 6, 9), nc = 3)
psolve(A)
}
\keyword{internal}
