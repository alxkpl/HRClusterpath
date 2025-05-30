% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization.R
\name{HR_Clusterpath}
\alias{HR_Clusterpath}
\title{Search of optimal lambda for the penalty in optimization}
\usage{
HR_Clusterpath(data, chi, lambda, it_max = 1000)
}
\arguments{
\item{data}{n x d matrix : the data.}

\item{chi}{a positive number : tuned parameter for the exponential weights}

\item{l_grid}{a numerix vector of psoitive number : the grid line for lambda.}

\item{include_zero}{Boolean : if FALSE (default) avoid computation for
non-penalization setting (i.e. lambda = 0).}
}
\value{
Returns the optimal optimization from get_clusters() with the best
lambda in the grid line.
}
\description{
Search of optimal lambda for the penalty in optimization
}
