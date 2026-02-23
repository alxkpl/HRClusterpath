#ifndef LIKELIHOOD_HPP
#define LIKELIHOOD_HPP

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

double Likelihood_penalised(
    Eigen::MatrixXd R,
    List clusters,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd P,
    const Eigen::MatrixXd W,
    double lambda
);

#endif
