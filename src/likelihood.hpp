#ifndef LIKELIHOOD_HPP
#define LIKELIHOOD_HPP

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

double Likelihood_penalised(
    Eigen::MatrixXd R_init,
    List clusters,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd P,
    const Eigen::MatrixXd W_cluster,
    const Eigen::MatrixXd W_lasso,
    double lambda,
    double mu,
    double eps_lasso
);

#endif
