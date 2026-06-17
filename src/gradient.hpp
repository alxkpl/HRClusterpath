#ifndef GRADIENT_HPP
#define GRADIENT_HPP

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

Eigen::MatrixXd Gradient_penalised(
    Eigen::MatrixXd R_matrix,
    List clusters,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd P,
    const Eigen::MatrixXd W_cc,
    const Eigen::MatrixXd W_lc,
    double lambda,
    double mu,
    double eps_lasso,
    int idx_m
);

#endif