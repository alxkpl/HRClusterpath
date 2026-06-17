#ifndef HESSIAN_HPP
#define HESSIAN_HPP

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

Eigen::MatrixXd Hessian(
    Eigen::MatrixXd R_matrix,
    List clusters,
    Eigen::MatrixXd P,
    Eigen::MatrixXd W_cc,
    Eigen::MatrixXd W_lc,
    double lambda,
    double mu,
    double eps_lasso,
    int idx_m
);

#endif