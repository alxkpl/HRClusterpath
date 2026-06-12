#ifndef HESSIAN_HPP
#define HESSIAN_HPP

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

Eigen::MatrixXd Hessian(
    Eigen::MatrixXd R,
    List clusters,
    Eigen::MatrixXd P,
    Eigen::MatrixXd tildeW,
    Eigen::MatrixXd tildeZ,
    double lambda,
    double mu,
    double eps_lasso,
    int m
);

#endif