#ifndef GRADIENT_HPP
#define GRADIENT_HPP

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

Eigen::MatrixXd Gradient_penalised(
    Eigen::MatrixXd R,
    List clusters,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd P,
    const Eigen::MatrixXd tildeW,
    const Eigen::MatrixXd tildeZ,
    double lambda,
    double mu,
    double eps_lasso,
    int m
);

#endif