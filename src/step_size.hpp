#ifndef STEP_SIZE_HPP
#define STEP_SIZE_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;      // to use Liste as Rcpp::List

double max_step(
    Eigen::MatrixXd R,
    List clusters,
    Eigen::MatrixXd gradient,
    int m
);

double Gradient_step_cpp(
    Eigen::MatrixXd R_init,
    List clusters,
    Eigen::MatrixXd Gamma,
    Eigen::MatrixXd P,
    Eigen::MatrixXd W,
    double lambda,
    Eigen::MatrixXd gradient,
    double tol,
    double lo,
    double hi
);

#endif