#ifndef MODEL_HPP
#define MODEL_HPP

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;     // for using List as Rcpp::List

Eigen::MatrixXd create_U(const int D_VARIABLE, List clusters);

Eigen::MatrixXd build_theta_cpp(
    const int D_VARIABLE,
    Eigen::MatrixXd R_matrix,
    List clusters
);

Eigen::VectorXd cluster_number(List clusters);

Eigen::MatrixXd clustered_weights(Eigen::MatrixXd weights, List clusters);

#endif