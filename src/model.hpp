#ifndef MODEL_HPP
#define MODEL_HPP

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;     // for using List as Rcpp::List

Eigen::MatrixXd create_U(List clusters);

Eigen::MatrixXd build_theta_cpp(
    Eigen::MatrixXd R,
    List clusters
);

Eigen::VectorXd cluster_number(List clusters);

Eigen::MatrixXd clustered_weights(Eigen::MatrixXd W, List clusters);

#endif