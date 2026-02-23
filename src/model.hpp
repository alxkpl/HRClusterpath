#ifndef MODEL_HPP
#define MODEL_HPP

#include <RcppEigen.h>

Eigen::MatrixXd create_U(List clusters);

Eigen::MatrixXd build_theta_cpp(
    Eigen::MatrixXd R,
    List clusters
);

#endif