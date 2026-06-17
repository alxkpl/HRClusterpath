#ifndef DISTANCE_MATRIX_HPP
#define DISTANCE_MATRIX_HPP

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

Eigen::MatrixXd distance_matrix(Eigen::MatrixXd R_matrix, List clusters);
double D_tilde2_r_term(Eigen::MatrixXd R_matrix, Eigen::VectorXd p_vector, int cluster_k, int cluster_l);

#endif
