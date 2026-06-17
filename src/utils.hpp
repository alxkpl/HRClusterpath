#ifndef UTILS_HPP
#define UTILS_HPP

#include <RcppEigen.h>

using namespace Rcpp;   // to use List as Rcpp::List

// Functions declaration
int min_indx_cpp(int k, int l);

int max_indx_cpp(int k, int l);

Eigen::VectorXi which_min_upper(Eigen::MatrixXd matrix);

Eigen::MatrixXd inverse(Eigen::MatrixXd A);

Eigen::MatrixXd non_singular_P(int dim);

List simple_list(int dim);

Eigen::MatrixXd E_matrix(int dim, int idx_1, int idx_2);

double abs_penalty(double value, double eps_smooth);

#endif