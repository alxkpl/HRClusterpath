#ifndef UTILS_HPP
#define UTILS_HPP

#include <RcppEigen.h>

using namespace Rcpp;   // to use List as Rcpp::List

// Functions declaration
int min_indx_cpp(int k, int l);

int max_indx_cpp(int k, int l);

Eigen::VectorXi which_min_upper(Eigen::MatrixXd mat);

Eigen::MatrixXd inverse(Eigen::MatrixXd A);

Eigen::MatrixXd non_singular_P(int d);

List simple_list(int d);

#endif