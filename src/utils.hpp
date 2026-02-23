#ifndef UTILS_HPP
#define UTILS_HPP

#include <RcppEigen.h>

// Functions declaration
int min_indx_cpp(int k, int l);

int max_indx_cpp(int k, int l);

Eigen::MatrixXd inverse(Eigen::MatrixXd A);

Eigen::MatrixXd non_singular_P(int d);

#endif