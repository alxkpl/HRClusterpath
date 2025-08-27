#ifndef DISTANCE_MATRIX_HPP
#define DISTANCE_MATRIX_HPP

#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix distance_matrix(NumericMatrix R, List clusters);
double D_tilde2_r_term(NumericMatrix R, NumericVector p, int k, int l);

#endif
