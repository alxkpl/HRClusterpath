#ifndef ADD_MATRICES_H
#define ADD_MATRICES_H

#include <Rcpp.h>
using namespace Rcpp;

// Déclaration de la fonction
void add_matrices_inplace(NumericMatrix& out, const NumericMatrix& A);

#endif