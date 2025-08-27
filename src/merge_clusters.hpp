#ifndef GRADIENT_DECENT_HPP
#define GRADIENT_DECENT_HPP

#include <Rcpp.h>
using namespace Rcpp;

// déclaration des fonctions
List merge_clusters_rcpp(NumericMatrix R, List clusters, double eps_f = 1e-1);

#endif