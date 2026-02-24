#ifndef GRADIENT_DECENT_HPP
#define GRADIENT_DECENT_HPP

#include <Rcpp.h>
using namespace Rcpp;

// déclaration des fonctions

void cluster_fusion(Eigen::MatrixXd& R, List& clusters, int k, int l);

#endif