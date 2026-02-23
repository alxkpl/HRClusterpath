#include <Rcpp.h>
#include <RcppEigen.h>
#include "distance_matrix.hpp"
#include "model.hpp"

using namespace Rcpp;     // to use List as Rcpp::List

double D_tilde2_r_term(Eigen::MatrixXd R, Eigen::VectorXd p, int k, int l) {
  
  int K = p.size();

  double result = 0.0;

  for (int j = 0; j < K; j++) {
    double mask = (j == k || j == l) ? 1.0 : 0.0;
    double diff = R(k, j) - R(l, j);
    result += (p[j] - mask) * diff * diff;
  }

  return result;
}

// [[Rcpp::export]]
Eigen::MatrixXd distance_matrix(Eigen::MatrixXd R, List clusters) {
    /* Compute the squared distance matrix of the clustered precision matrix
     *
     * Inputs :
     * R : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     *
     * Ouput :
     * A distance matrix
     */
    int K = clusters.size();

    Eigen::VectorXd p = cluster_number(clusters);

    Eigen::MatrixXd  out(K, K);

    for (int k = 0; k < K; k++) {
        for (int l = k; l < K; l++) {
                double tmp = D_tilde2_r_term(R, p, k, l);
                out(k, l) = tmp;
                out(l, k) = tmp;
        }
    }
    return out;
}
