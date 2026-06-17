#include <Rcpp.h>
#include <RcppEigen.h>
#include "distance_matrix.hpp"
#include "model.hpp"

using namespace Rcpp;

double D_tilde2_r_term(Eigen::MatrixXd R_matrix, Eigen::VectorXd p_vector, int cluster_k, int cluster_l) {
    /* Compute the squared distance between two clusters
     *
     * Inputs :
     * R : a matrix, the reduced matrix
     * p_vector : a vector, the size of each cluster
     *
     * Ouput :
     * A distance
     */
  // ---- INITIALIZATION ---- //
  int K_CLUSTER = p_vector.size();    // Number of cluster
  double distance = 0.0;              // Distance initialization

  // ---- COMPUTATION ---- //
  for (int cluster_j = 0; cluster_j < K_CLUSTER; cluster_j++) {
    // Ponderated distance which depends of the index values
    double mask = (cluster_j == cluster_k || cluster_j == cluster_l) ? 1.0 : 0.0;  // p_k - 1 if the index is k (or l)
    double diff = R_matrix(cluster_k, cluster_j) - R_matrix(cluster_l, cluster_j);

    distance += (p_vector[cluster_j] - mask) * diff * diff;
  }

  // ---- OUTPUT ---- //
  return distance;
}

// [[Rcpp::export(.distance_matrix)]]
Eigen::MatrixXd distance_matrix(Eigen::MatrixXd R_matrix, List clusters) {
    /* Compute the squared distance matrix of the clustered precision matrix
     *
     * Inputs :
     * R_matrix : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     *
     * Ouput :
     * A distance matrix
     */
    // ---- INITIALIZATION ---- //
    int K_CLUSTER = clusters.size();                        // Number of cluster
    Eigen::VectorXd p_vector = cluster_number(clusters);    // Vector with cluster's size
    Eigen::MatrixXd result(K_CLUSTER, K_CLUSTER);

    // ---- COMPUTATION ---- //
    for (int k = 0; k < K_CLUSTER; k++) {
        for (int l = k; l < K_CLUSTER; l++) {
                double distance = D_tilde2_r_term(R_matrix, p_vector, k, l);   // Compute the distance
                result(k, l) = distance;
                result(l, k) = distance;
        }
    }

    // ---- OUTPUT ---- //
    return result;
}
