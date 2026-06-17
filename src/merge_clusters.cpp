#include <RcppEigen.h>
#include "merge_clusters.hpp"
#include "utils.hpp"
#include "model.hpp"
#include "distance_matrix.hpp"

using namespace Rcpp;

void drop_column(Eigen::MatrixXd& R_matrix, int m) {
  /* Drop a column and a row from a matrix
   *
   * Inputs:
   * R_matrix : a matrix, the reduced matrix
   * m : an integer, the index of the column to drop
   *
   * Output :
   * Void, drop the column of the matrix inplace
   */
  // ---- INITIALIZATION ---- //
  const int K_CLUSTER = R_matrix.rows();

  // ---- COMPUTATION ---- //
  // To exit the function if the matrix is already of size 1 and avoid
  // memory issues
  if(K_CLUSTER == 1){
    R_matrix.resize(0, 0);
    return;
  }

  // Shift each with index larger than m one position upwards
  for (int k = 0; k < K_CLUSTER; k++) {
      for (int l = m; l < K_CLUSTER - 1; l++) {
          R_matrix(l, k) = R_matrix(l + 1, k);
      }
  }

  // Transpose and do the same for symmetry
  R_matrix.transposeInPlace();

  for (int k = 0; k < K_CLUSTER; k++) {
      for (int l = m; l < K_CLUSTER - 1; l++) {
          R_matrix(l, k) = R_matrix(l + 1, k);
      }
  }

  // Drop the last row and column to resize the matrix
  R_matrix.conservativeResize(K_CLUSTER - 1, K_CLUSTER - 1);
}

Eigen::VectorXd  merge_vector(Eigen::VectorXd vector_1, Eigen::VectorXd vector_2) {
  /* Merge two vectors into one
   *
   * Inputs:
   * vector_1 : a vector of size n
   * vector_2 : a vector of size m
   *
   * Output :
   * A vector of size n+m
   */
  // ---- INITIALIZATION ---- //
  int size_1 = vector_1.size();                // Size of the first vector
  int size_2 = vector_2.size();                // Size of the second vector
  Eigen::VectorXd merged(size_1 + size_2);     // Initialization of the output

  // ---- COMPUTATION ---- //
  for(int i = 0; i < (size_1 + size_2); i++) {
    if(i < size_1){
      // Fill with the first vector
      merged(i) = vector_1(i);
    } else {
      // Fill then with the oter
      merged(i) = vector_2(i - size_1);
    }
      
  }

  // ---- OUTPUT ---- //
  return merged;
}


void fuse_R(Eigen::MatrixXd& R_matrix, Eigen::VectorXd p_vector, int k, int l) {
  /* Update the column and size of R after merging step
   *
   * Inputs:
   * R_matrix : a matrix, the reduced matrix
   * p_vector : a vector, the size of each clusters
   * k : an integer, the first cluster to fuse
   * l : an integer, the second cluster to fuse
   *
   * Output:
   * Void
   */
  // ---- INITIALIZATION ---- //
  const int K_CLUSTER = p_vector.size();           // Size of the current R
  double size_total = p_vector[k] + p_vector[l];   // Size of the new cluster 
  double weight_k = p_vector[k] / size_total;      // Weight of cluster k
  double weight_l = p_vector[l] / size_total;      // Weight of cluster l

  // ---- COMPUTATION ---- //
  for(int i = 0; i < K_CLUSTER; i++) {
    if(i == k || i == l) continue;
    // Update the coefficient by taking the average
    double val = weight_k * R_matrix(k, i) + weight_l * R_matrix(l, i); 
    R_matrix(k, i) = val;
    R_matrix(i, k) = val;
  }

  R_matrix(k, k) = R_matrix(k, l);

  drop_column(R_matrix, l);   // Drop the l column
}


void cluster_fusion(Eigen::MatrixXd &R_matrix, List& clusters, int k, int l)
{
  /* Merging step of Clusterpath algorithm, update clusters list and R coefficients
   *
   * Inputs:
   * R_matrix : a matrix, the reduced matrix
   * clusters : a list of list, the list of clusters
   * k : an integer, the first cluster to fuse
   * l : an integer, the second cluster to fuse
   * 
   * Output :
   * Void
   */
  // ---- INITIALIZATION ---- //
  const int K_CLUSTER = clusters.size();
  Eigen::VectorXd p_vector = cluster_number(clusters);  // Vector of cluster's size

  // ---- COMPUTATION ---- //
  if (p_vector.size() == 2)
  {
    // Case for fusion to one cluster
    // Merge two last clusters
    clusters[0] = merge_vector(clusters[0], clusters[1]);
    clusters.erase(1);

    // Update R matrix
    R_matrix(0, 0) = (p_vector[0] * R_matrix(0, 0) + p_vector[1] * R_matrix(0, 1)) / (p_vector[0] + p_vector[1]);
    R_matrix.conservativeResize(K_CLUSTER - 1, K_CLUSTER - 1);
  }else{
    // Other cases
    // Cluster fusion step
    int f_indx = min_indx_cpp(k, l);
    int er_indx = max_indx_cpp(k, l);

    // Merge cluster k and l
    clusters[f_indx] = merge_vector(clusters[f_indx], clusters[er_indx]);
    // Delete the old cluster
    clusters.erase(er_indx);

    // Update R matrix step
    fuse_R(R_matrix, p_vector, k, l);
  }
}
