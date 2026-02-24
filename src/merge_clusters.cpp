#include <RcppEigen.h>
#include "merge_clusters.hpp"
#include "utils.hpp"
#include "model.hpp"
#include "distance_matrix.hpp"

using namespace Rcpp;


Eigen::VectorXd  merge_vector(Eigen::VectorXd a, Eigen::VectorXd b) {
  /* Merge two vectors into one
   *
   * Inputs:
   * a : a vector of size n
   * b : a vector of size m
   *
   * Output :
   * A vector of size n+m
   */
  int asize = a.size();
  int bsize = b.size();
  Eigen::VectorXd out(asize + bsize);

  for(int i = 0; i < (asize+bsize); i++) {
    if(i < asize){
      out(i) = a(i);
    } else {
      out(i) = b(i - asize);
    }
      
  }

  return out;
}


void fuse_R(Eigen::MatrixXd& R, Eigen::VectorXd p, int k, int l) {
  /* Update the column and size of R after merging step
   *
   * Inputs:
   * R : a matrix, the reduced matrix
   * p : a vector, the size of each clusters
   * k : an integer, the first cluster to fuse
   * l : an integer, the second cluster to fuse
   *
   * Output:
   * Void
   */
  // Initialization
  int K = p.size();   // Size of the current R

  Eigen::MatrixXd R_new(K - 1, K - 1);

  // The coefficients of the other clusters
  int row_idx = 0;
  for (int i = 0; i < K; i++) {
      if (i == l) continue; // To adjust the indices with no l row
      int col_idx = 0;
      for (int j = 0; j < K; j++) {
          if (j == l) continue; // To adjust the indices with no l column
          R_new(row_idx, col_idx) = R(i, j); // They do not change
          col_idx++;
      }
      row_idx++;
  }

  for (int i = 0; i < K - 1; i++) {
    if (i == k) continue; // skip the diagonal r_kk
    int old_i = (i >= l) ? i+1 : i; // adjust the indices without l coefficient
    // Update the coefficient by taking the average
    double val = ((p[k] * R(k, old_i) + p[l] * R(l, old_i)) / (p[k] + p[l])); 
    R_new(k, i) = val;
    R_new(i, k) = val;
  }
  // Update the value of r_kk
  R_new(k, k) = R(k, l);
  // Update the entire matrix
  R = R_new;
}


void cluster_fusion(Eigen::MatrixXd &R, List& clusters, int k, int l)
{
  /* Merging step of Clusterpath algorithm, update clusters list and R coefficients
   *
   * Inputs:
   * R : a matrix, the reduced matrix
   * clusters : a list of list, the list of clusters
   * k : an integer, the first cluster to fuse
   * l : an integer, the second cluster to fuse
   * 
   * Output :
   * Void
   */
  // Initialization
  int K = clusters.size();

  // Vector of cluster's size
  Eigen::VectorXd p = cluster_number(clusters);

  
  if (p.size() == 2)
  {
    // Case for fusion to one cluster
    R(0,0) = (p[0] * R(0, 0) + p[1] * R(0, 1)) / (p[0] + p[1]);
    R.conservativeResize(K - 1, K - 1);
    clusters = merge_vector(clusters[0], clusters[1]); 
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
    fuse_R(R, p, k, l);
  }
}
