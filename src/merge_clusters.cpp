#include "merge_clusters.hpp"
#include "distance_matrix.hpp"


IntegerVector which_min_upper(NumericMatrix mat) {
  int n = mat.nrow();
  int min_i = -1;
  int min_j = -1;
  double min_val = R_PosInf;

  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {  // diagonale supérieure
      if (mat(i, j) < min_val) {
        min_val = mat(i, j);
        min_i = i; 
        min_j = j;
      }
    }
  }

  return IntegerVector::create(min_i, min_j);
}

NumericVector merge_vector(NumericVector a, NumericVector b) {
  
  NumericVector out = a;

  // Copier les éléments de b
  for(int i = 0; i < b.size(); i++) {
      out.push_back(b[i]);
  }

  return out;
}

//' Function which merges clusters
//'
//' @param R K x K symmetric matrix.
//' @param clusters a list of vector : each vector gives the element of
//' a cluster.
//' @param eps_f a positive number : minimal tolerance for merging clusters
//'
//' @returns Returns, if merging, a list of the new clusters and the
//' corresponding R matrix, where the coefficient of the new clustered
//' is computing by averaging the coefficient of the two previous clusters.
//'
//' @keywords internal
//' @noRd
//[[Rcpp::export]]
List merge_clusters_rcpp(NumericMatrix R, List clusters, double eps_f)
{
  // Initialization
  int K = clusters.size();

  //Computation of the distance matrix
  NumericMatrix D = distance_matrix(R, clusters);

  // Search of the two potential clusters to merge
  IntegerVector index = which_min_upper(D);
  int k = index[0];
  int l = index[1];

  if (D(k, l) > eps_f)
  {
    return List::create(_["R"] = R, _["clusters"] = clusters); 
  }

  // Vector of cluster's size
  NumericVector p(K);
  for (int i = 0; i < K; i++)
  {
    List sub = clusters[i];
    p[i] = sub.size();
  }

  if (p.size() == 2)
  { 
    NumericMatrix R_new(1,1);
    R_new(0,0) =(p[0] * R(0, 0) + p[1] * R(0, 1)) / (p[0] + p[1]);
    return List::create(_["R"] = R_new,
        _["clusters"] = merge_vector(clusters[0], clusters[1])); 
  }

  // New clusters
  List new_clusters;

  for (int i = 0; i < K; i++) {
    if (i != l) { // enlever l
      if (i == k) 
      {
        NumericVector merged = clusters[k];
        NumericVector to_add = clusters[l];
        for (int j = 0; j < to_add.size(); j++) {
        merged.push_back(to_add[j]);
        }
        new_clusters.push_back(merged);
      } else 
      {
        new_clusters.push_back(clusters[i]);
      }
    }
  }

  // New R matrix
  NumericMatrix R_new(K - 1, K - 1);

  int row_idx = 0;
  for (int i = 0; i < K; i++) {
      if (i == l) continue;
      int col_idx = 0;
      for (int j = 0; j < K; j++) {
          if (j == l) continue;
          R_new(row_idx, col_idx) = R(i, j);
          col_idx++;
      }
      row_idx++;
  }

  for (int i = 0; i < K - 1; i++) {
    if (i == k) continue;
    int old_i = (i >= l) ? i+1 : i;
    double val = ((p[k] * R(k, old_i) + p[l] * R(l, old_i)) / (p[k] + p[l]));
    R_new(k, i) = val;
    R_new(i, k) = val;
  }

  R_new(k, k) = R(k, l);

  return List::create(
    _["R"] = R_new,
    _["clusters"] = new_clusters
  );

}