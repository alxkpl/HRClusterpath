#include <Rcpp.h>
using namespace Rcpp;


double D_tilde2_r_term(NumericMatrix R, NumericVector p, int k, int l) {
  
  int K = p.size();

  double result = 0.0;

  for (int j = 0; j < K; j++) {
    double mask = (j == k || j == l) ? 1.0 : 0.0;
    double diff = R(k, j) - R(l, j);
    result += (p[j] - mask) * diff * diff;
  }

  return result;
}

//[[Rcpp::export]]
NumericMatrix distance_matrix(NumericMatrix R, List clusters) {
    int K = clusters.size();

    NumericVector p(K);
    for (int i = 0; i < K; i++) {
        List sub = clusters[i];
        p[i] = sub.size();
    }

    NumericMatrix out(K, K);

    for (int k = 0; k < K; k++) {
        for (int l = k; l < K; l++) {
                double tmp = D_tilde2_r_term(R, p, k, l);
                out(k, l) = tmp;
                out(l, k) = tmp;
        }
    }
    return out;
}
