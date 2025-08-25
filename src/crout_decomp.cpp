#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix crout_decomposition_rcpp(NumericMatrix A, double tol = 1e-8) {
  int n = A.nrow();

  NumericVector d(n, 0.0);
  NumericMatrix L(n, n);

  // Initialisation dor d and the diagonal of L 
  d[0] = A(0, 0);

  for (int i = 0; i < n; i ++) {
      L(i, i) = 1;
  }


  for (int i = 1; i < n; i++) {
    // Non diagonal coefficient of L 
    for (int j = 0; j < i; j++) {
      double s = 0.0;
      for (int k = 0; k < i; k++) {
        s += L(i, k) * L(j, k) * d[k];
      }
      L(i, j) = (A(i, j) - s) / d[j];
    }

    // Diagonal coefficient of L
    double s2 = 0.0;
    for (int k = 0; k < n; k++) {
      s2 += d[k] * L(i, k) * L(i, k);
    }

    // Coefficient of the diagonal matrix D
    d[i] = A(i, i) - s2;
  }

  // Check for zero eigen values
  NumericMatrix D(n, n);
  for (int i = 0; i < n; i++) {
    if (d[i] > tol) {
      D(i, i) = std::sqrt(d[i]);
    } else {
      D(i, i) = 0.0;
    }
  }

  // results = L %*% D
  NumericMatrix out(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double s = 0.0;
      for (int k = 0; k < n; k++) {
        s += L(i, k) * D(k, j);
      }
      out(i, j) = s;
    }
  }

  return out;
}