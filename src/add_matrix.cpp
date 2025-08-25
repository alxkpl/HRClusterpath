#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void add_matrices_inplace(NumericMatrix& out, const NumericMatrix& A) {
  int nrow = out.nrow(), ncol = out.ncol();
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      out(i, j) += A(i,j);
    }
  }
}
