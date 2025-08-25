#include <Rcpp.h>
#include "add_matrix.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix penalty_grad_rcpp(int K, Function f) {
    NumericMatrix out(K, K);
    std::fill(out.begin(), out.end(), 0.0);
    
    for (int i=0; i < K - 1; i++) {
        for (int j = i + 1; j < K; j++) {
            add_matrices_inplace(out, as<NumericMatrix>(f(i + 1, j + 1)));
        }
    }

    return out;
}
