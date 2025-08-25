#include <Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
NumericMatrix distance_matrix(int n, Function f) {
    NumericMatrix out(n, n);
    for(int i = 0; i < n - 1; i++) {
        for(int j = i + 1; j < n; j++) {
            out(i, j) = as<double>(f(i + 1, j + 1));
        }
    }
    return out;
}

