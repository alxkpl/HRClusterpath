#include <Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
double s_optimal(double s, Function f) {
    double s_out = s;
    while (as<bool>(f(s))) {
        s_out *= 0.95;
    }
    return s_out;
}