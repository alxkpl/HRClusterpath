#include <Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
double s_optimal(double s, Function f) {
    double s_out = s;
    int cpt = 0;
    while ((as<bool>(f(s))) && (cpt < 100)) {
        s_out *= 0.95;
        cpt++;
    }
    return s_out;
}