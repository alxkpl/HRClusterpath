#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
int min_indx_cpp(int k, int l) {
    /* Compute the minimum between two numbers
     *
     * Inputs:
     * k : an integer
     * l : an integer
     *
     * Output:
     * Minimum
     */
    if (k < l) {
        return k;
    }
    return l;
}

// [[Rcpp::export]]
int max_indx_cpp(int k, int l) {
    /* Compute the maximum between two numbers
     *
     * Inputs:
     * k : an integer
     * l : an integer
     *
     * Output:
     * Maximum
     */
    if (k < l) {
        return l;
    }
    return k;
}
