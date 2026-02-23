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


// [[Rcpp::export]]
Eigen::MatrixXd inverse(Eigen::MatrixXd A) {
    /* Compute the inverse matrix
     * 
     * Input :
     * A : a matrix
     * 
     * Output :
     * The inverse of A
     */
    const int p = A.rows();
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(p, p);
    Eigen::MatrixXd x = A.colPivHouseholderQr().solve(I);
    return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd inverse(Eigen::MatrixXd A) {
    /* Compute the inverse matrix
     * 
     * Input :
     * A : a matrix
     * 
     * Output :
     * The inverse of A
     */
    const int p = A.rows();
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(p, p);
    Eigen::MatrixXd x = A.colPivHouseholderQr().solve(I);
    return x;
}
