#include <Rcpp.h>
#include <RcppEigen.h>
#include "utils.hpp"

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

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
Eigen::MatrixXd non_singular_P(int d) {
    /* Compute the non sigular svd vectors of the projection in <1_n>^\perp
     * 
     * Input :
     * d : an integer, the dimension
     * 
     * Output :
     * The matrix U without the singular vector
     */
    // Conctruction of the projection matrix Pi
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(d);
    Eigen::MatrixXd P = Eigen::MatrixXd::Identity(d, d) - (1.0 / d) * ones * ones.transpose();

    // Computation of the SVD decomposition
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(P, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Extraction of the non-sigular vectors
    return svd.matrixU().leftCols(d - 1);
}


// [[Rcpp::export]]
Eigen::VectorXd cluster_number(List clusters) {
    /* Compute the vector of cluster's sizes
     *
     * Inputs:
     * clusters : a list of list, the list of the clusters
     *
     * Output:
     * Vector
     */
    const int K = clusters.size();
    Eigen::VectorXd results(K);

    for(int k = 0; k < K; k++){
        List current_list = clusters[k];
        results(k) = current_list.size();
    }

    return results;
}
