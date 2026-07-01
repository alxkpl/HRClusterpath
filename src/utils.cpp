#include <Rcpp.h>
#include <RcppEigen.h>
#include "utils.hpp"

using namespace Rcpp;     // for using List as Rcpp::List
// [[Rcpp::depends(RcppEigen)]]


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
    // ---- OUTPUT ---- //
    if (k < l) {
        return k;
    }
    return l;
}


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
    // ---- OUTPUT ---- //
    if (k < l) {
        return l;
    }
    return k;
}


Eigen::VectorXi which_min_upper(Eigen::MatrixXd matrix) {
    /* Extract the indices where the value is minimal
     *
     * Input :
     * matrix : a matrix
     *
     * Output:
     * A vector of size 2
     */
    // ---- INITIALIZATION ---- //
    int n_size = matrix.rows();
    if (n_size <= 1) {
        return Eigen::VectorXi::Constant(2, 0);
    }
    int min_i = -1;
    int min_j = -1;
    double min_val = R_PosInf;      // Infinity initialization

    // ---- COMPUTATION ---- //
    for (int i = 0; i < n_size - 1; i++) {
        for (int j = i + 1; j < n_size; j++) {  // Upper diagonal matrix
            if (matrix(i, j) < min_val) {
                min_val = matrix(i, j);
                min_i = i; 
                min_j = j;
            }
        }
    }
    Eigen::VectorXi indx(2);
    indx(0) = min_i;
    indx(1) = min_j;

    // ---- OUTPUT ---- //
    return indx;
}


Eigen::MatrixXd inverse(Eigen::MatrixXd A) {
    /* Compute the inverse matrix
     * 
     * Input :
     * A : a matrix
     * 
     * Output :
     * The inverse of A
     */
    // ---- INITIALIZATION ---- //
    const int p_size = A.rows();

    // ---- COMPUTATION ---- //
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(p_size, p_size);
    Eigen::MatrixXd inv_matrix = A.colPivHouseholderQr().solve(I);

    // ---- OUTPUT ---- //
    return inv_matrix;
}


// [[Rcpp::export(.non_singular_P)]]
Eigen::MatrixXd non_singular_P(int dim) {
    /* Compute the non sigular svd vectors of the projection in <1_n>^\perp
     * 
     * Input :
     * dim : an integer, the dimension
     * 
     * Output :
     * The matrix U without the singular vector
     */
    // ---- INITIALIZATION ---- //
    // Conctruction of the projection matrix Pi
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(dim);
    Eigen::MatrixXd P = Eigen::MatrixXd::Identity(dim, dim) - (1.0 / dim) * ones * ones.transpose();

    // ---- COMPUTATION ---- //
    // Computation of the SVD decomposition
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(P, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // ---- OUTPUT ---- //
    // Extraction of the non-sigular vectors
    return svd.matrixU().leftCols(dim - 1);
}


List simple_list(int size){
    /* Build the list from 1 to d
     * 
     * Input :
     * size : an integer, the size of the list
     * 
     * Output :
     * The list (1, 2, 3, ..., size)
     */
    // ---- INITIALIZATION ---- //
    List results;

    // ---- COMPUTATION ---- //
    for(int i = 1; i<=size; i++){
        results.push_back(i);
    }

    // ---- OUTPUT ---- //
    return results;
}

Eigen::MatrixXd E_matrix(int dim, int idx_1, int idx_2){
    /* Compute the symmetric matrix of indicators for the indices k and l
     * 
     * Input :
     * dim : an integer, the dimension
     * k : an integer, the first index
     * l : an integer, the second index
     * 
     * Output :
     * The matrix E_kl
     */
    // ---- INITIALIZATION ---- //
    Eigen::MatrixXd results = Eigen::MatrixXd::Zero(dim, dim);

    // ---- COMPUTATION ---- //
    results(idx_1, idx_2) = 1;
    results(idx_2, idx_1) = 1;

    // ---- OUTPUT ---- //
    return results;
}

double abs_penalty(double value, double eps_smooth) {
    /* Compute the value of the lasso penalty for a single value
     * 
     * Input :
     * value : a double, the value
     * eps_smooth : a double, the smooth parameter for the absolute value
     * 
     * Output :
     * The value of the lasso penalty for value
     */
    // ---- OUTPUT ---- //
    if (value >= -eps_smooth && value <= eps_smooth) {
        // Smoothed version if the value is under the tolerance
        return value * value / (2 * eps_smooth) + eps_smooth / 2;
    } else {
        // Absolute value otherwise
        return std::fabs(value);
    }
}