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
    if (k < l) {
        return l;
    }
    return k;
}


Eigen::VectorXd which_min_upper(Eigen::MatrixXd mat) {
  /* Extract the indices where the value is minimal
   *
   * Input :
   * mat : a matrix
   *
   * Output:
   * A vector of size 2
   */
  int n = mat.rows();
  int min_i = -1;
  int min_j = -1;
  double min_val = R_PosInf;

  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {  // Upper diagonal matrix
      if (mat(i, j) < min_val) {
        min_val = mat(i, j);
        min_i = i; 
        min_j = j;
      }
    }
  }
  Eigen::VectorXd indx(2);
  indx(0) = min_i;
  indx(1) = min_j;
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


List simple_list(int d){
    List results;
    for(int i = 1; i<=d; i++){
        results.push_back(i);
    }
    return results;
}