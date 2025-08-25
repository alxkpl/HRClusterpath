#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd psolve_rcpp(const Eigen::MatrixXd& L) {
    Eigen::MatrixXd M = L.transpose() * L;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(M.rows(), M.cols());
    Eigen::MatrixXd Minv = M.colPivHouseholderQr().solve(I);
    return L * Minv * Minv * L.transpose();
}