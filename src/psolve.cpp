#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix psolve_rcpp(const NumericMatrix& L_) {
    Eigen::MatrixXd L = Rcpp::as<Eigen::MatrixXd>(L_);
    Eigen::MatrixXd M = L.transpose() * L;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(M.rows(), M.cols());
    Eigen::MatrixXd Minv = M.colPivHouseholderQr().solve(I);
    Eigen::MatrixXd out =  L * Minv * Minv * L.transpose();

    return Rcpp::wrap(out); // to cahnge type to NumericMatrix
}