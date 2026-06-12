#include <Rcpp.h>
#include <RcppEigen.h>
#include "likelihood.hpp"
#include "utils.hpp"
#include "model.hpp"
#include "distance_matrix.hpp"

using namespace Rcpp;     // for using List as Rcpp::List


// [[Rcpp::export(.Likelihood_HR)]]
double Likelihood_raw(
    Eigen::MatrixXd R,
    List clusters,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd P
) {
    /* Compute the likelihood without penalty
     *
     * Inputs :
     * R : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     * Gamma : a matrix, the fixed variogram
     * P : a matrix, computed with non_sigular_P in the right dimension
     * tildeW : a matrix, the cumulative weights per cluster
     *
     * Output :
     * A real number
     */
    Eigen::MatrixXd Theta = build_theta_cpp(R, clusters);
    Eigen::MatrixXd det_mat = P.transpose() * Theta * P;
    Eigen::MatrixXd tr_mat = Theta * Gamma;

    double det_value = det_mat.determinant();
    double tr_value = tr_mat.trace();

    return  - std::log(det_value) - 0.5 * tr_value;
}

//[[Rcpp::export(.Penalty)]]
double Penalty(
    Eigen::MatrixXd R,
    List clusters,
    const Eigen::MatrixXd W
) {
    /* Compute the value of the Clusterpath penalty
     *
     * Inputs
     * R : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     * W : a matrix, the matrix of weights
     *
     * Output :
     * A real number
     */
    Eigen::MatrixXd Theta = build_theta_cpp(R, clusters);
        const int d = Theta.rows();
    List list_var = simple_list(d);
    Eigen::MatrixXd inter_mat = (W.array() * distance_matrix(Theta, list_var).array()).matrix();
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(d);

    return (ones.transpose() * inter_mat * ones)(0) / 2.0;
}

double Lasso_penalty(
    Eigen::MatrixXd R,
    List clusters,
    const Eigen::MatrixXd Z,
    double epsilon
) {
    /* Compute the value of the lasso penalty
     *
     * Inputs
     * R : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     * Z : a matrix, the matrix of lasso weights
     * epsilon : a double, the smooth parameter for the absolute value
     *
     * Output :
     * A real number
     */
    Eigen::MatrixXd Theta = build_theta_cpp(R, clusters);
    const int d = Theta.rows();
    Eigen::MatrixXd inter_mat = Z;
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(d);

    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            if (i == j) continue;
            inter_mat(i, j) *= abs_penalty(Theta(i, j), epsilon);
        }
    }

    return (ones.transpose() * inter_mat * ones)(0);
}

// [[Rcpp::export(.Likelihood_penalised)]]
double Likelihood_penalised(
    Eigen::MatrixXd R,
    List clusters,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd P,
    const Eigen::MatrixXd W,
    const Eigen::MatrixXd Z,
    double lambda,
    double mu,
    double eps_lasso
) {
    /* Compute the penalised likelihood
     *
     * Inputs :
     * R : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     * Gamma : a matrix, the fixed variogram
     * P : a matrix, computed with non_sigular_P in the right dimension
     * W : a matrix, the matrix of weights
     * Z : a matrix, the lasso weights
     * lambda : a double, the regularization parameter
     * mu : a double, the lasso parameter
     * eps_lasso : a double, the smooth parameter for the absolute value

     *
     * Output :
     * A real number
     */
    double nllh = Likelihood_raw(R, clusters, Gamma, P);
    double penalty = Penalty(R, clusters, W);

    if(mu > 0){
        double lasso_pen = Lasso_penalty(R, clusters, Z, eps_lasso);
        return nllh + lambda * penalty + mu * lasso_pen;
    }

    return nllh + lambda * penalty;
}
