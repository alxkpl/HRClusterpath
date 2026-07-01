#include <Rcpp.h>
#include <RcppEigen.h>
#include "likelihood.hpp"
#include "utils.hpp"
#include "model.hpp"
#include "distance_matrix.hpp"

using namespace Rcpp;     // for using List as Rcpp::List


// [[Rcpp::export(.Likelihood_HR)]]
double Likelihood_raw(
    Eigen::MatrixXd R_matrix,
    List clusters,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd P
) {
    /* Compute the likelihood without penalty
     *
     * Inputs :
     * R_matrix : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     * Gamma : a matrix, the fixed variogram
     * P : a matrix, computed with non_sigular_P in the right dimension
     * tildeW : a matrix, the cumulative weights per cluster
     *
     * Output :
     * A real number
     */
    // ---- INITIALIZATION --- //
    const int D_VARIABLE = Gamma.rows();           // Number of variables

    // ---- COMPUTATION ---- //
    Eigen::MatrixXd Theta = build_theta_cpp(D_VARIABLE, R_matrix, clusters);   // Rebuild Theta from R and the list of clusters
    Eigen::MatrixXd det_mat = P.transpose() * Theta * P;    // Reduce dimension of Theta for the determinant
    Eigen::MatrixXd tr_mat = P * P.transpose() * Theta * P * P.transpose() * Gamma; // Matrix in the trace

    double det_value = det_mat.determinant();       // Compute determinant
    double tr_value = tr_mat.trace();               // Compute trace

    // ---- OUTPUT ---- //
    return  - std::log(det_value) - 0.5 * tr_value;
}

//[[Rcpp::export(.Penalty)]]
double Penalty(
    Eigen::MatrixXd R_matrix,
    List clusters,
    const Eigen::MatrixXd W_cluster
) {
    /* Compute the value of the Clusterpath penalty
     *
     * Inputs
     * R_matrix : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     * W_cluster : a matrix, the matrix of weights of the clusterpath penalty
     *
     * Output :
     * A real number
     */

    // ---- INITIALIZATION ---- //
    const int D_VARIABLE = W_cluster.rows();                           // Number of variables
    List list_var = simple_list(D_VARIABLE);
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(D_VARIABLE);

    // ---- COMPUTATION ---- //
    Eigen::MatrixXd Theta = build_theta_cpp(D_VARIABLE, R_matrix, clusters);   // Rebuild Theta from R and the list of clusters
    Eigen::MatrixXd inter_mat = (W_cluster.array() * distance_matrix(Theta, list_var).array()).matrix();    // Matrix of weights times the distance


    // ---- OUTPUT ---- //
    return (ones.transpose() * inter_mat * ones)(0) / 2.0;
}

double Lasso_penalty(
    Eigen::MatrixXd R_init,
    List clusters,
    const Eigen::MatrixXd W_lasso,
    double eps_smooth
) {
    /* Compute the value of the lasso penalty
     *
     * Inputs
     * R_init : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     * W_lasso : a matrix, the matrix of lasso weights
     * eps_smooth : a double, the smooth parameter for the absolute value
     *
     * Output :
     * A real number
     */
    // ---- INITIALIZATION ---- //
    const int D_VARIABLE = W_lasso.rows();                     // Number of variables
    Eigen::MatrixXd inter_mat = W_lasso;                       // For the output
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(D_VARIABLE);  // Constant vector

    // ---- COMPUTATION ---- //
    Eigen::MatrixXd Theta = build_theta_cpp(D_VARIABLE, R_init, clusters);  // Rebuild Theta from R and the list of clusters

    for (int i = 0; i < D_VARIABLE; i++) {
        for (int j = 0; j < D_VARIABLE; j++) {
            if (i == j) continue;
            inter_mat(i, j) *= abs_penalty(Theta(i, j), eps_smooth);
        }
    }

    // ---- OUTPUT ---- //
    return (ones.transpose() * inter_mat * ones)(0);
}

// [[Rcpp::export(.Likelihood_penalised)]]
double Likelihood_penalised(
    Eigen::MatrixXd R_init,
    List clusters,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd P,
    const Eigen::MatrixXd W_cluster,
    const Eigen::MatrixXd W_lasso,
    double lambda,
    double mu,
    double eps_lasso
) {
    /* Compute the penalised likelihood
     *
     * Inputs :
     * R_init : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     * Gamma : a matrix, the fixed variogram
     * P : a matrix, computed with non_sigular_P in the right dimension
     * W_cluster : a matrix, the matrix of weights
     * W_lasso : a matrix, the lasso weights
     * lambda : a double, the regularization parameter
     * mu : a double, the lasso parameter
     * eps_lasso : a double, the smooth parameter for the absolute value

     *
     * Output :
     * A real number
     */
    // ---- COMPUTATION ---- //
    double nllh = Likelihood_raw(R_init, clusters, Gamma, P);   // Negative log-kelihood
    double penalty = Penalty(R_init, clusters, W_cluster);      // Clusterpath penalty

    // ---- OUTPUT ---- //
    if(mu > 0){
        // With lasso penalty if mu > 0
        double lasso_pen = Lasso_penalty(R_init, clusters, W_lasso, eps_lasso);  // Lasso penalty
        return nllh + lambda * penalty + mu * lasso_pen;
    }

    return nllh + lambda * penalty;
}
