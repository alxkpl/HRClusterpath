#include <RcppEigen.h>
#include <Rcpp.h>
#include "utils.hpp"
#include "model.hpp"
#include "hessian.hpp"

using namespace Rcpp;     // to use List as Rcpp::List

Eigen::MatrixXd Gradient_base(
    Eigen::MatrixXd Theta,
    Eigen::MatrixXd Gamma,
    Eigen::MatrixXd P
) {
    /* Compute the gradient of the likelihood
     * 
     * Input :
     * Theta : a matrix, the precision matrix
     * Gamma : a matrix, the fixed variogram
     * P : a matrix, computed with non_sigular_P in the right dimension
     * 
     * Output :
     * General gradient
     */
    // ---- COMPUTATION --- //
    // Log-determinant part of the likelihood
    Eigen::MatrixXd dlog = - P * inverse(P.transpose() * Theta * P) * P.transpose();
    // Trace part of the likelihood
    Eigen::MatrixXd dtr = - 0.5 * P * P.transpose() * Gamma * P * P.transpose();

    // ---- OUTPUT ---- //
    return dlog + dtr;
}


Eigen::VectorXd Gradient_block(
    Eigen::MatrixXd R_matrix,
    List clusters,
    Eigen::MatrixXd Gamma,
    Eigen::MatrixXd P,
    int idx_m
) {
    /* Compute the column gradient in a block structure with only URU^t part
     * 
     * Input :
     * R_matrix : a matrix, the reduced matrix
     * clusters : a list of list, the list of the clusters
     * Gamma : a matrix, the fixed variogram
     * P : a matrix, computed with non_sigular_P in the right dimension
     * idx_m : an integer, the column of the gradient
     *
     * Output :
     * The column gradient
     */
    // ---- INITIALIZATION ---- //
    const int D_VARIABLE = Gamma.rows();
    const int K_CLUSTER = clusters.size();                                          // Number of clusters
    Eigen::MatrixXd U_matrix = create_U(D_VARIABLE, clusters);                                  // Matrix of clusters U
    Eigen::MatrixXd diag_matrix = Eigen::MatrixXd::Identity(K_CLUSTER, K_CLUSTER);  // Symmetry multiplicator
    diag_matrix(idx_m, idx_m) = 0.5;

    // ---- COMPUTATION ---- //
    Eigen::MatrixXd gradient = Gradient_base(build_theta_cpp(D_VARIABLE, R_matrix, clusters), Gamma, P);

    // ---- OUTPUT ---- //
    return 2 * diag_matrix * U_matrix.transpose() * gradient * U_matrix.col(idx_m);
}

Eigen::VectorXd correction_gradient(
    Eigen::MatrixXd R_matrix,
    List clusters,
    Eigen::MatrixXd Gamma,
    Eigen::MatrixXd P,
    int idx_m) {

    /* Compute the column gradient in a block structure with only the A part
     * 
     * Input :
     * R_matrix : a matrix, the reduced matrix
     * clusters : a list of list, the list of the clusters
     * Gamma : a matrix, the fixed variogram
     * P : a matrix, computed with non_sigular_P in the right dimension
     * idx_m : an integer, the column of the gradient
     *
     * Output :
     * The correction gradient for the block gradient descent
     */
    // ---- INITIALIZATION ---- //
    const int D_VARIABLE = Gamma.rows();                          // Number of variables
    const int K_CLUSTER = clusters.size();                        // Number of clusters
    const Eigen::MatrixXd U_matrix = create_U(D_VARIABLE, clusters);          // Matrix of clusters U
    const Eigen::VectorXd p_vector = cluster_number(clusters);    // Vector with cluster's size
    Eigen::MatrixXd M = P * inverse(P.transpose() * build_theta_cpp(D_VARIABLE, R_matrix, clusters) * P) * P.transpose();      // Gradient log-determinant 
    Eigen::MatrixXd Gamma_P = P * P.transpose() * Gamma * P * P.transpose();                                // Gradient trace
    Eigen::VectorXd result = Eigen::VectorXd::Zero(K_CLUSTER);      // Output initialization

    // ---- COMPUTATION ---- //
    // d/dr_m =  p_m * sum_(i \in C_k) (M_ii + 0.5 * Gamma_P_ii)
    result(idx_m) = p_vector(idx_m) * (U_matrix.col(idx_m).asDiagonal() * (M + 0.5 * Gamma_P)).trace();

    for(int k = 0; k < K_CLUSTER; k++) {
        if (k == idx_m) continue; // Already done for m

        // d/dr_k = pk * sum_(i \in C_m) (M_ii + 0.5 * Gamma_P_ii) + pm * sum_(i \in C_k) (M_ii + 0.5 * Gamma_P_ii) for k != m
        result(k) = p_vector(k) * (U_matrix.col(idx_m).asDiagonal() * (M + 0.5 * Gamma_P)).trace() + p_vector(idx_m) * (U_matrix.col(k).asDiagonal() * (M + 0.5 * Gamma_P)).trace();
    }

    // ---- OUTPUT ---- //
    return result;
}

Eigen::VectorXd penalty_gradient(
    Eigen::MatrixXd R_matrix,
    List clusters,
    Eigen::MatrixXd W_cc,
    int idx_m
) {
    /* Compute the block gradient of the penalty
     * 
     * Input :
     * R_matrix : a matrix, the reduced matrix
     * clusters : a list of list, the list of the clusters
     * W_cc : a matrix, the cumulative weights per cluster
     * idx_m : an integer, the column of the gradient
     * 
     * Output :
     * Gradient for the penalty
     */
    
    // ---- INITIALIZATION ---- //
    const int K_CLUSTER = R_matrix.rows();                            // Number of clusters
    const Eigen::VectorXd p_vector = cluster_number(clusters);        // Vector with cluster's size
    Eigen::VectorXd results = Eigen::VectorXd::Zero(K_CLUSTER); // Output initialization

    // ---- COMPUTATION ---- //
    for(int k = 0; k < K_CLUSTER; k++){
        if(k == idx_m) continue;
        // Gradient for dd_km / dr_mi
        for(int i = 0; i < K_CLUSTER; i++){
            if(i==k || i==idx_m) continue;
            results(i) += 2 * W_cc(min_indx_cpp(k, idx_m), max_indx_cpp(k, idx_m)) * p_vector(i) * (R_matrix(idx_m, i) - R_matrix(k, i));
        }

        results(idx_m) += 2 * W_cc(min_indx_cpp(k, idx_m), max_indx_cpp(k, idx_m)) * (p_vector(idx_m) - 1) * (R_matrix(idx_m, idx_m) - R_matrix(k, idx_m));
        results(k) += 2 * W_cc(min_indx_cpp(k, idx_m), max_indx_cpp(k, idx_m)) * ((p_vector(k) - 1) * (R_matrix(k, idx_m) - R_matrix(k, k)) + (p_vector(idx_m) - 1) * (R_matrix(k, idx_m) - R_matrix(idx_m, idx_m)));

        if(k == K_CLUSTER - 1) continue;
        // Gradient for dd_kl / dr_m.
        for(int l = k + 1; l < K_CLUSTER; l++) {
            if(l == idx_m) continue;

            results(k) += 2 * W_cc(min_indx_cpp(k, l), max_indx_cpp(k, l)) * p_vector(idx_m) * (R_matrix(k, idx_m) - R_matrix(l, idx_m));
            results(l) += 2 * W_cc(min_indx_cpp(k, l), max_indx_cpp(k, l)) * p_vector(idx_m) * (R_matrix(l, idx_m) - R_matrix(k, idx_m));
        }
    }

    // ---- OUTPUT ---- //
    return results;
}

Eigen::VectorXd lasso_gradient(
    Eigen::MatrixXd R_matrix,
    Eigen::MatrixXd W_lc,
    double eps_smooth,
    int idx_m
) {
    /* Compute the block gradient of the lasso penalty
     * 
     * Input :
     * R_matrix : a matrix, the reduced matrix
     * W_lc : a matrix, the cumulative lasso weights per cluster
     * epsilon : a double, the regularization parameter
     * m : an integer, the column of the gradient
     * 
     * Output :
     * Gradient
     */
    // ---- INITIALIZATION ---- //
    const int K_CLUSTER = R_matrix.rows();                // Number of clusters
    Eigen::VectorXd results = W_lc.col(idx_m);      // Output initialization 

    // ---- COMPUTATION ---- //
    for(int k = 0; k < K_CLUSTER; k++){
        // Alternative gradient if the value of the coefficient is under the smoothness threshold
        if(R_matrix(idx_m, k) >= - eps_smooth && R_matrix(idx_m, k) <= eps_smooth){
            results(k) = W_lc(idx_m, k) * W_lc(idx_m, k) / eps_smooth;
        } else {
            results(k) = (2 - (k == idx_m)) * W_lc(idx_m, k) * ((R_matrix(idx_m, k) > 0) - (R_matrix(idx_m, k) < 0));
        }
    }

    // ---- OUTPUT ---- //
    return results;
}

Eigen::MatrixXd Gradient_penalised(
    Eigen::MatrixXd R_matrix,
    List clusters,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd P,
    const Eigen::MatrixXd W_cc,
    const Eigen::MatrixXd W_lc,
    double lambda,
    double mu,
    double eps_lasso,
    int idx_m
){
    /* Compute the gradient matrix for a column/row of the penalized likelihood
     *
     * Input :
     * R : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     * Gamma : a matrix, the fixed variogram
     * P : a matrix, computed with non_sigular_P in the right dimension
     * W_cc : a matrix, the cumulative weights per cluster
     * W_lc : a matrix, the lasso weights per cluster
     * lambda : a double, the regularisation parameter
     * mu : a double, the lasso parameter
     * eps_lasso : a double, the smooth parameter for the absolute value
     * idx_m : an integer, the column of the gradient
     *
     * Output :
     * Block Gradient matrix
     */
    // ---- INITIALIZATION ---- //
    const int K_CLUSTER = clusters.size();                                   // Number of clusters
    Eigen::MatrixXd results = Eigen::MatrixXd::Zero(K_CLUSTER, K_CLUSTER);   // Output Initialization

    // ---- COMPUTATION ---- //
    Eigen::VectorXd d_llh = Gradient_block(R_matrix, clusters, Gamma, P, idx_m);  // Gradient of the likelihood
    d_llh += correction_gradient(R_matrix, clusters, Gamma, P, idx_m);            // Correction for the block gradient descent with A matrix
    Eigen::VectorXd d_pen = penalty_gradient(R_matrix, clusters, W_cc, idx_m);    // Gradient of the penalty 
    Eigen::VectorXd gradient = d_llh + lambda * d_pen;                            // Gradient in the row/column

    if (mu > 0) {
        // Add the gradient of the lasso penalty
        gradient += mu * lasso_gradient(R_matrix, W_lc, eps_lasso, idx_m);        
    }

    Eigen::MatrixXd hessian = Hessian(R_matrix, clusters, P, W_cc, W_lc, lambda, mu, eps_lasso, idx_m); // Hessian int the m-th column

    Eigen::VectorXd value = inverse(hessian) * gradient;   // Newton-Raphson direction for the gradient descent

    results.col(idx_m) = value;
    results.row(idx_m) = value;

    // ---- OUTPUT ---- //
    return results;
}