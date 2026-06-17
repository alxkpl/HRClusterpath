#include <RcppEigen.h>
#include <Rcpp.h>
#include "utils.hpp"
#include "model.hpp"

Eigen::MatrixXd Hessian_base(
    Eigen::MatrixXd R_matrix,
    List clusters,
    Eigen::MatrixXd P,
    int idx_m
) {
    /* Compute the Hessian of the likelihood
     * 
     * Input :
     * R_matrix : a matrix, the reduced matrix
     * clusters : a list of list, the list of the clusters
     * P : a matrix, computed with non_sigular_P in the right dimension
     * idx_m : an integer, the column of the gradient
     * 
     * Output :
     * Hessian of the likelihood
     */
    // ---- INITIALIZATION ---- //
    const int D_VARIABLE = P.rows();                                        // Number of variables
    const int K_CLUSTER = clusters.size();                                  // Number of clusters
    const Eigen::MatrixXd p_vector = cluster_number(clusters);              // Vector with cluster's size
    const Eigen::MatrixXd U_matrix = create_U(D_VARIABLE, clusters);                    // Matrix of clusters U
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(K_CLUSTER, K_CLUSTER);        // Output Initialization
    Eigen::MatrixXd D = Eigen::MatrixXd::Identity(K_CLUSTER, K_CLUSTER);    // Symmetry multiplicator
    D(idx_m, idx_m) = 0.5;

    // ---- COMPUTATION ---- //

    // Gradient of the log-determinant part
    Eigen::MatrixXd M = P * inverse(P.transpose() * build_theta_cpp(D_VARIABLE, R_matrix, clusters) * P) * P.transpose();

    for(int l = 0; l < K_CLUSTER; l++){
        if(l == idx_m) continue;
        Eigen::MatrixXd E_ml = E_matrix(K_CLUSTER, idx_m, l);
        Eigen::MatrixXd B_ml = - p_vector(idx_m) * U_matrix.col(l).asDiagonal() - p_vector(l) * U_matrix.col(idx_m).asDiagonal();
        Eigen::MatrixXd M_ml = - M * (B_ml + U_matrix * E_ml * U_matrix.transpose()) * M;
        Eigen::MatrixXd UMU = U_matrix.transpose() * M_ml * U_matrix;

        // Hessian for k,l != m
        for(int k = l; k < K_CLUSTER; k++){
            if (k == idx_m) continue;
            double value = p_vector(idx_m) * (M_ml * U_matrix.col(k).asDiagonal()).trace() +  p_vector(k) * (M_ml * U_matrix.col(idx_m).asDiagonal()).trace() - 2 * UMU(idx_m, k);
            H(l, k) = value;
            H(k, l) = value;
        }

        // Hessian for k = m, l != m
        double value_ml = p_vector(idx_m) * (M_ml * U_matrix.col(idx_m).asDiagonal()).trace() - UMU(idx_m, idx_m);
        H(idx_m, l) = value_ml;
        H(l, idx_m) = value_ml;
    }

    // Hessian for l = k = m
    Eigen::MatrixXd E_mm = E_matrix(K_CLUSTER, idx_m, idx_m);
    Eigen::MatrixXd B_mm = - p_vector(idx_m) * U_matrix.col(idx_m).asDiagonal();
    Eigen::MatrixXd N_mm = - M * (B_mm + U_matrix * E_mm * U_matrix.transpose()) * M;
    H(idx_m, idx_m) = p_vector(idx_m) * (N_mm * U_matrix.col(idx_m).asDiagonal()).trace() - (U_matrix.transpose() * N_mm * U_matrix)(idx_m, idx_m);

    // ---- OUTPUT ---- //
    return H;

}

Eigen::MatrixXd Hessian_penalty(
    Eigen::MatrixXd R_matrix,
    List clusters,
    Eigen::MatrixXd W_cc,
    int idx_m
) {
    /* Compute the Hessian of the penalty
     * 
     * Input :
     * R_matrix : a matrix, the reduced matrix
     * clusters : a list of list, the list of the clusters
     * W_cc : a matrix, the cumulative weights per cluster
     * idx_m : an integer, the column of the gradient
     * 
     * Output :
     * Hessian of the clusterpath penalty
     */
    // ---- INITIALIZATION ---- //
    const int K_CLUSTER = R_matrix.rows();                           // Number of clusters
    const Eigen::VectorXd p_vector = cluster_number(clusters);       // Vector with cluster's size
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(K_CLUSTER, K_CLUSTER); // Output Initialization
    
    // ---- COMPUTATION ---- //
    double sum_W = W_cc.col(idx_m).sum() - W_cc(idx_m, idx_m);

    // The value of H out diagonal with column/row != m
    H = - 2 * p_vector(idx_m) * W_cc;

    // The value of H in the m-th column and row (symmetric)
    H.col(idx_m) = - W_cc.col(idx_m) * 2 * (p_vector(idx_m) - 1);
    H.row(idx_m) = - W_cc.col(idx_m) * 2 * (p_vector(idx_m) - 1);

    // A different value for the diagonal element (m, m)
    H(idx_m, idx_m) =  2 * (p_vector(idx_m) - 1) * sum_W;

    // The others values are 0 otherwise
    for(int q = 0; q < K_CLUSTER; q ++){
        if (q == idx_m) continue;
        double sum_q = W_cc.col(q).sum() - W_cc(q, q) - W_cc(q, idx_m);
        H(q, q) = 2 * p_vector(q) * (sum_W - W_cc(idx_m, q)) + 2 * (p_vector(idx_m) + p_vector(q) - 2) * W_cc(q, idx_m) + 2 * p_vector(idx_m) * sum_q;
    }

    // ---- OUTPUT ---- //
    return H;
}

Eigen::MatrixXd Hessian_lasso(
    Eigen::MatrixXd R_matrix,
    Eigen::MatrixXd W_lc,
    double eps_smooth,
    int idx_m
) {
    /* Compute the Hessian of the lasso penalty
     * 
     * Input :
     * R_matrix : a matrix, the reduced matrix
     * W_lc : a matrix, the cumulative lasso weights per cluster
     * eps_smooth : a double, the smooth parameter for the absolute value
     * idx_m : an integer, the column of the gradient
     * 
     * Output :
     * Hessian of the lasso penalty
     */
    // ---- INITIALIZATION ---- //
    const int K_CLUSTER = R_matrix.rows();                              // Number of clusters
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(K_CLUSTER, K_CLUSTER);    // Output Initialization

    // ---- COMPUTATION ---- //
    for (int k = 0; k < K_CLUSTER; k++){
        // Non null if the coefficient is under the smoothness threshold
        if(R_matrix(idx_m, k) >= -eps_smooth && R_matrix(idx_m, k) <= eps_smooth) H(k, k) = W_lc(idx_m, k) / eps_smooth;
    }

    // ---- OUTPUT ---- //
    return H;
}

Eigen::MatrixXd Hessian(
    Eigen::MatrixXd R_matrix,
    List clusters,
    Eigen::MatrixXd P,
    Eigen::MatrixXd W_cc,
    Eigen::MatrixXd W_lc,
    double lambda,
    double mu,
    double eps_lasso,
    int idx_m
) {
    /* Compute the Hessian of the penalised likelihood
     * 
     * Input :
     * R_matrix : a matrix, the reduced matrix
     * clusters : a list of list, the list of the clusters
     * P : a matrix, computed with non_sigular_P in the right dimension
     * W_cc : a matrix, the cumulative weights per cluster
     * W_lc : a matrix, the cumulative lasso weights per cluster
     * lambda : a double, the regularisation parameter
     * mu : a double, the lasso parameter
     * eps_lasso : a double, the smooth parameter for the absolute value
     * idx_m : an integer, the column of the gradient
     * 
     * Output :
     * Hessian of the lasso penalty
     */
    // ---- COMPUTATION ---- //
    Eigen::MatrixXd H_base = Hessian_base(R_matrix, clusters, P, idx_m);           // Hessian of the likelihood
    Eigen::MatrixXd H_pen = Hessian_penalty(R_matrix, clusters, W_cc, idx_m);      // Hessian of the clusterpath penalty

    // ---- OUTPUT ---- //
    if (mu > 0) {
        // Add the lasso penalty if mu > 0
        Eigen::MatrixXd H_lasso = Hessian_lasso(R_matrix, W_lc, eps_lasso, idx_m);  // Hessians of the lasso penalty
        return H_base + lambda * H_pen + mu * H_lasso;
    } else {
        return H_base + lambda * H_pen;
    }
}