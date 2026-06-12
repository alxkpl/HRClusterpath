#include <RcppEigen.h>
#include <Rcpp.h>
#include "utils.hpp"
#include "model.hpp"

Eigen::MatrixXd Hessian_base(
    Eigen::MatrixXd R,
    List clusters,
    Eigen::MatrixXd P,
    int m
) {
    /* Compute the Hessian of the likelihood
     * 
     * Input :
     * Theta : a matrix, the precision matrix
     * Gamma : a matrix, the fixed variogram
     * P : a matrix, computed with non_sigular_P in the right dimension
     * m : an integer, the column of the gradient
     * 
     * Output :
     * Hessian
     */
    //Initialization
    const int K = clusters.size();
    Eigen::MatrixXd p = cluster_number(clusters);
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(K, K);
    Eigen::MatrixXd U = create_U(clusters);
    Eigen::MatrixXd D = Eigen::MatrixXd::Identity(K, K);
    D(m, m) = 0.5;

    // Gradient of the log-determinant part
    Eigen::MatrixXd M = P * inverse(P.transpose() * build_theta_cpp(R, clusters) * P) * P.transpose();
    // 5 * sum(diag(U[, 1]) * Mkl) + 5 * sum(diag(U[, 2]) * Mkl) - 2 * (t(U) %*% Mkl %*% U)[1, 2] dr_ml dr_ml l!=m

    for(int l = 0; l < K; l++){
        if(l == m) continue;
        Eigen::MatrixXd E_ml = E_matrix(K, m, l);
        Eigen::MatrixXd B_ml = - p(m) * U.col(l).asDiagonal() - p(l) * U.col(m).asDiagonal();
        Eigen::MatrixXd M_ml = - M * (B_ml + U * E_ml * U.transpose()) * M;
        Eigen::MatrixXd UMU = U.transpose() * M_ml * U;

        // Hessian for k,l != m
        for(int k = l; k < K; k++){
            if (k == m) continue;
            double value = p(m) * (M_ml * U.col(k).asDiagonal()).trace() +  p(k) * (M_ml * U.col(m).asDiagonal()).trace() - 2 * UMU(m, k);
            H(l, k) = value;
            H(k, l) = value;
        }

        // Hessian for k = m, l != m
        double value_ml = p(m) * (M_ml * U.col(m).asDiagonal()).trace() - UMU(m, m);
        H(m, l) = value_ml;
        H(l, m) = value_ml;
    }

    // Hessian for l = k = m
    Eigen::MatrixXd E_mm = E_matrix(K, m, m);
    Eigen::MatrixXd B_mm = - p(m) * U.col(m).asDiagonal();
    Eigen::MatrixXd N_mm = - M * (B_mm + U * E_mm * U.transpose()) * M;
    H(m, m) = p(m) * (N_mm * U.col(m).asDiagonal()).trace() - (U.transpose() * N_mm * U)(m, m);

    
    return H;

}

Eigen::MatrixXd Hessian_penalty(
    Eigen::MatrixXd R,
    List clusters,
    Eigen::MatrixXd tildeW,
    int m
) {
    /* Compute the Hessian of the penalty
     * 
     * Input :
     * R : a matrix, the reduced matrix
     * clusters : a list of list, the list of the clusters
     * tildeW : a matrix, the cumulative weights per cluster
     * m : an integer, the column of the gradient
     * 
     * Output :
     * Hessian
     */
    const int K = R.rows();
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(K, K);
    Eigen::VectorXd p = cluster_number(clusters);
    
    double sum_W = tildeW.col(m).sum() - tildeW(m, m);

    // The value of H out diagonal with column/row != m
    H = - 2 * p(m) * tildeW;

    // The value of H in the m-th column and row (symmetric)
    H.col(m) = - tildeW.col(m) * 2 * (p(m) - 1);
    H.row(m) = - tildeW.col(m) * 2 * (p(m) - 1);

    // A different value for the diagonal element (m, m)
    H(m, m) =  2 * (p(m) - 1) * sum_W;

    // The others values are 0 otherwise
    for(int q = 0; q < K; q ++){
        if (q == m) continue;
        double sum_q = tildeW.col(q).sum() - tildeW(q, q) - tildeW(q, m);
        H(q, q) = 2 * p(q) * (sum_W - tildeW(m, q)) + 2 * (p(m) + p(q) - 2) * tildeW(q, m) + 2 * p(m) * sum_q;
    }

    return H;
    
}

Eigen::MatrixXd Hessian_lasso(
    Eigen::MatrixXd R,
    Eigen::MatrixXd tildeZ,
    double epsilon,
    int m
) {
    /* Compute the Hessian of the lasso penalty
     * 
     * Input :
     * R : a matrix, the reduced matrix
     * tildeZ : a matrix, the cumulative lasso weights per cluster
     * epsilon : a double, the smooth parameter for the absolute value
     * m : an integer, the column of the gradient
     * 
     * Output :
     * Hessian
     */
    const int K = R.rows();
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(K, K);
    for (int k = 0; k < K; k++){
        if(R(m, k) >= -epsilon && R(m, k) <= epsilon) H(k, k) = tildeZ(m, k) / epsilon;
    }

    return H;
}

Eigen::MatrixXd Hessian(
    Eigen::MatrixXd R,
    List clusters,
    Eigen::MatrixXd P,
    Eigen::MatrixXd tildeW,
    Eigen::MatrixXd tildeZ,
    double lambda,
    double mu,
    double eps_lasso,
    int m
) {
    Eigen::MatrixXd H_base = Hessian_base(R, clusters, P, m);
    Eigen::MatrixXd H_pen = Hessian_penalty(R, clusters, tildeW, m);
    if (mu > 0) {
        Eigen::MatrixXd H_lasso = Hessian_lasso(R, tildeZ, eps_lasso, m);
        return H_base + lambda * H_pen + mu * H_lasso;
    } else {
        return H_base + lambda * H_pen;
    }
}