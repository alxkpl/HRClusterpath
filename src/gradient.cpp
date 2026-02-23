#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

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
     * Gradient
     */
    // Log-determinant part
    Eigen::MatrixXd dlog = - P * inverse(P.transpose() * Theta * P) * P.transpose();
    // Trace part
    Eigen::MatrixXd dtr = - 0.5 * P * P.transpose() * Gamma * P * P.transpose();
    return dlog + dtr;
}


Eigen::VectorXd Gradient_block(
    Eigen::MatrixXd R,
    List clusters,
    Eigen::MatrixXd Gamma,
    Eigen::MatrixXd P,
    int m
) {
    /* Compute the column gradient in a block structure
     * 
     * Input :
     * R : a matrix, the reduced matrix
     * clusters : a list of list, the list of the clusters
     * Gamma : a matrix, the fixed variogram
     * P : a matrix, computed with non_sigular_P in the right dimension
     * m : an integer, the column of the gradient
     *
     * Output :
     * The column gradient
     */
    const int K = clusters.size();
    Eigen::MatrixXd U = create_U(clusters);
    Eigen::MatrixXd D = Eigen::MatrixXd::Identity(K, K);
    D(m, m) = 0.5;

    return D * U.transpose() * Gradient_base(build_theta_cpp(R, clusters), Gamma, P) * U.col(m);
}


Eigen::VectorXd penalty_gradient(
    Eigen::MatrixXd R,
    List clusters,
    Eigen::MatrixXd tildeW,
    int m // careful it is indexed with minus 1 in C++
) {
    /* Compute the block gradient of the penalty
     * 
     * Input :
     * R : a matrix, the reduced matrix
     * clusters : a list of list, the list of the clusters
     * tildeW : a matrix, the cumulative weights per cluster
     * m : an integer, the column of the gradient
     * 
     * Output :
     * Gradient
     */
    int K = R.rows();
    Eigen::VectorXd p = cluster_number(clusters);
    Eigen::VectorXd results = Eigen::VectorXd::Zero(K);
    for(int k = 0; k < K; k++){
        if(k == m) continue;
        
        for(int i = 0; i < K; i++){
            if(i==k || i==m) continue;
            results(i) += 2 * tildeW(min_indx_cpp(k, m), max_indx_cpp(k, m)) * (R(m, i) - R(k, i));
        }

        results(m) += 2 * tildeW(min_indx_cpp(k, m), max_indx_cpp(k, m)) * (p(m) - 1) * (R(m, m) - R(k, m));
        results(k) += 2 * tildeW(min_indx_cpp(k, m), max_indx_cpp(k, m)) * (p(k) - 1) * (R(k, k) - R(k, m));

        if(k == K - 1) continue;
        for(int l = k + 1; l < K; l++) {
            if(l == m) continue;

            results(k) += 2 * tildeW(min_indx_cpp(k, l), max_indx_cpp(k, l)) * p(k) * (R(k, m) - R(l, m));
            results(l) += 2 * tildeW(min_indx_cpp(k, l), max_indx_cpp(k, l)) * p(l) * (R(l, m) - R(k, m));
        }
    }
    return results;
}


// [[Rcpp::export]]
Eigen::MatrixXd Gradient_penalised(
    Eigen::MatrixXd R,
    List clusters,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd P,
    const Eigen::MatrixXd tildeW,
    double lambda,
    int m
){
    /* Compute the gradient matrix for a column/row
     *
     * Input :
     * R : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     * Gamma : a matrix, the fixed variogram
     * P : a matrix, computed with non_sigular_P in the right dimension
     * tildeW : a matrix, the cumulative weights per cluster
     * m : an integer, the column of the gradient
     *
     * Output :
     * Gradient matrix
     */
    // Initialization
    int K = clusters.size();   // Dimension of the matrix (number of cluster)
    Eigen::MatrixXd results = Eigen::MatrixXd::Zero(K, K);   // Zeros everywhere


    Eigen::VectorXd d_llh = Gradient_block(R, clusters, Gamma, P, m);  // Gradient of the likelihood
    Eigen::VectorXd d_pen = penalty_gradient(R, clusters, tildeW, m);  // Gradient of the penalty
    Eigen::VectorXd value = d_llh + lambda * d_pen;     // Gradient values in the row/column

    results.col(m) = value;
    results.row(m) = value;

    return results;
}