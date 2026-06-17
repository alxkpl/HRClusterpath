#include <RcppEigen.h>
#include <Rcpp.h>
#include "model.hpp"

using namespace Rcpp;     // for using List as Rcpp::List

//[[Rcpp::export(.create_U)]]
Eigen::MatrixXd create_U(const int D_VARIABLE, List clusters) {
     /* Compute the cluster matrix
     * 
     * Input :
     * D_VARIABLE : an integer, the number of variables
     * clusters : a list of list, the list of the clusters
     * 
     * Output :
     * A matrix
     */
    // ---- INITIALIZATION ---- //
    int K_CLUSTER = clusters.size();
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(D_VARIABLE, K_CLUSTER);

    for(int k = 0; k < K_CLUSTER; k++){
        List current_list = clusters(k);
        for(int i = 0; i < current_list.size(); i++){
            // Indicate if the variable j belongs to cluster k
            int j = current_list[i];
            U(j - 1, k) = 1.0;      // The index is minus 1 in C++
        }
    }
    return U;
}


// [[Rcpp::export(.build_theta)]]
Eigen::MatrixXd build_theta_cpp(
    const int D_VARIABLE,
    Eigen::MatrixXd R_matrix,
    List clusters
) {
     /* Compute the precision matrix
     * 
     * Input :
     * D_VARIABLE : an integer, the number of variables
     * R : a matrix, the reduced matrix
     * clusters : a list of list, the list of the clusters
     * 
     * Output :
     * The associated precision matrix
     */
    
    // ---- INITIALIZATION ---- //
    Eigen::MatrixXd Theta(D_VARIABLE, D_VARIABLE);

    // ---- COMPUTATION ---- //
    // Adaptation if there is only one cluster
    if (R_matrix.rows() == 1) {
        int K;
        // If the format of the cluster list is different
        if (clusters.size() > 1) {
            K = clusters.size();
        }
        else {
            List one_list = clusters[0];
            K = one_list.size();
        }
        // In that case, there is only the value of R in the non diagonal
        Eigen::VectorXd ones = Eigen::VectorXd::Ones(K);
        Theta = R_matrix(0, 0) * (ones * ones.transpose());
    }
    else
    { 
        // Expression in the non diagonal comes from URU^t
        Eigen::MatrixXd U = create_U(D_VARIABLE, clusters);
        Theta = U * R_matrix * U.transpose();
    }

    // To get null row/column sum
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(D_VARIABLE);
    Eigen::VectorXd diag = Theta * ones;

    for(int i = 0; i < D_VARIABLE; i++){
        Theta(i, i) -= diag(i);
    }

    // ---- OUTPUT ---- //
    return Theta;
}


Eigen::VectorXd cluster_number(List clusters) {
    /* Compute the vector of cluster's sizes
     *
     * Inputs:
     * clusters : a list of list, the list of the clusters
     *
     * Output:
     * Vector
     */
    // ---- INITIALIZATION ---- //
    const int K_CLUSTER = clusters.size();      // Number of clusters
    Eigen::VectorXd results(K_CLUSTER);         // Output initialization

    // ---- COMPUTATION ---- //
    for(int k = 0; k < K_CLUSTER; k++){
        List current_list = clusters[k];
        results(k) = current_list.size();      // Size of cluster k
    }

    // ---- OUTPUT ---- //
    return results;
}


Eigen::MatrixXd clustered_weights(Eigen::MatrixXd weights, List clusters){
    /* Compute cumulative weights per cluster
     *
     * Inputs:
     * weights : a matrix, weights
     * clusters : a list of list, the list of the clusters
     *
     * Output:
     * The matrix of cumulative weights
     */
    // ---- INITIALIZATION ---- //
    int D_VARIABLE = weights.rows();                              // Number of variables
    Eigen::MatrixXd U = create_U(D_VARIABLE, clusters);     // Cluster matrix U

    // Null diagonal constraint to get real weight matrix
    for(int i = 0; i < D_VARIABLE; i++){
        weights(i, i) = 0;
    }

    // ---- OUTPUT ---- //
    return U.transpose() * weights * U;
}