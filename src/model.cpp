#include <RcppEigen.h>
#include "model.hpp"

Eigen::MatrixXd create_U(List clusters) {
     /* Compute the cluster matrix
     * 
     * Input :
     * clusters : a list of list, the list of the clusters
     * 
     * Output :
     * A matrix
     */
    // Extraction of the number of variable
    List all_ind;
    int K = clusters.size();
    for(int k =0; k < K; k++){
        List current_list = clusters(k);
        for(int i = 0; i < current_list.size(); i++){
            all_ind.push_back(current_list(i));
        }
    }
    int d = all_ind.size();

    // Building of the matrix of clusters U
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(d, K);

    for(int k = 0; k < K; k++){
        List current_list = clusters(k);
        for(int i = 0; i < current_list.size(); i++){
            // Indicate if the variable j belongs to cluster k
            int j = current_list[i];
            U(j - 1, k) = 1.0;      // The index is minus 1 in C++
        }
    }
    return U;
}



// [[Rcpp::export(build_theta)]]
Eigen::MatrixXd build_theta_cpp(
    Eigen::MatrixXd R,
    List clusters
) {
     /* Compute the precision matrix
     * 
     * Input :
     * R : a matrix, the reduced matrix
     * clusters : a list of list, the list of the clusters
     * 
     * Output :
     * The associated precision matrix
     */
    
    // Initialization
    Eigen::MatrixXd Theta;

    // Adaptation if there is only one cluster
    if (R.rows() == 1) {
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
        Theta = R(0, 0) * (ones * ones.transpose());
    }
    else
    { 
        // Expression in the non diagonal comes from URU^t
        Eigen::MatrixXd U = create_U(clusters);
        Theta = U * R * U.transpose();
    }

    // To get null row/column sum
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(Theta.rows());
    Eigen::VectorXd diag = Theta * ones;

    for(int i = 0; i < Theta.rows(); i++){
        Theta(i, i) -= diag(i);
    }
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
    const int K = clusters.size();
    Eigen::VectorXd results(K);

    for(int k = 0; k < K; k++){
        List current_list = clusters[k];
        results(k) = current_list.size();
    }

    return results;
}