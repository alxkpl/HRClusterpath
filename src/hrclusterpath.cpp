#include <Rcpp.h>
#include <RcppEigen.h>
#include "utils.hpp"
#include "model.hpp"
#include "merge_clusters.hpp"
#include "distance_matrix.hpp"
#include "likelihood.hpp"
#include "step_size.hpp"
#include "gradient.hpp"

using namespace Rcpp;


void Gradient_update(
    Eigen::MatrixXd& R,
    List clusters,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd P,
    const Eigen::MatrixXd tildeW,
    const Eigen::MatrixXd W,
    double tol_opt,
    double lambda,
    int m
){
    /* Gradient update for t block gradient descent. 
     *
     * Input :
     * R : a matrixx, the current reduced matrix
     * clusters : a list, the list of clusters (list of vector of indices)
     * Gamma : a matrix, the estimated variogram
     * P : a matrix, the Projection matrix
     * tildeW : a matrix, the weights inside a cluster
     * W : a matrix, the weights for each variable
     * tol_opt : a double, the tolerance for the optimal step size
     * lambda : a double, the regularization parameter
     * m : an integer, the index of the block to update
     *
     * Output :
     * Void, the function updates R in place
     */
    Eigen::MatrixXd gradient = Gradient_penalised(
        R, clusters, Gamma, P, tildeW, lambda, m
    );

    // Computation of the optimal step for the block gradient descent
    double m_step = max_step(R, clusters, gradient, m);     // Max step using inequalities
    if (m_step > 10) m_step = 10.0; // To improve code velocity
    double s = Gradient_step_cpp(R, clusters, Gamma, P, W, lambda, gradient, tol_opt, 0, m_step);

    R -= s * gradient;
}


//[[Rcpp::export(.HRClusterpath)]]
List HRClusterpath_unique(
    Eigen::MatrixXd R_init,
    List clusters_init,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd W,
    double lambda,
    double eps_f,
    double eps_conv,
    double tol_opt,
    int iter_max
){
    /* Clusterpath algorithm for one value of lambda.
     *
     * Input :
     * R_init : a matrix, the initial reduced matrix
     * clusters_init : a list, the initial clusters
     * Gamma : a matrix, the estimated variogram
     * W : a matrix, the weights for each variable
     * lambda : a double, the regularization parameter
     * eps_f : a double, the threshold for merging clusters
     * eps_conv : a double, the threshold for convergence
     * tol_opt : a double, the tolerance for optimal step size computation
     * iter_max : an integer, the maximum number of iterations
     * 
     * Output :
     * A list containing : 
     * - R : a matrix, the final reduced matrix
     * - clusters : a list, the final clusters
     */
    // Projection matrix
    int d = Gamma.rows();           // Number of variables
    Eigen::MatrixXd P = non_singular_P(d);


    // Initialization the results
    Eigen::MatrixXd R = R_init;
    List clusters = clusters_init;

    // Likelihood values for the convergence criteria
    double l_new = Likelihood_penalised(R, clusters, Gamma, P, W, lambda);
    double l_old = 2 * l_new;

    Eigen::MatrixXd tildeW = clustered_weights(W, clusters);  // Weights inside a cluster

    // Initial value for the while loop
    int count = 0;
    int K_max = R.cols();

    Rcpp::Rcerr << "Optimization Initialization : \n";
    while(std::abs((l_old / l_new) - 1.0) > eps_conv & count < iter_max) {
        Rcpp::Rcerr << "Step : ";
        Rcpp::Rcerr << count;
        Rcpp::Rcerr << " | Variation : ";
        Rcpp::Rcerr << std::abs((l_old / l_new) - 1.0); 
        Rcpp::Rcerr << "\n";

        int k = 0;
        l_old = l_new;
        while (k < K_max){
            // Gradient update
            Gradient_update(R, clusters, Gamma, P, tildeW, W, tol_opt, lambda, k);

            // Merging step
            Eigen::MatrixXd d_mat = distance_matrix(R, clusters);
            Eigen::VectorXi indx = which_min_upper(d_mat);
            double min_dist = d_mat(indx(0), indx(1));

            // Check if the smaller distance is under the thresold eps_f
            if(min_dist < eps_f) {
                // Merge clusters and update R coefficients
                cluster_fusion(R, clusters, indx(0), indx(1));
                // Change dimension size
                K_max = R.cols();
                tildeW = clustered_weights(W, clusters);
            }
            k++;
        }
        l_new = Likelihood_penalised(R, clusters, Gamma, P, W, lambda);
        count++;
    }

    Rcpp::Rcerr << "Final variation : ";
    Rcpp::Rcerr << std::abs((l_old / l_new) - 1.0);
    Rcpp::Rcerr << ".\n";
    Rcpp::Rcerr << "Optimization finished.\n";


    // Message for the user
    if (count == iter_max) {
        Rcpp::Rcerr << "Warning : Maximum number of iterations reached. \n";
    } else {
        Rcpp::Rcerr << "Convergence reached. \n";
    }

    return List::create(
        _["R"] = R,
        _["clusters"] = clusters
    );
}
