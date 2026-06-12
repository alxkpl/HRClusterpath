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
    const Eigen::MatrixXd tildeZ,
    const Eigen::MatrixXd Z,
    double tol_opt,
    double lambda,
    double mu,
    double eps_lasso,
    int m
){
    /* Gradient update for t block gradient descent. 
     *
     * Input :
     * R : a matrix, the current reduced matrix
     * clusters : a list, the list of clusters (list of vector of indices)
     * Gamma : a matrix, the estimated variogram
     * P : a matrix, the Projection matrix
     * tildeW : a matrix, the weights inside a cluster
     * W : a matrix, the weights for each variable
     * tildeZ : a matrix, lasso weights inside a cluster
     * Z : a matrix, lasso weights for each variable
     * tol_opt : a double, the tolerance for the optimal step size
     * lambda : a double, the regularization parameter
     * mu : a double, the lasso parameter
     * eps_lasso : a double, the smooth parameter for the absolute value
     * m : an integer, the index of the block to update
     *
     * Output :
     * Void, the function updates R in place
     */
    Eigen::MatrixXd gradient = Gradient_penalised(
        R, clusters, Gamma, P, tildeW, tildeZ, lambda, mu, eps_lasso, m
    );

    // Computation of the optimal step for the block gradient descent
    double m_step = max_step(R, clusters, gradient, m);     // Max step using inequalities
    if (m_step > 10) m_step = 10.0; // To improve code velocity
    double s = Gradient_step_cpp(R, clusters, Gamma, P, W, Z, lambda, mu, eps_lasso, gradient, tol_opt, 0, m_step);

    R -= s * gradient;
}


//[[Rcpp::export(.HRClusterpath)]]
List HRClusterpath_unique(
    Eigen::MatrixXd R_init,
    List clusters_init,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd W,
    const Eigen::MatrixXd Z,
    double lambda,
    double mu,
    double eps_lasso,
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
     * W : a matrix, the clusterpath weights for each variable
     * Z : a matrix, the lasso weights for each variable
     * lambda : a double, the regularization parameter
     * mu : a double, the lasso parameter
     * eps_lasso : a double, the smooth parameter for the absolute value
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
    double l_new = Likelihood_penalised(R, clusters, Gamma, P, W, Z, lambda, mu, eps_lasso);
    double l_old = 2 * l_new;

    Eigen::MatrixXd tildeW = clustered_weights(W, clusters);  // Weights inside a cluster
    Eigen::MatrixXd tildeZ = clustered_weights(Z, clusters);  // Lasso weights inside a cluster

    // Initial value for the while loop
    int count = 0;
    int K_max = R.cols();

    // Rcpp::Rcerr << "Optimization Initialization : \n";
    while(std::abs((l_old / l_new) - 1.0) > eps_conv & count < iter_max) {
        // Rcpp::Rcerr << "Step : " << count << " | Variation : " << std::abs((l_old / l_new) - 1.0) << "\n";

        int k = 0;
        l_old = l_new;
        while (k < K_max){
            // Merging step
            Eigen::MatrixXd d_mat = distance_matrix(R, clusters);
            Eigen::VectorXi indx = which_min_upper(d_mat);
            double min_dist = d_mat(indx(0), indx(1));

            // Check if the smaller distance is under the thresold eps_f
            if(min_dist < eps_f && K_max > 1) {
                // Merge clusters and update R coefficients
                cluster_fusion(R, clusters, indx(0), indx(1));
                // Change dimension size
                K_max = R.cols();
                tildeW = clustered_weights(W, clusters);
                tildeZ = clustered_weights(Z, clusters);
            } else {
                // Gradient update
                Gradient_update(R, clusters, Gamma, P, tildeW, W, tildeZ, Z, tol_opt, lambda, mu, eps_lasso, k);
            }
            k++;
        }
        l_new = Likelihood_penalised(R, clusters, Gamma, P, W, Z, lambda, mu, eps_lasso);
        count++;
    }

    // Check if some coefficients are under the smooth threshold
    for(int k = 0; k < K_max; k++){
        for(int l = k; l < K_max; l++){
            if (R(k, l) > - eps_lasso && R(k, l) < eps_lasso) {
                R(k, l) = 0;
                R(l, k) = 0;
            }
        }
    }

    // Rcpp::Rcerr << "Final variation : " << std::abs((l_old / l_new) - 1.0)  << ".\n" << "Optimization finished.\n";


    // Message for the user
    // if (count == iter_max) {
    //     Rcpp::Rcerr << "Warning : Maximum number of iterations reached. \n";
    // } else {
    //     Rcpp::Rcerr << "Convergence reached. \n";
    // }

    return List::create(
        _["R"] = wrap(R),
        _["clusters"] = wrap(clusters)
    );
}
