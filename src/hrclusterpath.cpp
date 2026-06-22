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
    double TOL_OPT,
    double lambda,
    double mu,
    double eps_lasso,
    int m
){
    /* Gradient update for the block gradient descent. 
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
     * TOL_OPT : a double, the tolerance for the optimal step size
     * lambda : a double, the regularization parameter
     * mu : a double, the lasso parameter
     * eps_lasso : a double, the smooth parameter for the absolute value
     * m : an integer, the index of the block to update
     *
     * Output :
     * Void, the function updates R inplace
     */
    // ---- COMPUTATION ---- //
    Eigen::MatrixXd gradient = Gradient_penalised(
        R, clusters, Gamma, P, tildeW, tildeZ, lambda, mu, eps_lasso, m
    );

    // Computation of the optimal step for the block gradient descent
    double m_step = max_step(R, clusters, gradient, m);     // Max step using inequalities
    if (m_step > 10) m_step = 10.0; // To improve code velocity
    double s = Gradient_step_cpp(R, clusters, Gamma, P, W, Z, lambda, mu, eps_lasso, gradient, TOL_OPT, 0, m_step);

    // Update the reduced matrix inplace 
    R -= s * gradient;
}


//[[Rcpp::export(.HRClusterpath)]]
List HRClusterpath_unique(
    Eigen::MatrixXd R_init,
    List clusters_init,
    const Eigen::MatrixXd Gamma,
    const Eigen::MatrixXd W_cluster,
    const Eigen::MatrixXd W_lasso,
    double lambda,
    double mu,
    double eps_lasso,
    double eps_f,
    double EPS_CONV,
    double TOL_OPT,
    int MAX_ITER
){
    /* Clusterpath algorithm for one value of lambda.
     *
     * Input :
     * R_init : a matrix, the initial reduced matrix
     * clusters_init : a list, the initial clusters
     * Gamma : a matrix, the estimated variogram
     * W_cluster : a matrix, the clusterpath weights for each variable
     * Z : a matrix, the lasso weights for each variable
     * lambda : a double, the regularization parameter
     * mu : a double, the lasso parameter
     * eps_lasso : a double, the smooth parameter for the absolute value
     * eps_f : a double, the threshold for merging clusters
     * EPS_CONV : a double, the threshold for convergence
     * TOL_OPT : a double, the tolerance for optimal step size computation
     * MAX_ITER : an integer, the maximum number of iterations
     * 
     * Output :
     * A list containing : 
     * - R : a matrix, the final reduced matrix
     * - clusters : a list, the final clusters
     */
    // ---- INITIALIZATION ---- //
    const int D_VARIABLE = Gamma.rows();             // Number of variables
    Eigen::MatrixXd P = non_singular_P(D_VARIABLE);  // Projection matrix

    // Initialization the results
    Eigen::MatrixXd R = R_init;             // First guess for the R matrix
    List clusters = clusters_init;          // Initial clusters

    // Likelihood values for the convergence criteria
    double l_new = Likelihood_penalised(R, clusters, Gamma, P, W_cluster, W_lasso, lambda, mu, eps_lasso);
    double l_old = 2 * l_new;

    Eigen::MatrixXd W_cc = clustered_weights(W_cluster, clusters);  // Clusterpath weights inside a cluster
    Eigen::MatrixXd W_lc = clustered_weights(W_lasso, clusters);    // Lasso weights inside a cluster

    // Initial value for the while loop
    int count = 0;                          // Counter for the gradient descent
    int K_max = R.cols();                   // Current number of cluster

    // ---- COMPUTATION ---- //
    while(std::abs((l_old / l_new) - 1.0) > EPS_CONV & count < MAX_ITER & std::isfinite(l_new)) {
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
                W_cc = clustered_weights(W_cluster, clusters);
                W_lc = clustered_weights(W_lasso, clusters);
            } else {
                // Gradient update
                Gradient_update(R, clusters, Gamma, P, W_cc, W_cluster, W_lc, W_lasso, TOL_OPT, lambda, mu, eps_lasso, k);
            }
            k++;
        }
        l_new = Likelihood_penalised(R, clusters, Gamma, P, W_cluster, W_lasso, lambda, mu, eps_lasso);
        count++;
    }

    // Check if some coefficients are under the smooth threshold
    if(mu > 0){
        for(int k = 0; k < K_max; k++){
            for(int l = k; l < K_max; l++){
                if (R(k, l) > - eps_lasso && R(k, l) < eps_lasso) {
                    R(k, l) = 0;
                    R(l, k) = 0;
                }
            }
        }
    }

    // ---- OUTPUT ---- //
    return List::create(
        _["R"] = wrap(R),
        _["clusters"] = wrap(clusters)
    );
}
