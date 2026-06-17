#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include "step_size.hpp"
#include "model.hpp"
#include "likelihood.hpp"

using namespace Rcpp;       // to use List as Rcpp::List


double max_step(
    Eigen::MatrixXd R_matrix,
    List clusters,
    Eigen::MatrixXd gradient,
    int idx_m
) {
    /* Compute the maximal step size for the computation of the optimal step
     *
     * Inputs :
     * R_matrix : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     * gradient : a matrix, a direction of the gradient
     * idx_m : an integer, the column of the gradient
     *
     * Output :
     * A positive number
     */
    // ---- INITIALIZATION ---- //
    const int K_CLUSTER = R_matrix.rows();                          // Number of clusters
    const Eigen::VectorXd p_vector = cluster_number(clusters);      // Vector with cluster's size 

    // ---- COMPUTATION ---- //
    double s = 1e7;
    for(int k = 0; k < K_CLUSTER; k++) {
        if (p_vector(k) == 1) continue;
        double numerator = (p_vector.transpose() * R_matrix.col(k))(0);
        if(k == idx_m){
            double grad_sum = (p_vector.transpose() * gradient.col(idx_m))(0);
            double value = (grad_sum < 0) ? numerator / grad_sum : s + 1;

            s = (s < value) ? s : value;
        } else {
            double value = (gradient(k, idx_m) < 0) ? numerator / (p_vector(idx_m) * gradient(k, idx_m)) : s + 1;
            s = (s < value) ? s : value;
        }
    }

    // ---- OUTPUT ---- //
    if(s < 0) return 1e-12;
    return s;
}


double Gradient_step_cpp(
    Eigen::MatrixXd R_init,
    List clusters,
    Eigen::MatrixXd Gamma,
    Eigen::MatrixXd P,
    Eigen::MatrixXd W_cluster,
    Eigen::MatrixXd W_lasso,
    double lambda,
    double mu,
    double eps_lasso,
    Eigen::MatrixXd gradient,
    double tol,
    double lo,
    double hi
) {
    /* Compute the optimal step size for the gradient descent using Golden-section search
     *
     * Inputs :
     * R_init : a matrix, the reduced matrix
     * clusters : a list of list, the list of clusters
     * Gamma : a matrix, the fixed variogram
     * P : a matrix, computed with non_sigular_P in the right dimension
     * W_cluster : a matrix, the matrix of weights
     * W_lasso : a matrix, the matrix of latent variables
     * lambda : double, the regularization parameter
     * mu : double, the Lasso regularization parameter
     * eps_lasso : double, the epsilon for the Lasso penalty
     * gradient : a matrix, a direction of the gradient
     * tol : positive number, the tolerance for convergence of the method
     * lo : a positive number, the value of the lower bound of the interval
     * hi : a positive number, the value of the lower bound of the interval
     *
     * Output :
     * A positive number
     */
    // ---- INITIALIZATION ---- //
    double phi = (1.0 + sqrt(5.0)) / 2.0;       // Golden ratio
    double a = lo;                              // Left bound initalization
    double b = hi;                              // Right bound initalization
    double c = b - (b - a) / phi;
    double d = a + (b - a) / phi;

    // Number of step according the tolerance and the interval length
    int n_steps = std::ceil(std::log(tol / (b - a)) / std::log(1 / phi));
    n_steps = std::max(n_steps, 2);

    // Value of the likelihood without gradient step
    double val_0 = Likelihood_penalised(R_init, clusters, Gamma, P, W_cluster, W_lasso, lambda, mu, eps_lasso);

    // ---- COMPUTATION ---- //
    // Values of the likelihood with the two gradient steps
    double val_c = Likelihood_penalised(R_init - c * gradient, clusters, Gamma, P, W_cluster, W_lasso, lambda, mu, eps_lasso);
    double val_d = Likelihood_penalised(R_init - d * gradient, clusters, Gamma, P, W_cluster, W_lasso, lambda, mu, eps_lasso);

    for(int i = 0; i < n_steps; i++){
        if(std::abs(a-b) < tol) break;      // Break the loop if the tolerance is reached

        if(std::abs(val_c - val_d) < 1e-10) break; // break the loop to avoid oscillation
        if(val_c < val_d) {
            // Update the interval where the minimum should be
            b = d;
        }
        else {
            // Update the interval where the minimum should be
            if(std::isnan(val_d) | std::isnan(val_c)) { // Reduce the interval to get value closer of zero which gives a valid matrix
                b = d;
            } else {
                a = c;
            }
        }

        // Update of the two gradient step
        c = b - (b - a) / phi;
        d = a + (b - a) / phi;

        // Update of the associated likelihood values
        val_c = Likelihood_penalised(R_init - c * gradient, clusters, Gamma, P, W_cluster, W_lasso, lambda, mu, eps_lasso);
        val_d = Likelihood_penalised(R_init - d * gradient, clusters, Gamma, P, W_cluster, W_lasso, lambda, mu, eps_lasso);

    }

    // Candidate for the optimal step
    double s = (a + b) / 2.0;

    double val_s = Likelihood_penalised(R_init - s * gradient, clusters, Gamma, P, W_cluster, W_lasso, lambda, mu, eps_lasso);

    // ---- OUTPUT ---- //
    // Checking if it does better than initial matrix, return 0.0 if it does not
    if (val_s > val_0) return 0.0;
    return s;
}