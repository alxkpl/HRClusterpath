#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include "step_size.hpp"
#include "model.hpp"
#include "likelihood.hpp"

using namespace Rcpp;       // to use List as Rcpp::List

// [[Rcpp::export]]
double max_step(
    Eigen::MatrixXd R,
    List clusters,
    Eigen::MatrixXd gradient,
    int m
) {
    const int K = R.rows();
    const Eigen::VectorXd p = cluster_number(clusters);

    double s = 1e7;
    for(int k = 0; k < K; k++) {
        // Need to find something else for clusters with one element
        if (p(k) == 1) continue;
        double numerator = (p.transpose() * R.col(k))(0);
        if(k == m){
            double grad_sum = (p.transpose() * gradient.col(m))(0);
            double value = (grad_sum < 0) ? numerator / grad_sum : s + 1;

            s = (s < value) ? s : value;
        } else {
            double value = (gradient(k, m) < 0) ? numerator / (p(m) * gradient(k, m)) : s + 1;
            s = (s < value) ? s : value;
        }
    }

    if(s < 0) return 1e-12;
    return s;
}

// [[Rcpp::export]]
double Gradient_step_cpp(
    Eigen::MatrixXd R_init,
    List clusters,
    Eigen::MatrixXd Gamma,
    Eigen::MatrixXd P,
    Eigen::MatrixXd W,
    double lambda,
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
     * W : a matrix, the matrix of weights
     * lambda : double, the regularization parameter
     * gradient : a matrix, a direction of the gradient
     * tol : positive number, the tolerance for convergence of the method
     * lo : a positive number, the value of the lower bound of the interval
     * hi : a positive number, the value of the lower bound of the interval
     *
     * Output :
     * A positive number
     */
    // Initialization of the method
    double phi = (1.0 + sqrt(5.0)) / 2.0;
    double a = lo;
    double b = hi;
    double c = b - (b - a) / phi;
    double d = a + (b - a) / phi;

    // Number of step according the tolerance and the interval length
    int n_steps = std::ceil(std::log(tol / (b - a)) / std::log(1 / phi));
    n_steps = std::max(n_steps, 2);

    // Value of the likelihood without gradient step
    double val_0 = Likelihood_penalised(R_init, clusters, Gamma, P, W, lambda);

    // Values of the likelihood with the two gradient steps
    double val_c = Likelihood_penalised(R_init - c * gradient, clusters, Gamma, P, W, lambda);
    double val_d = Likelihood_penalised(R_init - d * gradient, clusters, Gamma, P, W, lambda);

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
        val_c = Likelihood_penalised(R_init - c * gradient, clusters, Gamma, P, W, lambda);
        val_d = Likelihood_penalised(R_init - d * gradient, clusters, Gamma, P, W, lambda);

    }

    // Candidate for the optimal step
    double s = (a + b) / 2.0;

    double val_s = Likelihood_penalised(R_init - s * gradient, clusters, Gamma, P, W, lambda);

    // Checking if it does better than initial matrix, return 0.0 if it does not
    if (val_s > val_0) return 0.0;
    return s;
}