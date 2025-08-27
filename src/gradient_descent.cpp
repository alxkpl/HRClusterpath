#include "merge_clusters.hpp"
#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
List gradient_descent_rcpp(
    NumericMatrix R, List clusters, Function step,
    double lambda, int it_max, double eps_g, double eps_f
)
{
    int count = 0;
    double norm_grad = eps_g * eps_g + 1;
    NumericMatrix R_current = clone(R);
    List clusters_current = clone(clusters);

    while (count < it_max) {
        if (R_current.nrow() == 1 && R_current.ncol() == 1) break;
        if (norm_grad <= eps_g) break;
        norm_grad = 0;

        List grad_step = step(R_current, clusters_current, lambda);

        NumericMatrix grad = grad_step["gradient"];
        double step_size = as<double>(grad_step["step"]);

        // Update R 
        for (int i = 0; i < R_current.nrow(); i++) {
            for (int j = 0; j < R_current.ncol(); j++) {
                norm_grad += grad(i, j) * grad(i, j);
                R_current(i, j) -= step_size * grad(i, j);
            }
        }

        // Try to merge clusters via R function
        List res_merge = merge_clusters_rcpp(R_current, clusters_current, eps_f);
        
        NumericMatrix R_new = res_merge["R"];

        if (R_new.nrow() != R_current.nrow() || R_new.ncol() != R_current.ncol()) {
            R_current = R_new;
            clusters_current = res_merge["clusters"];
        }

        count++;
    }

    return List::create(
        _["R"] = R_current,
        _["clusters"] = clusters_current
    );
}