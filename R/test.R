# test.R

pacman::p_load(
  graphicalExtremes,
  CVXR
)

library(Rcpp)
library(RcppEigen)
sourceCpp("src/hrclusterpath.cpp")
#-----------------------------------------------------------------------------------------
#------------------------------------- SIMULATION -----------------------------------------
# Setup of the model and simulation
R <- matrix(c(1, -3, 0,
              -3, 2, -2,
              0, -2, 1), nc = 3)
clusters <- list(1:5, 6:10, 11:15)
#clusters <- list(1:3, 4:10, 11:15)
#clusters <- list(1:17, 18:34, 35:50)

# Construction of induced Theta and corresponding variogram gamma
Theta <- build_theta(R, clusters)
Gamma <- Theta2Gamma(Theta)
d <- ncol(Gamma)

# Simulation
set.seed(1234)
data <- rmstable(n = 10000, model = "HR", d = d, par = Gamma)
norm_inf <- apply(abs(data), 1, max)
quantile <- sort(norm_inf, decreasing = TRUE) |> as.numeric()

k <- floor(0.1 * nrow(data))
u <- quantile[k]

data_par <- data[norm_inf > u, ] / u

# Estimation of variogram
Gamma_est <- emp_vario(data_par)
Theta2 <- Gamma2Theta(Gamma_est)
zeta <- log(10000)**2
lambda <- 5


res <- HR_Clusterpath(data_par, 1, zeta)
res


lambda <- seq(0,10,1)
results <- HR_Clusterpath(data_par, lambda, zeta)

results


parallel::mclapply(lambda[1], function(l) HR_Clusterpath(data_par, zeta, l))
