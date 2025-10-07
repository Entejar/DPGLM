rm(list = ls())
# Load Libraries
# require(foreach)
# require(doParallel)
source("load.R")

# # Parallel Setup
# num_cores <- 4
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)
# 
# # Stop Parallel Backend
# stopCluster(cl)

# Simulation ---------------------------------------------------------------------
f0_true_kde <- readRDS("sim_truth/f0_kde_unif_0_1.rds")
spt <- f0_true_kde$spt
f0_kde  <- f0_true_kde$f0 

dx <- diff(spt)[1]
m0 <- sum(spt * f0_kde) * dx

#sum(f0_kde) * diff(spt)[1]

truth1 <- list(beta = c(1, 0), 
               f0 = matrix(c(spt, f0_kde), 
                           ncol = 2, 
                           byrow = FALSE))

truth2 <- list(beta = c(0.2, 0.7), 
               f0 = matrix(c(spt, f0_kde), 
                           ncol = 2, 
                           byrow = FALSE))


sim_data <- function(n, true_beta) {
  X     <- cbind(1, runif(n, -sqrt(12) / 4, sqrt(12) / 4))
  X[, -1] <- scale(X[, -1])
  beta <- true_beta
  mu <-  exp(X %*% beta) / (1 + exp(X %*% beta))
  theta <- gldrm:::getTheta(
    spt = spt,
    f0  = f0_kde,
    mu  = mu,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = NULL
  )$theta
  
  #theta <- rep(0, n) # for testing
  
  # Compute the matrix of exponentiated outer products
  exp_matrix <- exp(outer(theta, spt, "*"))
  
  dx <- diff(spt)[1]
  # Compute the dot product row-wise and apply the log
  btheta <- log(exp_matrix %*% f0_kde * dx)
  
  X       <- X %>% as.matrix()
  fY <- t(sapply(1:n, function(i) {
    exp(theta[i] * spt - btheta[i]) * f0_kde 
  }))
  
  # check if density y | x is valid
  #rowSums(fY) * diff(spt)[1] # should be 1
  #lines(spt, fY[2, ], type = "l", col = "red", lwd = 2, xlab = "y", ylab = "f(y | x)")
  
  Y     <- sapply(1:n, function(i) {
    sample(spt, 1, prob = fY[i, ])
  })
  return(data.frame(Y, theta = theta, X))
}


# --- Example1: Polynomial-tilt simulator via Bernstein basis (degree m) -----------------
softplus <- function(z) log1p(exp(z))

# Build Bernstein basis B_{m,k}(y) for y in [0,1]
bernstein_basis <- function(y, m) {
  # returns a matrix of size length(y) x (m+1), columns k=0..m
  k <- 0:m
  # choose(m,k) y^k (1-y)^(m-k)
  outer(y, k, function(yy, kk) dbinom(kk, size = m, prob = yy))
}

# c_k(x) = softplus(delta_{k,0} + delta_{k,1} x + delta_{k,2} x^2)
# delta is a (m+1) x 3 matrix of coefficients (for 1, x, x^2)
sim_data_polytilt <- function(n, m = 3,
                              delta = matrix(c(
                                # k=0    k=1    k=2    k=3   (rows)
                                # intercept, x, x^2  for each k (stack rows)
                                0.5,  0.5,  0.2,   # k=0
                                0.2,  0.8,  0.3,   # k=1
                                0.2, -0.2,  0.4,   # k=2
                                0.5,  0.0,  0.1    # k=3
                              ), nrow = m + 1, byrow = TRUE)) {
  
  stopifnot(ncol(delta) == 3, nrow(delta) == (m + 1))
  
  X       <- cbind(1, runif(n, -sqrt(12) / 4, sqrt(12) / 4))
  X[, -1] <- scale(X[, -1])
  xstd    <- X[, 2]
  
  # Bernstein basis on the support grid spt \in (0,1)
  B <- bernstein_basis(spt, m)  # length(spt) x (m+1)
  
  # Calculate c_k(x_i) for each i and k
  # For each k: c_k(x) = softplus(delta[k,1] + delta[k,2]*x + delta[k,3]*x^2)
  # We create an (n x (m+1)) matrix C where C[i,k+1] = c_k(x_i)
  C <- sapply(0:m, function(k) {
    z <- delta[k + 1, 1] + delta[k + 1, 2] * xstd + delta[k + 1, 3] * xstd^2
    softplus(z)
  })  # n x (m+1)
  
  # Weight at each grid point for each i: w_i(y_j) = sum_k C[i,k] * B_{m,k}(y_j)
  # Compute W as n x length(spt)
  W <- C %*% t(B)  # (n x (m+1)) %*% ((m+1) x length(spt)) = n x length(spt)
  
  # Get tilted densities: f_Y[i, j] \propto W[i, j] * f0_kde[j]
  fY_unnorm <- W * matrix(f0_kde, nrow = n, ncol = length(spt), byrow = TRUE)
  
  # Normalize rows to integrate to 1 on the grid
  dx      <- diff(spt)[1]
  Z_i     <- rowSums(fY_unnorm) * dx
  fY      <- fY_unnorm / Z_i
  
  # Sample Y on the grid using row-wise probabilities
  Y <- sapply(1:n, function(i) sample(spt, 1, prob = fY[i, ]))
  
  # True conditional mean 
  mu_true <- as.numeric(fY %*% (spt * dx))
  
  data.frame(Y = Y, mu_true = mu_true, X)
}

set.seed(2025)
data3 <- sim_data_polytilt(250, m = 3)
X <- as.matrix(data3[, c("X1", "X2")])

library(gldrm)
fit3 <- gldrm(Y ~ X2, data = data3, link = "logit")
truth3 <- list(beta = fit3$beta %>% as.numeric(),
               f0 = matrix(c(spt, f0_kde), 
                           ncol = 2, 
                           byrow = FALSE))

scale_factors <- c(1, 1) %>% sqrt()
m0 <- mean(data3$Y)

iter <- 500
fit3_dpglm <- fit_func(data3, iter, truth3, truth3$beta, 2.5, scale_factors, m0)
beta_est <- colMeans(fit3_dpglm$dpglm$beta[(iter/2):iter, ]) 

# Calculate E(y|x) using DPGLM 
meanY_x <- plogis(X %*% beta_est) 

# Check E(y|x) calibration
plot(meanY_x, data3$mu_true)
abline(0,1,col='red')
cor(data3$mu_true, meanY_x) 
mean((data3$mu_true - meanY_x)^2) # Empirical MSE on the dataset

# --- Example2: Quadratic-tilt simulator -----------------

# w(y|x) = softplus(a0 + a1*x + a2*x^2 + b1*y + b2*y^2)
sim_data_polytilt_quad <- function(n,
                                   a0 = 0.0, a1 = 1.0, a2 = 0.4,
                                   b1 = 1.0, b2 = -0.6) {
  X       <- cbind(1, runif(n, -sqrt(12) / 4, sqrt(12) / 4))
  X[, -1] <- scale(X[, -1])
  xstd    <- X[, 2]
  
  W <- outer(xstd, spt, function(xi, yj) {
    softplus(a0 + a1*xi + a2*xi^2 + b1*yj + b2*yj^2)
  })
  
  fY_unnorm <- W * matrix(f0_kde, nrow = n, ncol = length(spt), byrow = TRUE)
  dx        <- diff(spt)[1]
  Z_i       <- rowSums(fY_unnorm) * dx
  fY        <- fY_unnorm / Z_i
  
  Y       <- sapply(1:n, function(i) sample(spt, 1, prob = fY[i, ]))
  mu_true <- as.numeric(fY %*% (spt * dx))
  
  data.frame(Y = Y, mu_true = mu_true, X)
}

set.seed(2025)
data3q <- sim_data_polytilt_quad(250)

X <- as.matrix(data3q[, c("X1", "X2")])

library(gldrm)
fit3q <- gldrm(Y ~ X2, data = data3q, link = "logit")
truth3q <- list(beta = fit3q$beta %>% as.numeric(),
               f0 = matrix(c(spt, f0_kde), 
                           ncol = 2, 
                           byrow = FALSE))

scale_factors <- c(1, 1) %>% sqrt()
m0 <- mean(data3q$Y)

iter <- 500
fit3q_dpglm <- fit_func(data3q, iter, truth3q, truth3q$beta, 2.5, scale_factors, m0)
beta_est <- colMeans(fit3q_dpglm$dpglm$beta[(iter/2):iter, ]) 

# Calculate E(y|x) using DPGLM 
meanY_x <- plogis(X %*% beta_est) 

# Check E(y|x) calibration
plot(data3q$mu_true, meanY_x)
abline(0,1,col='red')
cor(data3q$mu_true, meanY_x) 
mean((data3q$mu_true - meanY_x)^2) # Empirical MSE on the dataset
