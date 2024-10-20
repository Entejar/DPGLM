
# Loading -----------

source("load.R")

# # Required Packages ----
# require(foreach)
# require(doParallel)
# 
# 
# # Parallel Setup ----
# num_cores <- 6
# cl        <- makeCluster(num_cores)
# registerDoParallel(cl)

# Load Data --------------------------------------------------------------------
hustadTD <- readRDS("/Users/entejar/Library/CloudStorage/Box-Box/DP-GLM/Data/hustadTD.rds")

# Subset TD Data: Multi Word, Ordered y --------------------------------------------------
hustadTDMW <- filter(hustadTD %>% as.data.frame(),
                     intelligibility_type == 'multiword') %>% arrange(mean_intelligibility)

# KDE
dens  <- density(
  hustadTDMW$mean_intelligibility,
  kernel = 'gaussian',
  bw = 'nrd0',
  adjust = 1,
  n = 1000,
  cut = 3,
  ext = 4
)
dat   <- data.frame(y = dens$x, f0 = dens$y)
indx  <- dat$y >= 0 & dat$y <= 1
true_spt <- dat$y[indx]
true_f0    <- dat$f0[indx]
true_mu0 <- sum(true_spt * true_f0) / sum(true_f0)

# Beta(a, b) Fit
betafit <- MASS::fitdistr(hustadTDMW$mean_intelligibility, dbeta, start = list(shape1 = 1, shape2 = 1), lower = c(0, 0))
a <- betafit$estimate[1]
b <- betafit$estimate[2]
f00 <- dbeta(true_spt, shape1 = a, shape2 = b)

# Simulation: Setting I

sim_I <- function(n) {
  X     <- cbind(1, matrix(c(rnorm(n, 1, 0.5), rnorm(n, 2, 1)), ncol = 2))
  X       <- X %>% as.matrix()
  btheta <- theta <- rep(0, n)
  fY <- t(sapply(1:n, function(i) {
    true_f0
  }))
  Y     <- sapply(1:n, function(i) {
    sample(true_spt, 1, prob = fY[i, ])
  })
  return(data.frame(Y, X) %>% arrange(Y))
}

sim_II <- function(n, p0, p1) {
  X     <- cbind(1, matrix(c(rnorm(n, 1, 0.5), rnorm(n, 2, 1)), ncol = 2))
  X       <- X %>% as.matrix()
  btheta <- theta <- rep(0, n)
  fY <- t(sapply(1:n, function(i) {
    f00
  }))
  Y     <- sapply(1:n, function(i) {
    u <- runif(1)
    if(u < p0){
      0
    } else if (u < p0 + p1){
      1
    } else {
      sample(true_spt, 1, prob = f00)
    }
  })
  return(data.frame(Y, X) %>% arrange(Y))
}

sim_III <- function(n, link) {
  X     <- cbind(1, matrix(c(rnorm(n, 1, 1), rnorm(n, 2, 1)), ncol = 2))
  X       <- X %>% as.matrix()
  beta <- c(-0.7, 0.2, -0.1)
  mu <- link_fun(X %*% beta, link)$mu
  out <- theta_sol(true_spt, true_f0, mu, NULL)
  theta <- out$theta
  btheta <- out$btht
  fY <- t(sapply(1:n, function(i) {
    exp(theta[i] * true_spt - btheta[i]) * true_f0
  }))
  Y     <- sapply(1:n, function(i) {
    sample(true_spt, 1, prob = fY[i, ])
  })
  return(data.frame(Y, X) %>% arrange(Y)) # ATTENTION: data frame is sorted by Y
}

# Unit Testing
n <- 400
dat <- sim_III(n, link = 'logit')

rho  <- 1
M <- 20
alpha <- 1
G0.dist <- 6
delta <- 2
kdist <- 6                                # ATTENTION: Choose K as 6 or 7?
sigma_theta <- 0.001
a00 <- 0
b00 <- 1
c0 <- 0.025
beta.sigma <- 1

tuning.params <- list(rho = rho,
                      M = M,
                      alpha = alpha,
                      G0.dist = G0.dist,
                      delta = delta,
                      k.dist = kdist,
                      sigma_theta = sigma_theta,
                      a00 = a00,
                      b00 = b00,
                      c0 = c0,
                      beta.sigma = beta.sigma)

y <- dat[, 1] %>% as.numeric()
X <- dat[, -1] %>% as.matrix()

gldrm_fit <- gldrm(y ~ X[, -1], link = 'logit')
gldrm_fit
set.beta <- TRUE
out100 <- dpglm(
  y = y[1:100],
  X = X[1:100, ],
  iter = 10,
  tuning.params = tuning.params,
  set.beta = set.beta
)
