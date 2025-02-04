#source("load.R")

# Simulation III ---------------------------------------------------------------------
f0_true_kde <- readRDS("sim_truth/f0_true_kde.rds")
spt <- f0_true_kde$spt
f0_kde  <- f0_true_kde$f0 / sum(f0_true_kde$f0)

sim_III <- function(n) {
  X     <- cbind(1, runif(n, -sqrt(12) / 4, sqrt(12) / 4))
  X[, -1] <- scale(X[, -1])
  beta <- c(-0.7, 0.2)
  mu <-  exp(X %*% beta) / (1 + exp(X %*% beta))
  theta <- gldrm:::getTheta(
    spt = spt,
    f0  = f0_kde,
    mu  = mu,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = NULL
  )$theta
  
 # theta <- rmvnorm(1, mean = theta, sigma = 2 * diag(n)) %>% as.numeric()
  
  btheta <- b_theta(theta, spt, f0_kde)
  X       <- X %>% as.matrix()
  fY <- t(sapply(1:n, function(i) {
    exp(theta[i] * spt - btheta[i]) * f0_kde
  }))
  Y     <- sapply(1:n, function(i) {
    sample(spt, 1, prob = fY[i, ])
  })
  return(data.frame(Y, theta = theta, X))
}

true_beta <- c(0, 0)

truth <- list(beta = true_beta, 
              f0 = matrix(c(spt, f0_kde), 
                          ncol = 2, 
                          byrow = FALSE))

# Simulate Data ---------------------------------------------------------------------
dat <- sim_III(100)
y <- dat[, 1]
mu0 <- mean(y)
true_theta <- dat[, 2]
X <- dat[, -c(1,2)]

# Tuning Parameters --------------------------------------------------------------------
rho <- 1
M <- 20
alpha <- 1
delta <- 2
sigma_theta <- 0.001# sqrt(2)
c0 <- 0.25
beta.sigma <- 3 

# Data Preparation ---------------------------------------------------------------------
X <- X %>% as.matrix()
y <- y %>% as.numeric()
n <- length(y)
p <- dim(X)[2]

# Link Function -----------------------------------------------------------------------
link <- 'logit'

## Initialization -----------------------------------------------------------------------
iter <- 500
beta_samples <- matrix(NA, nrow = iter, ncol = p)
theta_samples <- matrix(NA, nrow = iter, ncol = n)       
z_samples <- matrix(NA, nrow = iter, ncol = n)
crm2_samples <- crm1_samples   <- list()
init <- gldrm(y ~ X[, -1], link = link)
beta_samples[1, ] <- beta <-  truth$beta %>% as.numeric() 
z.tld <- z_tld <-  truth$f0[, 1] %>% as.numeric()                  
J.tld <- J_tld <- truth$f0[, 2] %>% as.numeric()                  
theta_samples[1, ] <- theta <- true_theta   
z_samples[1, ] <- z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta)

btheta <- b_theta(theta, z.tld, J.tld)
T.vec <- exp(btheta)
u <- rgamma(n, shape = 1, rate = 1)
count <- 0

for(itr in 1:iter){
  log_post_wrapper <- function(beta) {
    return(log_post_beta(beta, z.tld, J.tld, X, theta, sigma_theta, beta.sigma))
  }
  
  # Optimization to find the mode
  result <- optim(par = beta, fn = log_post_wrapper, 
                  control = list(fnscale = -1), hessian = TRUE)
  beta_mode <- result$par 
  beta_cov <- - solve(result$hessian)
  beta <- rmvnorm(n = 1, mean = beta_mode, sigma = beta_cov) %>% as.numeric()
  meanY_x <- as.numeric(expit(X %*% beta))
  beta
  
  # theta update ------------------------------------
  theta_tilde <- gldrm:::getTheta(
    spt = z.tld,
    f0  = J.tld,
    mu  = meanY_x,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = theta
  )$theta
  
  theta_bar <- as.vector(mvrnormArma(1, 
                                     theta_tilde + z * sigma_theta^2, 
                                     diag(length(z)) * sigma_theta^2))
  
  btheta_bar <- b_theta(theta_bar, z.tld, J.tld)
  btheta <- b_theta(theta, z.tld, J.tld)
  log_theta_ratio <- btheta - btheta_bar
  for(i in 1:n){
    if(log(runif(1)) < log_theta_ratio[i]){
      theta[i] <- theta_bar[i]
    }
  }
  
  # u update ----------------------------------------
  u_old <- u
  u <- update_u(u, z, theta, alpha, delta)
  hist(u_old -u)
  
  # CRM update --------------------------------------
  crm_star <- crm_sampler_fractionalY(M, u, zstar, nstar, RL, RJ, theta, alpha, meanY_x, 
                                      y, shape_a, shape_b)
  z.tld_star <- c(crm_star$RL, crm_star$zstar)
  J.tld_star <- c(crm_star$RJ, crm_star$Jstar)
  
  # MH step
  theta_tilde_star <- gldrm:::getTheta(
    spt = z.tld_star,
    f0  = J.tld_star,
    mu  = meanY_x,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = theta
  )$theta
  
  # log acceptance probability
  log_r <- sum(dnorm(theta, theta_tilde_star, sigma_theta, log = TRUE)) -
    sum(dnorm(theta, theta_tilde, sigma_theta, log = TRUE))
  
  if(log(runif(1)) < log_r){
    count <- count + 1
    RL <- crm_star$RL
    RJ <- crm_star$RJ
    zstar <- crm_star$zstar
    Jstar <- crm_star$Jstar
    z.tld <- c(RL, zstar)
    J.tld <- c(RJ, zstar)
  }
  
  print(round(c(beta, sd(theta_tilde_star), count / itr), 3))
  
  # z update ------------------------------------------------------------------
  ord <- order(z.tld)
  sorted_crm <- data.frame(z.tld = z.tld[ord], J.tld = J.tld[ord])
  z.tld <- sorted_crm[, 1]
  J.tld <- sorted_crm[, 2]
  z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta) 
  
  # zstar and nstar update ----------------------------------------------------
  resampled_z <- resample_zstar(z)
  zstar <- resampled_z$zstar
  nstar <- resampled_z$nstar
  
  # Storing MCMC simulations --------------------------------------------------
  z_samples[itr, ] <- z
  beta_samples[itr,] <- beta
  theta_samples[itr, ] <- theta
  crm1_samples[[itr]] <- list(RL = RL, RJ = RJ)
  crm2_samples[[itr]] <- list(zstar = zstar, Jstar = Jstar)
}

#hist(c(z_samples))

#hist(beta_samples[51:100, 2])



