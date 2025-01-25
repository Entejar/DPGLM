#source("load.R")

# Simulation III ---------------------------------------------------------------------
f0_true_kde <- readRDS("sim_truth/f0_true_kde.rds")
spt <- f0_true_kde$spt
f0_kde  <- f0_true_kde$f0 / sum(f0_true_kde$f0)
#mu0 <- sum(spt * f0_kde) / sum(f0_kde)

sim_III <- function(n) {
  X     <- cbind(1, runif(n, -sqrt(12) / 4, sqrt(12) / 4))
  X[, -1] <- scale(X[, -1])
  beta <- c(-0.7, 0.2)
  mu <-  exp(X %*% beta) / (1 + exp(X %*% beta))
  # theta <- solve_thetas(mu, spt, f0_kde)
  theta <- gldrm:::getTheta(
    spt = spt,
    f0  = f0_kde,
    mu  = mu,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = NULL
  )$theta
  
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

true_beta <- c(-0.7, 0.2)

truth <- list(beta = true_beta, 
              f0 = matrix(c(spt, f0_kde), 
                          ncol = 2, 
                          byrow = FALSE))

# Simulate Data ---------------------------------------------------------------------
dat <- sim_III(250)
y <- dat[, 1]
mu0 <- mean(y)
true_theta <- dat[, 2]
#hist(y)
X <- dat[, -c(1,2)]

# Tuning Parameters --------------------------------------------------------------------
rho <- 1
M <- 20
alpha <- 1
delta <- 2
sigma_theta <- 1.5 #sd(true_theta) / 2
c0 <- 0.25
beta.sigma <- 2.5 #max(abs(true_beta)) / 2

# Data Preparation ---------------------------------------------------------------------
X <- X %>% as.matrix()
y <- y %>% as.numeric()
mu_y <- mu0 <- mean(y)
sd2_y <- var(y)
temp <- mu_y * (1 - mu_y) / sd2_y - 1
shape_a <- mu_y * temp
shape_b <- (1 - mu_y) * temp
n <- length(y)
p <- dim(X)[2]

# Link Function -----------------------------------------------------------------------
link <- 'logit'

## Initialization -----------------------------------------------------------------------
iter <- 100
beta_samples <- matrix(NA, nrow = iter, ncol = p)
theta_samples <- matrix(NA, nrow = iter, ncol = n)       
z_samples <- matrix(NA, nrow = iter, ncol = n)
crm2_samples <- crm1_samples   <- list()
init <- gldrm(y ~ X[, -1], link = link)
beta_samples[1, ] <- beta <-  truth$beta %>% as.numeric() #init$beta %>% as.numeric()  
z_tld <-  truth$f0[, 1] %>% as.numeric()                  #init$spt %>% as.numeric() 
J_tld <- truth$f0[, 2] %>% as.numeric()                   #init$f0 %>% as.numeric() 

RL <- runif(M, min(z_tld), max(z_tld))
RJ <- stick_breaking_init(M, alpha)

z.tld <- c(RL, z_tld)
J.tld <- c(RJ, J_tld)
crm1_samples[[1]] <- list(RL = RL, RJ = RJ)
theta_samples[1, ] <- theta <- true_theta                 #init$theta                    
btheta <- b_theta(theta, z.tld, J.tld)
z_samples[1, ] <- z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta)
temp <- resample_zstar(z)
zstar <- temp$zstar
nstar <- temp$nstar
Jstar <- rgamma(nstar, 1)
crm2_samples[[1]] <- list(zstar = zstar, Jstar = Jstar)
z.tld <- c(RL, zstar)
J.tld <- c(RJ, Jstar)

T.vec <- exp(btheta)
u <- rgamma(n, shape = 1, rate = T.vec)
mubetaprior <- rep(0, p)                     #true_beta 
Sigbetaprior <- beta.sigma * diag(p)
Sig <- rho * init$varbeta
#itr <- 1
count <- 0

for(itr in 2:iter){
  #itr <- itr + 1
  # beta update -------------------------------------
  # beta <- beta_sampler_fractionalY(beta, Sig, theta, z.tld, J.tld, 
  #                                  mubetaprior, Sigbetaprior, sigma_theta) 
  # Wrapper for log_post_beta function in R
  log_post_wrapper <- function(beta) {
    # Call the C++ function, which returns the log posterior value for the given beta
    return(log_post_beta(beta, z.tld, J.tld, X, theta, sigma_theta, beta.sigma))
  }
  
  # Optimization to find the mode
  result <- optim(par = beta, fn = log_post_wrapper, 
                  control = list(fnscale = -1), hessian = TRUE)
  beta_mode <- result$par 
  beta_cov <- - solve(result$hessian)
  beta_tilde <- rmvnorm(n = 1, mean = beta_mode, sigma = beta_cov) %>% as.numeric()
  log_post_current <- log_post_wrapper(beta)
  log_post_proposed <- log_post_wrapper(beta_tilde)
  acceptance_ratio <- exp(log_post_proposed - log_post_current)
  if (runif(1) < acceptance_ratio) {
    beta <- beta_tilde
  }
  meanY_x <- as.numeric(expit(X %*% beta))
  
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
  #print(sum(log_theta_ratio))
  for(i in 1:n){
    if(log(runif(1)) < log_theta_ratio[i]){
      theta[i] <- theta_bar[i]
    }
  }
  # if(log(runif(1)) < sum(log_theta_ratio)){
  #   theta <- theta_bar
  # }
  
  #hist(theta - true_theta)
  
  
  # u update ----------------------------------------
  u_old <- u
  u <- u_sampler(u, z, theta, alpha, delta)
  #hist(u - u_old)
  #plot(u, u_old, pch = 19)
  #abline(a = 0, b = 1)
  
  # CRM update --------------------------------------
  crm_star <- crm_sampler_fractionalY(M, u, zstar, nstar, RL, RJ, theta, alpha, meanY_x, 
                                      y, shape_a, shape_b)
  z.tld_star <- c(crm_star$RL, crm_star$zstar)
  J.tld_star <- c(crm_star$RJ, crm_star$Jstar)
  
  #hist(z.tld_star, breaks = 50)
  #hist(J.tld_star, breaks = 50)
  
  # MH step
  theta_tilde_star <- gldrm:::getTheta(
    spt = z.tld_star,
    f0  = J.tld_star,
    mu  = meanY_x,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = theta
  )$theta
  
  
  #plot(theta_tilde_star, theta, type = 'l')
  #abline(a = 0, b = 1)
  # 
  # theta <- theta_tilde_star1
  # 
  # crm_star <- crm_sampler_fractionalY(M, u, zstar, nstar, RL, RJ, theta, alpha, meanY_x, 
  #                                     y, shape_a, shape_b)
  # z.tld_star <- c(crm_star$RL, crm_star$zstar)
  # J.tld_star <- c(crm_star$RJ, crm_star$Jstar)
  # 
  # #hist(z.tld_star, breaks = 50)
  # #hist(J.tld_star, breaks = 50)
  # 
  # # MH step
  # theta_tilde_star2 <- gldrm:::getTheta(
  #   spt = z.tld_star,
  #   f0  = J.tld_star,
  #   mu  = meanY_x,
  #   sampprobs  = NULL,
  #   ySptIndex  = NULL,
  #   thetaStart = theta
  # )$theta
  # 
  # #hist(theta_tilde - theta_tilde_star)
  #plot(theta_tilde_star, theta_tilde, type = 'l')
  #abline(a = 0, b = 1)
  # theta_tilde_star <- theta_tilde_star2
  # hist(theta_tilde)
  # mean(abs(theta_tilde_star - theta_tilde))
  # 
  #theta_tilde_star <- theta_tilde_star - mean((theta_tilde_star - theta_tilde))
  
  # log acceptance probability
  log_r <- sum(dnorm(theta, theta_tilde_star, sigma_theta, log = TRUE)) -
    sum(dnorm(theta, theta_tilde, sigma_theta, log = TRUE))
  
  # IS approximation
  #log_r <- sum(dnorm(theta, theta_tilde_star, sigma_theta, log = TRUE)) 
  
  if(log(runif(1)) < log_r){
    count <- count + 1
    RL <- crm_star$RL
    RJ <- crm_star$RJ
    zstar <- crm_star$zstar
    Jstar <- crm_star$Jstar
    z.tld <- c(RL, zstar)
    J.tld <- c(RJ, zstar)
    #sigma_theta <- 2 * sd(theta_tilde_star) 
  }
  
  
  print(round(c(acceptance_ratio, beta, sd(theta_tilde_star), log_r, count / itr), 3))
  
  # z update ------------------------------------------------------------------
  ord <- order(z.tld)
  sorted_crm <- data.frame(z.tld = z.tld[ord], J.tld = J.tld[ord])
  z.tld <- sorted_crm[, 1]
  J.tld <- sorted_crm[, 2]
  
  #hist(z.tld, breaks = 50)
  #hist(J.tld, breaks = 50)
  
  z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta) 
  #hist(z, breaks = 50)
  
  # zstar and nstar update ----------------------------------------------------
  resampled_z <- resample_zstar(z)
  zstar <- resampled_z$zstar
  nstar <- resampled_z$nstar
  
  #range(z.tld)
  #length(J.tld)
  
  # theta0 <- gldrm:::getTheta(
  #   spt = z.tld,
  #   f0  = J.tld,
  #   mu  = rep(mu0, length(z.tld)),
  #   sampprobs  = NULL,
  #   ySptIndex  = NULL,
  #   thetaStart = NULL
  # )$theta
  # 
  # J.tld <- exp(theta0) * J.tld / sum(exp(theta0) * J.tld)
  # 
  # Storing MCMC simulations --------------------------------------------------
  z_samples[itr, ] <- z
  beta_samples[itr,] <- beta
  theta_samples[itr, ] <- theta
  crm1_samples[[itr]] <- list(RL = RL, RJ = RJ)
  crm2_samples[[itr]] <- list(zstar = zstar, Jstar = Jstar)
}

#hist(c(z_samples))

#hist(beta_samples[51:100, 2])





