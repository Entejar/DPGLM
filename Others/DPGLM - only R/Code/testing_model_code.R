#source("load.R")

# Simulation III ---------------------------------------------------------------------
f0_true_kde <- readRDS("sim_truth/f0_true_kde.rds")
spt <- f0_true_kde$spt
f0_kde  <- f0_true_kde$f0 / sum(f0_true_kde$f0)
true_beta <- c(-0.7, 0.2)

true_beta <- c(1, 0)

sim_I <- function(n) {
  X     <- cbind(1, runif(n, -sqrt(12) / 4, sqrt(12) / 4))
  X[, -1] <- scale(X[, -1])
  beta <- true_beta
  mu <- as.numeric(expit(X %*% beta))
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
  # btheta <- theta <- rep(0, n)
  # fY <- t(sapply(1:n, function(i) {
  #   f0_kde
  # }))
  # Y     <- sapply(1:n, function(i) {
  #   sample(spt, 1, prob = fY[i, ])
  # })
  return(data.frame(Y, theta = theta, X))
}


sim_III <- function(n) {
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

truth <- list(beta = true_beta, 
              f0 = matrix(c(spt, f0_kde), 
                          ncol = 2, 
                          byrow = FALSE))

# Simulate Data ---------------------------------------------------------------------
dat <- sim_III(25)
dat <- sim_I(250)
y <- dat[, 1]
hist(y)
mu0 <- mean(y)
true_theta <- dat[, 2]
X <- dat[, -c(1,2)]

# Tuning Parameters --------------------------------------------------------------------
rho <- 1
M <- 20
alpha <- 1
delta <- 2
#sigma_theta <- 1.5
c0 <- 0.025
mu_beta <- c(0,0)
sigma_beta <- 1

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
init$beta
# init$spt
# init$f0
init$llik
# loglikelihood(y, X, init$beta, init$spt, init$f0)
# 
# loglik(z = y, X = X, beta = true_beta, atoms = truth$f0[, 1] %>% as.numeric() , jumps = truth$f0[, 2] %>% as.numeric())
# loglik(z = y, X = X, beta = as.numeric(init$beta), atoms = init$spt, jumps = init$f0)
beta_samples[1, ] <- beta <-  truth$beta %>% as.numeric()
z.tld <- z_tld <-  truth$f0[, 1] %>% as.numeric()
J.tld <- J_tld <- truth$f0[, 2] %>% as.numeric()

# beta_samples[1, ] <- beta <-  init$beta %>% as.numeric() 
# z.tld <- z_tld <-  init$spt %>% as.numeric()                  
# J.tld <- J_tld <- init$f0 %>% as.numeric()  

ord <- order(J_tld)[1:M]

RL <- z_tld[ord]  
RJ <- J_tld[ord]

theta_samples[1, ] <- theta <- true_theta   
z_samples[1, ] <- z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta)


btheta <- b_theta(theta, z.tld, J.tld)
T_vec <- exp(btheta)
#z <- y
u <- rgamma(n, shape = 1, rate = T_vec)
# u_mat <- matrix(NA, nrow = 100, ncol = n)
# u_mat[1, ] <- u
# for(itr in 2:100){
# u_mat[itr,] <- sampler_u(u, z, theta, alpha, delta)
# hist(u_mat[itr,] - u_mat[itr-1,])
# }
# 
# u <- colMeans(u_mat[76:100,])

resampled_z <- resample_zstar(z)
zstar <- resampled_z$zstar
nstar <- resampled_z$nstar
Jstar <- rgamma(n = length(nstar), shape = nstar, rate = 1)
count <- 0

out_mat <- matrix(NA, nrow = iter, ncol = 4)
colnames(out_mat) <- c("beta1", "beta2", "crm_change", "loglik")

for(itr in 1:iter){
  #itr <- itr + 1
  theta0 <- gldrm:::getTheta(
    spt = z.tld,
    f0  = J.tld / sum(J.tld),
    mu  = mu0,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = NULL
  )$theta
  Jtld_0 <- exp(theta0 * z.tld) * J.tld / sum(exp(theta0 * z.tld) * J.tld)
  lnlik <- loglik(z = z, X = X, beta = beta, atoms = z.tld, jumps = Jtld_0)
  #lnlik <- 0
  # log_post_wrapper <- function(beta) {
  #   return(log_post_beta(beta, z.tld, J.tld, X, theta, sigma_theta, sigma_beta))
  # }
  # 
  # Optimization to find the mode
  result <- optim(par = beta, fn = logpost_beta, z = z, X = X, 
                  atoms = z.tld, jumps  = J.tld, mu_beta = mu_beta,
                  sigma_beta = sigma_beta, 
                  control = list(fnscale = -1), hessian = TRUE)
  beta_mode <- result$par 
  beta_cov <- - solve(result$hessian)
  beta_ <- rmvnorm(n = 1, mean = beta_mode, sigma = beta_cov) %>% as.numeric()
  mean_z_ <- as.numeric(expit(X %*% beta))
  if(min(mean_z_) > min(z) && max(mean_z_) < max(z)){
    beta <- beta_
    meanY_x <- mean_z <- mean_z_
  }
  
  # theta update ------------------------------------
  theta_tilde <- gldrm:::getTheta(
    spt = z.tld,
    f0  = J.tld / sum(J.tld),
    mu  = meanY_x,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = theta
  )$theta
  
  theta <- theta_tilde
  
  # theta_bar <- as.vector(mvrnormArma(1, 
  #                                    theta_tilde + z * sigma_theta^2, 
  #                                    diag(length(z)) * sigma_theta^2))
  # 
  # btheta_bar <- b_theta(theta_bar, z.tld, J.tld)
  # btheta <- b_theta(theta, z.tld, J.tld)
  # log_theta_ratio <- btheta - btheta_bar
  # for(i in 1:n){
  #   if(log(runif(1)) < log_theta_ratio[i]){
  #     theta[i] <- theta_bar[i]
  #   }
  # }
  
  # u update ----------------------------------------
  u_old <- u
  #u <- update_u(u, z, theta, alpha, delta)
  u <- sampler_u(u, z, theta, alpha, delta)
  #u <- sampler_u_(u, z, theta, alpha, delta) 
  hist(u_old -u)
  
  # CRM update --------------------------------------
  crm_star <- crm_sampler_fractionalY(M, u, zstar, nstar, RL, RJ, theta, alpha, meanY_x, 
                                      y, 1, 1)
  z.tld_star <- c(crm_star$RL, crm_star$zstar)
  J.tld_star <- c(crm_star$RJ, crm_star$Jstar)
  
  #if(min(meanY_x) > min(z.tld_star) && max(meanY_x) < max(z.tld_star)){
  
  # MH step
  theta_tilde_star <- gldrm:::getTheta(
    spt = z.tld_star,
    f0  = J.tld_star / sum(J.tld_star),
    mu  = meanY_x,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = theta
  )$theta
  
  # log acceptance probability
  # log_r <- sum(dnorm(theta, theta_tilde_star, sigma_theta, log = TRUE)) -
  #   sum(dnorm(theta, theta_tilde, sigma_theta, log = TRUE))
  
  b1 <- b_theta(theta_tilde_star, z.tld_star, J.tld_star)
  b2 <- b_theta(theta_tilde, z.tld, J.tld)
  b3 <- b_theta(theta_tilde_star, z.tld, J.tld)
  b4 <- b_theta(theta_tilde, z.tld_star, J.tld_star)
  log_r <- log(exp(sum(2*(theta_tilde_star - theta_tilde)*z - b1 + b2 - b3 + b4)))
  
  if(log(runif(1)) < log_r){
    count <- count + 1
    RL <- crm_star$RL
    RJ <- crm_star$RJ
    zstar <- crm_star$zstar
    Jstar <- crm_star$Jstar
    z.tld <- c(RL, zstar)
    J.tld <- c(RJ, Jstar)
    theta <- theta_tilde_star
  }
  #}
  
  print(round(c(beta, count / itr, lnlik), 3))
  out_mat[itr, ] <- round(c(beta, count / itr, lnlik), 3)
  
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
  
  plot(out_mat[1:itr, 1], type = "l", xlab = "Iteration", ylab = "Beta1")
  abline(h = true_beta[1], col = "red")
  plot(out_mat[1:itr, 2], type = "l", xlab = "Iteration", ylab = "Beta2")
  abline(h = true_beta[2], col = "red")
  plot(out_mat[1:itr, 3], type = "l", xlab = "Iteration", ylab = "CRM Acceptance Rate")
  #plot(out_mat[1:itr, 4], type = "l", xlab = "Iteration", ylab = "Log Likelihood")
}

#out_mat
# n_iter <- 101
# par(mfrow = c(2, 2))
# plot(out_mat[1:n_iter, 1], type = "l", xlab = "Iteration", ylab = "Beta1")
# abline(h = true_beta[1], col = "red")
# plot(out_mat[1:n_iter, 2], type = "l", xlab = "Iteration", ylab = "Beta2")
# abline(h = true_beta[2], col = "red")
# plot(out_mat[1:n_iter, 3], type = "l", xlab = "Iteration", ylab = "CRM Acceptance Rate")
# plot(out_mat[1:n_iter, 4], type = "l", xlab = "Iteration", ylab = "Log Likelihood")

#hist(c(z_samples))

#hist(beta_samples[51:100, 2])

# n <- 100
# y <- rbeta(n, 1, 1)
# M <- 10
# alpha <- 1
# theta <- rep(0, n)
# temp <- resample_zstar(y)
# zstar <- temp$zstar
# nstar <- temp$nstar
# delta <- 2
# c0 <- 0.025
# u <- rgamma(n, shape = 1, rate = 1)
# z <- y
# 
# crm <- update_crm(M, u, zstar, nstar, theta, alpha) 
# z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta) 
# u <- u_sampler(u, z, theta, alpha, delta)
# hist(u)
# 

# hist(y)
# r <- 450
# r <- r +1 
# hist(rep(crm1_samples[[r]]$RL, ceiling(crm1_samples[[r]]$RJ / min(crm1_samples[[r]]$RJ))))

data <- dat[, -2]

fit <- brm(
  bf(Y ~ X2, phi ~ X2),
  data = data,             
  family = Beta(),           
  chains = 1, iter = 2000, warmup = 1000
)
fit
