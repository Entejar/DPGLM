
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
betafit <- MASS::fitdistr(hustadTDMW$mean_intelligibility, dbeta, start = list(shape1 = 1, shape2 = 1))
a <- betafit$estimate[1]
b <- betafit$estimate[2]
f00 <- dbeta(true_spt, shape1 = a, shape2 = b)

# Simulation: Setting I

sim_I <- function(n) {
  X     <- cbind(1, matrix(c(rnorm(n, 0, 1), rnorm(n, 1, 1)), ncol = 2))
  X[, -1] <- scale(X[, -1])
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
  X     <- cbind(1, matrix(c(runif(n, 1, 3), runif(n, -1, 1)), ncol = 2))
  X[, -1] <- scale(X[, -1])
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
  X     <- cbind(1, matrix(c(runif(n, 1, 3), runif(n, -1, 1)), ncol = 2))
  X[, -1] <- scale(X[, -1])
  X       <- X %>% as.matrix()
  beta <- c(0.1, 0.25, 0.75)
  mu <- link_fn(X %*% beta, link)$mu
  out <- sapply(1:n, function(i){
    theta_solver(true_spt, true_f0, mu[i], NULL)
  })
  theta <- out[2, ] %>% unlist()
  btheta <- out[3, ] %>% unlist()
  fY <- t(sapply(1:n, function(i) {
    exp(theta[i] * true_spt - btheta[i]) * true_f0
  }))
  Y     <- sapply(1:n, function(i) {
    sample(true_spt, 1, prob = fY[i, ])
  })
  return(data.frame(Y, X) %>% arrange(Y))
}

# Unit Testing
n <- 500
dat <- sim_III(n, link = 'logit')

# # For beta_0 initialization
# Ybar <- mean(dat$Y)
# beta0Init <- log(Ybar / (1 - Ybar))    # a good initialization for beta_0 = log(true_mu0 / (1 - true_mu0))
# 
# Y <- dat[, 1] - Ybar
# X <- matrix(1, n, 2)
# data <- data.frame(Y, X)
# lm(Y ~ X, data = data.frame(Y, X))

rho  <- 1
M <- 20
alpha <- 1
G0.dist <- 6
delta <- 2
kdist <- 6                             # ATTENTION: Choose K as 6 or 7?
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

y <- dat[, 1]
X <- dat[, -1]
out100 <- dpglm(y = y[1:100], 
                X = X[1:100,], 
                iter = 100, 
                tuning.params = tuning.params)

out500 <- dpglm(y, X, 'logit', 500)
out50 <- dpglm(y[1:200], X[1:200,], 'logit', 500)


out100_total <- out100
out500_total <- out500
out50_total <- out50          # until this point, it is working fine

burn <- 10
out100$beta <- out100$beta[-c(1:burn), ]
out500$beta <- out500$beta[-c(1:burn), ]
out50$beta <- out50$beta[-c(1:burn), ]
true_beta <- c(1.28, 0, 0)
plot(abs(out50$beta[, 1] - true_beta[1]), type = 'l')
plot(abs(out100$beta[, 1] - true_beta[1]), type = 'l')
plot(abs(out500$beta[, 1] - true_beta[1]), type = 'l')

plot(abs(out100$beta[, 2] - true_beta[2]), type = 'l')
plot(abs(out500$beta[, 2] - true_beta[2]), type = 'l')
plot(abs(out50$beta[, 2] - true_beta[2]), type = 'l')

plot(abs(out100$beta[, 3] - true_beta[3]), type = 'l')
plot(abs(out500$beta[, 3] - true_beta[3]), type = 'l')
plot(abs(out50$beta[, 3] - true_beta[3]), type = 'l')

r <- 100
plot(out25$crm[[r]]$z.tld, out25$crm[[r]]$J.tld, type = 'l')
plot(out50$crm[[r]]$z.tld, out50$crm[[r]]$J.tld, type = 'l')
plot(out100$crm[[r]]$z.tld, out100$crm[[r]]$J.tld, type = 'l')
plot(out500$crm[[r]]$z.tld, out500$crm[[r]]$J.tld, type = 'l')

# out100$z <- out100$z[-c(1:burn), ]
# out500$z <- out500$z[-c(1:burn), ]

# abs(colMeans(out100$z) - y[1:100]) %>% mean()
# abs(colMeans(out500$z) - y) %>% mean()

out50$theta <- out50$theta[-c(1:burn), ]
out100$theta <- out100$theta[-c(1:burn), ]
out500$theta <- out500$theta[-c(1:burn), ]

hist(out50$theta)
hist(out100$theta)
hist(out500$theta)

theta_est_50 <- mean(colMeans(abs(out50$theta)))
theta_est_100 <- mean(colMeans(abs(out100$theta)))
theta_est_500 <- mean(colMeans(abs(out500$theta)))

c(theta_est_50, theta_est_100, theta_est_500)

theta_sd_50 <- mean(apply(out50$theta, 2, sd))
theta_sd_100 <- mean(apply(out100$theta, 2, sd))
theta_sd_500 <- mean(apply(out500$theta, 2, sd))

c(theta_sd_50, theta_sd_100, theta_sd_500)

plot_grid <- seq(1, length(true_spt), 1)
yGrid <- true_spt[plot_grid]

mu0 <- sum(true_spt * true_f0) / sum(true_f0)

theta <- 0
true_pro <- exp(theta * true_spt) * true_f0
true_pk  <- true_pro / sum(true_pro)
F0_true <- f0_true <- numeric(length(yGrid))
for (i in 1:length(yGrid)) {
  indx <- which(true_spt >= (yGrid[i] - c0) & true_spt <= (yGrid[i] + c0))
  f0_true[i] <- sum(true_pk[indx]) * (1 / (2 * c0))
  indx2 <- which(true_spt <= yGrid[i] - c0)
  F0_true[i] <- sum(true_pk[indx2]) +
    sum(true_pk[indx] * (yGrid[i] - true_spt[indx] + c0)/ (2 * c0))
}

c0 <- 0.025

f_est_fn <- function(out) {
    itr_indx <- 75:100
    F0_est <- f0_est <- matrix(NA, nrow = length(yGrid), ncol = length(itr_indx))
    for (j in 1:length(itr_indx)) {
      itr <- itr_indx[itr_indx[j]]
      spt <- out$crm[[itr_indx[j]]]$z.tld
      f0 <- out$crm[[itr_indx[j]]]$J.tld
      tht <- theta_solver(spt, f0, mu0, NULL)$tht
      f0 <- exp(tht * spt) * f0 / sum(exp(tht * spt) * f0)
      theta <- 0
      for (i in 1:length(yGrid)) {
        pro <- exp(theta * spt) * f0  # which is same as f0 for the null case
        pk  <- pro / sum(pro)
        indx <- which(spt >= (yGrid[i] - c0) & spt <= (yGrid[i] + c0))
        f0_est[i, j] <- sum(pk[indx]) * (1 / (2 * c0))
        indx2 <- which(spt <= yGrid[i] - c0)
        F0_est[i, j] <- sum(pk[indx2]) +
          sum(pk[indx] * (yGrid[i] - spt[indx] + c0)/ (2 * c0))
      }

    }
    return(list(F0_est = F0_est, f0_est = f0_est))
  }


baseDensity_100 <- f_est_fn(out100)
baseDensity_500 <- f_est_fn(out500)
baseDensity_50 <- f_est_fn(out50)


# Comparing f_0 and F_0 with truths

plot(yGrid, F0_true, type = 'l', col = 'black', lwd = 1, ylim = c(0, 1))
lines(yGrid, rowMeans(baseDensity_50$F0_est), col = 'red', lwd = 1)
lines(yGrid, rowMeans(baseDensity_100$F0_est), col = 'red', lwd = 1)
lines(yGrid, rowMeans(baseDensity_500$F0_est), col = 'green', lwd = 1)

# Total Variation distance
total_variation_distance <- function(f1, f2) {
  differences <- abs(f1 - f2)
  return(sum(differences) / 2)
}

# Comparing F_0 with truth

TV_50 <- total_variation_distance(rowMeans(baseDensity_50$F0_est) %>% as.vector(), F0_true)
TV_100 <- total_variation_distance(rowMeans(baseDensity_100$F0_est)%>% as.vector(), F0_true)
TV_500 <- total_variation_distance(rowMeans(baseDensity_500$F0_est) %>% as.vector(), F0_true)

c(TV_50, TV_100, TV_500)

# # Comparing f_0 with truth
# TV_50_f0 <- total_variation_distance(rowMeans(baseDensity_25$f0_est) %>% as.vector(), f0_true)
# TV_100_f0 <- total_variation_distance(rowMeans(baseDensity_100$f0_est) %>% as.vector(), f0_true)
# TV_500_f0 <- total_variation_distance(rowMeans(baseDensity_500$f0_est) %>% as.vector(), f0_true)
#
# c(TV_50_f0, TV_100_f0, TV_500_f0)

# # Parallel Loop
# num_datasets <- 1
# p0 <- 0.1
# p1 <- 0.4
# sim_study <- list()
# nVec <- c(50, 100, 500)
# scenarios <- 1:3
# set <- expand.grid(nVec, scenarios)
# R <- nrow(set)
# seed_init <-  sample.int(.Machine$integer.max, 1)
#
# start_time <- Sys.time()
#
# sim_study <- foreach(dat_index = 1:num_datasets,
#                      .packages = c("gldrm", "tidyverse", "truncnorm", "mvtnorm", "ggplot2")) %:%
#   foreach(setting = 1:R) %dopar% {
#     seed_val <- seed_init + dat_index + setting
#     set.seed(seed_val)
#     if (setting %in% 1:3) {dat <- sim_I(set[setting, 1])}
#     if (setting %in% 4:6) {dat <- sim_II(set[setting, 1], p0, p1)}
#     if (setting %in% 7:9) {dat <- sim_III(set[setting, 1])}
#     X <- dat[, -1]
#     y <- dat[, 1]
#     mcmc <- dpglm(y, X, 'logit', 100)
#     out <- list(seed = seed_val,
#                 set = setting,
#                 data = dat,
#                 mcmc = mcmc
#     )
#     return(out)
#   }
#
#
# # # Save result
# # saveRDS(sim_study, file = "0707_dpglm_simulation.rds")
#
# Sys.time() - start_time
#
# #
# # start_time <- Sys.time()
# #
# # sim_study <- foreach(dat_index = 1:num_datasets,
# #                      .packages = c("gldrm", "tidyverse", "truncnorm", "mvtnorm", "ggplot2")) %:%
# #   foreach(setting = 1:R) %dopar% {
# #     seed_val <- seed_init + dat_index + setting
# #     set.seed(seed_val)
# #     if (setting %in% 1:3) {dat <- sim_I(set[setting, 1])}
# #     if (setting %in% 4:6) {dat <- sim_II(set[setting, 1], p0, p1)}
# #     if (setting %in% 7:9) {dat <- sim_III(set[setting, 1])}
# #     X <- dat[, -1]
# #     y <- dat[, 1]
# #     mcmc <- dpglm(y, X, 'logit', 100)
# #     out <- list(seed = seed_val,
# #                 set = setting,
# #                 data = dat,
# #                 mcmc = mcmc
# #     )
# #     return(out)
# #   }
# #
# #
# # # # Save result
# # # saveRDS(sim_study, file = "0707_dpglm_simulation.rds")
# #
# # Sys.time() - start_time
#
# # Stop parallel backend
# stopCluster(cl)
