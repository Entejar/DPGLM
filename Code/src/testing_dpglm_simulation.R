`
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
f0_true_kde <- data.frame(spt = true_spt, f0 = true_f0)
saveRDS(f0_true_kde, "data/f0_true_kde.rds")
true_mu0 <- sum(true_spt * true_f0) / sum(true_f0)

# Beta(a, b) Fit
betafit <- MASS::fitdistr(hustadTDMW$mean_intelligibility, dbeta, start = list(shape1 = 1, shape2 = 1), lower = c(0, 0))
a <- betafit$estimate[1]
b <- betafit$estimate[2]
f00 <- dbeta(true_spt, shape1 = a, shape2 = b)
f0_true_beta <- data.frame(spt = true_spt, f0 = f00)
saveRDS(f0_true_beta, "data/f0_true_beta.rds")


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
delta <- 2                               # ATTENTION: Choose K as 6 or 7?
sigma_theta <- 0.001
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


out200 <- dpglm(
  y = y[1:200],
  X = X[1:200, ],
  iter = 1000,
  tuning.params = tuning.params,
  set.beta = set.beta
)

out400 <- dpglm(
  y = y,
  X = X,
  iter = 1000,
  tuning.params = tuning.params,
  set.beta = set.beta
)

burn <- 500
save <- seq(burn+1, 1000, 10)

sqrt(mean(out100$z[save, ] - y[1:100])^2)
sqrt(mean(out200$z[save, ] - y[1:200])^2)
sqrt(mean(out400$z[save, ] - y[1:400])^2)

true_beta <- c(-0.7, 0.2, -0.1)

beta_summary <- function(out) {
  beta_est <- out$beta[-(1:burn), ]
  beta_est_mean <- colMeans(beta_est)
  beta_est_sd <- apply(beta_est, 2, sd)
  beta_est_err <- abs(beta_est_mean - true_beta)
  return(list(beta_est_mean = beta_est_mean, beta_est_sd = beta_est_sd, 
              beta_est_err = beta_est_err))
}

beta_summary(out100)
beta_summary(out200)
beta_summary(out400)

plot(abs(out100$beta[-(1:burn), 1] - true_beta[1]), type = 'l')
plot(abs(out200$beta[-(1:burn), 1] - true_beta[1]), type = 'l')
plot(abs(out400$beta[-(1:burn), 1] - true_beta[1]), type = 'l')

plot(abs(out100$beta[-(1:burn), 2] - true_beta[2]), type = 'l')
plot(abs(out200$beta[-(1:burn), 2] - true_beta[2]), type = 'l')
plot(abs(out400$beta[-(1:burn), 2] - true_beta[2]), type = 'l')

plot(abs(out100$beta[-(1:burn), 3] - true_beta[3]), type = 'l')
plot(abs(out200$beta[-(1:burn), 3] - true_beta[3]), type = 'l')
plot(abs(out400$beta[-(1:burn), 3] - true_beta[3]), type = 'l')

## Relevant for testing null case i.e., theta = 0 or constant (theta0) over x

# out100$z <- out100$z[-c(1:burn), ]
# out500$z <- out500$z[-c(1:burn), ]

# abs(colMeans(out100$z) - y[1:100]) %>% mean()
# abs(colMeans(out500$z) - y) %>% mean()

# out100$theta <- out100$theta[-c(1:burn), ]
# out500$theta <- out500$theta[-c(1:burn), ]
# 
# hist(out100$theta)
# hist(out500$theta)
# 
# theta_est_100 <- mean(colMeans(abs(out100$theta)))
# theta_est_500 <- mean(colMeans(abs(out500$theta)))
# 
# c(theta_est_100, theta_est_500)
# 
# theta_sd_100 <- mean(apply(out100$theta, 2, sd))
# theta_sd_500 <- mean(apply(out500$theta, 2, sd))
# 
# c(theta_sd_100, theta_sd_500)

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
  itr_indx <- save
  F0_est <- f0_est <- matrix(0, nrow = length(yGrid), ncol = length(itr_indx))
  for (j in 1:length(itr_indx)) {
    spt <- out$crm[[itr_indx[j]]]$z.tld
    f0 <- out$crm[[itr_indx[j]]]$J.tld
    theta <- theta_sol(spt, f0, mu0, NULL)$theta
    f0 <- exp(theta * spt) * f0 / sum(exp(theta * spt) * f0)
    theta <- 0                       # For null case or f0 / F0
    for (i in 1:length(yGrid)) {
      pro <- exp(theta * spt) * f0   # Which is same as f0 for the null case
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
baseDensity_200 <- f_est_fn(out200)
baseDensity_400 <- f_est_fn(out400)

# Comparing F_0 with truths
plot(yGrid, F0_true, type = 'l', col = 'black', lwd = 1, ylim = c(0, 1))
lines(yGrid, rowMeans(baseDensity_100$F0_est), col = 'red', lwd = 1)
lines(yGrid, rowMeans(baseDensity_400$F0_est), col = 'green', lwd = 1)

F0_df <- data.frame(y = yGrid,
                   F0 = c(F0_true, rowMeans(baseDensity_100$F0_est), rowMeans(baseDensity_200$F0_est),
                          rowMeans(baseDensity_400$F0_est)),
                   n = rep(c('True', '100', '200', '400'), each = length(yGrid)))

ggplot(F0_df, aes(x = y, y = F0, color = n)) +
  geom_line() +
  theme_minimal()

# KS Test stat comparison

ks100 <- ks.test(rowMeans(baseDensity_100$F0_est), F0_true, alternative = 'two.sided', simulate.p.value = TRUE)
ks200 <- ks.test(rowMeans(baseDensity_200$F0_est), F0_true, alternative = 'two.sided', simulate.p.value = TRUE)
ks400 <- ks.test(rowMeans(baseDensity_400$F0_est), F0_true, alternative = 'two.sided', simulate.p.value = TRUE)
c(ks100$statistic, ks200$statistic, ks400$statistic)  # KS Test Statistic

# # Parallel Loop
# num_datasets <- 1
# p0 <- 0.1
# p1 <- 0.4
# sim_study <- list()
# nVec <- c(100, 500)
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
# # saveRDS(sim_study, file = "dpglm_simulation.rds")
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
# # # saveRDS(sim_study, file = "dpglm_simulation.rds")
# #
# # Sys.time() - start_time
#
# # Stop parallel backend
# stopCluster(cl)
