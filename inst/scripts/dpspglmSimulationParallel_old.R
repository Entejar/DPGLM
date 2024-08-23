# Clear Environment ----
#rm(list = ls())

# Required Packages ----
require(foreach)
require(doParallel)
require(tidyverse)
require(gldrm)
require(mvtnorm)


# Parallel Setup ----
num_cores <- 6
cl        <- makeCluster(num_cores)
registerDoParallel(cl)


# Source Files  -------------------------------------------------
setwd("~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/src")
source("functions.R")

setwd("/Users/entejar/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/inst/scripts")
source("dpglm.R")

# Load Data --------------------------------------------------------------------
setwd("/Users/entejar/Library/CloudStorage/Box-Box/DP-GLM/Data")
hustadTD <- readRDS("hustadTD.rds")
setwd("/Users/entejar/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM")

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

sim_III <- function(n) {
  X     <- cbind(1, matrix(c(runif(n, 1, 3), runif(n, -1, 1)), ncol = 2))
  X[, -1] <- scale(X[, -1])
  X       <- X %>% as.matrix()
  beta <- c(0.1, 0.25, 0.75)
  link <- 'logit'
  mu <- link_fn(X %*% beta, link)$mu
  out <- theta_solver(true_spt, true_f0, mu, NULL)
  theta <- out$tht
  btheta <- out$btht
  fY <- t(sapply(1:n, function(i) {
    exp(theta[i] * true_spt - btheta[i]) * true_f0
  }))
  Y     <- sapply(1:n, function(i) {
    sample(true_spt, 1, prob = fY[i, ])
  })
  return(data.frame(Y, X) %>% arrange(Y))
}

# Unit Testing
dat <- sim_I(50)
y <- dat[, 1]
X <- dat[, -1]
out <- dpglm(y, X, 'logit', 100)

# Parallel Loop
num_datasets <- 200
p0 <- 0.1
p1 <- 0.4
sim_study <- list()
nVec <- c(25, 100, 250)
scenarios <- 1:3
set <- expand.grid(nVec, scenarios)
R <- nrow(set)
seed_init <-  sample.int(.Machine$integer.max, 1)


start_time <- Sys.time()

sim_study <- foreach(dat_index = 1:num_datasets,
                     .packages = c("gldrm", "tidyverse", "truncnorm", "mvtnorm", "ggplot2")) %:%
  foreach(setting = 1:R) %dopar% {
    seed_val <- seed_init + dat_index + setting
    set.seed(seed_val)
    if (setting %in% 1:3) {dat <- sim_I(set[setting, 1])}
    if (setting %in% 4:6) {dat <- sim_II(set[setting, 1], p0, p1)}
    if (setting %in% 7:9) {dat <- sim_III(set[setting, 1])}
    X <- dat[, -1]
    y <- dat[, 1]
    mcmc <- dpglm(y, X, 'logit', 2000)
    out <- list(seed = seed_val,
                set = setting,
                data = dat,
                mcmc = mcmc
    )
    return(out)
  }


# Save result
saveRDS(sim_study, file = "0707_dpglm_simulation.rds")

Sys.time() - start_time

# Stop parallel backend
stopCluster(cl)
