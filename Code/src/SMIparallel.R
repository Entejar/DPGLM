# Clear the workspace
rm(list = ls())

# Packages Required -------------------------------------------------------------
library(gldrm)
library(tidyverse)
library(readr)
library(truncnorm)
library(splines)
library(mvtnorm)
library(ggplot2)
library(patchwork)
library(splines2)
library(ReIns)
library(boxr)
library(foreach)
library(doParallel)

# Parallelization setup
num_cores <- 6

# Initialize parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)


# Source files [TO BE DELETED] -------------------------------------------------
setwd("~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/src")

source("G0_generator.R")
source("K_dist.R")
source("link_fn.R")
source("theta_solver.R")
source("theta_sampler.R")
source("u_sampler.R")
source("crm_sampler.R")
source("beta_sampler.R")
source("z_sampler_unifK.R")
source("resample_zstar.R")



## Load Data --------------------------------------------------------------------
setwd("/Users/entejar/Library/CloudStorage/Box-Box/DP-GLM/Data")

hustadTD <- readRDS("hustadTD.rds")

setwd("/Users/entejar/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM")

## Subset TD Data: Multi Word, ordered y --------------------------------------------------
hustadTDMW <- filter(hustadTD %>% as.data.frame(),
                     intelligibility_type == "multiword") %>% select(-intelligibility_type) %>% arrange(mean_intelligibility)



## Data: X, y & link function --------------------------------------------------
prbs <- c(0.33, 0.67)

Xstar <- ns(
  hustadTDMW$age_months,
  knots = quantile(hustadTDMW$age_months, probs = prbs),
  intercept = FALSE
)

X    <- cbind(rep(1, nrow(Xstar)), Xstar)
colnames(X) <- c("intercept", "ageBase1", "ageBase2", "ageBase3")
y    <- hustadTDMW$mean_intelligibility
mu0 <- mean(y)
X <- X %>% as.matrix()
y <- y %>% as.numeric()
n <- length(y)
p <- dim(X)[2]
l <- unique(y)
max_spt <- max(y)
link <- "logit"

# Tuning Parameters --------------------------------------------------------------------
rho  <- 1
M <- 20
alpha <- 1
G0.dist <- 6
delta <- 2
kdist <- 7                    # ATTENTION: Choose K as 6 or 7?
sigmaTheta <- rep(0.001, n)
#a0 <- -Inf
#b0 <- +Inf
a00 <- 0
b00 <- 1
c0 <- 0.025

# ATTENTION: whenever switch choices of K, G0: 1 <-> 5, make sure to comment out / uncomment truncLoc
# here and also in z_sampler and also change a0, b0, a00, b00 from 0/1 <-> - + Inf

init_fn <- function(iter, X, y, link, ){
  n <- length(y)
  p <- ncol(X)
  beta_samples <- matrix(NA, nrow = iter, ncol = p)
  theta_samples <- matrix(NA, nrow = iter, ncol = n)
  z_samples <- matrix(NA, nrow = iter, ncol = n)
  crm_samples   <- list()
  init <- gldrm(y ~ X[, -1], link = link)
  beta <- beta_samples[1, ] <- init$beta %>% as.numeric()

  # spt <- unique(round(y, 2))
  # J.tld <- f0 #rep(1 / l, l)
  # tht0 <- gldrm:::getTheta(
  #   spt       = spt,
  #   f0        = J.tld,
  #   mu        = mu0,
  #   sampprobs = NULL,
  #   ySptIndex = NULL)$theta
  # Jstar <- (J.tld * exp(tht0 * spt)) %>% `/` (sum(.))
  z.tld <- spt <- init$spt
  J.tld <- Jstar <- init$f0
  crm_samples[[1]] <- list(z.tld = spt, J.tld = Jstar)

  mu <- link_fn(X %*% beta, link)$mu
  out <- theta_solver(spt, J.tld, mu, NULL)
  tht <- out$tht
  btht <- out$btht
  bpr2 <- out$bpr2
  z_samples[1, ] <- z <- y
  temp <- resample_zstar(z)
  zstar <- temp$zstar
  nstar <- temp$nstar

  T.vec <- exp(btht)
  u <- rgamma(n, shape = 1, rate = T.vec)
  mubetaprior <-  init$beta
  Sigbetaprior <- init$varbeta
  Sig <- rho * init$varbeta
  mu0G <- -99    #0.5
  sigma0G <- -99 #sqrt(sum((init$spt - sum(init$spt * init$f0))^2 * init$f0) / sum(init$f0))
}

seed <- sample.int(.Machine$integer.max, 1)
set.seed(seed)

k1 <- sample.int(.Machine$integer.max, 1)
set.seed(500)
s_trn_ind <- sample(1:n, 100)
#f_trn_ind <- 1:n
k2 <- sample.int(.Machine$integer.max, 1)
set.seed(k2)
f_trn_ind <- sample(1:n, 100)
trn_ls <- list(s_trn_ind, f_trn_ind)
data <- list(X[trn_ls[[1]], ], y[trn_ls[[1]]], X[trn_ls[[2]], ], y[trn_ls[[2]]])

dataa <- model.matrix(~ numiadl + age + sex + iwr + netwc, data = dat)
iter_v <- c(10000, 10000)

ahead_study <- list()
trn_ind <- list()

start_time <- Sys.time()

ahead_study <- foreach(
  index = 1:2,
  .packages = c("gldrm", "extraDistr", "mvtnorm", "tidyr")
) %dopar% {
  X <- as.matrix(data[[(2 * index) - 1]])
  X[, c(2,4)] <- scale(X[, c(2,4)])
  y <- as.matrix(data[[2 * index]])
  mu0 <- mean(y)

  trn_ind <- trn_ls[[index]]
  dat_spglm <- dataa[trn_ind, ]
  dat_spglm[, c(3, 5)] <- scale(dat_spglm[, c(3, 5)])
  dat_spglm <- dat_spglm %>% as.data.frame()
  dat_spglm$sex <- factor(dat_spglm$sex)
  dat_spglm$netwc <- factor(dat_spglm$netwc)
  fit <- gldrm(numiadl ~ age + sex + iwr + netwc, link = "log", mu0 = mu0,
               data = dat_spglm)

  iter <- iter_v[index]
  seed_val <-  sample.int(.Machine$integer.max, 1)
  set.seed(seed_val)
  mcmc_samples <- dir_spglm(X, y, spt = 0:5, init = fit, rho = 1, iter = iter)

  out <- list(seed_val = seed_val, trn_ind = trn_ind,
              mcmc_samples = mcmc_samples)
  return(out)
}


# Save the results
saveRDS(ahead_study, file = "dir-spglm_ahead_study_1125_small_sample.rds")

Sys.time() - start_time

# Stop the parallel backend
stopCluster(cl)



