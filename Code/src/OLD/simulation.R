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


# Source files [TO BE DELETED] -------------------------------------------------
setwd("~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/src")

source("G0_generator.R")
source("K_dist.R")
source("link_fn.R")
source("theta_solver.R")
source("theta_sampler.R")
source("u_sampler.R")
source("crm_sampler.R")
#source("beta_sampler_separate.R")
source("z_sampler_unifK.R")
#source("Sigma_beta.R")
source("resample_zstar.R")



## Load Data --------------------------------------------------------------------
setwd("/Users/entejar/Library/CloudStorage/Box-Box/DP-GLM/Data")

hustadTD <- readRDS("hustadTD.rds")

setwd("/Users/entejar/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM")

## Subset TD Data: Multi Word, ordered y --------------------------------------------------
hustadTDMW <- filter(hustadTD %>% as.data.frame(),
                    intelligibility_type == "multiword") %>% select(-intelligibility_type) %>% arrange(mean_intelligibility)

# hustadTDMW <- filter(hustadTD %>% as.data.frame(),
#                      intelligibility_type == "single-word") %>% select(-intelligibility_type) %>% arrange(mean_intelligibility)

y    <- hustadTDMW$mean_intelligibility
dens <- density(y, kernel = "gaussian", bw = 'nrd0', adjust = 1, n = 1000, cut = 3, ext = 4)
#plot(dens)
dat <- data.frame(y = dens$x, fy = dens$y)
indx <- dat$y >= 0 & dat$y <= 1
yGrid <- dat$y[indx]
f0 <- dat$fy[indx]
n <- 250
y <- sort(sample(dat$y[indx], n, prob = dat$fy[indx], replace = T))
#plot(density(y))
X <- cbind(1, matrix(c(runif(n, 1, 3), runif(n, -1, 1)), ncol = 2))
X[, -1] <- scale(X[, -1])
X <- X %>% as.matrix()
#beta <- c(0.05, 0.3, 0.9)
link <- 'logit'
#mu <- link_fn(X %*% beta, link)$mu
y <- y %>% as.numeric()
mu0 <- mean(y)
n <- length(y)
p <- dim(X)[2]
l <- unique(y)
max_spt <- max(y)

# Tuning Parameters --------------------------------------------------------------------
rho  <- 1
M <- 20
alpha <- 1
G0.dist <- 6
delta <- 2
kdist <- 6                    # ATTENTION: Choose K as 6 or 7?
sigmaTheta <- rep(0.001, n)
#a0 <- -Inf
#b0 <- +Inf
a00 <- 0
b00 <- 1
c0 <- 0.025

# ATTENTION: whenever switch choices of K, G0: 1 <-> 5, make sure to comment out / uncomment truncLoc
# here and also in z_sampler and also change a0, b0, a00, b00 from 0/1 <-> - + Inf

## Initialization ----
iter <- 3000
beta_samples <- matrix(NA, nrow = iter, ncol = p)
theta_samples <- matrix(NA, nrow = iter, ncol = n)
z_samples <- matrix(NA, nrow = iter, ncol = n)
G_samples   <- list()
init <- gldrm(y ~ X[, -1], link = link)
plot(density(init$theta))
beta <- beta_samples[1,] <- init$beta %>% as.numeric()

# n <- 250
# X <- cbind(1, matrix(c(runif(n, 1, 3), runif(n, -1, 1)), ncol = 2))
# beta <- c(0.05, 0.3, 0.9)
# link <- 'logit'
# mu <- link_fn(X %*% beta, link)$mu
# yGrid <- seq(0, 1, length.out = 100)
# f0 <- dbeta(yGrid, 5, 2)
# # f0[1] <- 2
# # f0[length(yGrid)] <- 1
# # f0 <- f0 / (sum(f0[2:length((yGrid))]) + f0[1] + f0[length(yGrid)])
# muf0 <- sum(yGrid * f0)
# muf0
# round(f0, 3)
# #f0[2:(length(yGrid) - 1)] <- (1 - f0[1] - f0[length(yGrid)]) * f0[2:(length(yGrid) - 1)]
# out <- theta_solver(yGrid, f0, mu, NULL)
# theta <- out$tht
# btheta <- out$btht
#
# probmt <- t(sapply(1:n, function(i) exp(theta[i]*yGrid - btheta[i])*f0))
# sampleY <- function(prob){
#   t <- runif(1)
#   y <- ifelse(t < prob[1], 0, ifelse(t < prob[1] + prob[length(yGrid)], 1,
#               sample(yGrid[-c(1, length(yGrid))], 1, prob = prob[-c(1, length(yGrid))])))
#   return(y)
# }
# y <- sapply(1:n, function(i) sampleY(prob = probmt[i, ])) ## simulated y
# hist(y)
# round(probmt[1, ], 2)
# round(f0, 2)
# p <- ncol(X)
#
f0y <- function(y, spt, f0) {
  n <- length(y)
  f0_y <- numeric(n)

  for (i in 1:n) {
    f0_y[i] <- sum(f0[y[i] == spt])
  }

  return(f0_y)
}

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
G_samples[[1]] <- list(z.tld = spt, J.tld = Jstar)

mu <- link_fn(X %*% beta, link)$mu
out <- theta_solver(spt, J.tld, mu, NULL)
tht <- out$tht
btht <- out$btht
bpr2 <- out$bpr2
z_samples[1, ] <- z <- y
temp <- resample_zstar(z)
zstar <- temp$zstar
nstar <- temp$nstar
SigmaY <- sigma0 <- rep(-99, n) # is this the game changer? should not use init$seMu for sparse y's

#SigmaY <- rep(sigma0, n) #sqrt(init$bPrime2) WRONG!! init$seMu
T.vec <- exp(btht)
u <- rgamma(n, shape = 1, rate = T.vec)
mubetaprior <-  init$beta
Sigbetaprior <- init$varbeta
#Sig <- Sigma_beta(X, mu, link,  bpr2, rho)
Sig <- rho * init$varbeta
mu0 <- round(mu0, 6)
mu0G <- -99 #mu0
sigma0G <- -99 #sqrt(sum((init$spt - sum(init$spt * init$f0))^2 * init$f0) / sum(init$f0))

seed <- sample.int(.Machine$integer.max, 1)
set.seed(seed)
llik_ <- lpr <- lpost <- numeric(iter)
zMat <- thtMat <- numeric(0)


#itr <- 1
#itr <- itr + 1
tm <- Sys.time()
for(itr in 2:iter){
  #rho <- ifelse(itr < 500, 0.7, ifelse(itr < 1500, 0.5, 0.25))
  #Sig <- rho * init$varbeta
  n <- dim(X)[1]
  p <- dim(X)[2]
  l <- length(z.tld)

  # Beta update --------------------------------
  pr_bt <- mvtnorm::rmvnorm(1, mean = beta, sigma = Sig) %>% as.vector()
  pr_mu <- link_fn(X %*% pr_bt, link)$mu %>% as.numeric()
  if(sum(min(z.tld) <= pr_mu & pr_mu <= max(z.tld)) == n){

    pr_pbt <- dmvnorm(pr_bt, mean = mubetaprior, sigma = Sigbetaprior, log = T)
    cr_pbt <- dmvnorm(beta, mean = mubetaprior, sigma = Sigbetaprior, log = T)
    cr_qbt <- dmvnorm(beta, mean = pr_bt, sigma = Sig, log = T)
    pr_qbt <- dmvnorm(pr_bt, mean = beta, sigma = Sig, log = T)

    alp <- min(0, pr_pbt - cr_pbt + cr_qbt - pr_qbt)

    if((log(runif(1)) < alp)){ #} & sum(pr_bt < 0) == 0){
      beta <- pr_bt
    }
  }
  print(cat("done beta"))
  # Theta update ------------------------------------
  mu <- link_fn(X %*% beta, link)$mu
  out <- theta_solver(z.tld, J.tld / sum(J.tld), mu, NULL)
  tht.tld <- out$tht
  sigma2Theta <- (sigmaTheta)^2
  muTheta <- tht.tld + (z * sigma2Theta)
  tht <- rnorm(n, mean = muTheta, sd = sigmaTheta)

  # U update ----------------------------------------
  mu0GGG <- mu0G #truncLoc(mu0G, sigma0G, a00, b00)
  u <- u_sampler(u, tht, z, n, alpha, G0.dist, mu0GGG, sigma0G, delta, a00, b00)

  # CRM update --------------------------------------
  N <- 3001
  w <- -log(seq(from = exp(-1e-06), to = exp(-5e-04), length = N))
  z.G0 <- sort(G0_generator(G0.dist, mu0G, sigma0G, a00, b00))
  R <- length(z.G0)
  mt <- matrix(NA, nrow = n, ncol = R)
  for(i in 1:R){
    mt[, i] <- tht * z.G0[i] %>% as.numeric()
  }
  mx <- apply(mt, 1, max)
  #mx <- 0 # is mx adjustment needed? for numerical issue.. Ans: YES.
  psi <- function(z){
    z <- rep(z, length(tht))
    sum(u * exp(tht * z - mx))
  }
  psi.z <- sapply(z.G0, psi)

  # First Part: random jumps
  fnS <- sapply(w, function(s) mean(exp( - (1 + psi.z) * s ) / s))
  dw <- diff(w)
  h <- (fnS[-1] + fnS[-N])/2
  Nv <- c(rev(cumsum(rev(dw * h))), 0)
  #plot(w, Nv)
  Nv <- alpha * Nv
  xi <- cumsum(rexp(M))
  iNv <- N
  J <- sapply(seq(M), function(i) {
    while (iNv > 0 && Nv[iNv] < xi[i]) {
      iNv <- iNv - 1
    }
    w[iNv + 1]
  })

  # First Part: random locations
  z <- sapply(seq(M), function(m) {
    xi <- runif(1)
    temp <- exp( - (1 + psi.z) * J[m] )
    cutoff <- sum(temp) * xi
    z.G0[min(which(cumsum(temp) - cutoff > 0))]
  })

  # Second Part: random jumps [for fixed locations]
  psi.zstar <- colSums(matrix(rep(u, length(zstar)) * exp(rep(tht, length(zstar)) * zstar), nrow = n))
  Jstar <- rgamma(length(zstar), shape = nstar, rate = psi.zstar + 1)

  z.tld_temp <- c(z, zstar)
  if(min(mu) >= min(z.tld_temp) && max(mu) <= max(z.tld_temp)){
    z.tld <- z.tld_temp
    J.tld <- c(J, Jstar)
  } else{
    z.tld <- z.tld
    J.tld <- J.tld
  }

  temp1 <- min(y) %in% z.tld
  temp2 <- max(y) %in% z.tld
  if(temp1){
    z.tld <- z.tld
    J.tld <- J.tld} else{
      z.tld <- c(min(y), z.tld)
      J.tld <- c(min(J.tld), J.tld)
    }
  if(temp2){
    z.tld <- z.tld
    J.tld <- J.tld} else{
      z.tld <- c(z.tld, max(y))
      J.tld <- c(J.tld, min(J.tld))
    }

  # Z update ------------------------------------------------------------------
  if(length(z.tld) > n){
    SigmaYY <- c(SigmaY, rep(mean(SigmaY), length(z.tld) - n))
  } else{
    SigmaYY <- SigmaY
  }
  temp <- data.frame(z.tld = z.tld, J.tld = J.tld) %>% arrange(z.tld)
  z.tld <- temp[, 1]
  J.tld <- temp[, 2]
  z <- z_sampler(y, n, c0, kdist, z.tld, J.tld, tht)
  #hist(z)
  #
  # zstar and nstar update ----------------------------------------------------
  temp <- resample_zstar(z)    ## add resampling part in resample_zstar.R to reduce 'sticky clusters effect'
  zstar <- temp$zstar
  nstar <- temp$nstar
  print(cat("done f0 and all"))

  # Storing MCMC simulations --------------------------------------------------
  z_samples[itr, ] <- z
  beta_samples[itr,] <- beta
  theta_samples[itr, ] <- tht
  G_samples[[itr]] <- list(z.tld = z.tld, J.tld = J.tld)

  zMat <- cbind(zMat, z)
  thtMat <- cbind(thtMat, tht)

  if(itr %% 10 == 0){
    dat <- data.frame(z = y, age = X[, 2])
    zMean <- apply(zMat, 1, mean) # median?
    thtMean <- apply(thtMat, 1, mean)
    zMat <- thtMat <- numeric(0)
    temp <- resample_zstar(zMean)
    yhat <- temp$zstar
    nYhat <- temp$nstar

    cat("Finishing iteration = ", itr, ",",
        "estErr =", mean(abs(y - zMean)), ",",
        "time taken = ", Sys.time() - tm, "\n")

    # p0 <- ggplot(data = data.frame(z = c(z, zMean), grp = rep(c("1", "2"), each = n))) +
    #   geom_density(aes(x = z, col = grp), alpha = 0.5) + theme_bw()

    f0est <- rep(z.tld, ceiling(1000 * rep(J.tld / sum(J.tld))))

    p1 <- ggplot(data = data.frame(x = f0est)) +
      geom_density(aes(x = x)) +
      theme_bw()

    p2 <- ggplot(data = data.frame(y = y, yhat = zMean)) +
      geom_point(aes(x = y, y = yhat), col = "blue", size = 0.4) + theme_bw() +
      geom_abline(slope = 1, intercept = 0)

    p3 <- ggplot(data = data.frame(theta = thtMean)) +
      geom_density(aes(x = theta)) +
      theme_bw()

    p4 <- ggplot(data = data.frame(iteration = 1:itr,
                                   beta0 = beta_samples[1:itr, 1])) +
      geom_line(aes(x = iteration, y = beta0)) + theme_bw()

    p5 <- ggplot(data = data.frame(iteration = 1:itr,
                                   beta1 = beta_samples[1:itr, 2])) +
      geom_line(aes(x = iteration, y = beta1)) + theme_bw()

    p6 <- ggplot(data = data.frame(iteration = 1:itr,
                                   beta2 = beta_samples[1:itr, 3])) +
      geom_line(aes(x = iteration, y = beta2)) + theme_bw()

    p <- (p1 + p2 + p3)/ (p4 + p5 + p6)
    print(p)
    tm <- Sys.time()
  }
}


dp_spglm <- list(z = z_samples, beta = beta_samples, theta = theta_samples, G = G_samples)
saveRDS(dp_spglm, "0701RDS3Simulation.rds")
