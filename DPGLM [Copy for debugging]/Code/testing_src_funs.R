crm1_samples <- crm2_samples <- list()
z <- y

for(i in 1:25){
u <- u_sampler(u, z, theta, alpha, delta)

# CRM update --------------------------------------
crm_star <- crm_sampler_fractionalY(M, u, zstar, nstar, RL, RJ, theta, alpha, meanY_x, 
                                    y, shape_a, shape_b)
RL <- crm_star$RL
zstar <- crm_star$zstar
RJ <- crm_star$RJ
Jstar <- crm_star$Jstar

atoms <- c(RL, zstar)
jumps <- c(RJ, Jstar)

crm1_samples[[i]] <- matrix(c(atoms, jumps), ncol = 2)
}

library(dirichletprocess)

dp <- DirichletProcessBeta(z, 1)
dp <- Fit(dp, 25)

library(BNPdensity)

out <- BNPdensity::MixNRMI1(y, Alpha = 1, Kappa = 0.015, Gama = 0.5, 
                distr.k = "3", distr.p0 = "3")




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
  theta <- rep(0, n)
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
hist(y)
X <- dat[, -c(1,2)]

# Tuning Parameters --------------------------------------------------------------------
rho <- 1
M <- 20
alpha <- 1
delta <- 2
sigma_theta <- 1.5 #sd(true_theta) / 2 #0.001 
c0 <- 0.25
beta.sigma <- 1  #max(abs(true_beta)) / 2

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
iter <- 50
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
