# *************************************************************************** #
#              Parallel code for the null cases simulation study              #
# *************************************************************************** #

# Load Libraries ---------------------------------------------------------------
require(foreach)
require(doParallel)
source("load.R")

# Parallel Setup ----------------------------------------------------------------
num_cores <- 6
cl        <- makeCluster(num_cores)
registerDoParallel(cl)

# Data simulation functions: null cases -------------------------------------------
# Set I

# True f0 = KDE on y -----------------------------------------------------------
f0_true_kde <- readRDS("data/f0_true_kde.rds")
spt <- f0_true_kde$spt
f0_kde  <- f0_true_kde$f0
mu0 <- sum(spt * f0_kde) / sum(f0_kde)

sim_I <- function(n) {
  X     <- cbind(1, matrix(c(rnorm(n, 1, 0.5), rnorm(n, 2, 1)), ncol = 2))
  X       <- X %>% as.matrix()
  btheta <- theta <- rep(0, n)
  fY <- t(sapply(1:n, function(i) {
    f0_kde
  }))
  Y     <- sapply(1:n, function(i) {
    sample(spt, 1, prob = fY[i, ])
  })
  return(data.frame(Y, X) %>% arrange(Y))
}

# Set II

# True f0 = Beta on y ----------------------------------------------------------
f0_true_beta <- readRDS("data/f0_true_beta.rds")
spt <- f0_true_beta$spt
f0_beta  <- f0_true_beta$f0

sim_II <- function(n, p0, p1) {
  X     <- cbind(1, matrix(c(rnorm(n, 1, 0.5), rnorm(n, 2, 1)), ncol = 2))
  X       <- X %>% as.matrix()
  btheta <- theta <- rep(0, n)          # theta <- const = say 1?
  fY <- t(sapply(1:n, function(i) {
    f0_beta
  }))
  Y     <- sapply(1:n, function(i) {
    u <- runif(1)
    if (u < p0) {
      0
    } else if (u < p0 + p1) {
      1
    } else {
      sample(spt, 1, prob = f0_beta)
    }
  })
  return(data.frame(Y, X) %>% arrange(Y))
}


# ***************** Parallel Loop ********************************************* #
n_datasets <- 20
n_settings <- 2
n_tot <- n_datasets * n_settings
out_list <- list()
start_time <- Sys.time()
seed_init <-  sample.int(.Machine$integer.max, 1)
n <- 400
set.beta <- TRUE
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

tuning.params <- list(
  rho = rho,
  M = M,
  alpha = alpha,
  G0.dist = G0.dist,
  delta = delta,
  k.dist = kdist,
  sigma_theta = sigma_theta,
  a00 = a00,
  b00 = b00,
  c0 = c0,
  beta.sigma = beta.sigma
)


out_list <- foreach(
  indx = 1:n_tot,
  .packages = c("gldrm", "tidyverse", "mvtnorm", "DPGLM")
) %do% {
  seed_val <- seed_init + indx
  if (indx <= (n_tot / 2)) {
    dat <- sim_I(n)
  } else {
    dat <- sim_II(n, 0.1, 0.4)
  }
  y <- dat[, 1] %>% as.numeric()
  X <- dat[, -1] %>% as.matrix()
  
  out100 <- dpglm(
    y = y[1:100],
    X = X[1:100, ],
    iter = 2000,
    tuning.params = tuning.params,
    set.beta = set.beta
  )
  
  out200 <- dpglm(
    y = y[1:200],
    X = X[1:200, ],
    iter = 2000,
    tuning.params = tuning.params,
    set.beta = set.beta
  )
  
  out400 <- dpglm(
    y = y,
    X = X,
    iter = 2000,
    tuning.params = tuning.params,
    set.beta = set.beta
  )
  
  out <- list(out100, out200, out400)
  return(out)
}

Sys.time() - start_time

# Stop parallel backend
stopCluster(cl)

# Save result ------------------------------------------------------------------
saveRDS(out_list, file = "oct_20_null_cases_sim.rds")




# *************************************************************************** #