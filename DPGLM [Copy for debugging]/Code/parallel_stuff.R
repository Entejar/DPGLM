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
  beta <- c(1, 0, 0)
  mu <- as.numeric(expit(X %*% beta))
  out <- theta_sol(spt, f0_kde, mu, NULL)
  theta <- out$theta
  btheta <- out$btht
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
n_datasets <- 10
n_settings <- 2                           # ATTENTION: DO NOT CHANGE THIS VALUE!
n_tot <- n_datasets * n_settings
n_iter <- 500
out_list <- list()
seed_init <-  sample.int(.Machine$integer.max, 1)
n <- 400                                  # Make sure n is divisible by 4
set.beta <- FALSE                         # ATTENTION: By default, set.beta = FALSE . Change to TRUE?
rho  <- 1
M <- 20
alpha <- 1
delta <- 2                                # ATTENTION: K and G are fixed as unif(z-c0, z+c0) and unif(0,1) respectively. For other 
# choices, see settings 1,2,4 in `implementation notes` in Rmarkdown folder.
sigma_theta <- 0.001
c0 <- 0.025
beta.sigma <- 1

tuning.params <- list(
  rho = rho,
  M = M,
  alpha = alpha,
  delta = delta,
  sigma_theta = sigma_theta,
  c0 = c0,
  beta.sigma = beta.sigma
)

start_time <- Sys.time()

out_list <- foreach(
  indx = 1:n_tot,
  .packages = c("gldrm", "tidyverse", "mvtnorm", "DPGLM")
) %dopar% {
  tryCatch({
    seed_val <- seed_init + indx
    if (indx <= (n_tot / 2)) {
      dat <- sim_I(n)
    } else {
      dat <- sim_II(n, 0.1, 0.4)
    }
    y <- dat[, 1] %>% as.numeric()
    X <- dat[, -1] %>% as.matrix()
    
    out100 <- dpglm(
      y = y[1:(n/4)],
      X = X[1:(n/4), ],
      iter = n_iter,
      tuning.params = tuning.params,
      set.beta = set.beta
    )
    
    out200 <- dpglm(
      y = y[1:(n/2)],
      X = X[1:(n/2), ],
      iter = n_iter,
      tuning.params = tuning.params,
      set.beta = set.beta
    )
    
    out400 <- dpglm(
      y = y,
      X = X,
      iter = n_iter,
      tuning.params = tuning.params,
      set.beta = set.beta
    )
    
    out <- list(mcmc_samples = list(MS1 = out100, MS2 = out200, MS3 = out400), data = dat)
    return(out)
  }, error = function(e) {
    return(list(error = TRUE, message = as.character(e)))
  })
}

time_tot <- Sys.time() - start_time

# Stop parallel backend --------------------------------------------------------
stopCluster(cl)

# Save result ------------------------------------------------------------------
saveRDS(out_list, file = "cache/oct_24_null_cases_sim_10_replicates_run2.rds")




# *************************************************************************** #