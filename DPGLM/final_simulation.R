rm(list = ls())
# Load Libraries
require(foreach)
require(doParallel)
source("load.R")

# Parallel Setup
num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Simulation III ---------------------------------------------------------------------
f0_true_kde <- readRDS("sim_truth/f0_kde_unif_0_1.rds")
spt <- f0_true_kde$spt
f0_kde  <- f0_true_kde$f0 

dx <- diff(spt)[1]
m0 <- sum(spt * f0_kde) * dx

#sum(f0_kde) * diff(spt)[1]

truth1 <- list(beta = c(1, 0), 
               f0 = matrix(c(spt, f0_kde), 
                           ncol = 2, 
                           byrow = FALSE))

truth2 <- list(beta = c(0.2, 0.7), 
               f0 = matrix(c(spt, f0_kde), 
                           ncol = 2, 
                           byrow = FALSE))


sim_data <- function(n, true_beta) {
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
  
  #theta <- rep(0, n) # for testing
  
  # Compute the matrix of exponentiated outer products
  exp_matrix <- exp(outer(theta, spt, "*"))
  
  dx <- diff(spt)[1]
  # Compute the dot product row-wise and apply the log
  btheta <- log(exp_matrix %*% f0_kde * dx)
  
  X       <- X %>% as.matrix()
  fY <- t(sapply(1:n, function(i) {
    exp(theta[i] * spt - btheta[i]) * f0_kde 
  }))
  
  # check if density y | x is valid
  #rowSums(fY) * diff(spt)[1] # should be 1
  #lines(spt, fY[2, ], type = "l", col = "red", lwd = 2, xlab = "y", ylab = "f(y | x)")
  
  Y     <- sapply(1:n, function(i) {
    sample(spt, 1, prob = fY[i, ])
  })
  return(data.frame(Y, theta = theta, X))
}

out <- list()
n_datasets <- 250 
seed_init <- 838160135 #sample.int(.Machine$integer.max, 1)
iter <- 1000

start_time <- Sys.time()

out <- foreach(
  indx = 1:n_datasets,
  .packages = c("tidyr", "gldrm", "mvtnorm")
) %dopar% {
  tryCatch({
    # Set unique seed for each worker
    set.seed(seed_init + indx)  
    dat1 <- sim_data(250, truth1$beta)
    
    scale_factors <- c(1, 1) %>% sqrt()
    r <- 1:50
    fit1 <- fit_func(dat1[r,], iter, truth1, truth1$beta, 2.5, scale_factors, m0)
    r <- 1:100
    fit2 <- fit_func(dat1[r,], iter, truth1, truth1$beta, 1.7, scale_factors, m0)
    r <- 1:250
    fit3 <- fit_func(dat1[r,], iter, truth1, truth1$beta, 1, scale_factors, m0)
    out1 <- list(fit1 = fit1, fit2 = fit2, fit3 = fit3)
    
    dat2 <- sim_data(250, truth2$beta)
    
    scale_factors <- c(1, 1) %>% sqrt()
    r <- 1:50
    fit1 <- fit_func(dat2[r,], iter, truth2, truth2$beta, 2.5, scale_factors, m0)
    r <- 1:100
    fit2 <- fit_func(dat2[r,], iter, truth2, truth2$beta, 1.7, scale_factors, m0)
    r <- 1:250
    fit3 <- fit_func(dat2[r,], iter, truth2, truth2$beta, 1, scale_factors, m0)
    out2 <- list(fit1 = fit1, fit2 = fit2, fit3 = fit3)
    
    
    
    return(list(out1 = out1, out2 = out2))
  }, error = function(e) {
    list(error = TRUE, message = as.character(e))
  })
}

# Save Results
saveRDS(out, file = "cache/final_sim_dpglm.rds")

# Stop Parallel Backend
stopCluster(cl)

# Compute Execution Time
time_tot <- Sys.time() - start_time
print(time_tot)

