out <- readRDS("cache/sim_2_Feb13_beta_0.2_0.7.rds")

null_idx <- numeric(0)
for(i in 1:length(out)){
  if(!is.null(out[[i]]$error)){
    null_idx <- c(null_idx, i)
  } else {
    acc_prob <- c(out[[i]]$fit1$dpglm$crm_acceptance, out[[i]]$fit2$dpglm$crm_acceptance, 
                  out[[i]]$fit3$dpglm$crm_acceptance, out[[i]]$fit1_max$dpglm$crm_acceptance,
                  out[[i]]$fit2_max$dpglm$crm_acceptance, out[[i]]$fit3_max$dpglm$crm_acceptance)
    if(min(acc_prob) < 0.05){
      null_idx <- c(null_idx, i)
    }
  }
}
  
out <- out[-unique(null_idx)]
out1 <- out2 <- out3 <- list()

for(i in 1:length(out)){
  out1[[i]] <- out[[i]]$fit1
  out2[[i]] <- out[[i]]$fit2
  out3[[i]] <- out[[i]]$fit3
}

mu0 <- 0.7833362

b_prime_theta <- function(s_k, f_k, theta_i) {
  exp_term <- exp(theta_i * s_k)
  numerator <- sum(s_k * f_k * exp_term) # element-wise multiplication and sum
  denominator <- sum(f_k * exp_term)    # element-wise multiplication and sum
  return(numerator / denominator)
}

solve_beta <- function(b_prime, X, init_beta = NULL, tol = 1e-6, max_iter = 1000) {
  # Ensure inputs are matrices or vectors
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(init_beta)) {
    init_beta <- rep(0, p) # Initialize beta to zeros if not provided
  }
  
  # Define the objective function to minimize
  objective <- function(beta) {
    eta <- X %*% beta
    f_beta <- exp(eta) / (1 + exp(eta))
    sum((b_prime - f_beta)^2) # L2 norm squared
  }
  
  # Use optim to minimize the objective function
  res <- optim(
    par = init_beta,
    fn = objective,
    method = "BFGS", # Quasi-Newton optimization
    control = list(reltol = tol, maxit = max_iter)
  )
  
  # Return the optimized beta and additional information
  list(beta = res$par, value = res$value, convergence = res$convergence)
}

beta_processed <- function(out, mu0){
  beta <- out$dpglm$beta[j,]
  spt <- c(out$dpglm$crm1[[j]]$RL, out2[[i]]$dpglm$crm2[[j-1]]$zstar) 
  f0 <- c(out[[i]]$dpglm$crm1[[j]]$RJ, out2[[i]]$dpglm$crm2[[j]]$Jstar)
  theta0 <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu0, sampprobs = NULL, ySptIndex = NULL)$theta
  f0 <- exp(theta0 * spt) * f0 / sum(exp(theta0 * spt) * f0)
  theta <- out2[[i]]$dpglm$theta[j,] - theta0
  b_prime <- sapply(1:length(theta), function(j) b_prime_theta(spt, f0, theta[j]))
  beta <- solve_beta(b_prime, X = out2[[i]]$data[, -c(1,2)], init_beta = c(1,0))$beta
  beta1_mat[l,] <- beta
}

#i <- j <- 1
beta1 <- beta2 <- beta3 <- matrix(0, nrow = length(out1), ncol = 6)
for(i in 1:length(out1)){
  beta1_mat <- beta2_mat <- beta3_mat <- matrix(0, nrow = 500, ncol = 4)
  for(l in 1:500){
    j <- l + 500 
    
    beta <- out1[[i]]$dpglm$beta[j,]
    spt <- c(out1[[i]]$dpglm$crm[[j]]$z.tld) 
    f0 <- c(out1[[i]]$dpglm$crm[[j]]$J.tld)
    theta0 <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu0, sampprobs = NULL, ySptIndex = NULL)$theta
    f0 <- exp(theta0 * spt) * f0 / sum(exp(theta0 * spt) * f0)
    theta <- out1[[i]]$dpglm$theta[j,] - theta0
    b_prime <- sapply(1:length(theta), function(j) b_prime_theta(spt, f0, theta[j]))
    betastar <- solve_beta(b_prime, X = out1[[i]]$data[, -c(1,2)], init_beta = c(0.2,0.7))$beta
    beta1_mat[l,] <- c(beta, betastar)
    
    beta <- out2[[i]]$dpglm$beta[j,]
    spt <- c(out2[[i]]$dpglm$crm[[j]]$z.tld) 
    f0 <- c(out2[[i]]$dpglm$crm[[j]]$J.tld)
    theta0 <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu0, sampprobs = NULL, ySptIndex = NULL)$theta
    f0 <- exp(theta0 * spt) * f0 / sum(exp(theta0 * spt) * f0)
    theta <- out2[[i]]$dpglm$theta[j,] - theta0
    b_prime <- sapply(1:length(theta), function(j) b_prime_theta(spt, f0, theta[j]))
    betastar <- solve_beta(b_prime, X = out2[[i]]$data[, -c(1,2)], init_beta = c(0.2,0.7))$beta
    beta2_mat[l,] <- c(beta, betastar)
    
    beta <- out3[[i]]$dpglm$beta[j,]
    spt <- c(out3[[i]]$dpglm$crm[[j]]$z.tld) 
    f0 <- c(out3[[i]]$dpglm$crm[[j]]$J.tld)
    theta0 <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu0, sampprobs = NULL, ySptIndex = NULL)$theta
    f0 <- exp(theta0 * spt) * f0 / sum(exp(theta0 * spt) * f0)
    theta <- out3[[i]]$dpglm$theta[j,] - theta0
    b_prime <- sapply(1:length(theta), function(j) b_prime_theta(spt, f0, theta[j]))
    betastar <- solve_beta(b_prime, X = out3[[i]]$data[, -c(1,2)], init_beta = c(0.2,0.7))$beta
    beta3_mat[l,] <- c(beta,betastar)
    
    
  }
  
  X = out1[[i]]$data[, -c(1,2)]
  y = out1[[i]]$data[, 1]
  fit <- gldrm(y ~ X[, -1], link = "logit")
  beta1[i,] <- c(colMeans(beta1_mat), fit$beta)
  
  X = out2[[i]]$data[, -c(1,2)]
  y = out2[[i]]$data[, 1]
  fit <- gldrm(y ~ X[, -1], link = "logit")
  beta2[i,] <- c(colMeans(beta2_mat), fit$beta)
  
  X = out3[[i]]$data[, -c(1,2)]
  y = out3[[i]]$data[, 1]
  fit <- gldrm(y ~ X[, -1], link = "logit")
  beta3[i,] <- c(colMeans(beta3_mat), fit$beta)
  
  
}

beta1[, -c(3,4)] %>% colMeans()
beta2[, -c(3,4)] %>% colMeans()
beta3[, -c(3,4)] %>% colMeans()
          
beta1[, -c(3,4)] %>% apply(2, sd) # these are wrong!!
beta2[, -c(3,4)] %>% apply(2, sd)
beta3[, -c(3,4)] %>% apply(2, sd)

par(mfrow = c(3,2))

i <- 3

plot(out1[[i]]$dpglm$beta[,1], type = "l", col = "blue", xlab = "Iteration", ylab = "Beta0")
lines(out2[[i]]$dpglm$beta[,1], type = "l", col = "green")
lines(out3[[i]]$dpglm$beta[,1], type = "l", col = "red")
abline(h = 0.2, col = "black", lty = 2)

plot(out1[[i]]$dpglm$beta[,2], type = "l", col = "blue", xlab = "Iteration", ylab = "Beta1")
lines(out2[[i]]$dpglm$beta[,2], type = "l", col = "green")
lines(out3[[i]]$dpglm$beta[,2], type = "l", col = "red")
abline(h = 0.7, col = "black", lty = 2)

plot(out1[[i]]$dpglm$loglik, type = "l", col = "blue", xlab = "Iteration", ylab = "Loglik")
plot(out2[[i]]$dpglm$loglik, type = "l", col = "green", xlab = "Iteration", ylab = "Loglik")
plot(out3[[i]]$dpglm$loglik, type = "l", col = "red", xlab = "Iteration", ylab = "Loglik")

plot(out1[[i]]$dpglm$crm_acceptance, type = "l", col = "blue", xlab = "Iteration", ylab = "CRM acceptance")
lines(out2[[i]]$dpglm$crm_acceptance, type = "l", col = "green")
lines(out3[[i]]$dpglm$crm_acceptance, type = "l", col = "red")
abline(h = 0.2, col = "black", lty = 2)






f01 <- f02 <- f03 <- f01_max <- f02_max <- f03_max <- matrix(0, nrow = length(out), ncol = length(y_grid))
F01 <- F02 <- F03 <- F01_max <- F02_max <- F03_max <- matrix(0, nrow = length(out), ncol = length(y_grid))
for(i in 1:length(out)){
  temp <- get_f0(out1, i)
  F01[i,] <- colMeans(temp$F0)
  temp <- get_f0(out2, i)
  F02[i,] <- colMeans(temp$F0)
  temp <- get_f0(out3, i)
  F03[i,] <- colMeans(temp$F0)
  temp <- get_f0(out1_max, i)
  F01_max[i,] <- colMeans(temp$F0)
  temp <- get_f0(out2_max, i)
  F02_max[i,] <- colMeans(temp$F0)
  temp <- get_f0(out3_max, i)
  F03_max[i,] <- colMeans(temp$F0)
}

plot(y_grid, colMeans(F01), type = "l", col = "cyan")
lines(true_spt, cumsum(f0_kde), type = 'l', col = "black")
lines(y_grid, colMeans(F02), type = "l", col = "green")
lines(y_grid, colMeans(F03), type = "l", col = "blue")

plot(y_grid, colMeans(F01_max), type = "l", col = "cyan")
lines(true_spt, cumsum(f0_kde), type = 'l', col = "black")
lines(y_grid, colMeans(F02_max), type = "l", col = "green")
lines(y_grid, colMeans(F03_max), type = "l", col = "blue")
    
    
    
f0_1 <- get_f0(out1, i)
f0_2 <- get_f0(out2, i)
f0_3 <- get_f0(out3, i)
f0_1_ <- get_f0(out1_max, i)
f0_2_ <- get_f0(out2_max, i)
f0_3_ <- get_f0(out3_max, i)

plot(true_spt, f0_kde, type = 'l', col = "black")

plot(y_grid, colMeans(f0_1$f0), type = "l", col = "cyan")
lines(y_grid, colMeans(f0_2$f0), type = "l", col = "green")
lines(y_grid, colMeans(f0_3$f0), type = "l", col = "blue")

par(mfrow = c(1,2))

plot(y_grid, colMeans(f0_1$F0), type = "l", col = "cyan")
lines(y_grid, colMeans(f0_2$F0), type = "l", col = "green")
lines(y_grid, colMeans(f0_3$F0), type = "l", col = "blue")
lines(true_spt, cumsum(f0_kde), type = 'l', col = "black")

plot(y_grid, colMeans(f0_1_$F0), type = "l", col = "cyan")
lines(y_grid, colMeans(f0_2_$F0), type = "l", col = "green")
lines(y_grid, colMeans(f0_3_$F0), type = "l", col = "blue")
lines(true_spt, cumsum(f0_kde), type = 'l', col = "black")

