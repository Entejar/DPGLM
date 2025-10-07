# -----------------------------------------------------------------------------#
#                               Update z                                      #
# -----------------------------------------------------------------------------#

# Sampling latent variables z from their full conditional with uniform kernel K

z_sampler_unifK <- function(y, c0, crm.atoms, crm.jumps, tht, min_y, max_y) {
  n <- length(y)
  eps <- 1e-6
  lower <- crm.atoms - c0
  upper <- crm.atoms + c0
  lower[lower < min_y] <- min_y + eps  
  upper[upper > max_y] <- max_y - eps
  width <- upper - lower
  z <- numeric(n)
  for (i in 1:n) {
    indices <- y[i] >= lower & y[i] <= upper
    indx <- which(indices)
    log_prob <- (1 / width[indx]) + (tht[i] * crm.atoms[indx]) + log(crm.jumps[indx])
    prob <- exp(log_prob - max(log_prob))
    if(sum(indices) == 1){
      z[i] <- crm.atoms[indices]
    } else {
      z[i] <- sample(crm.atoms[indx], 1, prob = prob)
    }
  }
  return(z)
}

c0_silverman <- function(y) {
  n <- length(y)
  sd_y <- sd(y)
  iqr_y <- IQR(y)
  return(0.9 * min(sd_y, iqr_y / 1.34) * n^(-1/5))
}

z_sampler_triK <- function(y, c0, crm.atoms, crm.jumps, tht) {
  n <- length(y)
  z <- numeric(n)
  
  for (i in 1:n) {
    # Compute triangular kernel weights
    tri_weight <- pmax(0, 1 - abs((y[i] - crm.atoms) / c0))
    
    # Compute log probabilities
    log_prob <- log(tri_weight) + (tht[i] * crm.atoms) + log(crm.jumps)
    
    # Normalize and sample
    prob <- exp(log_prob - max(log_prob))
    prob <- prob / sum(prob)
    
    # Sample z[i] based on triangular kernel weighting
    z[i] <- sample(crm.atoms, 1, prob = prob)
  }
  
  return(z)
}

z_sampler_epanK <- function(y, c0, crm.atoms, crm.jumps, tht) {
  n <- length(y)
  z <- numeric(n)
  
  for (i in 1:n) {
    # Compute Epanechnikov kernel weights
    epan_weight <- (3/4) * (1 - ((y[i] - crm.atoms) / c0)^2)
    epan_weight[abs(y[i] - crm.atoms) > c0] <- 0  # Ensure weights are 0 outside the range
    
    # Compute log probabilities
    log_prob <- log(epan_weight) + (tht[i] * crm.atoms) + log(crm.jumps)
    
    # Normalize and sample
    prob <- exp(log_prob - max(log_prob))
    prob <- prob / sum(prob)  
    
    # Sample z[i] based on Epanechnikov kernel weighting
    z[i] <- sample(crm.atoms, 1, prob = prob)
  }
  
  return(z)
}


# -----------------------------------------------------------------------------#
#                               Resample z                                     #
# -----------------------------------------------------------------------------#

resample_zstar <- function(z){
  z_table <- table(z)
  zstar   <- as.numeric(names(z_table))
  nstar   <- as.numeric(z_table)
  ## Write code for resampling zstar to avoid the ‘sticky clusters effect’
  
  return(list(zstar = zstar, nstar = nstar))
}

logpost_beta <- function(beta, z, X, atoms, jumps, mu_beta, sigma_beta){
  n <- length(z)
  #mu <- exp(X %*% beta) / (1 + exp(X %*% beta))
  mu <- plogis(X %*% beta)
  theta <- gldrm:::getTheta(spt = atoms, f0 = jumps, mu = mu, 
                            ySptIndex = NULL, sampprobs = NULL)$theta
  btheta <- b_theta(theta, atoms, jumps)
  log_post <- sum(theta * z - btheta) - sum((beta - mu_beta)^2 / (2 * sigma_beta^2)) # others const in beta
  return(log_post)
}

loglik <- function(z, X, beta, atoms, jumps){
  n <- length(z)
  #mu <- exp(X %*% beta) / (1 + exp(X %*% beta))   
  mu <- plogis(X %*% beta)
  theta <- gldrm:::getTheta(spt = atoms, f0 = jumps, mu = mu, ySptIndex = NULL, 
                            sampprobs = NULL)$theta
  
  btheta <- b_theta(theta, atoms, jumps) 
  
  f0_z <- numeric(n)
  for (i in 1:n) {
    f0_z[i] <- sum(jumps[z[i] == atoms])
  }
  loglik <- sum(theta * z - btheta + log(f0_z))
  return(loglik)
}



log_post_u <- function(u, zstar, nstar, theta, alpha, min_y, max_y) {
  # Number of grid points for the continuous part
  R <- 250
  eps <- 1e-6
  
  # Construct a grid for integration over the continuous part (G_0, uniform on (0,1))
  #z_grid <- seq(eps, 1 - eps, length.out = R)
  z_grid <- seq(min_y + eps, max_y - eps, length.out = R)
  diff_z <- diff(z_grid)[1]  # uniform grid spacing
  
  # ---- Continuous part ----
  # For each grid point v in z_grid, compute in a stabilized manner:
  #    log(1 + sum_i u_i * exp(theta_i * v))
  cont_vals <- sapply(z_grid, function(v) {
    A <- log(u) + theta * v
    max_A <- max(A)
    S <- exp(max_A) * sum(exp(A - max_A))
    log1p(S)  # log(1+S) computed in a numerically stable way
  })
  
  # Riemann-sum approximation of the integral over the continuous part:
  integral_continuous <- sum(cont_vals) * diff_z
  
  # ---- Discrete part ----
  # For each unique atom in zstar, with multiplicity given by nstar,
  # compute: nstar * log(1 + sum_i u_i * exp(theta_i * zstar))
  disc_vals <- sapply(zstar, function(v) {
    A <- log(u) + theta * v
    max_A <- max(A)
    S <- exp(max_A) * sum(exp(A - max_A))
    log1p(S)
  })
  
  sum_discrete <- sum(disc_vals * nstar)
  
  # Combine the continuous and discrete contributions
  neg_log_post <- alpha * integral_continuous + sum_discrete
  
  return(-neg_log_post)
}

sampler_u <- function(u, zstar, nstar, theta, alpha, delta, min_y, max_y) {
  n <- length(u) # Get the length of u
  
  for (i in seq_len(n)) {
    u_star <- u
    # Sample from Gamma distribution: shape = delta, scale = u[i] / delta
    u_star[i] <- rgamma(1, shape = delta, scale = u[i] / delta)
    
    # Compute logQ_ratio: log(q(ui | ui_star)) - log(q(ui_star | ui))
    logQ_ratio <- dgamma(u[i], shape = delta, scale = u_star[i] / delta, log = TRUE) - 
      dgamma(u_star[i], shape = delta, scale = u[i] / delta, log = TRUE)
    
    # Compute logratio
    logratio <- log_post_u(u_star, zstar, nstar, theta, alpha, min_y, max_y) - 
      log_post_u(u, zstar, nstar, theta, alpha, min_y, max_y) + logQ_ratio
    
    # Metropolis-Hastings acceptance step
    if (log(runif(1)) < logratio) {
      u[i] <- u_star[i]
    } # else, u[i] remains unchanged
  }
  
  return(u)
}