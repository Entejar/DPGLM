# -----------------------------------------------------------------------------#
#                               Update beta                                    #
# -----------------------------------------------------------------------------#

# beta_sampler_fractionalY <- function(y, X, z.tld, J.tld, beta, Sig, mubetaprior,
#                                      Sigbetaprior, c0, sigma_theta){
#   zMN <- min(z.tld)
#   zMX <- max(z.tld)
#   B <- 1000     # number of samples for the Monte Carlo integration
#   n <- dim(X)[1]
#   p <- dim(X)[2]
#   crm <- matrix(c(z.tld, J.tld), nrow = 2, byrow = T)
#   
#   # proposed beta
#   pr_bt <- as.vector(mvrnormArma(1, beta, Sig)) 
#   # proposed mu
#   pr_mu <- as.numeric(expit(X %*% pr_bt))                # fractional Y
#   
#   # check if proposed mu is within bounds
#   if(min(pr_mu) >= zMN & max(pr_mu) <= zMX){
#     cr_theta <- gldrm:::getTheta(
#       spt = crm[1,],
#       f0  = crm[2,],
#       mu  = as.numeric(expit(X %*% beta)),
#       sampprobs  = NULL,
#       ySptIndex  = NULL,
#       thetaStart = NULL
#     )$theta
#     
#     pr_theta <- gldrm:::getTheta(
#       spt = crm[1,],
#       f0  = crm[2,],
#       mu  = pr_mu,
#       sampprobs  = NULL,
#       ySptIndex  = NULL,
#       thetaStart = cr_theta
#     )$theta
#     
#     # log likelihood (exact except MC approximation): change it to laplace approximation / numerical integration?
#     #pr_llik <- llik_beta_unifK(y, X, pr_theta, crm, c0, B, sigma_theta)
#     #cr_llik <- llik_beta_unifK(y, X, cr_theta, crm, c0, B, sigma_theta)
#     # # log likelihood (approximate)
#     pr_llik <- sum(pr_theta * y - b_theta(pr_theta, crm[1, ], crm[2, ]))
#     cr_llik <- sum(cr_theta * y - b_theta(cr_theta, crm[1, ], crm[2, ]))
#     # log prior
#     pr_logpr <- logpdf_mvnorm(pr_bt, mubetaprior, Sigbetaprior)
#     cr_logpr <- logpdf_mvnorm(beta, mubetaprior, Sigbetaprior)
#     # log q(beta | beta_proposed, Sig) = log q(beta_proposed | beta, Sig) 
#     # hence, cr_logq - pr_logq = 0 => omitted from the logratio
#     
#     
#     # log acceptance probability
#     logratio <- min(0, pr_llik - cr_llik + pr_logpr - cr_logpr) 
#     
#     if(log(runif(1)) < logratio){
#       beta <- pr_bt
#     }
#   }
#   
#   return(beta)
# }

# beta_sampler_fractionalY <- function(beta, Sig, theta, z.tld, J.tld, 
#                                      mubetaprior, Sigbetaprior, sigma_theta){
#   
#   # proposed beta
#   pr_bt <- as.vector(mvrnormArma(1, beta, Sig)) 
#   # proposed mu
#   pr_mu <- as.numeric(expit(X %*% pr_bt))                # fractional Y
#   
#   # check if proposed mu is within bounds
#   if(min(pr_mu) >= min(z.tld) & max(pr_mu) <= max(z.tld)){
#    theta_tilde <- gldrm:::getTheta(
#       spt = z.tld,
#       f0  = J.tld,
#       mu  = as.numeric(expit(X %*% beta)),
#       sampprobs  = NULL,
#       ySptIndex  = NULL,
#       thetaStart = NULL
#     )$theta
#     
#     pr_theta_tilde <- gldrm:::getTheta(
#       spt = z.tld,
#       f0  = J.tld,
#       mu  = pr_mu,
#       sampprobs  = NULL,
#       ySptIndex  = NULL,
#       thetaStart = theta_tilde
#     )$theta
#     
#     # # log likelihood 
#     pr_llik <- sum(dnorm(theta, pr_theta_tilde, sigma_theta, log = TRUE))
#     cr_llik <- sum(dnorm(theta, theta_tilde, sigma_theta, log = TRUE))
#     # log prior
#     pr_logpr <- logpdf_mvnorm(pr_bt, mubetaprior, Sigbetaprior)
#     cr_logpr <- logpdf_mvnorm(beta, mubetaprior, Sigbetaprior)
#     # log q(beta | beta_proposed, Sig) = log q(beta_proposed | beta, Sig) 
#     # hence, cr_logq - pr_logq = 0 => omitted from the logratio
#     
#     
#     # log acceptance probability
#     logratio <- min(0, pr_llik - cr_llik + pr_logpr - cr_logpr) 
#     
#     if(log(runif(1)) < logratio){
#       beta <- pr_bt
#     }
#   }
#   
#   return(beta)
# }

beta_sampler_fractionalY <- function(beta, Sig, theta, z.tld, J.tld, 
                                     mubetaprior, Sigbetaprior, sigma_theta){
  
  # proposed beta
  for(j in 1:length(beta)){
    pr_bt <- beta
    pr_bt[j] <- rnorm(1, mean = beta[j], sd = sqrt(Sig[j,j]))
  # proposed mu
  pr_mu <- as.numeric(expit(X %*% pr_bt))                # fractional Y
  
  # check if proposed mu is within bounds
  if(min(pr_mu) >= min(z.tld) & max(pr_mu) <= max(z.tld)){
    theta_tilde <- gldrm:::getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = as.numeric(expit(X %*% beta)),
      sampprobs  = NULL,
      ySptIndex  = NULL,
      thetaStart = NULL
    )$theta
    
    pr_theta_tilde <- gldrm:::getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = pr_mu,
      sampprobs  = NULL,
      ySptIndex  = NULL,
      thetaStart = theta_tilde
    )$theta
    
    # # log likelihood 
    pr_llik <- sum(dnorm(theta, pr_theta_tilde, sigma_theta, log = TRUE))
    cr_llik <- sum(dnorm(theta, theta_tilde, sigma_theta, log = TRUE))
    # log prior
    pr_logpr <- dnorm(pr_bt[j], mubetaprior[j], sqrt(Sigbetaprior[j,j]), log = TRUE)
    cr_logpr <- dnorm(beta[j], mubetaprior[j], sqrt(Sigbetaprior[j,j]), log = TRUE)
    # log q(beta | beta_proposed, Sig) = log q(beta_proposed | beta, Sig) 
    # hence, cr_logq - pr_logq = 0 => omitted from the logratio
    
    
    # log acceptance probability
    logratio <- min(0, pr_llik - cr_llik + pr_logpr - cr_logpr) 
    
    if(log(runif(1)) < logratio){
      beta[j] <- pr_bt[j]
    }
  }
  }
  
  return(beta)
}



# -----------------------------------------------------------------------------#
#                               Update z                                      #
# -----------------------------------------------------------------------------#

z_sampler_unifK <- function(y, c0, crm.atoms, crm.jumps, tht) {
  n <- length(y)
  lower <- crm.atoms - c0
  lower[lower < 0] <- 0
  width <- crm.atoms + c0 - lower
  z <- numeric(n)
  for (i in 1:n) {
    indices <- y[i] >= lower & y[i] <= crm.atoms + c0
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

# update_u <- function(u, z, theta, alpha){
#   u_grid <- seq(0, 10, 0.1)
#   for(i in 1:length(u)){
#     ustar <- u
#     logprob <- numeric(length(u_grid))
#     for(j in 1:length(u_grid)){
#     ustar[i] <- u_grid[j]
#     logprob[j] <-  log_posterior_u(ustar, z, theta, alpha)
#     }
#     prob <- exp(logprob - max(logprob))
#     u[i] <- sample(u_grid, 1, prob = prob)
#   }
#   return(u)
# }

logpost_beta <- function(beta, z, X, atoms, jumps, mu_beta, sigma_beta){
  n <- length(z)
  mu <- exp(X %*% beta) / (1 + exp(X %*% beta))
  theta <- gldrm:::getTheta(spt = atoms, f0 = jumps, mu = mu, 
                                  ySptIndex = NULL, sampprobs = NULL)$theta
  btheta <- b_theta(theta, atoms, jumps)
  log_post <- sum(theta * z - btheta) - sum((beta - mu_beta)^2) / (2 * sigma_beta^2) # others const in beta
  return(log_post)
}

loglik <- function(z, X, beta, atoms, jumps){
  n <- length(z)
  mu <- exp(X %*% beta) / (1 + exp(X %*% beta))
  theta <- gldrm:::getTheta(spt = atoms, f0 = jumps, mu = mu, ySptIndex = NULL, 
                            sampprobs = NULL)$theta
  btheta <- apply(exp(outer(theta, atoms, "*")), 1, function(j) log(sum(j * jumps)))
  f0_z <- numeric(n)
  for (i in 1:n) {
    f0_z[i] <- sum(jumps[z[i] == atoms])
  }
  loglik <- sum(theta * z - btheta + log(f0_z))
  return(loglik)
}

log_post_u <- function(u, z, theta, alpha) {
  R <- 250
  eps <- 1e-6
  
  # Construct zstar: concatenate linspace from (eps, 1-eps) with z
  zstar <- c(seq(eps, 1 - eps, length.out = R), z)
  
  # Create the exp_mat: equivalent to exp(diagmat(theta) * repmat(zstar, 1, length(theta)).t())
  exp_mat <- exp(outer(theta, zstar, `*`)) # Outer product with element-wise multiplication
  
  # Calculate u_exp_mat: equivalent to diagmat(u) * exp_mat
  u_exp_mat <- u * exp_mat
  
  # Sum across columns, and take log(1 + ...)
  log_1_plus_u_exp_mat_sum <- log(1 + colSums(u_exp_mat))
  
  # Compute the negative log posterior
  neg_log_post <- alpha * mean(log_1_plus_u_exp_mat_sum[1:R]) +
    sum(log_1_plus_u_exp_mat_sum[(R + 1):length(zstar)])
  
  return(-neg_log_post)
}

sampler_u <- function(u, z, theta, alpha, delta) {
  n <- length(u) # Get the length of u
  
  for (i in seq_len(n)) {
    u_star <- u
    # Sample from Gamma distribution: shape = delta, scale = u[i] / delta
    u_star[i] <- rgamma(1, shape = delta, scale = u[i] / delta)
    
    # Compute logQ_ratio: log(q(ui | ui_star)) - log(q(ui_star | ui))
    logQ_ratio <- dgamma(u[i], shape = delta, scale = u_star[i] / delta, log = TRUE) - 
      dgamma(u_star[i], shape = delta, scale = u[i] / delta, log = TRUE)
    
    # Compute logratio
    logratio <- log_post_u(u_star, z, theta, alpha) - 
      log_post_u(u, z, theta, alpha) + logQ_ratio
    
    # Metropolis-Hastings acceptance step
    if (log(runif(1)) < logratio) {
      u[i] <- u_star[i]
    } # else, u[i] remains unchanged
  }
  
  return(u)
}

sampler_u_ <- function(u, z, theta, alpha, delta) {
  n <- length(u) # Get the length of u
  u_star <- numeric(n)
  logQ_ratio <- 0
  for (i in seq_len(n)) {
    # Sample from Gamma distribution: shape = delta, scale = u[i] / delta
    u_star[i] <- rgamma(1, shape = delta, scale = u[i] / delta)
    
    # Compute logQ_ratio: log(q(ui | ui_star)) - log(q(ui_star | ui))
    logQ_ratio <- logQ_ratio + dgamma(u[i], shape = delta, scale = u_star[i] / delta, log = TRUE) - 
      dgamma(u_star[i], shape = delta, scale = u[i] / delta, log = TRUE)
  }
    
    # Compute logratio
    logratio <- log_post_u(u_star, z, theta, alpha) - 
      log_post_u(u, z, theta, alpha) + logQ_ratio
    
    # Metropolis-Hastings acceptance step
    if (log(runif(1)) < logratio) {
      u <- u_star
    } # else, u remains unchanged
  
  return(u)
}


loglikelihood <- function(y, X, beta, spt, f0, sampprobs = NULL) {
  # Ensure valid inputs
  spt <- as.vector(spt)
  f0 <- as.vector(f0)
  theta <- as.vector(theta)
  sptN <- length(spt)
  
  if (length(f0) != sptN) stop("spt and f0 must be vectors of equal length.")
  if (any(f0 < 0)) stop("f0 values cannot be negative.")
  
  # Map observations y to indices in spt
  ySptIndex <- match(y, spt)
  if (any(is.na(ySptIndex))) stop("Some values in y are not present in spt.")
  
  mu <- exp(X %*% beta) / (1 + exp(X %*% beta))
  theta <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu, 
                            ySptIndex = NULL, sampprobs = NULL)$theta
  
  # Compute unnormalized probabilities
  fUnstd <- f0 * exp(tcrossprod(spt, theta))
  
  # Normalize probabilities
  fTilt <- fUnstd / rep(colSums(fUnstd), each = sptN)
  
  # Adjust for sampling probabilities if provided
  if (!is.null(sampprobs)) {
    fTilt <- fTilt * t(sampprobs)
    fTilt <- fTilt / rep(colSums(fTilt), each = nrow(fTilt))
  }
  
  # Compute log-likelihood
  llik <- sum(log(fTilt[cbind(ySptIndex, seq_along(ySptIndex))]))
  
  return(llik)
}

idx <- numeric(n)
for (i in 1:n) {
  idx[i] <- spt[round(y[i], 4) == round(spt, 4)]
}
