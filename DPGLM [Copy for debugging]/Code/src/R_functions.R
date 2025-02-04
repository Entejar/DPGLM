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

log_posterior_u <- function(u, z, theta, alpha) {
  R <- 250
  eps <- 1e-6
  zstar <- c(seq(0 + eps, 1 - eps, length.out = R), z)
  
  exp_mat <- exp(outer(theta, zstar, "*"))  
  u_exp_mat <- u * exp_mat 
  
  log_1_plus_u_exp_mat_sum <- log(1 + colSums(u_exp_mat))  
  neg_log_post <- alpha * mean(log_1_plus_u_exp_mat_sum[1:R]) + 
    mean(log_1_plus_u_exp_mat_sum[(R + 1):length(log_1_plus_u_exp_mat_sum)])
  
  return(-neg_log_post)
}

logpdf_gamma <- function(x, shape, scale) {
  (shape - 1) * log(x) - x / scale - shape * log(scale) - lgamma(shape)
}

update_u <- function(u, z, theta, alpha, delta) {
  for (i in seq_along(u)) {
    u_star <- u
    u_star[i] <- rgamma(1, shape = delta, scale = u[i] / delta)
    
    logQ_ratio <- logpdf_gamma(u[i], delta, u_star[i] / delta) - 
      logpdf_gamma(u_star[i], delta, u[i] / delta)
    
    logratio <- log_posterior_u(u_star, z, theta, alpha) - 
      log_posterior_u(u, z, theta, alpha) + logQ_ratio
    
    if (log(runif(1)) < logratio) {
      u <- u_star
    }
  }
  return(u)
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


