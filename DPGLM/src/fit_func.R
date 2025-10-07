fit_func <- function(dat, iter, truth, mu_beta, sigma_beta, scale_factors, m0) {
  y <- dat[, 1]
  min_y <- 0
  max_y <- 1
  true_theta <- dat[, 2]
  X <- dat[, -c(1, 2)]
  n <- length(y)
  # Tuning Parameters --------------------------------------------------------------------
  rho <- 1
  M <- 20
  alpha <- 1
  delta <- 2
  c0 <- c0_silverman(y) / 4 #0.25
  # mu_beta <- truth$beta
  # sigma_beta <- 2.5
  # sigma_beta <- rho * sigma_beta
  
  # Data Preparation ---------------------------------------------------------------------
  X <- X %>% as.matrix()
  y <- y %>% as.numeric()
  n <- length(y)
  p <- dim(X)[2]
  
  # Link Function -----------------------------------------------------------------------
  link <- 'logit'
  
  ## Initialization -----------------------------------------------------------------------
  beta_samples <- matrix(NA, nrow = iter, ncol = p)
  theta_samples <- matrix(NA, nrow = iter, ncol = n)       
  u_samples <- z_samples <- matrix(NA, nrow = iter, ncol = n)
  crm_samples   <- list()
  lnlik_samples <- numeric(iter)
  gldrm_fit <- gldrm(y ~ X[, -1], link = link)
  beta_samples[1, ] <- beta <- truth$beta %>% as.numeric() #gldrm_fit$beta %>% as.numeric()
  beta_mode <- beta
  beta_cov <- gldrm_fit$varbeta
  #meanY_x <- exp(X %*% beta) / (1 + exp(X %*% beta))
  meanY_x <- plogis(X %*% beta)
  z.tld <- z_tld <-  truth$f0[, 1] %>% as.numeric() #gldrm_fit$spt %>% as.numeric() 
  J.tld <- J_tld <- truth$f0[, 2] %>% as.numeric() #gldrm_fit$f0 %>% as.numeric()   
  crm_samples[[1]] <- list(z.tld = z.tld, J.tld = J.tld)
  
  ord <- order(J_tld)[1:M]
  
  RL <- z_tld[ord]
  RJ <- J_tld[ord]
  theta_samples[1, ] <- theta <- true_theta #gldrm_fit$theta %>% as.numeric() #
  z_samples[1, ] <- z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y)
  #z_samples[1, ] <- z <- z_sampler_triK(y, c0, z.tld, J.tld, theta)
  #z_samples[1, ] <- z <- z_sampler_epanK(y, c0, z.tld, J.tld, theta)
  
  btheta <- b_theta(theta, z.tld, J.tld)
  
  T_vec <- exp(btheta)
  u_samples[1, ] <- u <- rgamma(n, shape = 1, rate = T_vec)
  
  resampled_z <- resample_zstar(z)
  zstar <- resampled_z$zstar
  nstar <- resampled_z$nstar
  Jstar <- rgamma(n = length(nstar), shape = nstar, rate = 1)
  
  count1 <- count2 <- 0
  #scale_factors <- c(1.5, 2.5) %>% sqrt() 
  #beta_cov <- 0.5 * gldrm_fit$varbeta
  
  theta0 <- gldrm:::getTheta(
    spt = z.tld,
    f0  = J.tld,
    mu  = m0,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = NULL
  )$theta
  #Jtld_0 <- exp(theta0 * z.tld) * J.tld / sum(exp(theta0 * z.tld) * J.tld)
  temp <- exp(theta0 * z.tld - max(theta0 * z.tld))
  Jtld_0 <- (temp * J.tld) / sum(temp * J.tld)
  lnlik_samples[1] <- lnlik <- loglik(z = z, X = X, beta = beta, atoms = z.tld, jumps = Jtld_0)
  #itr <- 1
  for(itr in 2:iter){
    #itr <- itr + 1
    # Optimization to find the mode
    result <- tryCatch({
      optim(par = beta, fn = logpost_beta, z = z, X = X,
            atoms = z.tld, jumps  = J.tld, mu_beta = mu_beta,
            sigma_beta = sigma_beta,
            control = list(fnscale = -1), hessian = TRUE)
    }, error = function(e) return(NULL))  # If error occurs, return NULL
    
    # Check if `optim()` failed
    if (!is.null(result) && is.finite(det(result$hessian))) {
      beta_mode_ <- result$par
      
      # Safely compute Hessian inverse
      beta_cov_ <- tryCatch({
        -solve(result$hessian)
      }, error = function(e) return(NULL))  # If inversion fails, return NULL
      
      # Only proceed if covariance matrix is valid and positive definite
      if (!is.null(beta_cov_)) {
        beta_cov_ <- diag(scale_factors) %*% beta_cov_ %*% diag(scale_factors)
        # Check positive semi-definiteness
        if (all(eigen(beta_cov_, symmetric = TRUE, only.values = TRUE)$values >= 
                -sqrt(.Machine$double.eps))) {
          
          # Optional: Force positive defineness if needed
          beta_cov_ <- Matrix::nearPD(beta_cov_)$mat %>% as.matrix()
          
          # Sample new beta from proposal distribution
          beta_ <- rmvnorm(n = 1, mean = beta_mode_, sigma = beta_cov_) %>% as.numeric()
          
          
          #mean_z_ <- exp(X %*% beta_) / (1 + exp(X %*% beta_))
          mean_z_ <- plogis(X %*% beta_)
          
          if (min(mean_z_) >= min(z.tld) && max(mean_z_) <= max(z.tld)) {
            
            # Compute log proposal values
            
            pr_logprop_beta <- dmvnorm(x = beta_, mean = beta_mode, sigma = beta_cov, log = TRUE)
            cr_logprop_beta <- dmvnorm(x = beta, mean = beta_mode_, sigma = beta_cov_, log = TRUE)
            
            # Compute log-posterior values
            cr_logpost_beta <- logpost_beta(beta = beta, z = z, X = X, atoms = z.tld, 
                                            jumps  = J.tld, mu_beta = mu_beta, 
                                            sigma_beta = sigma_beta)
            pr_logpost_beta <- logpost_beta(beta = beta_, z = z, X = X, atoms = z.tld, 
                                            jumps  = J.tld, mu_beta = mu_beta, 
                                            sigma_beta = sigma_beta)
            
            # Metropolis-Hastings acceptance step
            log_acc_prob <- pr_logpost_beta - cr_logpost_beta + cr_logprop_beta - pr_logprop_beta
            if (log(runif(1)) < log_acc_prob) {
              count1 <- count1 + 1
              beta <- beta_
              beta_mode <- beta_mode_
              beta_cov <- beta_cov_
              meanY_x <- mean_z <- mean_z_
            }
          }
        } 
      }
    }
    # If `optim()` fails, keep beta as it is and continue to next step.
    
    
    
    
    # theta update ------------------------------------
    theta_tilde <- gldrm:::getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = meanY_x,
      sampprobs  = NULL,
      ySptIndex  = NULL,
      thetaStart = theta
    )$theta
    
    theta <- theta_tilde
    
    # u update ----------------------------------------
    u <- sampler_u(u, zstar, nstar, theta, alpha, delta, min_y, max_y)
    
    # CRM update --------------------------------------
    crm_star <- crm_sampler(M, u, zstar, nstar, theta, alpha, min_y, max_y)
    z.tld_star <- c(crm_star$RL, crm_star$zstar)
    J.tld_star <- c(crm_star$RJ, crm_star$Jstar)
    
    if(min(meanY_x) >= min(z.tld_star) && max(meanY_x) <= max(z.tld_star)){
      # MH step
      theta_tilde_star <- gldrm:::getTheta(
        spt = z.tld_star,
        f0  = J.tld_star,
        mu  = meanY_x,
        sampprobs  = NULL,
        ySptIndex  = NULL,
        thetaStart = theta
      )$theta
      
      b1 <- b_theta(theta_tilde_star, z.tld_star, J.tld_star)
      b2 <- b_theta(theta_tilde, z.tld, J.tld)
      b3 <- b_theta(theta_tilde_star, z.tld, J.tld)
      b4 <- b_theta(theta_tilde, z.tld_star, J.tld_star)
      log_r <- log(exp(sum(2*(theta_tilde_star - theta_tilde)*z - b1 + b2 - b3 + b4)))
      
      if(log(runif(1)) < log_r){
        count2 <- count2 + 1
        RL <- crm_star$RL
        RJ <- crm_star$RJ
        zstar <- crm_star$zstar
        Jstar <- crm_star$Jstar
        z.tld <- c(RL, zstar)
        J.tld <- c(RJ, Jstar)
        theta <- theta_tilde_star
      }
    }
    

    # z update ------------------------------------------------------------------
    z <- z_sampler_unifK(y, c0, z.tld, J.tld, theta, min_y, max_y)
    
    # zstar and nstar update ----------------------------------------------------
    resampled_z <- resample_zstar(z)
    zstar <- resampled_z$zstar
    nstar <- resampled_z$nstar
    
    theta0 <- gldrm:::getTheta(
      spt = z.tld,
      f0  = J.tld,
      mu  = m0,
      sampprobs  = NULL,
      ySptIndex  = NULL,
      thetaStart = NULL
    )$theta
    
    temp <- exp(theta0 * z.tld - max(theta0 * z.tld))
    Jtld_0 <- (temp * J.tld) / sum(temp * J.tld)
    lnlik <- loglik(z = z, X = X, beta = beta, atoms = z.tld, jumps = Jtld_0)
    
    # Storing MCMC simulations --------------------------------------------------
    z_samples[itr, ] <- z
    u_samples[itr, ] <- u
    beta_samples[itr,] <- beta
    theta_samples[itr, ] <- theta
    crm_samples[[itr]] <- list(z.tld = z.tld, J.tld = J.tld)
    lnlik_samples[itr] <- lnlik
  }
  
  dpglm_fit <- list(z = z_samples, beta = beta_samples, 
                    crm = crm_samples, 
                    beta_acceptance = count1 / iter,
                    crm_acceptance = count2 / iter)
  
  
  out <- list(data = dat, dpglm = dpglm_fit)
  return(out)
}
