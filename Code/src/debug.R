dpglm <- function(y, X, iter, tuning.params, set.beta = TRUE){
  
  # Link Function -----------------------------------------------------------------------
  link <- 'logit'
  
  # Tuning Parameters --------------------------------------------------------------------
  rho <- tuning.params$rho
  M <- tuning.params$M
  alpha <- tuning.params$alpha
  delta <- tuning.params$delta
  sigma_theta <- tuning.params$sigma_theta
  c0 <- tuning.params$c0
  beta.sigma <- tuning.params$beta.sigma
  
  # Data Preparation ---------------------------------------------------------------------
  X <- X %>% as.matrix()
  y <- y %>% as.numeric()
  mu0 <- mean(y)
  n <- length(y)
  p <- dim(X)[2]
  
  ## Initialization -----------------------------------------------------------------------
  beta_samples <- matrix(NA, nrow = iter, ncol = p)
  theta_samples <- matrix(NA, nrow = iter, ncol = n)
  z_samples <- matrix(NA, nrow = iter, ncol = n)
  crm_samples   <- list()
  init <- gldrm(y ~ X[, -1], link = link)
  beta0 <- log(mu0 / (1 - mu0))  
  # if(set.beta == TRUE){
  #   beta_samples[1, ] <- beta <- c(beta0, rep(0, p-1))
  # } else {
  #   beta_samples[1, ] <- beta <- init$beta %>% as.numeric()
  # }
  beta_samples[1, ] <- beta <- c(1, 0, 0)
  # z.tld <- spt <- init$spt
  # J.tld <- Jstar <- init$f0
  z.tld <- spt
  J.tld <- Jstar <- f0_kde
  crm_samples[[1]] <- list(z.tld = spt, J.tld = Jstar)
  mu <- expit(X %*% beta)
  out <- theta_solver(locations = spt, jumps = J.tld, meanY_x = mu, thetastart = NULL)
  theta_samples[1, ] <- tht <- out$theta
  btht <- out$btht
  bpr2 <- out$bpr2
  z_samples[1, ] <- z <- y
  temp <- resample_zstar(z)
  zstar <- temp$zstar
  nstar <- temp$nstar
  
  T.vec <- exp(btht)
  u <- rgamma(n, shape = 1, rate = T.vec)
  mubetaprior <-  rep(0, p)
  Sigbetaprior <- beta.sigma * diag(p)
  Sig <- rho * init$varbeta
  itr <- 1
  for(itr in 2:iter){
    itr <- itr + 1
    # beta update -------------------------------------
    zMN <- min(z.tld)
    zMX <- max(z.tld)
    beta <- beta_sampler(y, X, z.tld, J.tld, zMN, zMX, beta, Sig, mubetaprior, Sigbetaprior, c0, sigma_theta) 
    meanY_x <- as.numeric(expit(X %*% beta))
    
    # theta_tilde update ------------------------------------
    theta_tilde <- theta_tilde_sampler(meanY_x, sigma_theta, z, locations = z.tld, jumps = J.tld, thetastart = tht)
    
    # u update ----------------------------------------
    u <- u_sampler(u, z, tht, alpha, delta)
    
    # CRM update --------------------------------------
    crm <- crm_sampler(M, u, zstar, nstar, tht, alpha, mu, y, z.tld, J.tld)
    z.tld <- crm$z_tld
    J.tld <- crm$J_tld
    
    # z update ------------------------------------------------------------------
    sorted_crm <- data.frame(z.tld = z.tld, J.tld = J.tld) %>% arrange(z.tld)
    z.tld <- sorted_crm[, 1]
    J.tld <- sorted_crm[, 2]
    z <- z_sampler(y, c0, tht, sorted_crm %>% as.matrix())
    
    # zstar and nstar update ----------------------------------------------------
    resampled_z <- resample_zstar(z)
    zstar <- resampled_z$zstar
    nstar <- resampled_z$nstar
    
    # Storing MCMC simulations --------------------------------------------------
    z_samples[itr, ] <- z
    beta_samples[itr,] <- beta
    theta_samples[itr, ] <- tht
    crm_samples[[itr]] <- list(z.tld = z.tld, J.tld = J.tld)
    
    # Diagnostics
    #par(mfrow = c(3, 1))
    # plot(beta_samples[1:itr, 1], type = "l", col = "blue", lwd = 2, xlab = "Iteration", ylab = "Beta0")
    # plot(beta_samples[1:itr, 2], type = "l", col = "blue", lwd = 2, xlab = "Iteration", ylab = "Beta1")
    # plot(beta_samples[1:itr, 3], type = "l", col = "blue", lwd = 2, xlab = "Iteration", ylab = "Beta2")
    hist(z.tld, breaks = 20, col = "blue", xlab = "Z", main = "CRM")
  
    
    }
  return(list(z = z_samples, beta = beta_samples, theta = theta_samples, crm = crm_samples))
}
