dpglm <- function(y, X, link, iter, tuning.params){

  # Tuning Parameters --------------------------------------------------------------------
  rho <- tuning.params$rho
  M <- tuning.params$M
  alpha <- tuning.params$alpha
  G0.dist <- tuning.params$G0.dist
  delta <- tuning.params$delta
  K.dist <- tuning.params$K.dist
  sigmaTheta <- tuning.params$sigmaTheta
  a00 <- tuning.params$a00
  b00 <- tuning.params$b00
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
  beta0 <- g(mu0)
  if(n <= 250){
    beta <- beta_samples[1,] <- c(beta0, rep(0, p - 1))
  } else{
    beta <- beta_samples[1,] <- c(beta0, init$beta[2:p] %>% as.numeric())
  }

  z.tld <- spt <- init$spt
  J.tld <- Jstar <- init$f0
  crm_samples[[1]] <- list(z.tld = spt, J.tld = Jstar)
  mu <- link_fn(X %*% beta, link)$mu
  out <- theta_solver(spt, J.tld, mu, NULL)
  theta_samples[1, ] <- tht <- out$tht
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
  mu0G <- -99
  sigma0G <- -99

  for(itr in 2:iter){
    # beta update -------------------------------------
    zMN <- min(z.tld)
    zMX <- max(z.tld)
    beta <- beta_sampler(y, X, z.tld, J.tld, zMN, zMX, beta, Sig, mubetaprior, Sigbetaprior, link, c0,
                         sigmaTheta[1])
    mu <- link_fn(X %*% beta, link)$mu
    # theta update ------------------------------------
    tht <- theta_sampler(X, beta, sigmaTheta, z, link, z.tld, J.tld)

    # u update ----------------------------------------
    u <- u_sampler(u, tht, z, n, alpha, G0.dist, mu0G, sigma0G, delta, a00, b00)

    # crm update --------------------------------------
    crm <- crm_sampler(M, u, zstar, nstar, tht, n, alpha, G0.dist, a00, b00, mu, z.tld, J.tld)
    z.tld <- crm$z.tld
    J.tld <- crm$J.tld

    # z update ------------------------------------------------------------------
    sorted_crm <- data.frame(z.tld = z.tld, J.tld = J.tld) %>% arrange(z.tld)
    z.tld <- sorted_crm[, 1]
    J.tld <- sorted_crm[, 2]
    z <- z_sampler(y, n, c0, kdist, z.tld, J.tld, tht)

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
    par(mfrow = c(3, 1))
    plot(beta_samples[1:itr, 1], type = "l", col = "blue", lwd = 2, xlab = "Iteration", ylab = "Beta0")
    plot(beta_samples[1:itr, 2], type = "l", col = "blue", lwd = 2, xlab = "Iteration", ylab = "Beta1")
    plot(beta_samples[1:itr, 3], type = "l", col = "blue", lwd = 2, xlab = "Iteration", ylab = "Beta2")
  }
  return(list(z = z_samples, beta = beta_samples, theta = theta_samples, crm = crm_samples))
}
