beta_sampler <- function(y, X, z.tld, J.tld, zMN, zMX, beta, Sig, mubetaprior,
                         Sigbetaprior, c0, sigma_theta){
  B <- 1000     # number of samples for the Monte Carlo integration
  n <- dim(X)[1]
  p <- dim(X)[2]
  crm <- matrix(c(z.tld, J.tld), nrow = 2, byrow = T)
  cr_theta <- theta_solver(locations = crm[1,], jumps = crm[2,], meanY_x = as.numeric(expit(X %*% beta)), thetastart = NULL)$theta
  
  # proposed beta
  pr_bt <- as.vector(mvrnormArma(1, beta, Sig)) 
  # proposed mu
  pr_mu <- as.numeric(expit(X %*% pr_bt))
  # proposed theta
  pr_theta <- theta_solver(locations = crm[1,], jumps = crm[2,], meanY_x = pr_mu, thetastart = cr_theta)$theta
  
  # check if proposed mu is within bounds
  if(sum(zMN <= pr_mu & pr_mu <= zMX) == n){
    # log likelihood
    pr_llik <- llik_beta(y, X, pr_theta, crm, c0, B, sigma_theta)
    cr_llik <- llik_beta(y, X, cr_theta, crm, c0, B, sigma_theta)
    # log prior
    pr_logpr <- logpdf_mvnorm(pr_bt, mubetaprior, Sigbetaprior)
    cr_logpr <- logpdf_mvnorm(beta, mubetaprior, Sigbetaprior)
    # log q(beta | beta_proposed, Sig) = log q(beta_proposed | beta, Sig) 
    # hence, cr_logq - pr_logq = 0 => omitted from the logratio
    
    
    # log acceptance probability
    logratio <- min(0, pr_llik - cr_llik + pr_logpr - cr_logpr) 

    if(log(runif(1)) < logratio){
      beta <- pr_bt
    }
  }
  
  return(beta)
}


