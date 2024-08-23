beta_sampler <- function(y, X, z.tld, J.tld, zMN, zMX, beta, Sig, mubetaprior,
                         Sigbetaprior, link, c0, sigma_theta){
  B <- 1000
  n <- dim(X)[1]
  p <- dim(X)[2]
  crm <- matrix(c(z.tld, J.tld), nrow = 2, byrow = T)
  # Beta update --------------------------------
  pr_bt <- mvtnorm::rmvnorm(1, mean = beta, sigma = Sig) %>% as.vector()
  pr_mu <- link_fn(X %*% pr_bt, link)$mu %>% as.numeric()
  pr_theta <- theta_solver(z.tld, J.tld, pr_mu, NULL)$tht
  cr_theta <- theta_solver(z.tld, J.tld, link_fn(X %*% beta, link)$mu, NULL)$tht
  if(sum(zMN <= pr_mu & pr_mu <= zMX) == n){
    pr_llik <- llik_beta(y, X, pr_theta, crm, c0, B, sigma_theta)
    cr_llik <- llik_beta(y, X, cr_theta, crm, c0, B, sigma_theta)
    pr_pbt <- dmvnorm(pr_bt, mean = mubetaprior, sigma = Sigbetaprior, log = T)
    cr_pbt <- dmvnorm(beta, mean = mubetaprior, sigma = Sigbetaprior, log = T)
    cr_qbt <- dmvnorm(beta, mean = pr_bt, sigma = Sig, log = T)
    pr_qbt <- dmvnorm(pr_bt, mean = beta, sigma = Sig, log = T)

    alp <- min(0, pr_llik - cr_llik + pr_pbt - cr_pbt + cr_qbt - pr_qbt)

    if((log(runif(1)) < alp)){
      beta <- pr_bt
    }
  }

  return(beta)
}


