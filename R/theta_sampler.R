theta_sampler <- function(X, beta, sigmaTheta, z, link, z.tld, J.tld){
  n <- nrow(X)
  mu <- link_fn(X %*% beta, link)$mu
  tht.tld <- theta_solver(z.tld, J.tld, # / sum(J.tld),
                          mu, NULL)$tht
  sigma2Theta <- (sigmaTheta)^2
  muTheta <- tht.tld + (z * sigma2Theta)
  tht <- rnorm(n, mean = muTheta, sd = sigmaTheta)
  return(tht)
}
