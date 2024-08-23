crm_sampler <-function(M, u, zstar, nstar, tht, n, alpha, G0.dist, a00, b00, mu, z.tld, J.tld)
{
  N <- 3001
  w <- -log(seq(from = exp(-1e-06), to = exp(-5e-04), length = N))
  z.G0 <- sort(G0_generator(G0.dist, mu0G, sigma0G, a00, b00))
  R <- length(z.G0)
  mt <- matrix(NA, nrow = n, ncol = R)
  for(i in 1:R){
    mt[, i] <- tht * z.G0[i] %>% as.numeric()
  }
  mx <- apply(mt, 1, max)
  psi <- function(z){
    z <- rep(z, length(tht))
    sum(u * exp(tht * z - mx))
  }
  psi.z <- sapply(z.G0, psi)

  # First Part: random jumps
  fnS <- sapply(w, function(s) mean(exp( - (1 + psi.z) * s ) / s))
  dw <- diff(w)
  h <- (fnS[-1] + fnS[-N])/2
  Nv <- c(rev(cumsum(rev(dw * h))), 0)
  Nv <- alpha * Nv
  xi <- cumsum(rexp(M))
  iNv <- N
  J <- sapply(seq(M), function(i) {
    while (iNv > 0 && Nv[iNv] < xi[i]) {
      iNv <- iNv - 1
    }
    w[iNv + 1]
  })

  # First Part: random locations
  z <- sapply(seq(M), function(m) {
    xi <- runif(1)
    temp <- exp( - (1 + psi.z) * J[m] )
    cutoff <- sum(temp) * xi
    z.G0[min(which(cumsum(temp) - cutoff > 0))]
  })

  # Second Part: random jumps [for fixed locations]
  psi.zstar <- colSums(matrix(rep(u, length(zstar)) * exp(rep(tht, length(zstar)) * zstar), nrow = n))
  Jstar <- rgamma(length(zstar), shape = nstar, rate = psi.zstar + 1)

  z.tld_temp <- c(z, zstar)
  if(min(mu) >= min(z.tld_temp) && max(mu) <= max(z.tld_temp)){
    z.tld <- z.tld_temp
    J.tld <- c(J, Jstar)
  } else{
    z.tld <- z.tld
    J.tld <- J.tld
  }

  temp1 <- min(y) %in% z.tld
  temp2 <- max(y) %in% z.tld
  if(temp1){
    z.tld <- z.tld
    J.tld <- J.tld} else{
      z.tld <- c(min(y), z.tld)
      J.tld <- c(min(J.tld), J.tld)
    }
  if(temp2){
    z.tld <- z.tld
    J.tld <- J.tld} else{
      z.tld <- c(z.tld, max(y))
      J.tld <- c(J.tld, min(J.tld))
    }

  return(list(z.tld = z.tld, J.tld = J.tld))
}

