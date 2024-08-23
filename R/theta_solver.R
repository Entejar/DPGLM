theta_solver <- function(z.tld, J.tld, mu, thtst){
  out <- gldrm:::getTheta(spt = z.tld, f0 = J.tld, mu = mu, sampprobs = NULL,
                          ySptIndex = NULL, thetaStart = thtst)
  tht <- out$theta
  bpr2 <- out$bPrime2
  btht <- apply(exp(outer(tht, z.tld, "*")), 1, function(row) log(sum(row * J.tld)))
  return(list(bpr2 = bpr2, tht = tht, btht = btht))
}




