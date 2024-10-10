theta_solver <- function(locations, jumps, meanY_x, thetastart) {
  out <- gldrm:::getTheta(
    spt = locations,
    f0  = jumps,
    mu  = meanY_x,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = thetastart)
  theta   <- out$theta
  bprime2 <- out$bPrime2
  btheta <- apply(exp(outer(theta, locations, "*")), 1, function(row)
    log(sum(row * jumps)))
  
  return(list(bpr2 = bprime2, 
              theta = theta, 
              btht = btheta))
}




