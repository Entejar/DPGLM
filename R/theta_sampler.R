## NOTE: theta_tilde = free parameters, theta = derived parameters

theta_tilde_sampler <- function(meanY_x,
                                sigma_theta,
                                z,
                                locations,
                                jumps) {
  theta <-  gldrm:::getTheta(
    spt = locations,
    f0  = jumps,
    mu  = meanY_x,
    sampprobs  = NULL,
    ySptIndex  = NULL,
    thetaStart = thetastart
  )$theta
  
  theta_tilde <- as.vector(mvrnormArma(1, theta + z * sigma_theta^2, 
                                       diag(length(z)) * sigma_theta^2))
  
  return(theta_tilde)
}
