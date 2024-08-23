## ----------------------------------------------------------------------------- ##
#     This is an old file for sampling u, which is now coded in C++, see u.cpp   #
## ----------------------------------------------------------------------------- ##



# u_sampler <-function(u_current, theta_current, z_current, 
#                      n, alpha, G0.dist, mu0, sigma0, delta, a00, b00)
#   {
#   z <- G0_generator(G0.dist, mu0, sigma0, a00, b00)
#   R <- length(z)
# 
#   u_proposed <- rgamma(n, shape = delta, rate = delta/u_current)
#   tiltA      <- exp(theta_current * matrix(rep(z, n), byrow = T, nrow = n, ncol = R))
#   A          <- mean(log((1 + crossprod(u_proposed, tiltA)) / (1 + crossprod(u_current, tiltA))))
#   tiltB      <- exp(theta_current * matrix(rep(z_current, n), byrow = T, nrow = n, ncol = n))
#   B          <- sum(log(1 + crossprod(u_proposed, tiltB)) / (1 + crossprod(u_current, tiltB)))
#   D          <- dgamma(u_current, shape = delta, rate = delta/u_proposed, log = T) %>% sum()
#   E          <- dgamma(u_proposed, shape = delta, rate = delta/u_current, log = T) %>% sum()
#   logratio   <- - alpha * A  - B + D - E
# 
#   p          <- min(0, logratio)
#   if(log(runif(1)) <= p){
#     u_current <- u_proposed
#   }
#   u_current
# }

