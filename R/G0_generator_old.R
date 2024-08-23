## ----------------------------------------------------------------------------- ##
#     This is an old file for sampling from G_0, which is now coded in C++      #
## ----------------------------------------------------------------------------- ##



# G0_generator <- function(G0.dist, mu0, sigma0, a00, b00){
#   R <- 25000
# 
#   if(G0.dist == 1){
#     z    <- rnorm(R, mean = mu0, sd = sigma0)
#   }
# 
#   if(G0.dist == 2){
#     z    <- rtruncnorm(R, mean = mu0, sd = sigma0, lower = a00)
#   }
# 
#   if(G0.dist == 3){
#     z    <- rgamma(R, shape = mu0^2 / sigma0^2, rate = mu0 / sigma0^2)
#   }
# 
#   if(G0.dist == 4){
#     nu <- (mu0 * (1 - mu0) / (sigma0^2) - 1)
#     a <- mu0 * nu
#     b <- (1 - mu0) * nu
#     z    <- rbeta(R, shape1 = a00, shape2 = b00)
#   }
# 
#   if(G0.dist == 5){
#     z    <- rtruncnorm(R, mean = mu0, sd = sigma0, a = a00, b = b00)
#   }
# 
#   if(G0.dist == 6){
#     z    <- runif(R, min = a00, max = b00)
#   }
# 
#   return(z)
# 
# }
# 

