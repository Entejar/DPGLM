## ----------------------------------------------------------------------------- ##
#     This is an old file for log(K(y | z)), which is now coded in C++           #
## ----------------------------------------------------------------------------- ##


# logK <- function(y, dist, z, sigma0, c0, mn, mx){
#   logK <- numeric(length(y))
#   
#   if(dist == "unif"){
#     a <- z - c0
#     b <- z + c0
# 
#     logK    <- dunif(y, min = a, max = b, log = T)
#   }
# 
#   if(dist == "restricted_unif"){
#     c00 <- apply(data.frame(rep(c0, length(z)), z - mn, mx - z), 1, min)
#     a <- z - c00
#     b <- z + c00
# 
#     logK    <- dunif(y, min = a, max = b, log = T)
#   }
# 
#   return(logK)
# }
