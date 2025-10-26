
## ----------------------------------------------------------------------------- ##
#     This is an old file for sampling z, which is now coded in C++, see z.cpp   #
## ----------------------------------------------------------------------------- ##



# z_sampler <- function(y, c0, K.dist, crm.atoms, crm.jumps, tht) {
#   n <- length(y)
#   mn <- min(y)
#   mx <- max(y)
#   z <- numeric(n)
#   # eps <- 0.000001        # for K.dist = 7
#   for (i in 1:n) {
#     indx <- crm.atoms > (y[i] - c0) & crm.atoms < (y[i] + c0)
#     #& crm.atoms > mn + eps &
#     #  crm.atoms < mx - eps   # last two AND statements are for Kdist = 7 but remove when Kdist = 6
#     log.prob <- logK(
#       y[i],
#       dist = kdist,
#       mu = crm.atoms[indx],
#       sigma = -99,
#       deci = T,
#       -99,
#       -99,
#       c0,
#       mn,
#       mx
#     ) +
#       log(crm.jumps[indx]) + (tht[i] * crm.atoms[indx])
#     prob <- exp(log.prob - max(log.prob))
#     z[i] <- sample(crm.atoms[indx], 1, prob = prob)
#   }
#   return(z)
# }
# 
