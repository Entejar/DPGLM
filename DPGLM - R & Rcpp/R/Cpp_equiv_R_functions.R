# --------------------------------------------------------------------------- #
#                          For beta update                                   #
# --------------------------------------------------------------------------- #

#' Compute b(theta) for discrete base measure
#'
#' Calculates \eqn{b(\theta) = \log \left( \int \exp(\theta z) dG_0(z) \right)} for a discrete base measure \eqn{G_0}
#' with support points \code{spt} and weights \code{f0}, using log-sum-exp trick for numerical stability.
#'
#' @param theta Numeric vector of parameter values.
#' @param spt Numeric vector of support points for the base measure.
#' @param f0 Numeric vector of weights or density values at each support point.
#'
#' @return Numeric vector of \eqn{b(\theta)} values, one for each \code{theta}.
#' @export

b_theta <- function(theta, spt, f0) {
  log_f0 <- log(f0) # log density values
  
  # Create the log weights matrix of dimensions length(theta) x length(spt)
  # (i, j)th element is theta[i]*spt[j] + log(f0[j])
  log_weights <- outer(theta, spt, "*") +
    matrix(log_f0, nrow = length(theta), ncol = length(spt), byrow = TRUE)
  
  # Use log-sum-exp trick to compute b(theta) for each theta
  result <- apply(log_weights, 1, function(x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
  })
  
  return(as.vector(result))  # Ensure the output is a vector
}

# --------------------------------------------------------------------------- #
#                               CRM update                                   #
# --------------------------------------------------------------------------- #

crm_sampler <- function(M, u, zstar, nstar, tht, alpha, min_y, max_y){
  N <- 3001
  R <- 3001
  eps <- 1e-6
  s <- -log(seq(exp(-eps), exp(-5e-4), length.out = N))
  
  # Sorted, ascending order, needed in RL
  #z <- seq(eps, 1 - eps, length.out = R)
  z <- seq(min_y + eps, max_y - eps, length.out = R)
  
  # Assume u is a vector and z is the grid where you want to compute psi_z.
  # Here we compute psi_z for each grid point z[j] by summing over u.
  
  # Compute a matrix of log-terms:
  log_terms <- outer(log(u), z, function(lu, z_val) lu + tht * z_val)
  
  # Now compute psi_z using the log-sum-exp trick along the u dimension (rows):
  psi_z <- apply(log_terms, 2, function(log_vec) {
    max_val <- max(log_vec)
    exp(max_val + log(sum(exp(log_vec - max_val))))
  })
  
  
  # Assume:
  #   psi_z: vector of length R computed as before over the z grid
  #   s: vector of length N representing the s grid
  #   eps: small positive number for numerical boundaries
  dz <- (1 - 2 * eps) / (R - 1)   # Grid spacing for the uniform measure
  
  # Preallocate the result vector
  fnS <- numeric(length(s))
  
  for (j in seq_along(s)) {
    # For a fixed s[j], the integrand at each z is:
    #    f(z, s[j]) = exp( -(1+psi_z) * s[j] ) / s[j]
    # Taking logs gives:
    #    log f(z, s[j]) = -(1+psi_z) * s[j] - log(s[j])
    log_vals <- -(1 + psi_z) * s[j] - log(s[j])
    
    # Use the log-sum-exp trick for the summation over z grid:
    max_val <- max(log_vals)
    log_sum_exp <- max_val + log(sum(exp(log_vals - max_val)))
    
    # Multiply by the grid spacing to approximate the integral:
    log_integral <- log(dz) + log_sum_exp
    
    # Store the stabilized value (exponentiating back)
    fnS[j] <- exp(log_integral)
  }
  
  
  ds <- diff(s)
  h <- (fnS[-N] + fnS[-1]) / 2
  Nv <- rev(cumsum(rev(ds * h)))
  Nv <- Nv * alpha
  Nv <- c(Nv, 0)
  
  
  #Generate random jumps RJ and random locations RL
  xi <- cumsum(rexp(M, rate = 1.0))
  RJ <- numeric(M)
  iNv <- N - 1
  
  for (i in seq_len(M)) {
    while (iNv > 0 && Nv[iNv] < xi[i]) {
      iNv <- iNv - 1
    }
    RJ[i] <- s[iNv + 1]
  }
  
  
  RL <- numeric(M)
  for (m in seq_len(M)) {
    xi_rl <- runif(1)
    # Compute log probabilities: log_temp = -(1 + psi_z) * RJ[m]
    log_temp <- -(1 + psi_z) * RJ[m]
    
    # Stabilize by subtracting the maximum log value
    max_log <- max(log_temp)
    
    # Convert to probabilities in a numerically stable way
    p <- exp(log_temp - max_log)
    p <- p / sum(p)
    
    # Compute the cumulative probabilities
    cumsum_p <- cumsum(p)
    
    # Select the index where the cumulative sum exceeds xi_rl
    RL[m] <- z[min(which(cumsum_p > xi_rl))]
  }
  
  
  # Second Part: random jumps [for fixed locations]
  # Suppose u is a vector and zstar is a vector.
  # We want to compute, for each fixed zstar[j]:
  #   psi_star[j] = sum(u_i * exp(tht * zstar[j]))
  # We'll compute it in a numerically stable way:
  
  # Create a matrix of log-terms: each element is log(u_i) + tht * zstar[j]
  log_terms <- outer(log(u), zstar, function(lu, z) lu + tht * z)
  
  # For each column (corresponding to a fixed zstar), use the log-sum-exp trick
  psi_star <- sapply(seq_along(zstar), function(j) {
    max_val <- max(log_terms[, j])
    exp(max_val + log(sum(exp(log_terms[, j] - max_val))))
  })
  
  
  Jstar <- rgamma(length(zstar), shape = nstar, rate = psi_star + 1)
  
  return(list(
    RL = RL,
    RJ = RJ,
    zstar = zstar,
    Jstar = Jstar
  ))
}
