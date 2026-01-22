# ------------------------------------------------------------------
# Compute observation-specific marginal CDFs F_{x_j}(y_j) using a
# CRM-induced mixture of truncated uniform kernels.
#
# The mixture weights are covariate-dependent and vary across
# observations through theta[j].
#
# Arguments:
#   y         : numeric vector of outcomes (length n)
#   crm.atoms : vector of CRM atoms (z_l)
#   crm.jumps : vector of CRM jump sizes (J_l)
#   theta   : numeric vector of covariate effects (length n)
#   c0        : half-width of the uniform kernel
#   min_y     : global lower bound of the outcome support
#   max_y     : global upper bound of the outcome support
#
# Returns:
#   Fy        : numeric vector of marginal CDF values F_{x_j}(y_j)
# ------------------------------------------------------------------
marginal_cdf_unif_kde <- function(y, crm.atoms, crm.jumps, theta, c0,
                         min_y, max_y) {
  
  eps <- 1e-6
  n <- length(y)
  
  ## 1. Construct truncated kernel bounds (shared across observations)
  lower <- crm.atoms - c0
  upper <- crm.atoms + c0
  
  lower[lower < min_y] <- min_y + eps
  upper[upper > max_y] <- max_y - eps
  
  width <- upper - lower
  
  ## 2. Evaluate CDF observation by observation
  Fy <- numeric(n)
  
  for (j in seq_len(n)) {
    
    ## Observation-specific mixture weights (log-scale)
    logw <- log(crm.jumps) + theta[j] * crm.atoms
    logw <- logw - max(logw)   # log-sum-exp stabilization
    w <- exp(logw)
    w <- w / sum(w)
    
    ## Mixture CDF at y_j
    yj <- y[j]
    Fy[j] <- sum(
      w * ifelse(
        yj < lower, 0,
        ifelse(
          yj > upper, 1,
          (yj - lower) / width
        )
      )
    )
  }
  
  Fy
}



# ------------------------------------------------------------------
# Compute the log-density of a Gaussian copula with an exchangeable
# (compound symmetry) correlation structure.
#
# The copula density is evaluated at a vector of marginal probabilities
# u = (u_1, ..., u_n), where u_j = F_x(y_ij).
#
# Arguments:
#   u   : numeric vector of marginal CDF values in (0, 1)
#   rho : dependence parameter (pairwise correlation)
#
# Returns:
#   logc : scalar value of log c_rho(u)
# ------------------------------------------------------------------
log_copula_gaussian <- function(u, rho) {
  n <- length(u)
  
  # avoid boundary
  eps <- 1e-10
  u <- pmin(pmax(u, eps), 1 - eps)
  
  ## Transform to latent Gaussian scale
  z <- qnorm(u)
  
  # Exchangeable correlation structure
  Sigma <- matrix(rho, n, n)
  diag(Sigma) <- 1
  
  # # AR(1) correlation
  # Sigma <- outer(1:n, 1:n, function(i, j) rho^abs(i - j))
  
  Sigma_inv <- solve(Sigma)
  logdet <- determinant(Sigma, logarithm = TRUE)$modulus
  
  ## Gaussian copula log-density
  # log c(u) = -1/2 log|Sigma| - 1/2 z' (Sigma^{-1} - I) z
  val <- -0.5 * logdet -
    0.5 * t(z) %*% (Sigma_inv - diag(n)) %*% z
  
  as.numeric(val)
}


# ------------------------------------------------------------------
# Compute observation-level Gaussian copula log-likelihood contributions.
#
# Marginal CDFs F_j = F_{x_j}(y_j) are computed internally via a
# CRM-induced mixture of truncated uniform kernels using
# marginal_cdf_unif_kde().
#
# Observations are coupled within groups using a Gaussian copula with
# exchangeable correlation parameter rho.
#
# For a group i with n_i observations, each observation j in that group
# receives an equal share of the group copula log-density:
#
#   log c_rho(F_{i1}, ..., F_{in_i}) / n_i
#
# Arguments:
#   y           : numeric vector of outcomes
#   group_index : vector indicating group membership for each y
#   rho         : copula dependence parameter
#   crm.atoms   : vector of CRM atoms (z_l)
#   crm.jumps   : vector of CRM jump sizes (J_l)
#   theta       : numeric vector of covariate effects (length = length(y))
#   c0          : half-width of the uniform kernel
#   min_y       : global lower bound of the outcome support
#   max_y       : global upper bound of the outcome support
#
# Returns:
#   logcop_obs  : numeric vector of observation-level copula
#                 log-likelihood contributions (same length as y)
# ------------------------------------------------------------------
log_copula_contribution_by_obs <- function(y, group_index, rho,
                                           crm.atoms, crm.jumps, theta, c0,
                                           min_y, max_y) {
  
  ## 1. Compute marginal CDFs internally
  Fy <- marginal_cdf_unif_kde(y = y,
                     crm.atoms = crm.atoms, crm.jumps = crm.jumps, theta = theta, c0 = c0,
                     min_y = min_y, max_y = max_y)
  
  ## 2. Initialize observation-level output
  logcop_obs <- numeric(length(y))
  
  ## 3. Compute copula contribution per group, then distribute to observations
  groups <- unique(group_index)
  
  for (g in groups) {
    idx <- which(group_index == g)
    u <- Fy[idx]
    n_i <- length(u)
    
    logc_group <- log_copula_gaussian(u, rho)
    
    # Assign equal share to each observation in the group
    logcop_obs[idx] <- logc_group / n_i
  }
  
  logcop_obs
}





# ------------------------------------------------------------------
# Log posterior kernel for the copula dependence parameter rho
# (simplified interface).
#
# This function evaluates the log posterior of rho under the
# simplified joint distribution
#
#   pi(rho | y, lambda, beta) ∝ ∏_i c(y_i, lambda, beta, rho) p(rho),
#
# where all marginal likelihood terms are handled elsewhere and the
# input 'logc_obs' already contains the log copula density contributions
# evaluated at the current value of rho.
#
# Arguments:
#   rho      : scalar copula dependence parameter, must lie in (0, 1)
#   logc_obs : numeric vector of log copula contributions whose sum
#              equals log ∏_i c(y_i, lambda, beta, rho)
#
# Returns:
#   Scalar log posterior value log pi(rho | ·), up to an additive constant
# ------------------------------------------------------------------
log_post_rho <- function(rho, logc_obs) {
  
  ## Enforce support of the copula parameter
  if (rho <= 0 || rho >= 1) return(-Inf)
  
  ## Copula log-likelihood contribution
  loglik <- sum(logc_obs)
  
  ## Beta(8, 2) prior on rho
  ## We probably need to make the prior on rho as an argument
  logprior <- dbeta(rho, shape1 = 8, shape2 = 2, log = TRUE)
  
  ## Log posterior kernel
  loglik + logprior
}
