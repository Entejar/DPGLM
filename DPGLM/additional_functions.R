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
unif_kde_cdf <- function(y, crm.atoms, crm.jumps, theta, c0,
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
# Compute per-group Gaussian copula log-likelihood contributions,
# with observation-specific marginal CDFs constructed internally.
#
# For each observation j, the marginal CDF F_{x_j}(y_j) is computed
# via a CRM-induced mixture of truncated uniform kernels using
# unif_kde_cdf(). These marginal probabilities are then coupled
# within groups using a Gaussian copula with exchangeable correlation.
#
# For each group i with n_i observations, the function returns:
#   (1 / n_i) * log c_rho(F_{i1}, ..., F_{in_i})
#
# Arguments:
#   y           : numeric vector of outcomes
#   group_index : vector indicating group membership for each y
#   rho         : copula dependence parameter
#   crm.atoms   : vector of CRM atoms (z_l)
#   crm.jumps   : vector of CRM jump sizes (J_l)
#   theta     : numeric vector of covariate effects (length = length(y))
#   c0          : half-width of the uniform kernel
#   min_y       : global lower bound of the outcome support
#   max_y       : global upper bound of the outcome support
#
# Returns:
#   out : named numeric vector of per-group copula log-likelihood
#         contributions, normalized by group size
# ------------------------------------------------------------------
copula_contribution_by_group <- function(y, group_index, rho,
                                         crm.atoms, crm.jumps, theta, c0,
                                         min_y, max_y) {
  ## 1. Compute marginal CDFs internally
  Fy <- unif_kde_cdf(y = y,
                     crm.atoms = crm.atoms, crm.jumps = crm.jumps, theta = theta, c0 = c0,
                     min_y = min_y, max_y = max_y)
  
  ## 2. Per-group copula contributions
  groups <- unique(group_index)
  out <- numeric(length(groups))
  names(out) <- groups
  
  for (g in groups) {
    idx <- which(group_index == g)
    u <- Fy[idx]
    n_i <- length(u)
    
    logc <- log_copula_gaussian(u, rho)
    out[as.character(g)] <- logc / n_i
  }
  
  out
}

