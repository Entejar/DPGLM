get_F0 <- function(out, save, yGrid, mu0, c0) {
  itr_indx <- save
  F0_est <- matrix(0, nrow = length(yGrid), ncol = length(itr_indx))
  for (j in 1:length(itr_indx)) {
    spt <- out$crm[[itr_indx[j]]]$z.tld
    f0 <- out$crm[[itr_indx[j]]]$J.tld
    f0 <- f0 / sum(f0)
    # theta <- theta_sol(spt, f0, mu0, NULL)$theta
    # f0 <- exp(theta * spt) * f0 / sum(exp(theta * spt) * f0)
    theta <- 0                       # For null case and point masses or f0 / F0
    for (i in 1:length(yGrid)) {
      pro <- exp(theta * spt) * f0   # Which is same as f0 for the null case
      pk  <- pro / sum(pro)
      indx <- which(spt >= (yGrid[i] - c0) & spt <= (yGrid[i] + c0))
      indx2 <- which(spt <= yGrid[i] - c0)
      F0_est[i, j] <- sum(pk[indx2]) +
        sum(pk[indx] * (yGrid[i] - spt[indx] + c0)/ (2 * c0))
    }
  }
  return(F0_est)
}