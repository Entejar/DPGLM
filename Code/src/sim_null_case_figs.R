# Load packages and helper functions
source("load.R")

# Truths: f0, F0
true <- readRDS("data/f0_true_kde.rds")
true_spt <- true$spt
true_f0 <- true$f0 / sum(true$f0)
mu0 <- sum(true_spt * true_f0)
true_F0 <- get_cdf(true_spt, true_f0)

# Load MCMC samples
out <- readRDS("cache/oct_24_null_case_sim_20_replicates.rds")
# remove NULL if any
out <- out[!sapply(out, is.null)]

# Plots!!
c0 <- 0.025
save <- seq(300, 500, 10)
plot_indx <- seq(1, length(true_spt), 10)
res <- plot_fn(out, true_F0[plot_indx], save, true_spt[plot_indx], mu0, c0)
res$p_F
tikzprint(res$p_F, "01_sim_null_cases_baseline_cdf", width = 10, height = 5)
res$ks_summary 

# Rough work and debugging
plot(true_spt, true_F0, type = "l", col = "red", lwd = 2, xlab = "spt", ylab = "f0", main = "True f0 vs. Estimated f0")

spt2 <- spt1 <- numeric(0)
for(i in 1:length(out)){
  for(j in save){
    spt1 <- c(spt1, out[[i]]$mcmc_samples$MS1$crm[[j]]$z.tld)
    spt2 <- c(spt2, out[[i]]$mcmc_samples$MS2$crm[[j]]$z.tld)
  }
}

hist(spt1, breaks = 100)
hist(spt2, breaks = 100)

summary(spt1)
mean(spt1 > 0.85)
mean(spt2 > 0.85)
summary(spt2)

# Notes: Yes, spt1 has very few obs > 0.85 (prob = 0.03), but spt2 has a lot of obs > 0.85 (prob = 0.38). Why this diff?
