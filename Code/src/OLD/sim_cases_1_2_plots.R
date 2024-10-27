library(KScorrect)
library(mistr)
library(tidyverse)
library(patchwork)
library(latex2exp)
source("load.R")

out1 <- readRDS("cache/oct_24_null_case_sim_20_replicates.rds")
out2 <- readRDS("cache/oct_24_point_masses_sim_20_replicates.rds")

# remove NULL if any

out1 <- out1[!sapply(out1, is.null)]
out2 <- out2[!sapply(out2, is.null)]

cdf <- calculate_cdf(true_spt, true_f0)

f0_true_kde <- readRDS("data/f0_true_kde.rds")
true_spt <- f0_true_kde$spt
true_f0 <- f0_true_kde$f0 / sum(f0_true_kde$f0)
mu0_1 <- sum(true_spt * true_f0) / sum(true_f0)
yGrid <- true_spt[seq(1, length(true_spt), 1)]
F0_true_1 <- calculate_cdf(true_spt, true_f0)

f0_true_beta <- readRDS("data/f0_true_beta.rds")
true_f0 <- f0_true_beta$f0 / sum(f0_true_beta$f0)
true_spt <- f0_true_beta$spt
#eps <- 0.00001
#yGrid2 <- c(seq(eps, min(yGrid), eps), yGrid, seq(max(yGrid), 1 - eps, eps)) 
F0_true_2 <- true_F0_fn(true_f0, true_spt, yGrid2, theta, c0)
p0 <- 0.1
p1 <- 0.4
F0_true_2 <- p0 + (F0_true_2 * (1 - p0 - p1))
F0_true_2 <- c(0.1, F0_true_2, 1)
true_f0 <- c(0.1, true_f0 * 0.5, 0.4)
true_spt <- c(0, true_spt, 1)
yGrid2 <- c(0, yGrid2, 1)
mu0_2 <- sum(true_spt * true_f0) / sum(true_f0)
F0_true_3 <- get_cdf(true_spt, true_f0)

plot(yGrid2, F0_true_2, type = "l", col = "blue", lwd = 2, xlab = "y", ylab = "f0")
lines(true_spt, F0_true_3, col = "red", lwd = 2)

F0_est_fn <- function(out, save, yGrid, mu0) {
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

F0_est_fn <- function(out, save, yGrid, mu0) {
  itr_indx <- save
  F0_est <- matrix(0, nrow = length(yGrid), ncol = length(itr_indx))
  for (j in 1:length(itr_indx)) {
    spt <- out$crm[[itr_indx[j]]]$z.tld
    f0 <- out$crm[[itr_indx[j]]]$J.tld
    f0 <- f0 / sum(f0)
    #theta <- theta_sol(spt, f0, mu0, NULL)$theta
    #f0 <- exp(theta * spt) * f0 / sum(exp(theta * spt) * f0)
    theta <- 0                       # For null case and point masses or f0 / F0
    for (i in 2:(length(yGrid)-1)) {
      pro <- exp(theta * spt) * f0   # Which is same as f0 for the null case
      pk  <- pro / sum(pro)
      indx <- which(spt >= (yGrid[i] - c0) & spt <= (yGrid[i] + c0))
      indx2 <- which(spt <= yGrid[i] - c0)
      F0_est[i, j] <- sum(pk[indx2]) +
        sum(pk[indx] * (yGrid[i] - spt[indx] + c0)/ (2 * c0))
    }
    indx <- which(spt >= (yGrid[1] - c0) & spt <= (yGrid[1] + c0))
    p0_hat <- sum(pk[indx] * (yGrid[1] - spt[indx] + c0)/ (2 * c0))
    F0_est[1, j] <- p0_hat
    indx <- which(spt >= (yGrid[length(yGrid)] - c0) & spt <= (yGrid[length(yGrid)] + c0))
    p1_hat <- sum(pk[indx] * (yGrid[length(yGrid)] - spt[indx] + c0)/ (2 * c0))
    F0_est[2:(length(yGrid)-1), j] <- p0_hat + (F0_est[2:(length(yGrid)-1), j] * (1 - p0_hat - p1_hat))
    F0_est[length(yGrid), j] <- 1
  }
  return(F0_est)
}

spt2 <- spt1 <- numeric(0)
for(i in 1:20){
  for(j in save){
  spt1 <- c(spt1, out1[[i]]$mcmc_samples$MS1$crm[[j]]$z.tld)
  spt2 <- c(spt2, out2[[i]]$mcmc_samples$MS1$crm[[j]]$z.tld)
  }
}

hist(spt1, breaks = 100)
hist(spt2, breaks = 100)

big_fn <- function(out, F0_true, save, yGrid, mu0){
err <- numeric(0)
basedensity_200 <- array(NA, dim = c(length(out), length(yGrid), length(save)))
for (i in 1:length(out)) {
  tryCatch({
    basedensity_200[i,,] <- F0_est_fn(out[[i]]$mcmc_samples$MS1, save, yGrid, mu0)
  }, error = function(e) {
    err <<- c(err, i)
    basedensity_200[i,,] <<- NA           # Fill NA if mu0 outside support of f_0: can happen for posterior samples
    print(paste("Error in iteration", i))
  })
}

# Function to process each non-NA slice
process_slice <- function(slice) {
  if(all(is.na(slice))) {
    return(NULL)  # Return NULL for completely NA slices
  } else {
    return(rowMeans(slice, na.rm = TRUE))
  }
}


df_list <- list()

for(i in 1:dim(basedensity_200)[1]) {
  slice_mean <- process_slice(basedensity_200[i,,])
  if(!is.null(slice_mean)) {
    df_list[[length(df_list) + 1]] <- data.frame(yGrid = yGrid, mean = slice_mean)
  }
}

# Combine all data frames
final_df <- do.call(rbind, df_list)

final_df <- cbind(final_df, indx = rep(1:length(df_list), each = length(yGrid)), 
                  F0_true = rep(F0_true, length(df_list)))

ks1 <- numeric(0)
t <- nrow(final_df) / length(yGrid)

for(i in 1:t){
  ks1 <- c(ks1, ks.test(final_df$mean[(i-1)*length(yGrid) + 1:i*length(yGrid)], final_df$F0_true[(i-1)*length(yGrid) + 1:i*length(yGrid)], alternative = 'two.sided', simulate.p.value = TRUE)$statistic)
}

# Plotting

p1 <- ggplot(final_df, aes(x = yGrid, y = mean, group = indx)) +
  geom_line(color = 'blue') +
  geom_line(aes(y = F0_true), color = 'red') +
  labs(x = "y", y = TeX("$F_0$")) +
  theme_bw() +
  annotate("text", x = 0.15, y = 0.8, label = "n = 200", color = "black") +
  annotate("rect", xmin = 0.05, xmax = .25, ymin = 0.75, ymax = 0.85, fill = "cyan", alpha = 0.2)



err <- numeric(0)
basedensity_400 <- array(NA, dim = c(length(out), length(yGrid), length(save)))
for (i in 1:length(out)) {
  tryCatch({
    basedensity_400[i,,] <- F0_est_fn(out[[i]]$mcmc_samples$MS2, save, yGrid, mu0)
  }, error = function(e) {
    err <<- c(err, i)
    basedensity_400[i,,] <<- NA           # Fill NA if mu0 outside support of f_0: can happen for posterior samples
    print(paste("Error in iteration", i))
  })
}

# Function to process each non-NA slice
process_slice <- function(slice) {
  if(all(is.na(slice))) {
    return(NULL)  # Return NULL for completely NA slices
  } else {
    return(rowMeans(slice, na.rm = TRUE))
  }
}


df_list <- list()
for(i in 1:dim(basedensity_400)[1]) {
  slice_mean <- process_slice(basedensity_400[i,,])
  if(!is.null(slice_mean)) {
    df_list[[length(df_list) + 1]] <- data.frame(yGrid = yGrid, mean = slice_mean)
  }
}

# Combine all data frames
final_df <- do.call(rbind, df_list)

final_df <- cbind(final_df, indx = rep(1:length(df_list), each = length(yGrid)), 
                  F0_true = rep(F0_true, length(df_list)))

ks2 <- numeric(0)
t <- nrow(final_df) / length(yGrid)

for(i in 1:t){
  ks2 <- c(ks2, ks.test(final_df$mean[(i-1)*length(yGrid) + 1:i*length(yGrid)], final_df$F0_true[(i-1)*length(yGrid) + 1:i*length(yGrid)], alternative = 'two.sided', simulate.p.value = TRUE)$statistic)
}

# Plotting

p2 <- ggplot(final_df, aes(x = yGrid, y = mean, group = indx)) +
  geom_line(color = 'blue') +
  geom_line(aes(y = F0_true), color = 'red') +
  labs(x = "y", y = TeX("$F_0$")) +
  theme_bw() +
  annotate("text", x = 0.15, y = 0.8, label = "n = 400", color = "black") +
  annotate("rect", xmin = 0.05, xmax = .25, ymin = 0.75, ymax = 0.85, fill = "cyan", alpha = 0.2)

p_F <- p1 + p2 

# Calculate KS statistics

ks <- data.frame(ks = c(ks1, ks2), group = rep(c(200, 400), c(length(ks1), length(ks2))))

ks_summary <- ks %>% group_by(group) %>% summarise(mean = mean(ks), sd = sd(ks))

p_ks <- ggplot(ks, aes(x = factor(group), y = ks)) +
  geom_boxplot() +
  theme_bw()

return(list(p_F = p_F, p_ks = p_ks, ks_summary = ks_summary))
}

save <- seq(300, 500, 10)
res1 <- big_fn(out1, F0_true_1, save, yGrid, mu0_1)
res2 <- big_fn(out2, F0_true_2, save, yGrid2, mu0_2)
res1$p_F
res2$p_F
re3 <- big_fn(out2, F0_true_2, save, yGrid2, mu0_2)
res2$p_F
re3$p_F
tikzprint(res1$p_F, "01_sim_null_cases_baseline_cdf", width = 10, height = 5)
tikzprint(res2$p_F, "02_sim_point_masses_baseline_cdf", width = 10, height = 5)

res1$ks_summary 
res2$ks_summary
res1$p_F
