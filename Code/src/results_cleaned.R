out_c  <- readRDS("~/Documents/GitHub/DPGLM/Code/cache/sim_1_2_Feb21_beta_1_0_minus_0.7_0.2.rds")
out_d  <- readRDS("~/Documents/GitHub/DPGLM/Code/cache/sim_1_2_Feb14_beta_1_0_minus_0.7_0.2.rds")
out_cd <- list()
for(i in 1:length(out_d)){
  out_cd[[36+i]] <- out_d[[i]]
}

out <- list()
for(i in 1:length(out_cd)){
  out[[i]] <- out_cd[[i]]$fit_
}

null_idx <- numeric(0)
for(i in 1:length(out)){
  if(is.null(out[[i]]$fit1$dpglm$beta)){
    null_idx <- c(null_idx, i)
  } else {
    acc_prob <- c(out[[i]]$fit1$dpglm$crm_acceptance, out[[i]]$fit2$dpglm$crm_acceptance, 
                  out[[i]]$fit3$dpglm$crm_acceptance, out[[i]]$fit1_max$dpglm$crm_acceptance,
                  out[[i]]$fit2_max$dpglm$crm_acceptance, out[[i]]$fit3_max$dpglm$crm_acceptance)
    if(min(acc_prob) < 0.05){
      null_idx <- c(null_idx, i)
    }
  }
}

if(length(null_idx) >= 1){
  out <- out[-unique(null_idx)]
}

out1 <- out2 <- out3 <- out1_max <- out2_max <- out3_max <- list()

for(i in 1:length(out)){
  out1[[i]] <- out[[i]]$fit1
  out2[[i]] <- out[[i]]$fit2
  out3[[i]] <- out[[i]]$fit3
  out1_max[[i]] <- out[[i]]$fit1_max
  out2_max[[i]] <- out[[i]]$fit2_max
  out3_max[[i]] <- out[[i]]$fit3_max
}

mu0 <- 0.7833362
beta1 <- beta2 <- beta3 <- matrix(0, nrow = length(out1), ncol = 8)
spt1 <- spt2 <- spt3 <- f01 <- f02 <- f03 <- numeric(0)
for(i in 1:length(out1)){
  beta1_mat <- beta2_mat <- beta3_mat <- matrix(0, nrow = 250, ncol = 2)
  for(l in 1:250){
    j <- l + 250 
    
    beta <- out1[[i]]$dpglm$beta[j,]
    spt <- c(out1[[i]]$dpglm$crm[[j]]$z.tld) 
    f0 <- c(out1[[i]]$dpglm$crm[[j]]$J.tld)
    theta0 <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu0, sampprobs = NULL, ySptIndex = NULL)$theta
    f0 <- exp(theta0 * spt) * f0 / sum(exp(theta0 * spt) * f0)
    theta <- out1[[i]]$dpglm$theta[j,] - theta0
    beta1_mat[l,] <- beta
    spt1 <- c(spt1, spt) 
    f01 <- c(f01, f0)
    
    beta <- out2[[i]]$dpglm$beta[j,]
    spt <- c(out2[[i]]$dpglm$crm[[j]]$z.tld) 
    f0 <- c(out2[[i]]$dpglm$crm[[j]]$J.tld)
    theta0 <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu0, sampprobs = NULL, ySptIndex = NULL)$theta
    f0 <- exp(theta0 * spt) * f0 / sum(exp(theta0 * spt) * f0)
    theta <- out2[[i]]$dpglm$theta[j,] - theta0
    beta2_mat[l,] <- beta
    spt2 <- c(spt2, spt) 
    f02 <- c(f02, f0)
    
    beta <- out3[[i]]$dpglm$beta[j,]
    spt <- c(out3[[i]]$dpglm$crm[[j]]$z.tld) 
    f0 <- c(out3[[i]]$dpglm$crm[[j]]$J.tld)
    theta0 <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu0, sampprobs = NULL, ySptIndex = NULL)$theta
    f0 <- exp(theta0 * spt) * f0 / sum(exp(theta0 * spt) * f0)
    theta <- out3[[i]]$dpglm$theta[j,] - theta0
    beta3_mat[l,] <- beta
    spt3 <- c(spt3, spt) 
    f03 <- c(f03, f0)
    
  }
  
  X = out1[[i]]$data[, -c(1,2)]
  y = out1[[i]]$data[, 1]
  fit <- gldrm(y ~ X[, -1], link = "logit")
  beta1[i,] <- c(colMeans(beta1_mat), apply(beta1_mat, 2, sd),
                 fit$beta, fit$seBeta)
  
  X = out2[[i]]$data[, -c(1,2)]
  y = out2[[i]]$data[, 1]
  fit <- gldrm(y ~ X[, -1], link = "logit")
  beta2[i,] <- c(colMeans(beta2_mat), apply(beta2_mat, 2, sd),
                 fit$beta, fit$seBeta)
  
  X = out3[[i]]$data[, -c(1,2)]
  y = out3[[i]]$data[, 1]
  fit <- gldrm(y ~ X[, -1], link = "logit")
  beta3[i,] <- c(colMeans(beta3_mat), apply(beta3_mat, 2, sd),
                 fit$beta, fit$seBeta)
  
  
}

beta1[, ] %>% colMeans()
beta2[, ] %>% colMeans()
beta3[, ] %>% colMeans()

f0_true_kde <- readRDS("sim_truth/f0_true_kde.rds")
true_spt <- f0_true_kde$spt
f0_kde <- f0_true_kde$f0 / sum(f0_true_kde$f0)
y_grid <- seq(0, 1, 0.001)

get_f0 <- function(out, i){
  F0_mat <- f0_mat <- matrix(NA, nrow = 250, ncol = length(y_grid))
  c0 <- 0.025
  for(k in 1:250){
    j <- 250 + k
    spt <- c(out[[i]]$dpglm$crm[[j]]$z.tld) 
    f0 <- c(out[[i]]$dpglm$crm[[j]]$J.tld)
    theta0 <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu0, sampprobs = NULL, ySptIndex = NULL)$theta
    f0 <- exp(theta0 * spt) * f0 / sum(exp(theta0 * spt) * f0)
    theta <- 0
    lower <- spt - c0
    upper <- spt + c0
    lower[lower < 0] <- 0
    upper[upper > 1] <- 1
    width <- upper - lower
    for(l in 1:length(y_grid)){
      pro <- exp(theta * spt) * f0
      pk  <- pro / sum(pro)
      indx <- y_grid[l] >= lower & y_grid[l] <= upper
      indx2 <- y_grid[l] > upper
      f0_mat[k, l] <- sum(pk[indx] / width[indx])
      F0_mat[k, l] <- sum(pk[indx] * (y_grid[l] - lower[indx])/ width[indx]) + sum(pk[indx2])
    }
  }
  return(list(f0 = f0_mat, F0 = F0_mat))
}

# Initialize matrices for CDF estimates and confidence intervals
F01 <- F02 <- F03 <- F01_max <- F02_max <- F03_max <- matrix(0, nrow = length(out), ncol = length(y_grid))
F01L <- F02L <- F03L <- F01_maxL <- F02_maxL <- F03_maxL <- matrix(0, nrow = length(out), ncol = length(y_grid))
F01U <- F02U <- F03U <- F01_maxU <- F02_maxU <- F03_maxU <- matrix(0, nrow = length(out), ncol = length(y_grid))

# Loop through each sample
for(i in 1:length(out)){
  # Compute the function values for each dataset
  temp <- get_f0(out1, i)
  F01[i,] <- colMeans(temp$F0)
  F01L[i,] <- apply(temp$F0, 2, quantile, probs = 0.025)  # Lower bound
  F01U[i,] <- apply(temp$F0, 2, quantile, probs = 0.975)  # Upper bound
  
  temp <- get_f0(out2, i)
  F02[i,] <- colMeans(temp$F0)
  F02L[i,] <- apply(temp$F0, 2, quantile, probs = 0.025)
  F02U[i,] <- apply(temp$F0, 2, quantile, probs = 0.975)
  
  temp <- get_f0(out3, i)
  F03[i,] <- colMeans(temp$F0)
  F03L[i,] <- apply(temp$F0, 2, quantile, probs = 0.025)
  F03U[i,] <- apply(temp$F0, 2, quantile, probs = 0.975)
  
  temp <- get_f0(out1_max, i)
  F01_max[i,] <- colMeans(temp$F0)
  F01_maxL[i,] <- apply(temp$F0, 2, quantile, probs = 0.025)
  F01_maxU[i,] <- apply(temp$F0, 2, quantile, probs = 0.975)
  
  temp <- get_f0(out2_max, i)
  F02_max[i,] <- colMeans(temp$F0)
  F02_maxL[i,] <- apply(temp$F0, 2, quantile, probs = 0.025)
  F02_maxU[i,] <- apply(temp$F0, 2, quantile, probs = 0.975)
  
  temp <- get_f0(out3_max, i)
  F03_max[i,] <- colMeans(temp$F0)
  F03_maxL[i,] <- apply(temp$F0, 2, quantile, probs = 0.025)
  F03_maxU[i,] <- apply(temp$F0, 2, quantile, probs = 0.975)
}

# Plot the first set of CDF estimates with confidence intervals
plot(y_grid, colMeans(F01), type = "l", col = "cyan", ylim = c(0, 1), ylab = "CDF", xlab = "y", main = "CDF Estimates")
lines(true_spt, cumsum(f0_kde), type = 'l', col = "black")
lines(y_grid, colMeans(F02), type = "l", col = "green")
lines(y_grid, colMeans(F03), type = "l", col = "blue")

# Add confidence intervals as shaded regions
polygon(c(y_grid, rev(y_grid)), c(colMeans(F01L), rev(colMeans(F01U))), col = rgb(0,1,1,0.2), border = NA)
polygon(c(y_grid, rev(y_grid)), c(colMeans(F02L), rev(colMeans(F02U))), col = rgb(0,1,0,0.2), border = NA)
polygon(c(y_grid, rev(y_grid)), c(colMeans(F03L), rev(colMeans(F03U))), col = rgb(0,0,1,0.2), border = NA)

# Second plot: Max CDF estimates
plot(y_grid, colMeans(F01_max), type = "l", col = "cyan", ylim = c(0, 1), ylab = "CDF", xlab = "y", main = "Max CDF Estimates")
lines(true_spt, cumsum(f0_kde), type = 'l', col = "black")
lines(y_grid, colMeans(F02_max), type = "l", col = "green")
lines(y_grid, colMeans(F03_max), type = "l", col = "blue")

# Add confidence intervals as shaded regions
polygon(c(y_grid, rev(y_grid)), c(colMeans(F01_maxL), rev(colMeans(F01_maxU))), col = rgb(0,1,1,0.2), border = NA)
polygon(c(y_grid, rev(y_grid)), c(colMeans(F02_maxL), rev(colMeans(F02_maxU))), col = rgb(0,1,0,0.2), border = NA)
polygon(c(y_grid, rev(y_grid)), c(colMeans(F03_maxL), rev(colMeans(F03_maxU))), col = rgb(0,0,1,0.2), border = NA)

library(gldrm)

# Load true distribution
f0_true_kde <- readRDS("sim_truth/f0_true_kde.rds")
true_spt <- f0_true_kde$spt
f0_kde <- f0_true_kde$f0 / sum(f0_true_kde$f0)

# Function to compute f0 and F0
y_grid <- seq(0, 1, 0.001)
c0 <- 0.025

generate_f0 <- function(out, i) {
  F0_mat <- f0_mat <- matrix(NA, nrow = 250, ncol = length(y_grid))
  for (k in 1:250) {
    j <- 250 + k
    spt <- out[[i]]$dpglm$crm[[j]]$z.tld
    f0 <- out[[i]]$dpglm$crm[[j]]$J.tld
    theta0 <- gldrm:::getTheta(spt = spt, f0 = f0, mu = mu0, sampprobs = NULL, ySptIndex = NULL)$theta
    f0 <- exp(theta0 * spt) * f0 / sum(exp(theta0 * spt) * f0)
    lower <- pmax(spt - c0, 0)
    upper <- pmin(spt + c0, 1)
    width <- upper - lower
    for (l in seq_along(y_grid)) {
      pro <- exp(theta0 * spt) * f0
      pk  <- pro / sum(pro)
      indx <- y_grid[l] >= lower & y_grid[l] <= upper
      indx2 <- y_grid[l] > upper
      f0_mat[k, l] <- sum(pk[indx] / width[indx])
      F0_mat[k, l] <- sum(pk[indx] * (y_grid[l] - lower[indx]) / width[indx]) + sum(pk[indx2])
    }
  }
  list(f0 = f0_mat, F0 = F0_mat)
}

library(ggplot2)

# Compute empirical CDF and quantiles
ecdf_vals <- cumsum(f0_kde) / sum(f0_kde)
y_percentiles <- approx(ecdf_vals, true_spt, xout = c(0.1, 0.25, 0.5, 0.75, 0.9))$y

# Function to compute exceedance probabilities
compute_exceedance_prob <- function(F0_matrix) {
  apply(F0_matrix, 1, function(row) 1 - approx(y_grid, row, y_percentiles, rule = 2)$y)
}

# Function to compute coverage probability
compute_coverage <- function(estimates, true_values) {
  apply(estimates, 1, function(est) mean(true_values >= est[2] & true_values <= est[3]))
}

# Function to extract F0 estimates and UQ for a given dataset
generate_f0 <- function(out, i) {
  temp <- get_f0(out, i)
  list(
    meanF0 = colMeans(temp$F0),
    lowerF0 = apply(temp$F0, 2, quantile, probs = 0.025),
    upperF0 = apply(temp$F0, 2, quantile, probs = 0.975),
    exceed_prob = compute_exceedance_prob(temp$F0)
  )
}

# List of all situations
situations <- list("F01" = out1, "F02" = out2, "F03" = out3,
                   "F01_max" = out1_max, "F02_max" = out2_max, "F03_max" = out3_max)

# Compute exceedance probabilities and UQ
results <- list()
for (name in names(situations)) {
  results[[name]] <- list()
  for (i in seq_along(situations[[name]])) {
    results[[name]][[i]] <- generate_f0(situations[[name]], i)
  }
}

# Prepare data for plotting
plot_data <- data.frame()
for (name in names(results)) {
  all_iterations <- do.call(rbind, lapply(results[[name]], function(res) {
    data.frame(
      y_grid = y_grid,
      meanF0 = res$meanF0,
      lowerF0 = res$lowerF0,
      upperF0 = res$upperF0
    )
  }))
  
  avg_results <- all_iterations %>%
    group_by(y_grid) %>%
    summarise(
      meanF0 = mean(meanF0),
      lowerF0 = mean(lowerF0),
      upperF0 = mean(upperF0)
    )
  
  avg_results$method <- name
  plot_data <- rbind(plot_data, avg_results)
}

plot_data$method <- factor(plot_data$method, levels = names(results))

# Plot estimated CDFs with confidence intervals
ggplot() +
  geom_line(data = plot_data, mapping = aes(x = y_grid, y = meanF0, 
  ), color = "blue") +
  geom_ribbon(data = plot_data, mapping = aes(x = y_grid, ymin = lowerF0, ymax = upperF0), 
              alpha = 0.2) +
  geom_line(data = data.frame(y_grid = true_spt, meanF0 = cumsum(f0_kde)),
            aes(x = y_grid, y = meanF0), color = "black", linetype = "dashed") +
  facet_wrap(. ~ method) +
  labs(x = "y", y = "F0(y)") +
  theme_minimal() +
  theme(legend.position = "right")

exceed_data <- data.frame()
for (name in names(results)) {
  all_iterations <- do.call(rbind, lapply(results[[name]], function(res) {
    t(res$exceed_prob)
  }))
  
  avg_results <- c(colMeans(all_iterations), name)
  
  exceed_data <- rbind(exceed_data, avg_results)
}

true_exceed_prob <- 1 - approx(true_spt, cumsum(f0_kde), y_percentiles, rule = 2)$y
exceed_data <- rbind(exceed_data, c(true_exceed_prob, "truth"))
colnames(exceed_data) <- c("10%", "25%", "50%", "75%", "90%", "Situation")
exceed_data

# Calculate coverage for each case
coverage <- list()
for (name in names(results)) {
  exceed_matrix <- do.call(rbind, lapply(results[[name]], function(res) {
    temp <- apply(t(res$exceed_prob), 2, quantile, probs = c(0.025, 0.975))
    c(true_exceed_prob >= temp[1,] &  true_exceed_prob <= temp[2,], temp[2,] - temp[1,])
  }))
  coverage[[name]] <- colMeans(exceed_matrix)
}

coverage
