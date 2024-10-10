require(KScorrect)
library(mistr)

setwd("~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/src")
source("quantile_fns_mix_trunc_norm.R")

setwd("/Users/entejar/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM")

out <- readRDS("0621_mcmc_samples_without_monotonicity0150pm.rds")
beta_samples <- out$beta
G_samples <- out$G
theta_samples <- dp_spglm$theta
age <- seq(min(unique(hustadTDMW$age_months)), max(unique(hustadTDMW$age_months)), 0.5)
prbs <- c(0.33, 0.67)
Xtest <- cbind(1, ns(age,
                     knots = quantile(age, probs = prbs),
                     intercept = FALSE))
colnames(Xtest) <- c("intercept", "ageKnot1", "ageKnot2", "ageKnot3")

itr_indx <- seq(2000, 3000, 10)
perScale <- 100
a <- 0.000005
resultDFF <- list()
for(itr in itr_indx){
  spt <- G_samples[[itr]]$z.tld
  f0 <- G_samples[[itr]]$J.tld
  f0 <- f0 / sum(f0)
  #theta0 <- theta_solver(spt, f0, mu0, NULL)$tht                  # identifiability constraint?
  #f0 <- f0 * exp(theta0 * spt) / sum(f0 * exp(theta0 * spt))
  #temp <- data.frame(spt = spt, f0 = f0) %>% arrange(spt)
  #spt <- temp[, 1]
  #f0 <- temp[, 2]
  l <- length(spt)
  beta <- beta_samples[itr, ]
  mu <- link_fn(Xtest %*% beta, link)$mu
  out <- theta_solver(spt, f0, mu, NULL)
  theta <- out$tht
  sigma0 <- SigmaY
  #sigma0 <- rep(max(SigmaY), length(sigma0))
  fTiltMat <- matrix(NA, nrow = nrow(Xtest), ncol = l)
  q5 <- q10 <- q25 <- q50 <- q75 <- q90 <- q95 <- numeric(nrow(Xtest))
  for(i in 1:nrow(Xtest)){
    pro <- (theta[i] * spt) +  log(f0)
    fTiltMat[i, ] <- exp(pro)
    # fTiltMat[i, ] <- fTiltMat[i, ] / (pnorm(b0 + 0.05, spt, rep(sigma0[i], length(spt))) -
    #                                    pnorm(a0 - 0.05, spt, rep(sigma0[i], length(spt))))
    q5[i] <-  KScorrect::qmixnorm(.05, mean = spt, sd = rep(sigma0[i], length(spt)), pro = fTiltMat[i, ], expand = a)
    q10[i] <- KScorrect::qmixnorm(.10, mean = spt, sd = rep(sigma0[i], length(spt)), pro = fTiltMat[i, ], expand = a)
    q25[i] <- KScorrect::qmixnorm(.25, mean = spt, sd = rep(sigma0[i], length(spt)), pro = fTiltMat[i, ], expand = a)
    q50[i] <- KScorrect::qmixnorm(0.5, mean = spt, sd = rep(sigma0[i], length(spt)), pro = fTiltMat[i, ], expand = a)
    q75[i] <- KScorrect::qmixnorm(0.75, mean = spt, sd = rep(sigma0[i], length(spt)), pro = fTiltMat[i, ], expand = a)
    q90[i] <- KScorrect::qmixnorm(0.90, mean = spt, sd = rep(sigma0[i], length(spt)), pro = fTiltMat[i, ], expand = a)
    q95[i] <- KScorrect::qmixnorm(0.95, mean = spt, sd = rep(sigma0[i], length(spt)), pro = fTiltMat[i, ], expand = a)
  }

  resultDFF[[itr]] <- data.frame(q5 = ifelse(perScale*q5 > 100, 100, perScale*q5),
                                q10 = ifelse(perScale*q10 > 100, 100, perScale*q10),
                                q25 = ifelse(perScale*q25 > 100, 100, perScale*q25),
                                q50 = ifelse(perScale*q50 > 100, 100, perScale*q50),
                                q75 = ifelse(perScale*q75 > 100, 100, perScale*q75),
                                q90 = ifelse(perScale*q90 > 100, 100, perScale*q90),
                                q95 = ifelse(perScale*q95 > 100, 100, perScale*q95))
  print(itr)

}

resultDF <- resultDFF[sapply(resultDFF, function(x) !is.null(x))]
itr_indx <- 1:length(resultDF)

hustadTD_temp <- filter(hustadTD %>% as.data.frame(),
                        intelligibility_type == "multiword")
ooo <- order(hustadTD_temp$age)
hustadTD_temp <- hustadTD_temp[ooo,]
meanAGE <- mean(hustadTD_temp$age_months)
sdAGE <- sd(hustadTD_temp$age_months)

hustadTD_dat <- data.frame(age_months = hustadTD_temp$age_months,
                           mean_intelligibility = hustadTD_temp$mean_intelligibility)


age_months <- age

K <- nrow(Xtest)
q5mat <- matrix(NA, nrow = K, ncol = length(itr_indx))
for(i in 1:K){
  for(j in 1:length(itr_indx)){
    q5mat[i, j] <- (resultDF[[itr_indx[j]]])[i, 1]
  }
}
#q5mat <- q5mat[, -unique(ceiling(which(q5mat == "NaN") / K))]
q5_mean <- rowMeans(q5mat)
q5_low <- apply(q5mat, 1, quantile, prob = 0.05)
q5_high <- apply(q5mat, 1, quantile, prob = 0.95)

q10mat <- matrix(NA, nrow = K, ncol = length(itr_indx))
for(i in 1:K){
  for(j in 1:length(itr_indx)){
    q10mat[i, j] <- (resultDF[[itr_indx[j]]])[i, 2]
  }
}
#q10mat <- q10mat[, -unique(ceiling(which(q10mat == "NaN") / K))]
q10_mean <- rowMeans(q10mat)
q10_low <- apply(q10mat, 1, quantile, prob = 0.05)
q10_high <- apply(q10mat, 1, quantile, prob = 0.95)

q25mat <- matrix(NA, nrow = K, ncol = length(itr_indx))
for(i in 1:K){
  for(j in 1:length(itr_indx)){
    q25mat[i, j] <- (resultDF[[itr_indx[j]]])[i, 3]
  }
}
#q25mat <- q25mat[, -unique(ceiling(which(q25mat == "NaN") / K))]
q25_mean <- rowMeans(q25mat)
q25_low <- apply(q25mat, 1, quantile, prob = 0.05)
q25_high <- apply(q25mat, 1, quantile, prob = 0.95)

q50mat <- matrix(NA, nrow = K, ncol = length(itr_indx))
for(i in 1:K){
  for(j in 1:length(itr_indx)){
    q50mat[i, j] <- (resultDF[[itr_indx[j]]])[i, 4]
  }
}
#q50mat <- q50mat[, -unique(ceiling(which(q50mat == "NaN") / K))]
q50_mean <- rowMeans(q50mat)
q50_low <- apply(q50mat, 1, quantile, prob = 0.05)
q50_high <- apply(q50mat, 1, quantile, prob = 0.95)

q75mat <- matrix(NA, nrow = K, ncol = length(itr_indx))
for(i in 1:K){
  for(j in 1:length(itr_indx)){
    q75mat[i, j] <- (resultDF[[itr_indx[j]]])[i, 5]
  }
}
#q50mat <- q50mat[, -unique(ceiling(which(q50mat == "NaN") / K))]
q75_mean <- rowMeans(q75mat)
q75_low <- apply(q75mat, 1, quantile, prob = 0.05)
q75_high <- apply(q75mat, 1, quantile, prob = 0.95)

q90mat <- matrix(NA, nrow = K, ncol = length(itr_indx))
for(i in 1:K){
  for(j in 1:length(itr_indx)){
    q90mat[i, j] <- (resultDF[[itr_indx[j]]])[i, 6]
  }
}
#q50mat <- q50mat[, -unique(ceiling(which(q50mat == "NaN") / K))]
q90_mean <- rowMeans(q90mat)
q90_low <- apply(q90mat, 1, quantile, prob = 0.05)
q90_high <- apply(q90mat, 1, quantile, prob = 0.95)

q95mat <- matrix(NA, nrow = K, ncol = length(itr_indx))
for(i in 1:K){
  for(j in 1:length(itr_indx)){
    q95mat[i, j] <- (resultDF[[itr_indx[j]]])[i, 7]
  }
}
#q50mat <- q50mat[, -unique(ceiling(which(q50mat == "NaN") / K))]
q95_mean <- rowMeans(q95mat)
q95_low <- apply(q95mat, 1, quantile, prob = 0.05)
q95_high <- apply(q95mat, 1, quantile, prob = 0.95)

low <- c(q5_low, q10_low, q25_low, q50_low, q75_low, q90_low, q95_low)
high <- c(q5_high, q10_high, q25_high, q50_high, q75_high, q90_high, q95_high)
qmean <- c(q5_mean, q10_mean, q25_mean, q50_mean, q75_mean, q90_mean, q95_mean)
data.frame(qmean, low, high)

resultDF3 <- data.frame(age = age_months, q = qmean,
                        low = low, high = high, q_type = rep(c("1", "2", "3", "4", "5", "6", "7"), each = nrow(Xtest)))

pal <- ggthemes::colorblind_pal()(8)



p <- ggplot() +
  geom_point(data = hustadTD_dat,
             aes(x = age_months, y = perScale*mean_intelligibility), col = pal[2], alpha = 0.7, size = 0.8) +
  geom_line(data = resultDF3, aes(x = age, y = q, col = q_type), lwd = 0.6) +
  geom_ribbon(data = resultDF3, aes(x = age, ymin = low, ymax = high,
                                    fill = q_type), alpha = 0.3) +
  scale_color_brewer(name = "Quantile",
                     labels = c("1" = "5%", "2" = "10%", "3" = "25%", "4" = "50%", "5" = "75%", "6" = "90%", "7" = "95%"),
                     type = "seq",
                     palette = 1,
                     direction = -1,
                     aesthetics = "colour") +
  scale_fill_brewer(name = "Quantile",
                     labels = c("1" = "5%", "2" = "10%", "3" = "25%", "4" = "50%", "5" = "75%", "6" = "90%", "7" = "95%"),
                     type = "seq",
                     palette = 1,
                     direction = -1, aesthetics = "fill")  +
  theme_bw() +
  labs(x = "age (months)", y = "intelligibility (in %)")

p

ggsave("~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/Plots/plots/quantilesRDS3.pdf", plot = p)

plot_data <- list(resultDF3 = resultDF3, hustadTD_dat = hustadTD_dat)
saveRDS(plot_data, "~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/Plots/data/quantilesRDS3.rds")


