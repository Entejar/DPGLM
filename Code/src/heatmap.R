require(KScorrect)
library(mistr)
library(hrbrthemes)
library(viridis)

setwd("~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/src")
source("quantile_fns_mix_trunc_norm.R")

setwd("/Users/entejar/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM")

out <- readRDS("0627RDS3.rds")
beta_samples <- out$beta
G_samples <- out$G
theta_samples <- out$theta
age <- seq(min(unique(hustadTDMW$age_months)), max(unique(hustadTDMW$age_months)), 2)
prbs <- c(0.33, 0.67)
Xtest <- cbind(1, ns(age,
                     knots = quantile(age, probs = prbs),
                     intercept = FALSE))
colnames(Xtest) <- c("intercept", "ageKnot1", "ageKnot2", "ageKnot3")

itr_indx <- seq(1000, 3000, 20)
perScale <- 100
c0 <- 0.025

resultDFF <- list()
count <- 0
yGrid <- seq(min(unique(hustadTDMW$mean_intelligibility)),
             max(unique(hustadTDMW$mean_intelligibility)), 0.05)
for(itr in itr_indx){
  count <- count + 1
  spt <- G_samples[[itr]]$z.tld
  f0 <- G_samples[[itr]]$J.tld
  # temp <- data.frame(spt = spt, f0 = f0) %>% arrange(spt)
  # spt <- temp[, 1]
  # f0 <- temp[, 2]
  l <- length(spt)
  beta <- beta_samples[itr, ]
  mu <- link_fn(Xtest %*% beta, link)$mu
  out <- theta_solver(spt, f0, mu, NULL)
  theta <- out$tht
  p_mat <- matrix(NA, nrow = nrow(Xtest), ncol = length(yGrid))
  for(i in 1:nrow(Xtest)){
    for(j in 1:length(yGrid)){
    pro <- exp(theta[i] * spt) * f0
    pk  <- pro / sum(pro)
    indx <- spt >= (yGrid[j] - c0) & spt <= (yGrid[j] + c0)
    p_mat[i, j] <- sum(pk[indx]) * (1 / (2*c0))
    }
  }
  resultDFF[[count]] <- p_mat
  print(itr)
}

resultDF <- list()
count <- 0
for(i in 1:length(resultDFF)){
  if(sum(is.na(resultDFF[[i]])) == 0){
    count <- count + 1
    resultDF[[count]] <- resultDFF[[i]]
  }
}

element_wise_mean <- function(matrix_list) {
  array_3d <- simplify2array(matrix_list)
  result <- apply(array_3d, c(1, 2), mean, na.rm = TRUE)
  return(result)
}

XY <- expand.grid(age, 100 * yGrid)
XYZ <- data.frame(XY, Z = c(Reduce('+', resultDF[1:length(resultDF)]) / length(resultDF)))
#element_wise_mean(resultDF)))
colnames(XYZ) <- c("X", "Y", "Z")

pal <- ggthemes::colorblind_pal()(8)

#XYZ <- readRDS("~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/Plots/data/heatmapRDS3.rds")

p <- ggplot(XYZ, aes(X, Y, fill= Z)) +
  geom_tile(color = "grey") +
  scale_fill_distiller(palette = "Blues",
                       direction = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Age", y = "Intelligibility (in %)")

deep_color <- rgb(0, 0.5, 0)
p <- ggplot(XYZ, aes(X, Y, fill= Z)) +
  geom_tile(color = "grey") +
  scale_fill_gradient(low = "white", high = "darkblue") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Age (in months)", y = "Speech Intelligibility (in %)")

p

ggsave("~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/Plots/plots/heatmapRDSSingleWord.pdf", plot = p)

plot_data <- XYZ
saveRDS(plot_data, "~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/Plots/data/heatmapRDS3SingleWord.rds")



