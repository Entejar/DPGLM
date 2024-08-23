# Fig: survival fn + CI
require(KScorrect)
library(mistr)
library(hrbrthemes)
library(viridis)
library(latex2exp)

setwd("~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/src")
source("quantile_fns_mix_trunc_norm.R")

setwd("/Users/entejar/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM")

out <- readRDS("0630RDS3SingleWord.rds")
beta_samples <- out$beta
G_samples <- out$G
theta_samples <- out$theta
age <- seq(min(unique(hustadTDMW$age_months)), max(unique(hustadTDMW$age_months)), 10)
prbs <- c(0.33, 0.67)
Xtest <- cbind(1, ns(age,
                     knots = quantile(age, probs = prbs),
                     intercept = FALSE))
colnames(Xtest) <- c("intercept", "ageKnot1", "ageKnot2", "ageKnot3")

itr_indx <- seq(1000, 3000, 1)
perScale <- 100
c0 <- 0.025

resultDFF <- list()

count <- 0
yGrid <- seq(min(unique(hustadTDMW$mean_intelligibility)),
             max(unique(hustadTDMW$mean_intelligibility)), 0.005)

tailsum <- function(x) {
  rev(cumsum(rev(x)))
}

for(itr in itr_indx){
  count <- count + 1
  spt <- G_samples[[itr]]$z.tld
  f0 <- G_samples[[itr]]$J.tld
  l <- length(spt)
  beta <- beta_samples[itr, ]
  mu <- link_fn(Xtest %*% beta, link)$mu
  out <- theta_solver(spt, f0, mu, NULL)
  theta <- out$tht
  p_mat <- matrix(NA, nrow = nrow(Xtest), ncol = length(yGrid))
  survivalDF <- matrix(0, nrow = nrow(Xtest), ncol = length(yGrid))
  for(i in 1:nrow(Xtest)){
    for(j in 1:length(yGrid)){
      pro <- exp(theta[i] * spt) * f0
      pk  <- pro / sum(pro)
      indx <- spt >= (yGrid[j] - c0) & spt <= (yGrid[j] + c0)
      p_mat[i, j] <- sum(pk[indx]) * (1 / (2*c0))
    }
    p_mat[i, ] <- p_mat[i, ] / sum(p_mat[i, ])
    survivalDF[i, ] <- sum(p_mat[i, ]) - c(0, cumsum(p_mat[i, 1:(length(yGrid)-1)]))
  }
  resultDFF[[count]] <- survivalDF
  print(itr)
}

get_data_for_row <- function(row_index) {
  data <- sapply(resultDFF, function(iter) iter[row_index, ])
  colnames(data) <- names(resultDFF)
  return(data)
}

# Get the number of rows in the matrices
num_rows <- nrow(resultDFF[[1]])

# Extract data for each row
row_data <- lapply(1:num_rows, get_data_for_row)

low <- high <- est <- numeric(0)
for(i in 1:length(age)){
  low <- c(low, apply(row_data[[i]], 1, quantile, prob = 0.05))
  high <- c(high, apply(row_data[[i]], 1, quantile, prob = 0.95))
  est <- c(est, apply(row_data[[i]], 1, mean))
}

plotDF <- data.frame(y = rep((yGrid * 100), length(age)),
                      low = low, high = high, est = est, Age = factor(rep(age, each = length(yGrid))))


#plotDF <- readRDS("~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/Plots/data/exceedanceRDS3.rds")

p <- ggplot(data = plotDF) +
  geom_ribbon(aes(x = y, ymin = low, ymax = high, fill = Age), alpha = 0.2) +
  geom_line(aes(x = y, y = est, col = Age), lwd = 1) +
  theme_minimal() +
  labs(
    x = "Speech Intelligibility (in %)",
    y = TeX("$P(Y \\geq Intelligibility | X = Age)$"),
    color = "Age (in months)",
    fill = "Age (in months)"
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_viridis_d(option = "plasma", end = 0.9) +
  scale_fill_viridis_d(option = "plasma", end = 0.9) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    #plot.title = element_text(size = 16, face = "bold"),
    #plot.subtitle = element_text(size = 12, color = "gray50"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 5, 0),
    panel.grid.minor = element_blank()
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 2), nrow = 1, title.position = "left"),
    fill = guide_legend(nrow = 1, title.position = "left")
  )

p

ggsave("~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/Plots/plots/exceedanceRDS3SingleWord.pdf", plot = p)

plot_data <- plotDF
saveRDS(plot_data, "~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/Plots/data/exceedanceRDS3SingleWord.rds")
