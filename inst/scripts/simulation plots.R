require(KScorrect)
library(mistr)
library(tidyverse)
library(patchwork)
library(latex2exp)

setwd("/Users/entejar/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM")

out <- readRDS("0720_dpglm_simulation_all.rds")

# remove out[[]][[]] if NULL
null_idx <- numeric(0)
for(dat_idx in 1:100){
  for(setting in 1:9){
    if(is.null(out[[dat_idx]][[setting]])){
      null_idx <- c(null_idx, dat_idx)
    }
  }
}

out <- out[-unique(null_idx)]

mcmc_out <- lapply(out, function(sublist) {
  lapply(sublist, function(inner_sublist) inner_sublist$mcmc)
})

dat_out <- lapply(out, function(sublist) {
  lapply(sublist, function(inner_sublist) inner_sublist$data)
})

beta_samples <- lapply(mcmc_out, function(sublist) {
  lapply(sublist, function(inner_sublist) inner_sublist$beta)
})

theta_samples <- lapply(mcmc_out, function(sublist) {
  lapply(sublist, function(inner_sublist) inner_sublist$theta)
})

z_samples <- lapply(mcmc_out, function(sublist) {
  lapply(sublist, function(inner_sublist) inner_sublist$z)
})

crm_samples <- lapply(mcmc_out, function(sublist) {
  lapply(sublist, function(inner_sublist) inner_sublist$crm)
})

checkfn <- function(data_index, setting){
  y <- dat_out[[data_index]][[setting]][, 1]
  yest <- apply(z_samples[[data_index]][[setting]], 2, mean)
  mean((y - yest)^2)
}

check1 <- check2 <- check3 <- numeric(length(dat_out))
for(i in 1:length(dat_out)){
    check1[i] <- checkfn(i, 7)
    check2[i] <- checkfn(i, 8)
    check3[i] <- checkfn(i, 9)
}

data.frame(check1, check2, check3) %>% colMeans()  # choose n grid as (25, 100, 250)

setwd("/Users/entejar/Library/CloudStorage/Box-Box/DP-GLM/Data")
hustadTD <- readRDS("hustadTD.rds")

setwd("/Users/entejar/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM")

# Subset TD Data: Multi Word, Ordered y --------------------------------------------------
hustadTDMW <- filter(hustadTD %>% as.data.frame(),
                     intelligibility_type == 'multiword') %>% arrange(mean_intelligibility)

# KDE
dens  <- density(
  hustadTDMW$mean_intelligibility,
  kernel = 'gaussian',
  bw = 'nrd0',
  adjust = 1,
  n = 1000,
  cut = 3,
  ext = 4
)
dat   <- data.frame(y = dens$x, f0 = dens$y)
indx  <- dat$y >= 0 & dat$y <= 1

##===========================================================================##
#               Setting I (Null case): simulation truth                       #                                    #
##===========================================================================##

R <- 1
true_spt   <- (dat$y[dat$y >= 0 &
                       dat$y <= 1])[seq(1, length(dat$y[dat$y >= 0 &
                                                          dat$y <= 1]), R)]
true_f0    <- (dat$f0[dat$y >= 0 &
                        dat$y <= 1])[seq(1, length(dat$y[dat$y >= 0 &
                                                           dat$y <= 1]), R)]

true_f0 <- true_f0 / sum(true_f0)
mu0 <- sum(true_spt * true_f0) / sum(true_f0)

yGrid <- seq(0.05, 0.95, 0.05) #true_spt
c0 <- 0.025

theta <- 0
# true_pro <- exp(theta * true_spt) * true_f0
true_pk  <- true_f0 #true_pro / sum(true_pro)
f0_true <- numeric(length(yGrid))
for(i in 1:length(yGrid)){
  indx <- true_spt >= (yGrid[i] - c0) & true_spt <= (yGrid[i] + c0)
  f0_true[i] <- sum(true_pk[indx]) * (1 / (2*c0))
}

f0_true <- f0_true / sum(f0_true)
F0_true <- cumsum(f0_true)

itr_indx <- seq(500, 1000, 10)

f0_est_fn <- function(setting){
f0_est_list <- list()
for (dat_indx in 1:length(dat_out)) {
  f0_est <- matrix(NA, nrow = length(yGrid), ncol = length(itr_indx))
  for (j in 1:length(itr_indx)) {
    itr <- itr_indx[j]
    spt <- crm_samples[[dat_indx]][[setting]][[itr]]$z.tld
    f0 <- crm_samples[[dat_indx]][[setting]][[itr]]$J.tld
    tht <- theta_solver(spt, f0, mu0, NULL)$tht
    f0 <- exp(tht * spt) * f0 / sum(exp(tht * spt) * f0)
    theta <- 0
    for (i in 1:length(yGrid)) {
      pro <- exp(theta * spt) * f0
      pk  <- pro / sum(pro)
      indx <- spt >= (yGrid[i] - c0) & spt <= (yGrid[i] + c0)
      f0_est[i, j] <- sum(pk[indx]) * (1 / (2 * c0))
    }
  }
  f0_est_list[[dat_indx]] <- f0_est
}

# Mean est, 5% and 95% quantiles
fn1 <- function(x) {
  x <- x / sum(x)
  cumsum(x)
}

F0_est_list <- lapply(f0_est_list, function(x) apply(x, 2, fn1))

F0_est_mean <- lapply(F0_est_list, function(x) rowMeans(x))
F0_est_rmse <- lapply(F0_est_mean, function(x) abs(F0_true - x) / F0_true)
F0_est_rmse2 <- lapply(F0_est_mean, function(x) sqrt(mean((F0_true - x)^2))) %>% unlist()

F0_est_mean <- reduce(F0_est_mean, `+`) / length(F0_est_mean)


F0_est_q <- lapply(F0_est_list, function(x) apply(x, 1, quantile, prob = c(0.025, 0.975)))
F0_est_q <- reduce(F0_est_q, `+`) / length(F0_est_q)
F0_est_q5 <- F0_est_q[1, ]
F0_est_q95 <- F0_est_q[2, ]

fn2 <- function(x) {
  x <- x / sum(x)
  x
}

f0_est_list <- lapply(f0_est_list, function(x) apply(x, 2, fn2))
f0_est_mean <- lapply(f0_est_list, function(x) rowMeans(x))
f0_est_rmse <- lapply(f0_est_mean, function(x) abs(f0_true - x) / f0_true)
f0_est_rmse2 <- lapply(f0_est_mean, function(x) sqrt(mean((f0_true - x)^2))) %>% unlist()
f0_est_mean <- reduce(f0_est_mean, `+`) / length(f0_est_mean)

f0_est_q <- lapply(f0_est_list, function(x) apply(x, 1, quantile, prob = c(0.05, 0.95)))
f0_est_q <- reduce(f0_est_q, `+`) / length(f0_est_q)
f0_est_q5 <- f0_est_q[1, ]
f0_est_q95 <- f0_est_q[2, ]


# Create a data frame

f0_data_all <- data.frame(
  y = yGrid,
  true = f0_true,
  est = f0_est_mean,
  q5 = f0_est_q5,
  q95 = f0_est_q95
)

F0_data_all <- data.frame(
  y = yGrid,
  true = F0_true,
  est = F0_est_mean,
  q5 = F0_est_q5,
  q95 = F0_est_q95
)

return(list(f0 = f0_data_all,
            F0 = F0_data_all,
            mseF0 = F0_est_rmse,
            mseF0_2 = F0_est_rmse2,
            msef0 = f0_est_rmse,
            msef0_2 = f0_est_rmse2
            ))
}


beta_est <- lapply(beta_samples, function(sublist) {
  lapply(sublist, function(inner_sublist) apply(inner_sublist, 2, mean))
})

beta_df <- data.frame(beta = beta_est %>% unlist(),
                      j = rep(1:3, length(beta_est) * 9),
                      set = rep(rep(1:9, each = 3), length(beta_est)))

ggplot(beta_df) +
  geom_boxplot(aes(x = factor(j), y = beta, fill = factor(set)),
               position = position_dodge(width = 0.9)) +
  theme_bw() +
  labs(x = "j", y = TeX("Estimate of $\\beta_j$")) +
  scale_fill_manual(
    name = "Setting",
    labels = c(1:9),
    values = ggthemes::colorblind_pal()(9)
  ) +
  theme(legend.position = "right")

# Quantitative summary
beta_df %>% group_by(j, set) %>% summarise(
  mean = mean(beta),
) %>% spread(set, mean)

Y <- dat_out[[1]][[9]][, 1]
X <- dat_out[[1]][[9]][, -1]

gldrm(Y ~ X[, -1], link = 'logit')
gldrm(Y ~ X2 + X3, data = dat_out[[1]][[6]], link = 'logit')


# Create the plot
out1 <- f0_est_fn(1)
out2 <- f0_est_fn(2)
out3 <- f0_est_fn(3)


f0_df <- cbind(rbind(out1$f0, out2$f0, out3$f0),
                 sample_size = factor(rep(
                   c("n = 25", "n = 100", "n = 250"),
                   each = length(yGrid))))
f0_df$sample_size <- factor(f0_df$sample_size, levels = c("n = 25", "n = 100", "n = 250"))

pal <- ggthemes::colorblind_pal()(8)

p1 <- ggplot(f0_df) +
  geom_line(aes(y = true, x = y, color = "Truth")) +
  geom_line(aes(y = est, x = y, color = sample_size,
                group = sample_size))  +
  scale_color_manual(
    name = "",
    labels = c(
      "Truth" = "Truth",
      "n = 25" = "N = 25",
      "n = 100" = "N = 100",
      "n = 250" = "N = 250"
    ),
    values = c(
      "Truth" = pal[4],
      "n = 25" = pal[3],
      "n = 100" = pal[8],
      "n = 250" = pal[6]
    )
  ) +
  theme_bw() +
  labs(x = "y", y = "f_0") +
  theme(legend.position = "top")



F0_df <- cbind(rbind(out1$F0, out2$F0, out3$F0),
                 sample_size = factor(rep(
                   c("n = 25", "n = 100", "n = 250"),
                   each = length(yGrid))))
F0_df$sample_size <- factor(F0_df$sample_size, levels = c("n = 25", "n = 100", "n = 250"))

p2 <- ggplot(F0_df) +
  geom_line(aes(y = true, x = y, color = "Truth")) +
  geom_line(aes(y = est, x = y, color = sample_size, group = sample_size)) +
  geom_ribbon(
    aes(
      ymin = q5,
      ymax = q95,
      x = y,
      fill = sample_size,
      group = sample_size
    ),
    alpha = 0.2,
    show.legend = F
  ) +
  scale_color_manual(
    name = "",
    labels = c(
      "Truth" = "Truth",
      "n = 25" = "N = 25",
      "n = 100" = "N = 100",
      "n = 250" = "N = 250"
    ),
    values = c(
      "Truth" = pal[4],
      "n = 25" = pal[3],
      "n = 100" = pal[8],
      "n = 250" = pal[6]
    )
  ) +
  theme_bw() +
  labs(x = "y", y = "F_0") +
  theme(legend.position = "top") +
  guides(fill = guide_legend(title = "")) +
  scale_fill_manual(name = "",
                    values = c("n = 25" = pal[3],
                               "n = 100" = pal[8],
                               "n = 250" = pal[6]))

p1 + p2

# some quantitative summary
f0_df %>% group_by(sample_size) %>% summarise(
  mse = sqrt(mean((true - est)^2)),
  mse_q5 = sqrt(mean((true - q5)^2)),
  mse_q95 = sqrt(mean((true - q95)^2))
)

F0_df %>% group_by(sample_size) %>% summarise(
  mse = sqrt(mean((true - est)^2)),
  mse_q5 = sqrt(mean((true - q5)^2)),
  mse_q95 = sqrt(mean((true - q95)^2))
)

p2

plot_grid <- (1:length(yGrid))[seq(2, 18, 2)]
yGrid_plot <- yGrid[plot_grid]

df <- data.frame(mse = c(lapply(out1$msef0, function(x) x[plot_grid]) %>% unlist(),
                         lapply(out2$msef0, function(x) x[plot_grid]) %>% unlist(),
                         lapply(out3$msef0, function(x) x[plot_grid]) %>% unlist()),
                 y = rep(yGrid_plot, 3 * length(out1$msef0)),
                 group = factor(rep(c("n = 25", "n = 100", "n = 250"),
                                    each = length(yGrid_plot) * length(out1$msef0))))
df$group <- factor(df$group, levels = c("n = 25", "n = 100", "n = 250"))

ggplot() +
  geom_boxplot(data = df, aes(x = factor(round(y, 3)), y = mse, col = group),
               outliers = F) +
  theme_bw() +
  labs(x = "y", y = TeX("RMSE in estimating $f_0(y)$")) +
  theme(legend.position = "top") +
  scale_color_manual(
    name = "sample size (n)",
    labels = c("n = 25" = "25", "n = 100" = "100", "n = 250" = "250"),
    values = c("n = 25" = pal[3], "n = 100" = pal[4], "n = 250" = pal[8]))

df <- data.frame(mse = c(lapply(out1$mseF0, function(x) x[plot_grid]) %>% unlist(),
                         lapply(out2$mseF0, function(x) x[plot_grid]) %>% unlist(),
                         lapply(out3$mseF0, function(x) x[plot_grid]) %>% unlist()),
                 y = rep(yGrid_plot, 3 * length(out1$mseF0)),
                 group = factor(rep(c("n = 25", "n = 100", "n = 250"),
                                    each = length(yGrid_plot) * length(out1$mseF0))))
df$group <- factor(df$group, levels = c("n = 25", "n = 100", "n = 250"))

ggplot() +
  geom_boxplot(data = df, aes(x = factor(round(y, 3)), y = mse, col = group),
               outliers = F) +
  theme_bw() +
  labs(x = "y", y = TeX("RMSE in estimating $F_0(y)$")) +
  theme(legend.position = "top") +
  scale_color_manual(
    name = "sample size (n)",
    labels = c("n = 25" = "25", "n = 100" = "100", "n = 250" = "250"),
    values = c("n = 25" = pal[3], "n = 100" = pal[4], "n = 250" = pal[8]))


yGrid <- seq(0.05, 0.95, 0.05)
coverage_fn <- function(setting, yGrid){
  f0_est_list <- list()
  for (dat_indx in 1:length(dat_out)) {
    f0_est <- matrix(NA, nrow = length(yGrid), ncol = length(itr_indx))
    for (j in 1:length(itr_indx)) {
      itr <- itr_indx[j]
      spt <- crm_samples[[dat_indx]][[setting]][[itr]]$z.tld
      f0 <- crm_samples[[dat_indx]][[setting]][[itr]]$J.tld
      tht <- theta_solver(spt, f0, mu0, NULL)$tht
      f0 <- exp(tht * spt) * f0 / sum(exp(tht * spt) * f0)
      theta <- 0
      for (i in 1:length(yGrid)) {
        pro <- exp(theta * spt) * f0
        pk  <- pro / sum(pro)
        indx <- spt >= (yGrid[i] - c0) & spt <= (yGrid[i] + c0)
        f0_est[i, j] <- sum(pk[indx]) * (1 / (2 * c0))
      }
    }
    f0_est_list[[dat_indx]] <- f0_est
  }

  # Mean est, 5% and 95% quantiles
  fn1 <- function(x) {
    x <- x / sum(x)
    cumsum(x)
  }

  fn2 <- function(x) {
    x <- x / sum(x)
    x
  }

  f0_est_list <- lapply(f0_est_list, function(x) apply(x, 2, fn2))

  F0_est_list <- lapply(f0_est_list, function(x) apply(x, 2, fn1))

  F0_est_q <- lapply(F0_est_list, function(x) apply(x, 1, quantile, prob = c(0.025, 0.975)))
  coverageF0 <- lapply(F0_est_q, function(x) {
    (x[2, ] >= F0_true & x[1, ] <= F0_true)
  })

  f0_est_q <- lapply(f0_est_list, function(x) apply(x, 1, quantile, prob = c(0.025, 0.975)))
  coveragef0 <- lapply(f0_est_q, function(x) {
    (x[2, ] >= f0_true & x[1, ] <= f0_true)
  })


  return(list(coveragef0 = coveragef0,
              CIf0 = f0_est_q,
              coverageF0 = coverageF0,
              CIF0 = F0_est_q))
}

out1_ <- coverage_fn(1, yGrid)
out2_ <- coverage_fn(2, yGrid)
out3_ <- coverage_fn(3, yGrid)

df <- data.frame(coverage = c(reduce(out1_$coveragef0, `+`) / length(out1_$coveragef0),
                              reduce(out2_$coveragef0, `+`) / length(out2_$coveragef0),
                              reduce(out3_$coveragef0, `+`) / length(out3_$coveragef0)),
                 y = rep(yGrid, 3),
                 group = factor(rep(c("n = 25", "n = 100", "n = 250"), each = length(yGrid)))
)
df$group <- factor(df$group, levels = c("n = 25", "n = 100", "n = 250"))

df <- df %>% filter(y %in% c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

# instead of geom_bar use geom_point
ggplot() +
  geom_point(data = df, aes(x = factor(round(y, 3)), y = coverage, col = group),
             position = position_dodge(width = 0.1)) +
  theme_classic() +
  labs(x = "y", y = TeX("Coverage of 95% UI for $f_0(y)$")) +
  theme(legend.position = "top") +
  scale_color_manual(
    name = "sample size (n)",
    labels = c("n = 25" = "25", "n = 100" = "100", "n = 250" = "250"),
    values = c("n = 25" = pal[3], "n = 100" = pal[8], "n = 250" = pal[6])
  ) +
  geom_hline(yintercept = 0.95, linetype = "dashed", col = "black") +
  annotate(
    "text",
    x = 4,
    y = 0.952,
    label = "95% coverage line",
    col = "red"
  )



df2 <- data.frame(coverage = c(reduce(out1_$coverageF0, `+`) / length(out1_$coverageF0),
                               reduce(out2_$coverageF0, `+`) / length(out2_$coverageF0),
                               reduce(out3_$coverageF0, `+`) / length(out3_$coverageF0)),
                  y = rep(yGrid, 3),
                  group = factor(rep(c("n = 25", "n = 100", "n = 250"), each = length(yGrid)))
)

df2$group <- factor(df2$group, levels = c("n = 25", "n = 100", "n = 250"))
df2 <- df2 %>% filter(y %in% c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

ggplot() +
  geom_point(data = df2, aes(x = factor(round(y, 3)), y = coverage, col = group),
             position = position_dodge(width = 0.2)) +
  theme_classic() +
  labs(x = "y", y = TeX("Coverage of 95% UI for $F_0(y)$")) +
  theme(legend.position = "top") +
  scale_color_manual(
    name = "sample size (n)",
    labels = c("n = 25" = "25", "n = 100" = "100", "n = 250" = "250"),
    values = c("n = 25" = pal[3], "n = 100" = pal[8], "n = 250" = pal[6])
  ) +
  geom_hline(yintercept = 0.95, linetype = "dashed", col = "black") +
  annotate(
    "text",
    x = 4,
    y = 0.951,
    label = "95% coverage line",
    col = "red"
  )

# CI length comparison for f0; divide by the max CI length
df <- data.frame(CI_length = c(reduce(out1_$CIf0, `+`) / length(out1_$CIf0),
                               reduce(out2_$CIf0, `+`) / length(out2_$CIf0),
                               reduce(out3_$CIf0, `+`) / length(out3_$CIf0)),
                  y = rep(yGrid, 3),
                  group = factor(rep(c("n = 25", "n = 100", "n = 250"), each = length(yGrid)))
)

df <- df %>% group_by(y) %>% mutate(CI_length = CI_length / max(CI_length))


df$group <- factor(df$group, levels = c("n = 25", "n = 100", "n = 250"))
df <- df %>% filter(y %in% c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

ggplot() +
  geom_boxplot(data = df, aes(x = factor(round(y, 3)), y = CI_length, col = group),
               outliers = F) +
  theme_bw() +
  labs(x = "y", y = TeX("CI length for $f_0(y)$")) +
  theme(legend.position = "top") +
  scale_color_manual(
    name = "sample size (n)",
    labels = c("n = 25" = "25", "n = 100" = "100", "n = 250" = "250"),
    values = c("n = 25" = pal[3], "n = 100" = pal[4], "n = 250" = pal[8])
  )
