# Source Files  -------------------------------------------------
setwd("~/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/src")
source("functions.R")

setwd("/Users/entejar/Documents/PhD Thesis Works/Paul and Peter/Project 2/DPGLM/inst/scripts")
source("dpglm.R")

# Load Data --------------------------------------------------------------------
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
true_spt <- dat$y[indx]
true_f0    <- dat$f0[indx]

# theta = 0 and hence, btheta = 0
# f(Y | X) = f(Y) = f0
sim_I <- function(n, p) {
  X <- matrix(runif(n * p), ncol = p)
  Y <- sample(true_spt, n, replace = TRUE, prob = true_f0)
  out <- data.frame(Y, X)
  return(out)
}

B <- 250
BetaII <- BetaI <- matrix(NA, nrow = B, ncol = 4)
for(b in 1:B){
  sim_dat <- sim_I(100, 1)

  fit1 <- glm(Y ~ 1 + ., data = sim_dat[1:25, ], family = gaussian(link = 'logit'))
  fit2 <- glm(Y ~ 1 + ., data = sim_dat,  family = gaussian(link = 'logit'))
  coef(fit1) %>% as.numeric() -> BetaI[b, c(1, 3)]
  coef(fit2) %>% as.numeric() -> BetaI[b, c(2, 4)]

  fit1 <- gldrm(Y ~ ., data = sim_dat[1:25, ], link = 'logit')
  fit2 <- gldrm(Y ~ ., data = sim_dat,  link = 'logit')
  fit1$beta %>% as.numeric() -> BetaII[b, c(1, 3)]
  fit2$beta %>% as.numeric() -> BetaII[b, c(2, 4)]

  Y <- sim_dat$Y
  X <- cbind(1, sim_dat$X)
  out1 <- dpglm(Y[1:25], X[1:25,], 'logit', 1000)
  out2 <- dpglm(Y, X, 'logit', 1000)
}
r <- 801:1000
colMeans(out1$beta[r, ])
apply(out1$beta[r, ], 2, sd)
colMeans(out2$beta[r, ])
apply(out2$beta[r, ], 2, sd)
plot(out1$beta[, 1], type = 'l', col = 'red')
plot(out2$beta[, 1], type = 'l', col = 'red')
plot(out1$beta[, 2], type = 'l', col = 'blue')
plot(out2$beta[, 2], type = 'l', col = 'blue')

colMeans(BetaI)
apply(BetaI, 2, sd)
colMeans(BetaII)
apply(BetaII, 2, sd)

X <- sim_dat[, -1]
y <- sim_dat[, 1]
out1 <- dpglm(y, X, 'logit', 100)
out2 <- dpglm(y, X, 'logit', 100)
