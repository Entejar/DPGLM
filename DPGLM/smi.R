rm(list = ls())
# Load Libraries
require(foreach)
require(doParallel)
source("load.R")
library(splines)
source("src/fit_func_smi_sw.R")
source("src/fit_func_smi_mw1.R")
source("src/fit_func_smi_mw2.R")

# Parallel Setup
num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)

hustadTDMW <- readRDS("data/hustadTDMW.rds")
hustadTDSW <- readRDS("data/hustadTDSW.rds")

## Data: X, y & link function --------------------------------------------------
link <- "logit"
prbs <- c(0.33, 0.67)

Xstar <- ns(
  hustadTDMW$age_months,
  knots = quantile(hustadTDMW$age_months, probs = prbs),
  intercept = FALSE
)

X    <- cbind(rep(1, nrow(Xstar)), Xstar)
colnames(X) <- c("intercept", "ageBase1", "ageBase2", "ageBase3")
y    <- hustadTDMW$mean_intelligibility
X <- X %>% as.matrix()
y <- y %>% as.numeric()

data_mw <- data.frame(y, X)

Xstar <- ns(
  hustadTDSW$age_months,
  knots = quantile(hustadTDSW$age_months, probs = prbs),
  intercept = FALSE
)

X    <- cbind(rep(1, nrow(Xstar)), Xstar)
colnames(X) <- c("intercept", "ageBase1", "ageBase2", "ageBase3")
y    <- hustadTDSW$mean_intelligibility
X <- X %>% as.matrix()
y <- y %>% as.numeric()

data_sw <- data.frame(y, X)


out_ <- list()
n_settings <- 8
iter <- 3000

start_time <- Sys.time()

out_ <- foreach(
  set = 1:n_settings,
  .packages = c("tidyr", "gldrm", "mvtnorm")
) %dopar% {
  tryCatch({
    if(set == 1){
      dat <- data_mw
      y <- dat[, 1]
      X <- dat[, -1]
      X <- X %>% as.matrix()
      y <- y %>% as.numeric()
      c0 <- c0_silverman(y) / 4 
      gldrm_fit <- gldrm(y ~ X[, -1], link = link)
      mu_beta <- gldrm_fit$beta
      sigma_beta <- gldrm_fit$seBeta
      out <- fit_func_smi_mw1(y, X, iter, c0, mu_beta, sigma_beta) 
    }
    
    if(set == 2){
      dat <- data_mw
      y <- dat[, 1]
      X <- dat[, -1]
      X <- X %>% as.matrix()
      y <- y %>% as.numeric()
      c0 <- c0_silverman(y) / 4 
      gldrm_fit <- gldrm(y ~ X[, -1], link = link)
      mu_beta <- gldrm_fit$beta
      sigma_beta <- gldrm_fit$seBeta
      out <- fit_func_smi_mw2(y, X, iter, c0, mu_beta, sigma_beta) 
    }
    
    if(set == 3){
      dat <- data_mw
      y <- dat[, 1]
      X <- dat[, -1]
      X <- X %>% as.matrix()
      y <- y %>% as.numeric()
      c0 <- 0.025
      gldrm_fit <- gldrm(y ~ X[, -1], link = link)
      mu_beta <- gldrm_fit$beta
      sigma_beta <- gldrm_fit$seBeta
      out <- fit_func_smi_mw1(y, X, iter, c0, mu_beta, sigma_beta) 
    }
    
    
    if(set == 4){
      dat <- data_mw
      y <- dat[, 1]
      X <- dat[, -1]
      X <- X %>% as.matrix()
      y <- y %>% as.numeric()
      c0 <- 0.025 
      gldrm_fit <- gldrm(y ~ X[, -1], link = link)
      mu_beta <- gldrm_fit$beta
      sigma_beta <- gldrm_fit$seBeta
      out <- fit_func_smi_mw2(y, X, iter, c0, mu_beta, sigma_beta) 
    }
    
    if(set == 5){
      dat <- data_sw
      y <- dat[, 1]
      X <- dat[, -1]
      X <- X %>% as.matrix()
      y <- y %>% as.numeric()
      c0 <- c0_silverman(y) / 4 
      gldrm_fit <- gldrm(y ~ X[, -1], link = link)
      mu_beta <- gldrm_fit$beta
      sigma_beta <- gldrm_fit$seBeta
      out <- fit_func_smi_mw1(y, X, iter, c0, mu_beta, sigma_beta) 
    }
    
    if(set == 6){
      dat <- data_sw
      y <- dat[, 1]
      X <- dat[, -1]
      X <- X %>% as.matrix()
      y <- y %>% as.numeric()
      c0 <- c0_silverman(y) / 4 
      gldrm_fit <- gldrm(y ~ X[, -1], link = link)
      mu_beta <- gldrm_fit$beta
      sigma_beta <- gldrm_fit$seBeta
      out <- fit_func_smi_mw2(y, X, iter, c0, mu_beta, sigma_beta) 
    }
    
    if(set == 7){
      dat <- data_sw
      y <- dat[, 1]
      X <- dat[, -1]
      X <- X %>% as.matrix()
      y <- y %>% as.numeric()
      c0 <- 0.025 
      gldrm_fit <- gldrm(y ~ X[, -1], link = link)
      mu_beta <- gldrm_fit$beta
      sigma_beta <- gldrm_fit$seBeta
      out <- fit_func_smi_mw1(y, X, iter, c0, mu_beta, sigma_beta) 
    }
    
    if(set == 8){
      dat <- data_sw
      y <- dat[, 1]
      X <- dat[, -1]
      X <- X %>% as.matrix()
      y <- y %>% as.numeric()
      c0 <- 0.025 
      gldrm_fit <- gldrm(y ~ X[, -1], link = link)
      mu_beta <- gldrm_fit$beta
      sigma_beta <- gldrm_fit$seBeta
      out <- fit_func_smi_mw2(y, X, iter, c0, mu_beta, sigma_beta) 
    }

    return(out)
  }, error = function(e) {
    list(error = TRUE, message = as.character(e))
  })
}

# Save Results
saveRDS(out_, file = "cache/final_smi_data_mar20.rds")

# Stop Parallel Backend
stopCluster(cl)

# Compute Execution Time
time_tot <- Sys.time() - start_time
print(time_tot)

