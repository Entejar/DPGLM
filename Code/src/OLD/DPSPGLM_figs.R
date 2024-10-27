library(truncnorm)
library(ggplot2)
library(plot3D)
library(patchwork)
library(latex2exp)
library(tidyverse)

# Fig: Posterior Levy intensity at an MCMC iteration
z <- seq(0.25, 1, 0.005)
s <- seq(1.420549e-05, 5e-04, length.out = 100)
zs <- expand.grid(z = z, s = s)

#theta <- tht # init$theta
psi <- function(z){
  z <- rep(z, length(theta))
  sum(u * exp(theta * z - mx))
}

mu0 <- mean(y)
sigma0 <- sd(y)

nuZS <- function(zs){
  nu_vec <- numeric(0)
  for(i in 1:dim(zs)[1]){
    psi.z <- psi(zs[i, 1])
    dz <- dtruncnorm(zs[i, 1], mean = mu0, sd = sigma0, a = 0, b = 1)
    temp <- exp( - (1 + psi.z) * zs[i, 2] ) / zs[i, 2]
    nu_vec[i] <- temp * dz
  }
  nu_vec
}
nustar <- nuZS(zs)
df <- cbind(zs, nustar = nustar)

fzs <- matrix(nustar, nrow = length(z), ncol = length(s), byrow = F)
persp3D(z, s, fzs, theta=30, phi=20, axes=T, scale=2, box=T, nticks=3, ticktype = "simple",
        xlab="z", ylab="s", zlab="nustar", contour = F, shade = 0.3)

ggplot() +
  geom_contour(data = df, mapping = aes(x = z, y = s, z = nustar))

## Fig: Quantile trajectory estimates (posterior mean) with adding uncertainties / spaghettis
age <- seq(-1.653570, 2.281983, 0.2)
age <- ns(age, df = 3)
Xtest <- matrix(cbind(1, age), nrow = nrow(age), ncol = 4)
sex <- c(rep(1, nrow(Xtest)), rep(0, nrow(Xtest)), rep(1, nrow(Xtest)), rep(0, nrow(Xtest)))
intelligibility_type <- c(rep(1, nrow(Xtest)), rep(1, nrow(Xtest)), rep(0, nrow(Xtest)), rep(0, nrow(Xtest)))
Xtest <- cbind(rbind(Xtest, Xtest, Xtest, Xtest), sex, intelligibility_type)
colnames(Xtest) <- c("intercept", "ageKnot1", "ageKnot2", "ageKnot3", "sex", "intelligibility_type")
itr_indx <- seq(1000, 1500, 1)
perScale <- 100
resultDF <- list()
for(itr in itr_indx){
spt <- G_samples[[itr]]$z.tld
f0 <- G_samples[[itr]]$J.tld
theta0 <- theta_solver(spt, f0, mu0, NULL)$tht
f0 <- f0 * exp(theta0 * spt) / sum(f0 * exp(theta0 * spt))
l <- length(spt)
beta <- beta_samples[itr, ]
#beta <- swIntellBs$beta %>% as.numeric()
mu <- link_fn(Xtest %*% beta, link)$mu
out <- theta_solver(spt, f0, mu, NULL) 
theta <- out$tht
fTiltMat <- matrix(NA, nrow = nrow(Xtest), ncol = l)
for(i in 1:nrow(Xtest)){
  fTiltMat[i, ] <- exp(theta[i] * spt) * f0 / sum(exp(theta[i] * spt) * f0)
}

resultDF[[itr]] <- data.frame(q5 = perScale*qGeneric(.05, fTiltMat, spt), 
                       q10 = perScale*qGeneric(.10, fTiltMat, spt),
                       q25 = perScale*qGeneric(.25, fTiltMat, spt),
                       q50 = perScale*qGeneric(.50, fTiltMat, spt))
}

resultDF1 <- Reduce(`+`, resultDF[itr_indx]) / length(itr_indx)
hustadTD_temp <- data.frame(hustad[hustad$group=="TD",])
ooo <- order(hustadTD_temp$age)
hustadTD_temp <- hustadTD_temp[ooo,]
meanAGE <- mean(hustadTD_temp$age_months)
sdAGE <- sd(hustadTD_temp$age_months)

hustadTD_plt <- data.frame(age_months = hustadTD_temp$age_months, mean_intelligibility = hustadTD_temp$mean_intelligibility)

temp <- seq(-1.653570, 2.281983, 0.2)
age_months <- temp * sdAGE + meanAGE
resultDF2 <- data.frame(age = rep(rep(age_months, 4), 4), q = c(resultDF1$q5, resultDF1$q10, 
                                                 resultDF1$q25, resultDF1$q50), 
                        q_type = rep(rep(c("1", "2", "3", "4"), each = nrow(Xtest)), each = 4),
                        subgrp = rep(rep(c("1", "2", "3", "4"), each = length(temp)), 4))

resultDF2 <- data.frame(age = rep(age_months, 4), q = c(resultDF1$q5, resultDF1$q10, 
                                                                resultDF1$q25, resultDF1$q50), 
                        q_type = rep(c("1", "2", "3", "4"), each = nrow(Xtest)))

pal <- ggthemes::colorblind_pal()(8)

ggplot() +
  geom_point(data = hustadTD_plt, aes(x = age_months, y = 100 * mean_intelligibility), col = pal[4], alpha = 0.3, size = 1) +
  geom_line(data = resultDF2[resultDF2$q_type == "4", ], mapping = aes(x = age, y = q, col = subgrp), lwd = 0.8) +
  scale_color_manual(name  = "quantile",
                     values = pal[c(2, 3, 6, 7)],
                     labels = c("1" = "5%", "2" = "10%", "3" = "25%", "4" = "50%")) +
  theme_bw() +
  labs(x = "age (months)", y = "intelligibility (in %)")

q5mat <- matrix(NA, nrow = length(age), ncol = length(itr_indx))
for(i in 1:length(age)){
  for(j in 1:length(itr_indx)){
    q5mat[i, j] <- (resultDF[[itr_indx[j]]])[i, 1]
  }
}

q5_low <- apply(q5mat, 1, quantile, prob = 0.05)
q5_high <- apply(q5mat, 1, quantile, prob = 0.95)

q10mat <- matrix(NA, nrow = length(age), ncol = length(itr_indx))
for(i in 1:length(age)){
  for(j in 1:length(itr_indx)){
    q10mat[i, j] <- (resultDF[[itr_indx[j]]])[i, 2]
  }
}

q10_low <- apply(q10mat, 1, quantile, prob = 0.05)
q10_high <- apply(q10mat, 1, quantile, prob = 0.95)

q25mat <- matrix(NA, nrow = length(age), ncol = length(itr_indx))
for(i in 1:length(age)){
  for(j in 1:length(itr_indx)){
    q25mat[i, j] <- (resultDF[[itr_indx[j]]])[i, 3]
  }
}

q25_low <- apply(q25mat, 1, quantile, prob = 0.05)
q25_high <- apply(q25mat, 1, quantile, prob = 0.95)

q50mat <- matrix(NA, nrow = length(age), ncol = length(itr_indx))
for(i in 1:length(age)){
  for(j in 1:length(itr_indx)){
    q50mat[i, j] <- (resultDF[[itr_indx[j]]])[i, 4]
  }
}

q50_low <- apply(q50mat, 1, quantile, prob = 0.05)
q50_high <- apply(q50mat, 1, quantile, prob = 0.95)

q50mat1 <- q50mat2 <- q50mat3 <- q50mat4 <- matrix(NA, nrow = nrow(Xtest) / 4, ncol = length(itr_indx))
l1 <- l2 <- l3 <- l4 <- 0
for(i in 1:nrow(Xtest)){
  if(Xtest[i,5] == 1 & Xtest[i,6] == 1){
    l1 <- l1 + 1
    for(j in 1:length(itr_indx)){
      q50mat1[l1, j] <- (resultDF[[itr_indx[j]]])[i, 4]
    }
  }
  if(Xtest[i,5] == 1 & Xtest[i,6] == 0){
    l2 <- l2 + 1
    for(j in 1:length(itr_indx)){
      q50mat2[l2, j] <- (resultDF[[itr_indx[j]]])[i, 4]
    }
  }
  if(Xtest[i,5] == 0 & Xtest[i,6] == 1){
    l3 <- l3 + 1
    for(j in 1:length(itr_indx)){
      q50mat3[l3, j] <- (resultDF[[itr_indx[j]]])[i, 4]
    }
  }
  if(Xtest[i,5] == 0 & Xtest[i,6] == 0){
    l4 <- l4 + 1
    for(j in 1:length(itr_indx)){
      q50mat4[l4, j] <- (resultDF[[itr_indx[j]]])[i, 4]
    }
  }
}


low <- high <- list()
q50mat <- list(q50mat1, q50mat2, q50mat3, q50mat4)
for(i in 1:4){
  low[[i]] <- apply(q50mat[[i]], 1, quantile, prob = 0.05)
  high[[i]] <- apply(q50mat[[i]], 1, quantile, prob = 0.95)
}

resultDF4 <- data.frame(age = rep(age_months, 4), q = resultDF1$q50, low = unlist(low), 
                        high = unlist(high), subgrp = rep(c("1", "2", "3", "4"), each = length(age_months)))

resultDF3 <- cbind(resultDF2, low = c(q5_low, q10_low, q25_low, q50_low), 
                   high = c(q5_high, q10_high, q25_high, q50_high))


  
ggplot(data = resultDF4) +
  #geom_point(aes(x = age, y = intelli), col = pal[4], alpha = 0.3, size = 1) +
  geom_line(aes(x = age, y = q, col = subgrp), lwd = 0.8) +
  #geom_ribbon(aes(x = age, ymin = low, ymax = high, fill = subgrp), alpha = 0.2) +
  scale_color_manual(name = "median line", 
                     values = pal[c(2, 4, 6, 8)], 
                     labels = c("1" = "male, single-word", "2" = "male, multi-word", 
                                "3" = "female, single-word", "4" = "female, multi-word")) +
  # scale_fill_manual(name = "median line", 
  #                   values = pal[c(2, 4, 6, 8)], 
  #                   labels = c("1" = "male, single-word", "2" = "male, multi-word", 
  #                              "3" = "female, single-word", "4" = "female, multi-word")) +
  theme_bw() +
  labs(x = "age (months)", y = "intelligibility (in %)") +
  theme(legend.position = "bottom")

ggplot() +
  geom_point(data = hustadTD_plt, aes(x = age_months, y = 100 *mean_intelligibility), col = pal[4], alpha = 0.3, size = 1) +
  geom_line(data = resultDF3, aes(x = age, y = q, col = q_type), lwd = 0.8) +
  geom_ribbon(data = resultDF3, aes(x = age, ymin = low, ymax = high, fill = q_type), alpha = 0.2) +
  scale_color_manual(name = "Quantile", 
                     values = pal[c(2, 4, 6, 8)], 
                     labels = c("1" = "5%", "2" = "10%", "3" = "25%", "4" = "50%")) +
  scale_fill_manual(name = "Quantile", 
                    values = pal[c(2, 4, 6, 8)], 
                    labels = c("1" = "5%", "2" = "10%", "3" = "25%", "4" = "50%")) +
  theme_bw() +
  labs(x = "age (months)", y = "intelligibility (in %)")



#Q: tilt f0 to have mean mu0 or not, quantile plots are same! Why?
  
# Fig: survival fn + CI
age <-  c(1, 0)
yGrid <- seq(0, 1, 0.001)
itr_indx <- seq(50, 235, 1)
survivalDF <- matrix(0, nrow = length(itr_indx), ncol = length(yGrid))
for(i in 1:length(itr_indx)){
  itr <- itr_indx[i]
  spt <- G_samples[[itr]]$z.tld
  f0 <- G_samples[[itr]]$J.tld
  theta0 <- theta_solver(spt, f0, mu0, NULL)$tht
  f0 <- f0 * exp(theta0 * spt) / sum(f0 * exp(theta0 * spt))
  l <- length(spt)
  beta <- beta_samples[itr, ]
  mu <- link_fn(age %*% beta, link)$mu
  out <- theta_solver(spt, f0, mu, NULL) 
  theta <- out$tht
  fTiltMat <- matrix(0, nrow = length(yGrid), ncol = l)
  for(j in 1:length(yGrid)){
    fTiltMat[j, ] <- exp(theta * spt) * f0 / sum(exp(theta* spt) * f0)
  }
  survivalDF[i, ] <- 1 - pGeneric(yGrid, fTiltMat, spt)
}

low <- apply(survivalDF, 2, quantile, prob = 0.05)
high <- apply(survivalDF, 2, quantile, prob = 0.95)
est <- apply(survivalDF, 2, mean)
plot3DF <- data.frame(y = (yGrid * 100), low = low, high = high, est = est)

ggplot(data = plot3DF) +
  geom_line(aes(x = y, y = est), lwd = 0.8, col = pal[4]) +
  geom_ribbon(aes(x = y, ymin = low, ymax = high), fill = pal[4], alpha = 0.2) +
  theme_bw() +
  labs(x = "intelligibility (in %)", y = TeX("1 - $F_{\\widetilde{\\mu}}$")) +
  theme(axis.title.y = element_text(size = 13))


# Fig: Exceedance plot (est + CIs)
age <- quantile(X[, 2], probs = c(0.1, 0.5, 0.9)) %>% as.numeric()
Xtest <- matrix(c(rep(1, length(age)), age), nrow = length(age), ncol = 2)
yGrid <- seq(0, 1, 0.001)
itr_indx <- seq(50, 235, 1)
survivalDF <- list()
for(i in 1:length(itr_indx)){
  itr <- itr_indx[i]
  spt <- G_samples[[itr]]$z.tld
  f0 <- G_samples[[itr]]$J.tld
  theta0 <- theta_solver(spt, f0, mu0, NULL)$tht
  f0 <- f0 * exp(theta0 * spt) / sum(f0 * exp(theta0 * spt))
  l <- length(spt)
  beta <- beta_samples[itr, ]
  mu <- link_fn(Xtest %*% beta, link)$mu
  out <- theta_solver(spt, f0, mu, NULL) 
  theta <- out$tht
  DF <- matrix(0, nrow = length(age), ncol = length(yGrid))
  for(k in 1:length(age)){
    fTiltMat <- matrix(0, nrow = length(yGrid), ncol = l)
  for(j in 1:length(yGrid)){
    fTiltMat[j, ] <- exp(theta[k] * spt) * f0 / sum(exp(theta[k] * spt) * f0)
  }
    DF[k, ] <- 1 - pGeneric(yGrid, fTiltMat, spt)
  }
  survivalDF[[i]] <- DF
}

get_data_for_row <- function(row_index) {
  data <- sapply(survivalDF, function(iter) iter[row_index, ])
  colnames(data) <- names(survivalDF)
  return(data)
}

# Get the number of rows in the matrices
num_rows <- nrow(survivalDF[[1]])

# Extract data for each row
row_data <- lapply(1:num_rows, get_data_for_row)

low <- high <- est <- numeric(0)
for(i in 1:length(age)){
  low <- c(low, apply(row_data[[i]], 1, quantile, prob = 0.05))
  high <- c(high, apply(row_data[[i]], 1, quantile, prob = 0.95))
  est <- c(est, apply(row_data[[i]], 1, mean))
}

age <- age * sdAGE + meanAGE
plot5DF <- data.frame(y = rep((yGrid * 100), length(age)), low = low, high = high, est = est, age = factor(rep(age, each = length(yGrid))))

ggplot(data = plot5DF) +
  geom_line(aes(x = y, y = est, col = age), lwd = 0.8) +
  geom_ribbon(aes(x = y, ymin = low, ymax = high, fill = age), alpha = 0.2) +
  theme_bw() +
  labs(x = "intelligibility (in %)", y = TeX("$P(y \\geq intelligibility | x = age)$")) +
  theme(axis.title.y = element_text(size = 13))

# Fig: Tail fn at an MCMC iteration
plot6DF <- data.frame(v = w, Nv = Nv)
plot6DF1 <- data.frame(xi = xi[1:3], J = J[1:3], Nv = temp[1:3])
ggplot() +
  geom_line(data = plot6DF, mapping = aes(x = v, y = Nv), col = pal[4], lwd = 0.8) +
  geom_segment(data = plot6DF1, mapping = aes(x = 0, y = xi, xend = J, yend = Nv), linetype = "dashed", col = pal[6], lwd = 0.6) + 
  geom_segment(data = plot6DF1, mapping = aes(x = J, y = 0, xend = J, yend = Nv), linetype = "dashed", col = pal[6], lwd = 0.6) +
  theme_bw() +
  labs(x = TeX("$v$"), y = TeX("$N(v)$")) +
  annotate("text", x = J[1], y = -0.2, label = TeX("$J_1$"), col = "red") +
  annotate("text", x = J[2], y = -0.2, label = TeX("$J_2$"), col = "red") +
  annotate("text", x = J[3], y = -0.2, label = TeX("$J_3$"), col = "red") +
  annotate("text", y = temp[1], x = -0.00001, label = TeX("$\\xi_1$"), col = "red") +
  annotate("text", y = temp[2], x = -0.00001, label = TeX("$\\xi_2$"), col = "red") +
  annotate("text", y = temp[3], x = -0.00001, label = TeX("$\\xi_3$"), col = "red")


  
# Fig: Heatmap 

# Fig: beta comparison
