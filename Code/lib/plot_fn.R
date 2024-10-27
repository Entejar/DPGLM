plot_fn <- function(out, F0_true, save, yGrid, mu0, c0){
  err <- numeric(0)
  basedensity_200 <- array(NA, dim = c(length(out), length(yGrid), length(save)))
  for (i in 1:length(out)) {
    tryCatch({
      basedensity_200[i,,] <- get_F0(out[[i]]$mcmc_samples$MS1, save, yGrid, mu0, c0)
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
      basedensity_400[i,,] <- get_F0(out[[i]]$mcmc_samples$MS2, save, yGrid, mu0, c0)
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