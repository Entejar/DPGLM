get_cdf <- function(y, fy) {
  # Ensure y is sorted in ascending order
  order <- order(y)
  y <- y[order]
  fy <- fy[order]
  
  # Calculate grid spacing (assuming equidistant grid)
  dy <- diff(y)[1]
  
  # Normalize fy if not already normalized
  fy_normalized <- fy / (sum(fy) * dy)
  
  # Calculate CDF using cumulative sum
  cdf <- cumsum(fy_normalized) * dy
  
  # Ensure CDF starts at 0 and ends at 1
  cdf <- cdf - cdf[1]
  cdf <- cdf / cdf[length(cdf)]
  
  return(cdf)
}