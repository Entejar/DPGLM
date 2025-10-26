resample_zstar <- function(z){
  z_table <- table(z)
  zstar   <- as.numeric(names(z_table))
  nstar   <- as.numeric(z_table)
  ## write code for resampling zstar to avoid the ‘sticky clusters effect’
  
  return(list(zstar = zstar, nstar = nstar))
}


