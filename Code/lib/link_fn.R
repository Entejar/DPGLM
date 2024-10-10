link_fn <- function(v, link){
  # for logit link
  if(link == "logit")
  {
    mu <- exp(v) / (1 + exp(v))
    gpr <- 1 / (v * (1 - v))
  }
  # for log link
  if(link == "log"){
  mu <- exp(v)
  gpr <- 1 / v
  }
  # for identity link
  if(link == "identity"){
    mu <- v
    gpr <- 1
  }

  return(list(mu = mu, gpr = gpr))
}
