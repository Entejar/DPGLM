g <- function(mu, link){
  if(link == 'logit'){
    logit(mu)
  }
  if(link == 'log'){
    log(mu)
  } else {
    stop('Link function not included. Add it to g.R')
  }
}