#include <RcppArmadillo.h>
#include "u.h"

// [[Rcpp::export]]
arma::vec u_sampler(const arma::vec& u_old, const arma::vec& z, const arma::vec& theta, 
                    const double& alpha, const double& delta){
  arma::vec u_star = arma::zeros<arma::vec>(u_old.size());
  double logQ_ratio = 0.0;
  for(int i = 0; i < u_old.size(); ++i){
    u_star(i) = R::rgamma(delta, u_old(i) / delta); // u_proposal := u_star ~ Gamma(shape = delta, scale = u_old / delta)  
    logQ_ratio += logpdf_gamma(u_old(i), delta, u_star(i) / delta) - 
      logpdf_gamma(u_star(i), delta, u_old(i) / delta); // log(q(u_old | u_star)) - log(q(u_star | u_old))
  }
  double logratio = log_posterior_u(u_star, z, theta, alpha) -
    log_posterior_u(u_old, z, theta, alpha) + logQ_ratio;
  if(log(R::runif(0, 1)) < logratio){
    return u_star;
  } else {
    return u_old;
  }
}


