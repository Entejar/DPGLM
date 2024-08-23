#ifndef HELPERS_H
#define HELPERS_H

#include <RcppArmadillo.h>
#include "functions.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// -------------------------- Log-posterior of latent u -------------------------- //

double log_posterior_u(const arma::vec& u, const arma::vec& z, const arma::vec& theta,
                       const double& alpha);

#endif
