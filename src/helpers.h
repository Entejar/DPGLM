#ifndef HELPERS_H
#define HELPERS_H

#include <RcppArmadillo.h>
#include "functions.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// -------------------------- Log-posterior of latent u -------------------------- //

double log_posterior_u(const arma::vec& u, const arma::vec& z, const arma::vec& theta,
                       const double& alpha);


// -------------------------- For CRM update -------------------------- //
double psi(double z, arma::vec u, arma::vec tht);
double psi_zstar(double zstar, arma::vec u, arma::vec tht);

#endif
