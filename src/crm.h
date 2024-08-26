#ifndef CRM_H
#define CRM_H

#include <RcppArmadillo.h>
#include "functions.h"
#include "helpers.h"

List crm_sampler(const int& M, const arma::vec& u, const arma::vec& zstar, 
                 const arma::vec& nstar, const arma::vec& tht, const double& alpha, 
                 const arma::vec& mu, const arma::vec& y, arma::vec z_tld, 
                 arma::vec J_tld);
#endif