#ifndef CPP_FUNCTIONS_H
#define CPP_FUNCTIONS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


static double LN_2PI = 1.83787706640934548356065947281;

double logit(double x);

double expit(double x);

double cdf_normal(double x, double mean, double sd);

// ------------------------- log-pdfs  -----------------------------

double logpdf_beta(double x, double a, double b);

double logpdf_gamma(double x, double a, double b);

double logpdf_normal(double x, double mean, double sd);

double logpdf_truncnormal(double x, double mean, double sd, double lower, double upper);

double logpdf_unif(double x, double lower, double upper);

// ------------------- random variate generation -------------------

arma::vec rmvnorm(const arma::vec& mean, const arma::mat& Precision);

int rcategorical(const arma::vec& probs);

arma::mat mvrnormArma(const int& R, const arma::vec& mu, const arma::mat& sigma);

// -------------------------- For beta update -------------------------- //
arma::vec b_theta(const arma::vec& theta, const arma::vec& spt, const arma::vec& f0);


// -------------------------- For CRM update -------------------------- //
double psi(double z, arma::vec u, arma::vec tht);
double psi_zstar(double zstar, arma::vec u, arma::vec tht);

arma::vec stick_breaking_init(int K, double alpha);       // For initializing first part of CRM

List crm_sampler_fractionalY(const int& M, const arma::vec& u, arma::vec zstar, 
                             arma::vec nstar, arma::vec RL, 
                             arma::vec RJ, const arma::vec& tht, const double& alpha, 
                             const arma::vec& mu, const arma::vec& y,
                             const float& shape1, const float& shape2);

// -------------------------- For u update -------------------------- //

double log_posterior_u(const arma::vec& u, const arma::vec& z, const arma::vec& theta,
                       const double& alpha);

arma::vec u_sampler(arma::vec& u, const arma::vec& theta, const arma::vec& z,
                    const double& alpha, const double& delta);

#endif