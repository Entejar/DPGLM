#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

static double LN_2PI = 1.83787706640934548356065947281;

double logit(double x);

double expit(double x);

double cdf_normal(double x, double mean, double sd);

// ------------------------- log-pdf  -----------------------------

double logpdf_beta(double x, double a, double b);

double logpdf_gamma(double x, double a, double b);

double logpdf_normal(double x, double mean, double sd);

double logpdf_truncnormal(double x, double mean, double sd, double lower, double upper);

double logpdf_unif(double x, double lower, double upper);

// ------------------- random variate generation -------------------

arma::vec rmvnorm(const arma::vec& mean, const arma::mat& Precision);

int rcategorical(const arma::vec& probs);

arma::mat mvrnormArma(const int& R, const arma::vec& mu, const arma::mat& sigma);

#endif