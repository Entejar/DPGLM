#include "functions.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double logit(double x) {
  return log(x) - log(1.0-x);
}

double expit(double x) {
  return 1.0 / (1.0 + exp(-x));
}

double cdf_normal(double x, double mean, double sd) {
  return 0.5 * (1.0 + erf((x - mean) / (sd * sqrt(2.0))));
}

// -------------------- Log-pdf functions --------------------

double logpdf_beta(double x, double shape1, double shape2) {
  return (shape1-1.0) * log(x) + (shape2-1.0) * log(1 - x) - Rf_lbeta(shape1,shape2);
}

// [[Rcpp::export]]
double logpdf_unif(double x, double lower, double upper){
    return -log(upper - lower);
}

double logpdf_gamma(double x, double shape, double scale) {
  return (shape - 1.0) * log(x) - x / scale - shape * log(scale) - R::lgammafn(shape);

}

double logpdf_normal(double x, double mean, double sd){
  return -0.5 * log(2.0 * M_PI) - log(sd) - 0.5 * pow((x - mean) / sd, 2.0); // M_PI is \pi = 3.14..
}


double logpdf_truncnormal(double x, double mean, double sd, double lower, double upper){
  return logpdf_normal(x, mean, sd) - log(cdf_normal(upper, mean, sd) - cdf_normal(lower, mean, sd));
}

// -------------------- Random variate generators --------------------

// [[Rcpp::export]]
arma::vec rmvnorm(const arma::vec& mean, const arma::mat& Precision) {
  arma::vec z = arma::zeros<arma::vec>(mean.size());
  for(int i = 0; i < mean.size(); i++) {
    z(i) = norm_rand();
  }
  arma::mat Sigma = inv_sympd(Precision);
  arma::mat L = chol(Sigma, "lower");
  arma::vec h = mean + L * z;
  /* arma::mat R = chol(Precision); */
  /* arma::vec h = solve(R,z) + mean; */
  return h;
}

// [[Rcpp::export]]
int rcategorical(const arma::vec& probs) {
  double U = R::unif_rand();
  double cdf = probs(0);
  int k = 0;
  while(U > cdf) {
    k++;
    cdf += probs(k);
  }
  return k;
}

// [[Rcpp::export]]
arma::mat mvrnormArma(const int& R, const arma::vec& mu, const arma::mat& sigma) {
  int n = mu.size();
  arma::mat Y = arma::randn(R, n);
  return arma::repmat(mu, 1, R).t() + Y * arma::chol(sigma);
}
