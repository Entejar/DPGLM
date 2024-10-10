#include <RcppArmadillo.h>
#include "helpers.h"

// -------------------------- Log-posterior of latent u -------------------------- //

// [[Rcpp::export]]
double log_posterior_u(const arma::vec& u, const arma::vec& z, const arma::vec& theta,
                       const double& alpha) {
  int R = 25000;
  double eps = 1e-6;
  arma::vec zstar = arma::join_cols(arma::linspace<arma::vec>(0+eps, 1-eps, R), z); // R + length(z) = R + n vector
  arma::mat exp_mat = arma::exp(arma::diagmat(theta) * arma::repmat(zstar, 1, theta.size()).t());  // n x (R + n) matrix
  arma::mat u_exp_mat = arma::diagmat(u) * exp_mat;  // n x (R + n) matrix
  arma::vec log_1_plus_u_exp_mat_sum = arma::log(1 + arma::sum(u_exp_mat, 0).t());  // (R + n) x 1 vector
  double neg_log_post = alpha * mean(log_1_plus_u_exp_mat_sum.head(R)) +
                               mean(log_1_plus_u_exp_mat_sum.tail(z.size()));
  return - neg_log_post;
}

// ------------------------------- For CRM update ------------------------------ //
//[[Rcpp::export]]
double psi(double z, arma::vec u, arma::vec tht) {
  double mx = arma::max(tht * z);
  return arma::accu(u % arma::exp(tht * z - mx));
}

//[[Rcpp::export]]
double psi_zstar(double zstar, arma::vec u, arma::vec tht) {
  return arma::accu(u % arma::exp(tht * zstar));
}


