#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

// Function to compute b(theta_i)
// [[Rcpp::export]]
double b_theta(const arma::vec& s_k, const arma::vec& f_k, double theta_i) {
  arma::vec exp_term = arma::exp(theta_i * s_k);
  double sum = arma::dot(f_k, exp_term);
  return std::log(sum);
}

// Function to compute b'(theta_i)
// [[Rcpp::export]]
double b_prime_theta(const arma::vec& s_k, const arma::vec& f_k, double theta_i) {
  arma::vec exp_term = arma::exp(theta_i * s_k);
  double numerator = arma::dot(s_k % f_k, exp_term);
  double denominator = arma::dot(f_k, exp_term);
  return numerator / denominator;
}

// Function to compute L(theta_i)
// [[Rcpp::export]]
double L_theta(const arma::vec& s_k, const arma::vec& f_k, double theta_i,
               double c_0, double y) {
  arma::vec exp_term = arma::exp(theta_i * s_k); // exp(theta_i * s_k)
  double sum = arma::dot(f_k, exp_term); // sum(f_k * exp(theta_i * s_k))
  arma::vec p_z = exp_term * f_k / sum; // p(z = s_k) = exp(theta_i * s_k) / sum(f_k * exp(theta_i * s_k))
  arma::uvec indices = arma::find(s_k - c_0 < y && y < s_k + c_0); // indices of s_k such that s_k - c_0 < y < s_k + c_0
  return arma::sum(p_z(indices)) / (2 * c_0); // sum_{s_k} (K(y | z = s_k) * p(z = s_k | theta_i,..))
}


// Function to compute logL(beta)
// [[Rcpp::export]]
double llik_beta(const arma::vec& y, const arma::mat& x, const arma::vec& theta,
              const arma::mat& crm, double c_0, int B, double sigma_theta) {
  double llik = 0.0;
  for(int i = 0; i < y.size(); ++i) {
    double L_theta_i = 0.0;
    for(int b = 0; b < B; ++b) {
      double t_theta_i = theta(i) + R::rnorm(0, sigma_theta);
      L_theta_i += L_theta(crm.row(0).t(), crm.row(1).t(), t_theta_i, c_0, y(i));
    }
    llik += std::log(L_theta_i / B);
  }
  return llik;
}
