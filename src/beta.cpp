#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

// Function to compute b(theta_i)
double b_theta(const arma::vec& s_k, const arma::vec& f_k, double theta_i) {
  arma::vec exp_term = arma::exp(theta_i * s_k);
  double sum = arma::accu(f_k % exp_term);
  return log(sum); // log(sum(f_k * exp(theta_i * s_k)))
}

// Function to compute b(theta)
// [[Rcpp::export]]
arma::vec b_theta(const arma::vec& s_k, const arma::vec& f_k, 
                  const arma::vec& theta) {
  arma::vec b_theta_vec(theta.size());
  for(int i = 0; i < theta.size(); ++i) {
    b_theta_vec(i) = b_theta(s_k, f_k, theta(i));
  }
  return b_theta_vec;
}

// Function to compute b'(theta_i)
// [[Rcpp::export]]
double b_prime_theta(const arma::vec& s_k, const arma::vec& f_k, double theta_i) {
  arma::vec exp_term = arma::exp(theta_i * s_k);
  double numerator = arma::accu(s_k % f_k % exp_term); // sum(s_k * f_k * exp(theta_i * s_k))
  double denominator = arma::accu(f_k % exp_term); // sum(f_k * exp(theta_i * s_k))
  return numerator / denominator;
}

// Function to compute L(theta_i)
// [[Rcpp::export]]
double L_theta(const arma::vec& s_k, const arma::vec& f_k, double theta_i,
               double c_0, double y) {
  arma::vec exp_term = arma::exp(theta_i * s_k); // exp(theta_i * s_k)
  double sum = arma::accu(f_k % exp_term); // sum(f_k * exp(theta_i * s_k)
  arma::vec p_z = f_k % exp_term / sum; // p(z = s_k) = f_k * exp(theta_i * s_k) / sum(f_k * exp(theta_i * s_k))
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
    llik += log(L_theta_i / B);
  }
  return llik;
}
