// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <random>
#include <vector>
#include <algorithm>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec stick_breaking_init(int K, double alpha) {
  // Initialize the remaining stick length
  double stick_left = 1.0;

  // Initialize vectors for stick-breaking proportions and mixture weights
  arma::vec v(K);
  arma::vec pi(K);

  // Stick-breaking process
  for (int k = 0; k < K; k++) {
    v[k] = R::rbeta(1, alpha);

    // Calculate mixture weight
    pi[k] = stick_left * v[k];

    // Update remaining stick length
    stick_left *= (1 - v[k]);
  }

  // Return mixture weights
  return pi;
}

// [[Rcpp::export]]
arma::vec stick_breaking_update(int K, double alpha, arma::uvec z) {
  // Initialize the remaining stick length
  double stick_left = 1.0;

  // Initialize vectors for stick-breaking proportions and mixture weights
  arma::vec v(K);
  arma::vec pi(K);

  // Stick-breaking process
  for (int k = 0; k < K; k++) {
    arma::uvec cluster_indices = arma::find(z == k);
    arma::uvec cluster_greater = arma::find(z > k);
    int n_k = cluster_indices.n_elem;
    int n_k_ge = cluster_greater.n_elem;
    v[k] = R::rbeta(1 + n_k, alpha + n_k_ge);

    // Calculate mixture weight
    pi[k] = stick_left * v[k];

    // Update remaining stick length
    stick_left *= (1 - v[k]);
  }

  // Return mixture weights
  return pi;
}

// [[Rcpp::export]]
arma::vec sample_cluster_assignments(arma::vec data, arma::vec pi, arma::vec mu, arma::vec sigma2) {
  int n = data.size();
  int K = pi.size();

  // Random number generator
  std::random_device rd;
  std::mt19937 gen(rd());

  // Initialize vector for cluster assignments
  arma::vec z(n);

  // Assign each data point to a cluster
  for (int i = 0; i < n; i++) {
    // Calculate log probabilities for each cluster
    arma::vec log_probs(K);
    for (int k = 0; k < K; k++) {
      log_probs[k] = std::log(pi[k]) + R::dnorm(data[i], mu[k], std::sqrt(sigma2[k]), true);
    }

    // Convert log probabilities to probabilities
    arma::vec probs = exp(log_probs - max(log_probs));
    probs = probs / arma::sum(probs);

    // Sample cluster assignment
    std::discrete_distribution<> d(probs.begin(), probs.end());
    z[i] = d(gen);
  }

  return z;
}

// [[Rcpp::export]]
Rcpp::List sample_cluster_parameters(arma::vec data, arma::uvec z, arma::vec sigma2,
                                     arma::vec mu, int K, double mu0, double sigma02,
                                     double a_sigma, double b_sigma) {
  // Update parameters for each cluster
  for (int k = 0; k < K; k++) {
    // Find data points assigned to current cluster and calculate n_k
    arma::uvec cluster_indices = arma::find(z == k);
    int n_k = cluster_indices.n_elem;

    if (n_k > 0) {
      // Extract data for current cluster
      arma::vec cluster_data = data.elem(cluster_indices);

      // Cluster mean and SS
      double ybar_k = arma::mean(cluster_data);
      double ss_k = arma::sum(arma::pow(cluster_data - ybar_k, 2.0));

      // Update mu using full conditional
      double sigma2_kn = 1.0 / (1.0 / sigma02 + n_k / sigma2[k]);
      double mu_kn = sigma2_kn * (mu0 / sigma02  + n_k * ybar_k / sigma2[k]);
      mu[k] = R::rnorm(mu_kn, std::sqrt(sigma2_kn));

      // Update sigma2 using full conditional
      double a_kn = a_sigma + 0.5 * n_k;
      double b_kn = b_sigma + 0.5 * (ss_k + n_k * std::pow(ybar_k - mu[k], 2));
      sigma2[k] = 1 / R::rgamma(a_kn, 1 / b_kn);
    } else {
      // Sample from prior for empty clusters
      mu[k] = R::rnorm(mu0, std::sqrt(sigma02));
      sigma2[k] = 1 / R::rgamma(a_sigma, 1 / b_sigma);
    }
  }

  // Return updated mu and sigma2
  return Rcpp::List::create(
    Rcpp::Named("mu") = mu,
    Rcpp::Named("sigma2") = sigma2
  );
}


// [[Rcpp::export]]
double update_alpha(double alpha, int n, int K, double a, double b) {
  // Based on Escobar and West [1995]
  // a, b: shape and rate parameter for gamma prior on alpha. Escobar and West suggest a = 2, b = 4

  // Sample eta
  double eta = R::rbeta(alpha + 1, n);

  // Calculate mixture probs for Gamma
  double pi_eta = (a + K - 1) / (n * (b - std::log(eta)) + a + K - 1);

  // Sample from the mixture
  if (R::runif(0, 1) < pi_eta) {
    // Sample from the first Gamma component
    alpha = R::rgamma(a + K, 1 / (b - std::log(eta)));
  } else {
    // Sample from the second Gamma component
    alpha = R::rgamma(a + K - 1, 1 / (b - std::log(eta)));
  }

  return alpha;
}

// [[Rcpp::export]]
double update_mu0(arma::vec mu, double sigma02) {
  // Update mu0 based on cluster means
  // p(mu0) \propto 1. Q: should we instead take N(0, 1) for centered data?

  int K = mu.n_elem;
  double mu_bar = arma::mean(mu);
  return R::rnorm(mu_bar, std::sqrt(sigma02 / K));
}

// [[Rcpp::export]]
double update_sigma02(const arma::vec& mu, double mu0, double a0, double b0) {
  // sigma02 ~ Inverse-Gamma(a0, b0)

  int K = mu.n_elem;

  // Calculate updated shape parameter
  double a0_K = a0 + K / 2.0;

  // Calculate updated rate parameter
  double b0_K = b0 + 0.5 * arma::sum(arma::square(mu - mu0));

  // Sample from Inverse-Gamma distribution
  double sigma02 = 1.0 / R::rgamma(a0_K, 1.0 / b0_K);      // R::rgamma uses shape and scale

  return sigma02;
}


double log_posterior_a_sigma(double a_sigma, double shape, double rate, const arma::vec& sigma2, double b_sigma) {
  // Prior on a_sigma: Gamma(shape, rate)
  int K = sigma2.n_elem;
  double log_post = (shape - 1) * std::log(a_sigma) - rate * a_sigma  // Log of Gamma prior
  + K * (a_sigma * std::log(b_sigma) - lgamma(a_sigma))  // Log-likelihood
    - a_sigma * arma::sum(arma::log(sigma2));  // Log-likelihood
    return log_post;
}

// [[Rcpp::export]]
double slice_sampler_a_sigma(double a_sigma, double shape, double rate, const arma::vec& sigma2, double b_sigma, double w = 0.1, int m = 20) {
  // Step 1: Evaluate the log-posterior at the current value
  double log_z = log_posterior_a_sigma(a_sigma, shape, rate, sigma2, b_sigma) - R::exp_rand();

  // Step 2: Create a horizontal interval (L, R) around the current value a_sigma
  double L = a_sigma - w * R::runif(0, 1);
  if(L < 0.0) {L = 0.0;}  // Ensure that L is non-negative as a_sigma > 0
  double R = L + w;

  // Step 3: Expand the interval until it contains the slice or
  // reaches the maximum steps J for L and Q for R
  int J = floor(m * R::runif(0, 1));
  int Q = (m - 1) - J;

  while (J > 0 && log_z < log_posterior_a_sigma(L, shape, rate, sigma2, b_sigma)) {
    L -= w;
    J--;
  }

  while (Q > 0 && log_z < log_posterior_a_sigma(R, shape, rate, sigma2, b_sigma)) {
    R += w;
    Q--;
  }

  // Step 4: Sample uniformly from the interval and continuously shrink the
  // interval until arrive at the desired set := {a_sigma_val: f(a_sigma_val) >=
  // exp(log_z) = z ~ unif(0, f(a_sigma))} or its subset
  if(L < 0.0) {L = 0.0;}  // Ensure that L is non-negative as a_sigma > 0
  double a_sigma_new = R::runif(L, R);
  while (log_z > log_posterior_a_sigma(a_sigma_new, shape, rate, sigma2, b_sigma)) {
    if (a_sigma_new < a_sigma & a_sigma_new > 0.0) {
      L = a_sigma_new;
    } else if (a_sigma_new > a_sigma){
      R = a_sigma_new;
    }
    a_sigma_new = R::runif(L, R);
  }

  return a_sigma_new;  // f(a_sigma_new) >= exp(log_z) = z ~ unif(0, f(a_sigma))
}

// [[Rcpp::export]]
double update_b_sigma(double d0, double e0, arma::vec sigma2, double a_sigma) {
  // b_sigma ~ Gamma(d0, e0)

  int K = sigma2.n_elem;

  // Update shape parameter
  double shape = d0 + K * a_sigma;   // + K * a_sigma; This was a mistake in the previous code

  // Update rate parameter
  double rate = e0 + arma::sum(1.0 / sigma2);

  // Sample from Gamma distribution
  double b_sigma = R::rgamma(shape, 1.0 / rate);  // R::rgamma uses shape and scale

  return b_sigma;
}

// Update the main sampling function
// [[Rcpp::export]]
Rcpp::List dirichlet_process_mixture_sampler(arma::vec data, int n_iter, double alpha,
                                             double mu0, double sigma02, double a, double b, double a0, double b0,
                                             double c0, double d0, double e0, double f0,
                                             int max_clusters) {
  int n = data.n_elem;

  // Initialize
  arma::vec pi_temp = stick_breaking_init(max_clusters, alpha);
  arma::vec pi = stick_breaking_init(max_clusters, alpha);
  arma::vec mu(max_clusters, arma::fill::value(mu0));
  arma::vec sigma2_temp(max_clusters, arma::fill::value(sigma02));
  arma::vec sigma2(max_clusters, arma::fill::value(sigma02));
  arma::uvec z = arma::conv_to<arma::uvec>::from(sample_cluster_assignments(data, pi, mu, sigma2));

  double a_sigma = c0 / f0;
  double b_sigma = d0 / e0;
  // Storage for results
  Rcpp::List results(n_iter);

  // Main sampling loop
  for (int iter = 0; iter < n_iter; iter++) {
    // Update cluster assignments
    z = arma::conv_to<arma::uvec>::from(sample_cluster_assignments(data, pi, mu, sigma2));

    // Update cluster parameters
    Rcpp::List cluster_params = sample_cluster_parameters(data, z, sigma2, mu, max_clusters, mu0, sigma02, a_sigma, b_sigma);
    mu = Rcpp::as<arma::vec>(cluster_params["mu"]);
    sigma2_temp = Rcpp::as<arma::vec>(cluster_params["sigma2"]);

    // if (sigma2_temp.has_nan()) {
    //   Rcpp::Rcerr << "NaN detected in sigma2 at iteration " << iter << std::endl;
    //   sigma2_temp = sigma2;
    // }

    sigma2 = sigma2_temp;

    // Update hyper-parameters
    double alpha_temp = update_alpha(alpha, n, max_clusters, a, b);
    if(alpha_temp > 0.0) {
      alpha = alpha_temp;
    }

    mu0 = update_mu0(mu, sigma02);

    double sigma02_temp = update_sigma02(mu, mu0, a0, b0);
    if(sigma02_temp > 0.0) {
      sigma02 = sigma02_temp;
    }

    double a_sigma_temp = slice_sampler_a_sigma(a_sigma, c0, f0, sigma2, b_sigma);
    if(a_sigma_temp > 0.0) {
      a_sigma = a_sigma_temp;
    }

    double b_sigma_temp = update_b_sigma(d0, e0, sigma2, a_sigma);
    if(b_sigma_temp > 0.0) {
      b_sigma = b_sigma_temp;
    }


    // Update stick-breaking weights
    pi_temp = stick_breaking_update(max_clusters, alpha, z);
    // if (pi_temp.has_nan()) {
    //   Rcpp::Rcerr << "NaN detected in pi at iteration " << iter << std::endl;
    //   pi_temp = pi;
    // }

    pi = pi_temp;

    // Store results
    results[iter] = Rcpp::List::create(
      Rcpp::Named("z") = z,
      Rcpp::Named("pi") = pi,
      Rcpp::Named("mu") = mu,
      Rcpp::Named("sigma2") = sigma2,
      Rcpp::Named("alpha") = alpha,
      Rcpp::Named("mu0") = mu0,
      Rcpp::Named("sigma02") = sigma02,
      Rcpp::Named("a_sigma") = a_sigma,
      Rcpp::Named("b_sigma") = b_sigma
    );
  }

  return results;
}
