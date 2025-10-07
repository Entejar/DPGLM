#include <Rcpp.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

// f = b_prime (in our notation)

// Function to compute f(theta_i)
double compute_f(double theta, const std::vector<double>& z, const std::vector<double>& J) {
  int M = z.size();
  double numerator = 0.0;
  double denominator = 0.0;
  
  for (int m = 0; m < M; ++m) {
    double exp_term = std::exp(theta * z[m]);
    numerator += z[m] * exp_term * J[m];
    denominator += exp_term * J[m];
  }
  
  return numerator / denominator;
}

// Derivative of f(theta_i) with respect to theta_i (First derivative)
double compute_f_prime(double theta, const std::vector<double>& z, const std::vector<double>& J) {
  int M = z.size();
  double numerator = 0.0;
  double denominator = 0.0;
  double f = compute_f(theta, z, J); // f(theta)
  
  for (int m = 0; m < M; ++m) {
    double exp_term = std::exp(theta * z[m]);
    numerator += z[m] * z[m] * exp_term * J[m];
    denominator += exp_term * J[m];
  }
  
  return numerator / denominator - f * f;
}

// Second derivative of f(theta_i) with respect to theta_i (Hessian)
double compute_f_double_prime(double theta, const std::vector<double>& z, const std::vector<double>& J) {
  int M = z.size();
  double numerator = 0.0;
  double denominator = 0.0;
  double f_prime = compute_f_prime(theta, z, J);  // First derivative
  double f = compute_f(theta, z, J);              // f(theta)
  
  for (int m = 0; m < M; ++m) {
    double exp_term = std::exp(theta * z[m]);
    numerator += z[m] * z[m] * z[m] * z[m] * exp_term * J[m];
    denominator += exp_term * J[m];
  }
  
  // Hessian
  return numerator / denominator - 2 * f_prime * f + f * f;
}

// Newton-Raphson solver for f(theta_i) - mu_i = 0
double solve_theta(double mu, const std::vector<double>& z, const std::vector<double>& J, 
                   double tol = 1e-8, int max_iter = 100) {
  double theta = 0.0; // Initial guess
  for (int iter = 0; iter < max_iter; ++iter) {
    double f = compute_f(theta, z, J);
    double f_prime = compute_f_prime(theta, z, J);
    double step = (f - mu) / f_prime;
    
    theta -= step;
    
    if (std::fabs(step) < tol) {
      return theta; // Converged
    }
  }
  
  Rcpp::stop("Newton-Raphson did not converge for mu = " + std::to_string(mu));
  return theta; // Should not reach here
}

// [[Rcpp::export]]
NumericVector solve_thetas(const NumericVector& mu, const NumericVector& z, 
                           const NumericVector& J, double tol = 1e-8, int max_iter = 100) {
  int n = mu.size();
  std::vector<double> z_vec(z.begin(), z.end());
  std::vector<double> J_vec(J.begin(), J.end());
  NumericVector thetas(n);
  
  for (int i = 0; i < n; ++i) {
    // Solve for theta_i
    thetas[i] = solve_theta(mu[i], z_vec, J_vec, tol, max_iter);
  }
  
  return thetas;
}

// [[Rcpp::export]]
double log_post_beta(const NumericVector& beta, const NumericVector& z, 
                     const NumericVector& J, const NumericMatrix& X, 
                     const NumericVector& theta,
                     const double& sigma_theta,
                     const double& sigma_beta) {
  int n = X.nrow();          // Number of rows (observations)
  int p = X.ncol();          // Number of columns (features)
  double tol = 1e-8;         // Tolerance for solving theta
  int max_iter = 100;        // Maximum iterations for solving theta
  double log_post = 0.0;     // Log-posterior value
  NumericVector mu(n);       // Store mu values
  
  // Convert z and J to std::vector for solve_theta
  std::vector<double> z_vec(z.begin(), z.end());
  std::vector<double> J_vec(J.begin(), J.end());
  
  // Compute mu = exp(X %*% beta) / (1 + exp(X %*% beta))
  for (int i = 0; i < n; ++i) {
    double linear_pred = 0.0;
    for (int j = 0; j < p; ++j) {
      linear_pred += X(i, j) * beta[j];
    }
    mu[i] = std::exp(linear_pred) / (1.0 + std::exp(linear_pred));
  }
  
  // Compute the log-posterior
  for (int i = 0; i < n; ++i) {
    double theta_tilde = solve_theta(mu[i], z_vec, J_vec, tol, max_iter); // Solve for theta_i
    log_post -= std::pow(theta[i] - theta_tilde, 2) / (2.0 * sigma_theta * sigma_theta);
  }
  
  for (int j = 0; j < p; ++j) {
    log_post -= std::pow(beta[j], 2) / (2.0 * sigma_beta * sigma_beta); // Prior on beta
  }
  
  return log_post;
}

// Compute the gradient (first derivative of log-posterior with respect to beta)
// [[Rcpp::export]]
NumericVector log_post_grad(const NumericVector& beta, const NumericVector& z, 
                            const NumericVector& J, const NumericMatrix& X, 
                            const NumericVector& theta,
                            const double& sigma_theta,
                            const double& sigma_beta) {
  int n = X.nrow();
  int p = X.ncol();
  double tol = 1e-8;
  int max_iter = 100;
  NumericVector grad(p); // Gradient vector
  
  // Convert z and J to std::vector for solve_theta
  std::vector<double> z_vec(z.begin(), z.end());
  std::vector<double> J_vec(J.begin(), J.end());
  
  // Compute mu
  NumericVector mu(n);
  for (int i = 0; i < n; ++i) {
    double linear_pred = 0.0;
    for (int j = 0; j < p; ++j) {
      linear_pred += X(i, j) * beta[j];
    }
    mu[i] = std::exp(linear_pred) / (1.0 + std::exp(linear_pred));
  }
  
  // Compute the gradient of log-posterior
  for (int j = 0; j < p; ++j) {
    double grad_j = 0.0;
    
    for (int i = 0; i < n; ++i) {
      double theta_tilde = solve_theta(mu[i], z_vec, J_vec, tol, max_iter);
      grad_j += (theta[i] - theta_tilde) * X(i, j) / (sigma_theta * sigma_theta);
    }
    
    grad_j -= beta[j] / (sigma_beta * sigma_beta); // Prior gradient
    grad[j] = grad_j;
  }
  
  return grad;
}


// Compute Hessian matrix for the log-posterior
// [[Rcpp::export]]
NumericMatrix compute_hessian(const NumericVector& beta, const NumericVector& z, 
                         const NumericVector& J, const NumericMatrix& X, 
                         const NumericVector& theta,
                         const double& sigma_theta,
                         const double& sigma_beta) {
  int p = beta.size();
  NumericMatrix hessian(p, p);
  
  std::vector<double> z_vec(z.begin(), z.end());
  std::vector<double> J_vec(J.begin(), J.end());
  
  // Compute the second derivatives of the log-posterior
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      double hessian_ij = 0.0;
      for (int k = 0; k < X.nrow(); ++k) {
        double mu_k = 0.0;
        for (int l = 0; l < X.ncol(); ++l) {
          mu_k += X(k, l) * beta[l];
        }
        mu_k = std::exp(mu_k) / (1 + std::exp(mu_k));
        
        double theta_tilde = solve_theta(mu_k, z_vec, J_vec, 1e-8, 100);
        
        double second_derivative = compute_f_double_prime(theta_tilde, z_vec, J_vec);
        
        hessian_ij += X(k, i) * X(k, j) * second_derivative;
      }
      hessian(i, j) = -hessian_ij / (2 * sigma_theta * sigma_theta);
    }
  }
  
  return hessian;
}

// BFGS Optimization to find the mode of the posterior
// [[Rcpp::export]]
NumericVector optimize_beta(const NumericVector& beta_init, const NumericVector& z, 
                            const NumericVector& J, const NumericMatrix& X, 
                            const NumericVector& theta, 
                            const double& sigma_theta, const double& sigma_beta, 
                            int max_iter = 1000, double tol = 1e-8) {
  NumericVector beta = beta_init;
  
  for (int iter = 0; iter < max_iter; ++iter) {
    // Compute the gradient of log-posterior
    NumericVector grad = log_post_grad(beta, z, J, X, theta, sigma_theta, sigma_beta);
    
    // Update beta using BFGS (or a simple gradient step)
    for (int j = 0; j < beta.size(); ++j) {
      beta[j] -= 0.1 * grad[j];  // Learning rate of 0.1, adjust if necessary
    }
    
    // Check for convergence
    if (sum(abs(grad)) < tol) {
      break;
    }
  }
  
  return beta;
}

