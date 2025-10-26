#include "Cpp_functions.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// --------------------------------------------------------------------------- //
//                          Basic functions                                   //
// --------------------------------------------------------------------------- // 

//[[Rcpp::export]]
double logit(double x) {
  return log(x) - log(1.0-x);
}

double expit(double x) {
  return 1.0 / (1.0 + exp(-x));
}

//[[Rcpp::export]]
arma::vec expit(const arma::vec& x) {
  return 1.0 / (1.0 + exp(-x));
}

double cdf_normal(double x, double mean, double sd) {
  return 0.5 * (1.0 + erf((x - mean) / (sd * sqrt(2.0))));
}


// --------------------------------------------------------------------------- //
//                          Log-pdf functions                                  //
// --------------------------------------------------------------------------- // 

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

//[[Rcpp::export]]
double logpdf_mvnorm(const arma::vec& x, const arma::vec& mean, const arma::mat& sigma) {
  return -0.5 * (x.size() * log(2.0 * M_PI) + 
                 log(det(sigma)) + 
                 as_scalar((x - mean).t() * inv_sympd(sigma) * (x - mean))); 
}


double logpdf_truncnormal(double x, double mean, double sd, double lower, double upper){
  return logpdf_normal(x, mean, sd) - log(cdf_normal(upper, mean, sd) - cdf_normal(lower, mean, sd));
}



// --------------------------------------------------------------------------- //
//                         Random variate generators                          //
// --------------------------------------------------------------------------- // 

// // [[Rcpp::export]]
// arma::vec rmvnorm(const arma::vec& mean, const arma::mat& Precision) {
//   arma::vec z = arma::zeros<arma::vec>(mean.size());
//   for(int i = 0; i < mean.size(); i++) {
//     z(i) = norm_rand();
//   }
//   arma::mat Sigma = inv_sympd(Precision);
//   arma::mat L = chol(Sigma, "lower");
//   arma::vec h = mean + L * z;
//   /* arma::mat R = chol(Precision); */
//   /* arma::vec h = solve(R,z) + mean; */
//   return h;
// }

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





// --------------------------------------------------------------------------- //
//                          For beta update                                   //
// --------------------------------------------------------------------------- //                          

// Function to compute b(theta_i)
double b_theta_i(const arma::vec& s_k, const arma::vec& f_k, double theta_i) {
  arma::vec exp_term = arma::exp(theta_i * s_k);
  double sum = arma::accu(f_k % exp_term);
  return log(sum); // log(sum(f_k * exp(theta_i * s_k)))
}

// [[Rcpp::export]]
arma::vec b_theta(const arma::vec& theta, const arma::vec& spt, const arma::vec& f0) {
  int n = theta.n_elem;
  int m = spt.n_elem;
  
  // Matrix of exponentiated outer product of theta and spt
  arma::mat exp_matrix(n, m);
  for (int i = 0; i < n; ++i) {
    exp_matrix.row(i) = exp(theta(i) * spt.t());
  }
  
  // Apply the sum with f0 and take log
  arma::vec result(n);
  for (int i = 0; i < n; ++i) {
    result(i) = log(arma::dot(exp_matrix.row(i), f0));
  }
  
  return result;
}


// // Function to compute b(theta)
// // [[Rcpp::export]]
// arma::vec b_theta(const arma::vec& s_k, const arma::vec& f_k, 
//                   const arma::vec& theta) {
//   arma::vec b_theta_vec(theta.size());
//   for(int i = 0; i < theta.size(); ++i) {
//     b_theta_vec(i) = b_theta_i(s_k, f_k, theta(i));
//   }
//   return b_theta_vec;
// }

// Function to compute b'(theta_i)
// [[Rcpp::export]]
double b_prime_theta(const arma::vec& s_k, const arma::vec& f_k, double theta_i) {
  arma::vec exp_term = arma::exp(theta_i * s_k);
  double numerator = arma::accu(s_k % f_k % exp_term); // sum(s_k * f_k * exp(theta_i * s_k))
  double denominator = arma::accu(f_k % exp_term); // sum(f_k * exp(theta_i * s_k))
  return numerator / denominator;
}

// // For uniform kernel

// Function to compute L(theta_i)
// [[Rcpp::export]]
double L_theta_unifK(const arma::vec& s_k, const arma::vec& f_k, double theta_i,
                     double c_0, double y) {
  arma::vec exp_term = arma::exp(theta_i * s_k); // exp(theta_i * s_k)
  double sum = arma::accu(f_k % exp_term); // sum(f_k * exp(theta_i * s_k)
  arma::vec p_z = f_k % exp_term / sum; // p(z = s_k) = f_k * exp(theta_i * s_k) / sum(f_k * exp(theta_i * s_k))
  // arma::uvec indices = arma::find(s_k - c_0 < y && y < s_k + c_0); // indices of s_k such that s_k - c_0 < y < s_k + c_0
  // return arma::sum(p_z(indices)) / (2 * c_0); // sum_{s_k} (K(y | z = s_k) * p(z = s_k | theta_i,..))
  
  arma::vec lower_bound = arma::max(s_k - c_0, arma::zeros(s_k.n_elem)); // max(0, s_k - c_0)
  arma::vec width = s_k + c_0 - lower_bound;                             // width of the kernel's support
  arma::uvec indices = arma::find(lower_bound <= y && y <= s_k + c_0);    // Valid indices for y in support
  
  arma::vec K_y_z;
  K_y_z = arma::zeros(s_k.n_elem);              // Initialize K(y | z) to zero
  K_y_z(indices) = 1.0 / width(indices);                                // Adjusted uniform kernel value
  
  // Compute the weighted sum of K(y | z) * p(z = s_k)
  return arma::accu(K_y_z % p_z);
}


// Function to compute logL(beta)
// [[Rcpp::export]]
double llik_beta_unifK(const arma::vec& y, const arma::mat& x, const arma::vec& theta,
                       const arma::mat& crm, double c_0, int B, double sigma_theta) {
  double llik = 0.0;
  for(int i = 0; i < y.size(); ++i) {
    double L_theta_i = 0.0;
    for(int b = 0; b < B; ++b) {
      double t_theta_i = theta(i) + R::rnorm(0, sigma_theta);
      L_theta_i += L_theta_unifK(crm.row(0).t(), crm.row(1).t(), t_theta_i, c_0, y(i));
    }
    llik += log(L_theta_i / B);
  }
  return llik;
}

// // For exponential kernel

// Function to compute L(theta_i)
// [[Rcpp::export]]
double L_theta_expK(const arma::vec& s_k, const arma::vec& f_k, double theta_i,
                    double c_0, double y) {
  arma::vec exp_term = arma::exp(theta_i * s_k); // exp(theta_i * s_k)
  double sum = arma::accu(f_k % exp_term); // sum(f_k * exp(theta_i * s_k)
  arma::vec p_z = f_k % exp_term / sum; // p(z = s_k) = f_k * exp(theta_i * s_k) / sum(f_k * exp(theta_i * s_k))
  arma::vec diff = arma::abs(y - s_k);      // |y - z|
  arma::vec K_y_z;
  K_y_z = arma::exp(-diff / c_0) / c_0;    // Exponential kernel: exp(-|y - z| / c_0) / c_0
  
  // Compute the weighted sum of K(y | z) * p(z = s_k)
  return arma::accu(K_y_z % p_z);
}

// Function to compute logL(beta)
// [[Rcpp::export]]
double llik_beta_expK(const arma::vec& y, const arma::mat& x, const arma::vec& theta,
                      const arma::mat& crm, double c_0, int B, double sigma_theta) {
  double llik = 0.0;
  for(int i = 0; i < y.size(); ++i) {
    double L_theta_i = 0.0;
    for(int b = 0; b < B; ++b) {
      double t_theta_i = theta(i) + R::rnorm(0, sigma_theta);
      L_theta_i += L_theta_expK(crm.row(0).t(), crm.row(1).t(), t_theta_i, c_0, y(i));
    }
    llik += log(L_theta_i / B);
  }
  return llik;
}




// --------------------------------------------------------------------------- //
//                               CRM update                                   //
// --------------------------------------------------------------------------- // 

//[[Rcpp::export]]
double psi(double z, arma::vec u, arma::vec tht) {
  double mx = arma::max(tht * z);
  return arma::accu(u % arma::exp(tht * z - mx));
}

//[[Rcpp::export]]
double psi_zstar(double zstar, arma::vec u, arma::vec tht) {
  return arma::accu(u % arma::exp(tht * zstar));
}


// For initializing first part of CRM
//[[Rcpp::export]]
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
List update_crm(const int& M, const arma::vec& u, const arma::vec zstar, 
                             arma::vec nstar, const arma::vec& theta, 
                             const double& alpha) {
  int n = theta.n_elem;
  int N = 3001;
  int R = 1000;
  double eps = 1e-6;
  arma::vec s = -arma::log(arma::linspace(exp(-1e-6), exp(-5e-4), N));
  arma::vec z = arma::linspace(eps, 1-eps, R);
  
  // Compute psi(z; u, tht)
  arma::vec psi_z(R);
  for (int i = 0; i < R; ++i) {
    psi_z(i) = psi(z(i), u, theta);
  }
  
  // First Part: random locations and jumps
  arma::vec fnS(N);
  for (int i = 0; i < N; ++i) {
    fnS(i) = arma::mean(arma::exp(-(1 + psi_z) * s(i)) / s(i));
  }
  
  arma::vec ds = arma::diff(s);
  arma::vec h = (fnS.subvec(0, N-2) + fnS.subvec(1, N-1)) / 2;
  arma::vec Nv = arma::reverse(arma::cumsum(arma::reverse(ds % h)));
  Nv = arma::join_cols(Nv, arma::zeros(1));
  Nv *= alpha;
  
  arma::vec xi = arma::cumsum(as<arma::vec>(Rcpp::rexp(M, 1.0)));
  arma::vec RJ(M); // RJ = Random jumps, M = Finite truncation for CRM 
  int iNv = N - 1;
  
  for (int i = 0; i < M; ++i) {
    while (iNv > 0 && Nv(iNv) < xi(i)) {
      --iNv;
    }
    RJ(i) = s(iNv + 1);
  }

  arma::vec RL(M); // RL = Random locations
  for (int m = 0; m < M; ++m) {
    double xi = arma::randu();
    arma::vec temp = arma::exp(-(1 + psi_z) * RJ(m));
    double cutoff = arma::accu(temp) * xi;
    RL(m) = z(arma::as_scalar(arma::find(arma::cumsum(temp) - cutoff > 0, 1, "first")));
  }
  
  // Second Part: random jumps [for fixed locations]
  arma::vec psi_star(zstar.n_elem);
  for (int i = 0; i < zstar.n_elem; ++i) {
    psi_star(i) = psi_zstar(zstar(i), u, theta);
  }
  
  arma::vec Jstar(zstar.n_elem);
  for (int i = 0; i < zstar.n_elem; ++i) {
    Jstar(i) = R::rgamma(nstar(i), 1 / (psi_star(i) + 1));
  }
  
  return List::create(Named("RL") = RL,
                      Named("RJ") = RJ,
                      Named("zstar") = zstar, 
                      Named("Jstar") = Jstar
  );
}


// [[Rcpp::export]]
List crm_sampler_fractionalY(const int& M, const arma::vec& u, arma::vec zstar,
                             arma::vec nstar, arma::vec RL,
                             arma::vec RJ, const arma::vec& tht, const double& alpha,
                             const arma::vec& mu, const arma::vec& y,
                             const float& shape1, const float& shape2) {

  int n = tht.n_elem;
  int N = 3001;
  int R = 3001;
  double eps = 1e-6;
  arma::vec s = -arma::log(arma::linspace(exp(-1e-6), exp(-5e-4), N));

  // Sorted, ascending order, needed in RL
  // arma::vec z = arma::sort(as<arma::vec>(Rcpp::runif(R, 0, 1))); // Centering dist: Uniform(0, 1)
  // arma::vec z = arma::sort(as<arma::vec>(Rcpp::rbeta(R, shape1, shape2))); // Centering dist: Beta(shape1, shape2)
  arma::vec z = arma::linspace(eps, 1-eps, R);

  // Compute psi(z; u, tht)
  arma::vec psi_z(R);
  for (int i = 0; i < R; ++i) {
    psi_z(i) = psi(z(i), u, tht);
  }

  // First Part: random locations and jumps
  arma::vec fnS(N);
  for (int i = 0; i < N; ++i) {
    fnS(i) = arma::mean(arma::exp(-(1 + psi_z) * s(i)) / s(i));
  }

  arma::vec ds = arma::diff(s);
  arma::vec h = (fnS.subvec(0, N-2) + fnS.subvec(1, N-1)) / 2;
  arma::vec Nv = arma::reverse(arma::cumsum(arma::reverse(ds % h)));
  Nv = arma::join_cols(Nv, arma::zeros(1));
  Nv *= alpha;

  arma::vec xi = arma::cumsum(as<arma::vec>(Rcpp::rexp(M, 1.0)));
  // arma::vec RJ(M); // RJ = Random jumps, M = Finite truncation for CRM
  int iNv = N - 1;

  for (int i = 0; i < M; ++i) {
    while (iNv > 0 && Nv(iNv) < xi(i)) {
      --iNv;
    }
    RJ(i) = s(iNv + 1);
  }

  // arma::vec RL(M); // RL = Random locations
  for (int m = 0; m < M; ++m) {
    double xi = arma::randu();
    arma::vec temp = arma::exp(-(1 + psi_z) * RJ(m));
    double cutoff = arma::accu(temp) * xi;
    RL(m) = z(arma::as_scalar(arma::find(arma::cumsum(temp) - cutoff > 0, 1, "first")));
  }

  // Second Part: random jumps [for fixed locations]
  arma::vec psi_star(zstar.n_elem);
  for (int i = 0; i < zstar.n_elem; ++i) {
    psi_star(i) = psi_zstar(zstar(i), u, tht);
  }

  arma::vec Jstar(zstar.n_elem);
  for (int i = 0; i < zstar.n_elem; ++i) {
    Jstar(i) = R::rgamma(nstar(i), 1 / (psi_star(i) + 1)); // shape, scale
  }

  // arma::vec z_tld_temp = arma::join_cols(RL, zstar);
  // if (arma::min(mu) >= arma::min(z_tld_temp) && arma::max(mu) <= arma::max(z_tld_temp)) {
  //   RL = RL;
  //   RJ = RJ;
  //   zstar = zstar;
  //   Jstar = Jstar;
  // }

  return List::create(Named("RL") = RL,
                      Named("RJ") = RJ,
                      Named("zstar") = zstar,
                      Named("Jstar") = Jstar
  );
}




// --------------------------------------------------------------------------- //
//                                  u update                                   //
// --------------------------------------------------------------------------- // 

// [[Rcpp::export]]
double log_posterior_u(const arma::vec& u, const arma::vec& z, const arma::vec& theta,
                       const double& alpha) {
  int R = 250;
  double eps = 1e-6;
  arma::vec zstar = arma::join_cols(arma::linspace<arma::vec>(0+eps, 1-eps, R), z); // R + length(z) = R + n vector
  arma::mat exp_mat = arma::exp(arma::diagmat(theta) * arma::repmat(zstar, 1, theta.size()).t());  // n x (R + n) matrix
  arma::mat u_exp_mat = arma::diagmat(u) * exp_mat;  // n x (R + n) matrix
  arma::vec log_1_plus_u_exp_mat_sum = arma::log(1 + arma::sum(u_exp_mat, 0).t());  // (R + n) x 1 vector
  double neg_log_post = alpha * mean(log_1_plus_u_exp_mat_sum.head(R)) +
    mean(log_1_plus_u_exp_mat_sum.tail(z.size()));
  return - neg_log_post;
}

// // [[Rcpp::export]]
// arma::vec u_sampler(const arma::vec& u, const arma::vec& z, const arma::vec& theta,
//                     const double& alpha, const double& delta){
//   arma::vec u_star = arma::zeros<arma::vec>(u.size());
//   double logQ_ratio = 0.0;
//   for(int i = 0; i < u.size(); ++i){
//     u_star(i) = R::rgamma(delta, u(i) / delta); // u_proposal := u_star ~ Gamma(shape = delta, scale = u / delta)
//     logQ_ratio += logpdf_gamma(u(i), delta, u_star(i) / delta) -
//       logpdf_gamma(u_star(i), delta, u(i) / delta); // log(q(u | u_star)) - log(q(u_star | u))
//   }
//   double logratio = log_posterior_u(u_star, z, theta, alpha) -
//     log_posterior_u(u, z, theta, alpha) + logQ_ratio;
//   if(log(R::runif(0, 1)) < logratio){
//     return u_star;
//   } else {
//     return u;
//   }
// }

// [[Rcpp::export]]
arma::vec u_sampler(arma::vec& u, const arma::vec& z, const arma::vec& theta,
                    const double& alpha, const double& delta){
  for(int i = 0; i < u.size(); ++i){
    arma::vec u_star = u;
    u_star(i) = R::rgamma(delta, u(i) / delta); // u_proposal := u_star ~ Gamma(shape = delta, scale = u / delta)
    double logQ_ratio = logpdf_gamma(u(i), delta, u_star(i) / delta) -
      logpdf_gamma(u_star(i), delta, u(i) / delta); // log(q(u | u_star)) - log(q(u_star | u))
    double logratio = log_posterior_u(u_star, z, theta, alpha) -
    log_posterior_u(u, z, theta, alpha) + logQ_ratio;
  if(log(R::runif(0, 1)) < logratio){
    u(i) = u_star(i);
  } else {
      u(i) = u(i);
    }
  }
  return u;
}