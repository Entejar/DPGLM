#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// --------------------------------------------------------------------------- //
//                          For beta update                                   //
// --------------------------------------------------------------------------- //                          

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


// --------------------------------------------------------------------------- //
//                               CRM update                                   //
// --------------------------------------------------------------------------- // 

//[[Rcpp::export]]
double psi(double z, arma::vec u, arma::vec tht) {
  double mx = 0; // arma::max(tht * z);
  return arma::accu(u % arma::exp(tht * z - mx));
}

//[[Rcpp::export]]
double psi_zstar(double zstar, arma::vec u, arma::vec tht) {
  double mx = 0; // arma::max(tht * zstar); 
  return arma::accu(u % arma::exp(tht * zstar - mx));
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

  return List::create(Named("RL") = RL,
                      Named("RJ") = RJ,
                      Named("zstar") = zstar,
                      Named("Jstar") = Jstar
  );
}


