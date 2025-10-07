#include <RcppArmadillo.h>
#include "crm.h"

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
List crm_sampler(const int& M, const arma::vec& u, const arma::vec& zstar, 
                 const arma::vec& nstar, const arma::vec& tht, const double& alpha, 
                 const arma::vec& mu, const arma::vec& y, arma::vec z_tld, 
                 arma::vec J_tld) {
  
  int n = tht.n_elem;
  int N = 3001;
  int R = 25000;
  arma::vec s = -arma::log(arma::linspace(exp(-1e-6), exp(-5e-4), N));
  arma::vec z = arma::sort(as<arma::vec>(Rcpp::runif(R, 0, 1))); // From centering distribution
                                                                 // Sorted, ascending order, needed in RL
  
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
    psi_star(i) = psi_zstar(zstar(i), u, tht);
  }
  
  arma::vec RJ_fixedL(zstar.n_elem);
  for (int i = 0; i < zstar.n_elem; ++i) {
    RJ_fixedL(i) = R::rgamma(nstar(i), 1 / (psi_star(i) + 1));
  }
  
  arma::vec z_tld_temp = arma::join_cols(RL, zstar);
  if (arma::min(mu) >= arma::min(z_tld_temp) && arma::max(mu) <= arma::max(z_tld_temp)) {
    z_tld = z_tld_temp;
    J_tld = arma::join_cols(RJ, RJ_fixedL);
  }
  
  bool temp1 = arma::min(y) == arma::min(z_tld);
  bool temp2 = arma::max(y) == arma::max(z_tld);
  
  if (!temp1) {
    z_tld = arma::join_cols(arma::vec({arma::min(y)}), z_tld);
    J_tld = arma::join_cols(arma::vec({arma::min(J_tld)}), J_tld);
  }
  if (!temp2) {
    z_tld = arma::join_cols(z_tld, arma::vec({arma::max(y)}));
    J_tld = arma::join_cols(J_tld, arma::vec({arma::min(J_tld)}));
  }
  
  return List::create(Named("z_tld") = z_tld, Named("J_tld") = J_tld);
}