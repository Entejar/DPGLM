#include <RcppArmadillo.h>
#include "functions.h"

 
// [[Rcpp::export]]
arma::vec z_sampler(const arma::vec& y, const double& c0, const arma::vec& theta, 
                    const arma::mat& crm_mat) {
  arma::vec z = arma::zeros<arma::vec>(y.size());
  arma::uvec sort_indices = arma::sort_index(crm_mat.row(0).t());
  
  arma::mat crm_sorted = crm_mat.cols(sort_indices); // Sort by locations, Todo: need to implement the 
                                                     // clever starting point for finding indices. For  
                                                     // that, data := (X,y) to be arranged by ascending y
  
  arma::vec locations = crm_sorted.row(0).t();
  arma::vec jumps = crm_sorted.row(1).t();
  for(int i = 0; i < y.size(); ++i) {
    arma::uvec indices = arma::find(locations - c0 < y(i) && y(i) < locations + c0);
    if (indices.size() == 0) { // if no indices are found
      z(i) = y(i);
      continue; // skip the rest of the loop
    }
    arma::vec logprob = arma::zeros<arma::vec>(indices.size());
    for(int j = 0; j < indices.size(); ++j) {
      logprob(j) = logpdf_unif(y(i), locations(indices(j)) - c0, locations(indices(j)) + c0) +
        theta(i) * locations(indices(j)) + log(jumps(indices(j)));
    }
    logprob -= arma::max(logprob);
    arma::vec prob = arma::exp(logprob);
    prob /= arma::sum(prob);
    int k = rcategorical(prob);
    z(i) = locations(indices(k));
  }
  return z;
}