#ifndef U_H
#define U_H

#include <RcppArmadillo.h>
#include "functions.h"
#include "helpers.h"

arma::vec u_sampler(const arma::vec& u_old, const arma::vec& theta, const arma::vec& z,
                    const double& alpha, const double& delta);

#endif