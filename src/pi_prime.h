#ifndef PI_PRIME_H
#define PI_PRIME_H

#include "blas.h"

double pi_prime(fveci yTxg_iter, fvec &xTx, fvec &beta, double pi, double sigma2, int V){

  double *xTxp, *betap;
  xTxp = &xTx[0];
  betap = &beta[0];
  double out = log(pi)  - 1.0 / (2.0 * sigma2) * (quad_form(V, betap, xTxp) -
                   2.0 * std::inner_product(yTxg_iter, yTxg_iter + V, beta.begin(), 0.0));

  return out;
}

#endif // PI_PRIME_H
