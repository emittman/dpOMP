#ifndef PI_PRIME_H
#define PI_PRIME_H

#include "blas.h"

double pi_prime(fveci yTx_g_iter, fvec &xTx, fveci beta_k_iter, double pi, double sigma2, int V){

  double out = log(pi)  - 1.0 / (2.0 * sigma2) * (quad_form(V, &(*beta_k_iter), &xTx[0]) -
                   2.0 * std::inner_product(yTx_g_iter, yTx_g_iter + V, beta_k_iter, 0.0));

  return out;
}

#endif // PI_PRIME_H
