#ifndef PI_PRIME_H
#define PI_PRIME_H

#include <Rmath.h>
#include <numeric>
#include "types.h"
#include "blas.h"

double pi_prime(fveci xTy_g_iter, fvec &xTx, fveci beta_k_iter, double pi, double sigma2, int V){

  double out = log(pi)  - 1.0 / (2.0 * sigma2) * (quad_form(V, &(*beta_k_iter), &xTx[0]) -
                   2.0 * std::inner_product(xTy_g_iter, xTy_g_iter + V, beta_k_iter, 0.0));

  return out;
}

#endif // PI_PRIME_H
