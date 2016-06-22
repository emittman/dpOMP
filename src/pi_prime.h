#ifndef PI_PRIME_H
#define PI_PRIME_H

#include <Rmath.h>
#include <numeric>
#include "types.h"
#include "blas.h"

double pi_prime(double *xTy_g_ptr, fvec &xTx, double *beta_k_ptr, double pi, double sigma2, int V){

  double out = log(pi)  - 1.0 / (2.0 * sigma2) * (quad_form(V, beta_k_ptr, &(xTx[0])) -
                   2.0 * inner_prod_vec(V, beta_k_ptr, xTy_g_ptr));

  return out;
}

#endif // PI_PRIME_H
