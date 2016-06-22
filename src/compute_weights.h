#ifndef COMPUTE_WEIGHTS_H
#define COMPUTE_WEIGHTS_H

#include "pi_prime.h"
/*
 * Compute the weights required for sampling Z.
 * Advice on subsetting from http://stackoverflow.com/questions/421573/best-way-to-extract-a-subvector-from-a-vector
 */

void compute_weights(fvec &weights, fvec &xTy, fvec &xTx, fvec &beta, fvec &pi, double sigma2, int G, int K, int V){
#pragma omp parallel for
  for(int g=0; g<G; g++){
    
    double *xTy_g_ptr = &(xTy[V*g]);

    for(int k=0; k<K; k++) {

      double *beta_k_ptr = &(beta[V*k]);
      double pi_k = pi[k];
      weights[g*K + k] = pi_prime(xTy_g_ptr, xTx, beta_k_ptr, pi_k, sigma2, V);

    }
  }
}


#endif //COMPUTE_WEIGHTS_H