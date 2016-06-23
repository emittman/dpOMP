#ifndef COMPUTE_WEIGHTS_H
#define COMPUTE_WEIGHTS_H

#include "pi_prime.h"
#include "chain.h"
/*
 * Compute the weights required for sampling Z.
 * Advice on subsetting from http://stackoverflow.com/questions/421573/best-way-to-extract-a-subvector-from-a-vector
 */
void compute_weights(chain_t &chain){
#pragma omp parallel for
  for(int g=0; g<chain.G; g++){
    
    double *xTy_g_ptr = &(chain.xTy[chain.V*g]);

    for(int k=0; k<chain.K; k++) {

      double *beta_k_ptr = &(chain.beta[chain.V*k]);
      double pi_k = chain.pi[k];
      chain.weights[g*chain.K + k] = pi_prime(xTy_g_ptr, chain.xTx, beta_k_ptr, pi_k, chain.sigma2, chain.V);

    }
  }
}


#endif //COMPUTE_WEIGHTS_H