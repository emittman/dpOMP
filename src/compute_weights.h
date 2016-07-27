#ifndef COMPUTE_WEIGHTS_H
#define COMPUTE_WEIGHTS_H

#include "pi_prime.h"
#include "chain.h"
/*
 * Compute the weights required for sampling Z.
 * Advice on subsetting from http://stackoverflow.com/questions/421573/best-way-to-extract-a-subvector-from-a-vector
 */
void compute_weights(chain_t &chain){
  int G = chain.G;
  int K = chain.K;
  int V = chain.V;
//#pragma omp parallel for
  for(int g=0; g<G; g++){
    
    double yTy_g = chain.yTy[g];
    double *xTy_g_ptr = &(chain.xTy[V*g]);

    for(int k=0; k<K; k++) {

      double *beta_k_ptr = &(chain.beta[V*k]);
      double pi_k = chain.pi[k];
      double sigma2_k = chain.sigma2[k];
      chain.weights[g*K + k] = pi_prime(yTy_g, xTy_g_ptr, chain.xTx, beta_k_ptr, pi_k, sigma2_k, V);

    }
  }
}


#endif //COMPUTE_WEIGHTS_H