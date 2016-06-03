#ifndef COMPUTE_WEIGHTS_H
#define COMPUTE_WEIGHTS_H

/*
 * Compute the weights required for sampling Z.
 * Advice on subsetting from http://stackoverflow.com/questions/421573/best-way-to-extract-a-subvector-from-a-vector
 */

void compute_weights(fvec &weights, fvec &yTx, fvec &xTx, fvec &beta, fvec &pi, int G, int K, int n, int V){
#pragma omp parallel for
  for(int g=0; g<G; g++){

    fveci yTx_g_iter = yTx.begin() + V*g;

    for(int k=0; k<K; k++) {
      fveci beta_k_iter = beta.begin() + V*k;
      double pi_k = pi[k];
      weights[g*K + k] = pi_prime(yTx_g_iter, xTx, beta_k_iter, pi_k, 1.0, V);
    }
  }
}


#endif //COMPUTE_WEIGHTS_H