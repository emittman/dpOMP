#include "blas.h"
#include "cluster_stats.h"

void draw_theta(uvec &z, fvec &yTx, fvec &xTx, fvec &Gk, int G, int K, int V, int n){
fvec IGscale(K);
#pragma omp parallel for //private k?
  for(int k = 0; k<K; k++){
    fvec beta_hat(V);
    fvec chol_S(V*V);
    cluster_stats(k, yTx, xTx, G, V, n, z, Gk.begin() + k, beta_hat, chol_S, IGscale.begin() + k);
  }
  
}
