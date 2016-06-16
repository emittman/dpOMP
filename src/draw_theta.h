#ifndef DRAW_THETA_H
#define DRAW_THETA_H

#include <Rmath.h>
#include "blas.h"
#include "cluster_stats.h"
#include "sample_wrappers.h"

void draw_theta(fvec &beta, double *sigma2, uvec &z, double *yTy, fvec &xTy, fvec &xTx, fvec &Gk, int G, int K, int V, int n){
  fvec IGscale(K);
  //private k?  
  // #pragma omp parallel for
  for(int k = 0; k<K; k++){
    fvec beta_hat(V);
    fvec chol_S(V*V);
    fvec beta_raw(V);
    cluster_stats(k, xTy, xTx, G, V, n, z, Gk.begin() + k, beta_hat, chol_S, IGscale.begin() + k);
    for(int v=0; v<V; v++)
      beta_raw[v] = rnorm(0, 1); //draw beta_raw
    
//     Rprintf("\nbeta raw:\n");
//     print_mat(beta_raw, 1, V);
    
    multiply_lowertri_vec(V, chol_S, beta_raw);
//     Rprintf("\nbeta raw (scaled):\n");
//     print_mat(beta_raw, 1, V);
   
    linear_comb_vec(V, sqrt(*sigma2), beta_raw, beta_hat);
//     Rprintf("\nbeta draw:\n");
//     print_mat(beta_hat, 1, V);
    
    
    for(int v=0; v<V; v++)
      beta[k*V + v] = beta_hat[v];
  }
  double scale = std::accumulate(IGscale.begin(), IGscale.end(), 0.0);
  *sigma2 = rinvgamma(G*V*n/2, (*yTy + scale)/2);
  Rprintf("IG shape = %d\n", G*V*n/2);
  Rprintf("IG scale = %lf\n",(*yTy +scale)/2);
  Rprintf("sigma2 = %lf\n", *sigma2);
}

#endif