#ifndef DRAW_THETA_H
#define DRAW_THETA_H

#include <Rmath.h>
#include "blas.h"
#include "cluster_stats.h"
#include "sample_wrappers.h"
#include "chain.h"

void draw_theta(chain_t &chain){
  fvec IGscale(chain.K);
  int V = chain.V;
//   Rprintf("just inside draw_theta, xTx:\n");
//   print_mat(xTx, V, V);
  //private k?  
  #pragma omp parallel for
  for(int k = 0; k<chain.K; k++){
    fvec beta_hat(V);
    fvec chol_S(V*V);
    fvec beta_raw(V);
    cluster_stats(k, chain.xTy, chain.xTx, chain.G, chain.V, chain.N, chain.z, chain.Gk.begin() + k, beta_hat, chol_S, IGscale.begin() + k);
    for(int v=0; v<V; v++)
      beta_raw[v] = rnorm(0, 1); //draw beta_raw ~ind. N(0,1)
    
//     Rprintf("\beta hat: \n");
//     print_mat(beta_hat, 1, V);
    
//     Rprintf("\nbeta raw:\n");
//     print_mat(beta_raw, 1, V);
    
    multiply_lowertri_vec(V, chol_S, beta_raw); //scale beta_raw by S_lower
//     Rprintf("\nbeta raw (scaled):\n");
//     print_mat(beta_raw, 1, V);
   
    linear_comb_vec(V, sqrt(chain.sigma2), beta_raw, beta_hat); // translate beta_raw by beta_hat
//     Rprintf("\nbeta draw:\n");
//     print_mat(beta_hat, 1, V);
    
    
    for(int v=0; v<V; v++)
      chain.beta[k*V + v] = beta_hat[v];
  }
//   Rprintf("compnents of IGscale:\n");
//   print_mat(IGscale, 1, K);
  double scale = std::accumulate(IGscale.begin(), IGscale.end(), 0.0);
  chain.sigma2 = rinvgamma(chain.G*V*chain.N/2, (chain.yTy - scale)/2);
//   Rprintf("IG shape = %d\n", G*V*n/2);
//   Rprintf("IG scale = %lf\n",(*yTy - scale)/2);
}

#endif