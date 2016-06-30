#ifndef DRAW_THETA_H
#define DRAW_THETA_H

#include <Rmath.h>
#include "blas.h"
#include "cluster_stats.h"
#include "sample_wrappers.h"
#include "chain.h"

void draw_theta(chain_t &chain){
  int G = chain.G;
  int K = chain.K;
  int V = chain.V;
  int N = chain.N;
  fvec IGscale(K);
  
  //private k?  
  #pragma omp parallel for
  for(int k = 0; k<K; k++){
    
    fveci Gkk = chain.Gk.begin() + k;
    fvec xTyk(V);
    fvec beta_hat(V);
    fvec beta_raw(V);
    fvec chol_S(V*V);
    fveci IGscalek = IGscale.begin() + k;
    
    //compute Gkk and xTyk
    cluster_sums(k, chain.xTy, G, V, N, chain.z, Gkk, xTyk);
    
    //     Rprintf("xTy:\n");
    //     print_mat(xTy, G, V);
    
    //     Rprintf("indices:\n");
    //     print_mat(indic, 1, G);
    
    // S_inv = (X^T * X + sigma2/lambda2 * I)
    // at this point, lambda2 == 1.0
    construct_precision_mat(chain.xTx, chol_S, *Gkk, chain.sigma2, 1.0, V, N);
    //     Rprintf("prec:\n");
    //     print_mat(chol_S, V, V);
    //
    std::copy(xTyk.begin(), xTyk.end(), beta_hat.begin());
    
    // solve normal equation X^Ty = S_inv * beta_hat, convert S_inv to cholesky factor
    solve_normaleq_symm_mat(V, &(chol_S[0]), &(beta_hat[0])); 
    //     Rprintf("betahat:\n");
    //     print_mat(beta_hat, 1, V);
    
    // invert S_inv so that chol_S now contains cholesky factor of (XTX + sigma2/lambda2 I)^{-1}
    invert_lower_tri(V, &(chol_S[0])); 
    // Rprintf("var:\n");
    // print_mat(chol_S, V, V);
    for(int v=0; v<V; v++)
      beta_raw[v] = rnorm(0, 1);
    
//     Rprintf("\nbeta raw:\n");
//     print_mat(beta_raw, 1, V);

//scale beta_raw by S_lower -> beta_raw
    multiply_lowertri_vec(V, chol_S, beta_raw); 
//     Rprintf("\nbeta raw (scaled):\n");
//     print_mat(beta_raw, 1, V);

// translate (scaled)beta_raw by beta_hat -> beta_hat
    linear_comb_vec(V, sqrt(chain.sigma2), beta_raw, beta_hat); 
//     Rprintf("\nbeta draw:\n");
//     print_mat(beta_hat, 1, V);

    for(int v=0; v<V; v++)
      chain.beta[k*V + v] = beta_hat[v];
  
  // compute contribution from cluster K to IGscale
  increment_IGscale(IGscalek, xTyk, beta_hat, *Gkk, V);
}
  //     Rprintf("IG:\n");
  //     Rprintf("%lf\n",*IGscale);
//   Rprintf("compnents of IGscale:\n");
//   print_mat(IGscale, 1, K);
  double scale = std::accumulate(IGscale.begin(), IGscale.end(), 0.0);
  chain.sigma2 = rinvgamma(chain.G*V*chain.N/2, (chain.yTy - scale)/2);
//   Rprintf("IG shape = %d\n", G*V*n/2);
//   Rprintf("IG scale = %lf\n",(*yTy - scale)/2);
}

#endif