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

  //private k?  
  #pragma omp parallel for
  for(int k = 0; k<K; k++){
    // Rprintf("cluster %d:",k);
    fveci Gkk = chain.Gk.begin() + k;
    double yTyk;
    fvec xTyk(V);
    fvec beta_hat(V);
    fvec beta_raw(V);
    fvec chol_S(V*V);
    double IGscalek = 0.0;
    double sigma2k = chain.sigma2[k];
    
    //compute Gkk and xTyk
    cluster_sums(k, chain.yTy, chain.xTy, G, V, chain.z, Gkk, &yTyk, xTyk);

    if(*Gkk > 0){
      // S_inv = (X^T * X + sigma2/lambda2 * I)
      construct_precision_mat(chain.xTx, chol_S, *Gkk, sigma2k, chain.lambda2, V);
      std::copy(xTyk.begin(), xTyk.end(), beta_hat.begin());
    
      // solve normal equation X^Ty = S_inv * beta_hat, convert S_inv to cholesky factor
      solve_normaleq_symm_mat(V, &(chol_S[0]), &(beta_hat[0])); 

      // invert S_inv so that chol_S now contains cholesky factor of (XTX + sigma2/lambda2 I)^{-1}
      invert_lower_tri(V, &(chol_S[0])); 

      for(int v=0; v<V; v++)
        beta_raw[v] = rnorm(0, 1);

      //scale beta_raw by S_lower -> beta_raw
      multiply_lowertri_vec(V, chol_S, beta_raw); 

      // translate (scaled)beta_raw by beta_hat -> beta_hat
      linear_comb_vec(V, sqrt(sigma2k), beta_raw, beta_hat); 

      for(int v=0; v<V; v++)
        chain.beta[k*V + v] = beta_hat[v];
  
      // compute IGscale
      calculate_IGscale(&IGscalek, yTyk, xTyk, chain.xTx, beta_hat, *Gkk, V);
      chain.sigma2[k] = rinvgamma(*Gkk*V*N/2 + chain.a, IGscalek/2 + chain.b);

    } else {
      for(int v=0; v<V; v++){
        chain.beta[k*V + v] = rnorm(0, 1) * sqrt(chain.lambda2);
      }
      
      chain.sigma2[k] = rinvgamma(chain.a, chain.b);
    }
     
  }
//   double scale = std::accumulate(IGscale.begin(), IGscale.end(), 0.0);
//   chain.sigma2 = rinvgamma(chain.G*V*chain.N/2, (chain.yTy - scale)/2);
//   Rprintf("IG shape = %d\n", G*V*n/2);
//   Rprintf("IG scale = %lf\n",(*yTy - scale)/2);
}

#endif