#ifndef CLUSTER_STATS_H
#define CLUSTER_STATS_H

#include<algorithm>
#include<numeric>
#include<functional>
#include "types.h"
#include "blas.h"
#include "debug_print.h"

class is_equal_to {
public:
  is_equal_to(unsigned int n) : n(n) {}
  unsigned int operator()(unsigned int i) const { return (double)(i==n); }
private:
  const unsigned int n;
};

//separate fn for indices
void get_k_indices(const uvec &z, fvec &indic, int k){
  std::transform(z.begin(), z.end(), indic.begin(), is_equal_to(k)); // pick out genes in group k
}

void construct_precision_mat(const fvec &xTx, fvec &chol_S, double Gk, double lambda, int V, int n){
  std::transform(xTx.begin(), xTx.end(), chol_S.begin(), std::bind2nd(std::multiplies<double>(), Gk*n));
  for(int i=0; i<V; i++)
    chol_S[i*V+i] += lambda;
}

void increment_IGscale(fveci IGscale, const fvec &chol_Sinv, fvec beta_hat){
  //pass beta_hat by copy since dtrmv will modify argument
  //multiply beta_hat by cholesky factor
  //increment IGscale by product * 1/2
  int n = beta_hat.size();
  multiply_lowertri_vec(n, chol_Sinv, beta_hat);
  *IGscale = ddot(n, &(beta_hat[0]), &(beta_hat[0]));
}

double get_cluster_size(fvec &indices){
  double result = std::accumulate(indices.begin(), indices.end(), 0); // num. genes alloc. to k
  return result;
}

void cluster_stats(int k, fvec &yTx, const fvec &xTx, int G, int V, int n,
                   const uvec &z, fveci Gkk,
                   fvec &beta_hat, fvec &chol_S, fveci IGscale){
    fvec indic(G);
    get_k_indices(z, indic, k);
    
    *Gkk = get_cluster_size(indic);
    
//     Rprintf("yTx:\n");
//     print_fmat(yTx, V, G);
//     
//     Rprintf("indices:\n");
//     print_fmat(indic, 1, G);
    
    multiply_mat_vec(G, V, yTx, indic, beta_hat, 1); // sum yTx for g: z_g==k, last arg means trans = "T"
//     Rprintf("yTx:\n");
//     print_fmat(beta_hat, 1, V);
    
    construct_precision_mat(xTx, chol_S, *Gkk, 1.0, V, n); // S_inv = (X^T * X + lambda * I)
//     Rprintf("chol_prec:\n");
//     print_fmat(chol_S, V, V);
//       
    
    solve_normaleq_symm_mat(V, &(chol_S[0]), &(beta_hat[0])); // solve normal equation X^Ty = S_inv * beta_hat, convert S_inv to cholesky factor
//     Rprintf("betahat:\n");
//     print_fmat(beta_hat, 1, V);
    
    increment_IGscale(IGscale, chol_S, beta_hat);
    
    invert_lower_tri(V, &(chol_S[0])); // invert S_inv 
  
}

#endif