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

void construct_precision_mat(const fvec &xTx, fvec &chol_S, double Gk, double sigma2, double lambda2, int V, int n){
  std::transform(xTx.begin(), xTx.end(), chol_S.begin(), std::bind2nd(std::multiplies<double>(), Gk));
  for(int i=0; i<V; i++)
    chol_S[i*V+i] += sigma2/lambda2;
}

void increment_IGscale(fveci IGscale, const fvec &chol_Sinv, fvec beta_hat){
  //pass beta_hat by copy since dtrmv will modify argument
  //multiply beta_hat by cholesky factor
  //increment IGscale by product * 1/2
  int n = beta_hat.size();
  multiply_lowertri_vec(n, chol_Sinv, beta_hat);
  *IGscale = inner_prod_vec(n, &(beta_hat[0]), &(beta_hat[0]));
}

double get_cluster_size(fvec &indices){
  double result = std::accumulate(indices.begin(), indices.end(), 0); // num. genes alloc. to k
  return result;
}

void cluster_sums(int k, fvec &xTy, int G, int V, int n, const uvec &z, fveci Gkk, fvec &xTyk){
    fvec indic(G);
    get_k_indices(z, indic, k);
    *Gkk = get_cluster_size(indic);
    // sum xTy for g: z_g==k, last arg means trans = "F"    
    multiply_mat_vec(V, G, xTy, indic, xTyk, 0); 
    //     Rprintf("xTyk:\n");
    //     print_mat(beta_hat, 1, V);
}

#endif