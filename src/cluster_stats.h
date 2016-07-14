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

void construct_precision_mat(const fvec &xTx, fvec &chol_S, double Gk, double sigma2, double lambda2, int V){
  
  std::transform(xTx.begin(), xTx.end(), chol_S.begin(), std::bind2nd(std::multiplies<double>(), Gk));
  for(int i=0; i<V; i++)
    chol_S[i*V+i] += sigma2/lambda2;
}

void calculate_IGscale(double *IGscale, double yTyk, fvec &xTyk, fvec &xTx, fvec &beta, double Gk, int V){
  
  double out = 0.0;
  if(Gk != 0) {
    out = yTyk;
    out = out - 2.0 * inner_prod_vec(V, &(xTyk[0]), &(beta[0]));
    out = out + Gk * quad_form(V, &(beta[0]), &(xTx[0]));
  }
  
  *IGscale = out;
}

double get_cluster_size(fvec &indices){
  
  double result = std::accumulate(indices.begin(), indices.end(), 0); // num. genes alloc. to k
  return result;
}

void cluster_sums(int k, fvec &yTy, fvec &xTy, int G, int V, const uvec &z, fveci Gkk, double *yTyk, fvec &xTyk){
  
    fvec indic(G);
    get_k_indices(z, indic, k);
    *Gkk = get_cluster_size(indic);
    *yTyk = inner_prod_vec(G, &(indic[0]), &(yTy[0]));
    // sum xTy for g: z_g==k, last arg means trans = "F"    
    multiply_mat_vec(V, G, xTy, indic, xTyk, 0); 
    //     Rprintf("xTyk:\n");
    //     print_mat(beta_hat, 1, V);
}

#endif