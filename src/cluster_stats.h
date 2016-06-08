#ifndef CLUSTER_STATS_H
#define CLUSTER_STATS_H

#include<algorithm>
#include<numeric>
#include "types.h"
#include "blas.h"

class is_equal_to {
public:
  is_equal_to(unsigned int n) : n(n) {}
  unsigned int operator()(unsigned int i) const { return (double)(i==n); }
private:
  const unsigned int n;
};

//separate fn for indices
void get_k_indices(uvec &z, fvec &indic, int k){
  transform(z.begin(), z.end(), indic.begin(), is_equal_to(k)); // pick out genes in group k
}


void cluster_stats(uvec &z, fvec &yTx, fvec &Gk, int G, int K, int V){
#pragma omp parallel for //private k?
  for(int k = 0; k<K; k++){
    fvec indic(G);
    get_k_indices(z, indic, k);
    
    Gk[k] = std::accumulate(indic.begin(), indic.end(), 0); // num. genes alloc. to k
    
    fvec sum_yTx(V);
    dgemv(V, G, 1, &(yTx[0]), &(indic[0]), 1, &(sum_yTx[0])); // sum yTx for g: z_g==k
//     for(int i=0; i<V; i++)
//       printf("sum_yTx[%d,%d]=%lf\t", k, i, sum_yTx[i]);
//     printf("\n\n");
  

  }
  
}

#endif