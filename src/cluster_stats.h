#ifndef CLUSTER_STATS_H
#define CLUSTER_STATS_H

#include "types.h"
#include<algorithm>
#include<numeric>

class is_equal_to {
public:
  is_equal_to(unsigned int n) : n(n) {}
  unsigned int operator()(unsigned int i) const { return i==n; }
private:
  const unsigned int n;
};

void cluster_stats(uvec &z, ivec& Gk, int G, int K){
  uvec indic(G);
  for(int k = 0; k<K; k++){
    transform(z.begin(), z.end(), indic.begin(), 
              is_equal_to(k));
    Gk[k] = std::accumulate(indic.begin(), indic.end(), 0);
  }
  
}

#endif