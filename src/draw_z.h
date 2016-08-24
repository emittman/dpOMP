#ifndef DRAW_Z_H
#define DRAW_Z_H
#include "rcategorical.h"
#include "chain.h"
#include "types.h"
/*
  * Compute the weights required for sampling Z.
  * Advice on subsetting from http://stackoverflow.com/questions/421573/best-way-to-extract-a-subvector-from-a-vector
  */
void draw_z(chain_t &chain){
  int K = chain.K;
  // #pragma omp parallel for
  for(int g=0; g<chain.G; g++){
    
    fveci weights_iter = chain.weights.begin() + K*g;
    chain.z[g] = rcategorical(weights_iter, K);
    
  }
}


#endif //DRAW_Z_H