#ifndef DRAW_Z_H
#define DRAW_Z_H
#include "rcategorical.h"
#include "chain.h"
/*
  * Compute the weights required for sampling Z.
  * Advice on subsetting from http://stackoverflow.com/questions/421573/best-way-to-extract-a-subvector-from-a-vector
  */
void draw_z(chain_t &chain){
#pragma omp parallel for
  for(int g=0; g<chain.G; g++){
    
    fveci weights_iter = chain.weights.begin() + chain.K*g;
    chain.z[g] = rcategorical(weights_iter, chain.K);
    
  }
}


#endif //DRAW_Z_H