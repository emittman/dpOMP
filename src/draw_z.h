#ifndef DRAW_Z_H
#define DRAW_Z_H
#include "rcategorical.h"
/*
  * Compute the weights required for sampling Z.
  * Advice on subsetting from http://stackoverflow.com/questions/421573/best-way-to-extract-a-subvector-from-a-vector
  */

void draw_z(fvec &weights, uvec &z, int G, int K){
#pragma omp parallel for
  for(int g=0; g<G; g++){
    
    fveci weights_iter = weights.begin() + K*g;
    z[g] = rcategorical(weights_iter, K);
    
  }
}


#endif //DRAW_Z_H