#ifndef DEBUG_PRINT_H
#define DEBUG_PRINT_H

#include<R.h>
#include "types.h"

void print_mat(fvec &m, int M, int N){
  for(int i=0; i<M; i++){
    for(int j=0; j<N; j++)
      Rprintf("%.2lf\t",m[i*N + j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

void print_mat(uvec &m, int M, int N){
  for(int i=0; i<M; i++){
    for(int j=0; j<N; j++)
      Rprintf("%d\t",m[i*N + j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

void print_mat(const fvec &m, int M, int N){
  for(int i=0; i<M; i++){
    for(int j=0; j<N; j++)
      Rprintf("%.2lf\t",m[i*N + j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

#endif