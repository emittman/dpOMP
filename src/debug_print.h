#include<R.h>
#include "types.h"

void print_fmat(fvec &m, int M, int N){
  for(int i=0; i<M; i++){
    for(int j=0; j<N; j++)
      Rprintf("%.2lf\t",m[i*N + j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}