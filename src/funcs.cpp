#include <vector>
#include <iostream>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

typedef std::vector<float> fvec;
typedef std::vector<int> ivec;


extern "C" SEXP n_squaresR(SEXP n){

  int nn = INTEGER(n)[0];

  ivec vec(nn);
  
  #pragma omp parallel for
  for(int i=0; i<nn; i++){
    vec[i] = i * i;
    //std::cout <<"tid: " << i <<"\n";
  }
  
  SEXP result = PROTECT(allocVector(INTSXP, nn));
  for(int i=0; i<nn; i++)
    INTEGER(result)[i] = vec[i];
  
  UNPROTECT(1);
  
  return result;
}
