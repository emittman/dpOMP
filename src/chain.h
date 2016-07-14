#ifndef CHAIN_H
#define CHAIN_H

#include <Rmath.h>
#include "types.h"
#include "debug_print.h"

typedef struct {

  fvec yTy;
  fvec xTy;
  fvec xTx;
  fvec beta;
  fvec pi;
  fvec sigma2;
  double lambda2;
  double alpha;
  double a;
  double b;
  int G;
  int K;
  int V;
  int N;
  uvec z;
  fvec weights;
  fvec Gk;

} chain_t;

chain_t construct_chain(double *yTy_p, double *xTy_p, double *xTx_p, double lambda2, double alpha, double a, double b, int GG, int KK, int VV, int NN){
  
  chain_t chain;
  chain.yTy = fvec(yTy_p, yTy_p + GG);
  chain.xTy = fvec(xTy_p, xTy_p + GG*VV);
  chain.xTx = fvec(xTx_p, xTx_p + VV*VV);
  chain.beta = fvec(VV*KK);
  chain.pi = fvec(KK);
  chain.sigma2 = fvec(KK);
  chain.lambda2 = lambda2;
  chain.alpha = alpha;
  chain.a = a;
  chain.b = b;
  chain.G = GG;
  chain.K = KK;
  chain.V = VV;
  chain.N = NN;
  chain.z = uvec(GG);
  chain.weights = fvec(GG*KK);
  chain.Gk = fvec(KK);
  
  return chain;
}

void initialize_chain(chain_t &chain){
  // Fill with default/"agnostic" starting values
  for(int i=0; i<chain.K; i++){
    chain.pi[i] = 1.0/(double)chain.K;
    chain.sigma2[i] = 1.0;
  }
  
  for(int i=0; i<chain.K*chain.V; i++){
    chain.beta[i] = rnorm(0,chain.lambda2);
  }
  
}

void print_chain_state(chain_t &chain){
  
  Rprintf("z:\n");
  print_mat(chain.z, 1, chain.G);
  Rprintf("Gk:\n");
  print_mat(chain.Gk, 1, chain.K);
  Rprintf("beta:\n");
  print_mat(chain.beta, chain.V, chain.K);
  Rprintf("pi:\n");
  print_mat(chain.pi, 1, chain.K);
  Rprintf("sigma2:\n");
  print_mat(chain.sigma2, 1, chain.K);
  Rprintf("lambda2:\n %lf \n", chain.lambda2);
  Rprintf("alpha:\n %lf \n", chain.alpha);
  Rprintf("a:\n %lf \n", chain.a);
  Rprintf("b:\n %lf \n", chain.b);
}

#endif // CHAIN_H