#include "includes.h"

#define Ralloc_Real(len)  PROTECT(allocVector(REALSXP, (len)))
#define Ralloc_List(len)  PROTECT(allocVector(VECSXP, (len)))

extern "C" SEXP dpgmmR(SEXP data, SEXP design, SEXP G, SEXP K, SEXP N, SEXP V){
  int GG, KK, NN, VV, data_len, design_len, beta_len;
  GG = INTEGER(G)[0];
  KK = INTEGER(K)[0];
  NN = INTEGER(N)[0];
  VV = INTEGER(V)[0];
  
  data_len = GG*NN*VV;
  design_len = VV*VV;
  beta_len = KK*VV;

  double *data_p, *design_p;
  data_p = NUMERIC_POINTER(data);
  design_p = NUMERIC_POINTER(design);
  
  fvec data_in(data_p, data_p + data_len);
  fvec design_in(design_p, design_p + design_len);
  
  fvec beta_h(KK*VV);
  fvec pi_h(KK);
  
  SEXP beta_out = Ralloc_Real(beta_len);
  SEXP pi_out = Ralloc_Real(KK);
  SEXP list_out = Ralloc_List(2);
  
  for(int i=0; i<beta_len; i++)
    REAL(beta_out)[i] = beta_h[i];
  
  for(int i=0; i<KK; i++)
    REAL(pi_out)[i] = pi_h[i];
  
  SET_VECTOR_ELT(list_out, 0, beta_out);
  SET_VECTOR_ELT(list_out, 1, pi_out);
  
  UNPROTECT(3);
  
  return list_out;
}

extern "C" SEXP pi_primeR(SEXP yTx, SEXP xTx, SEXP beta, SEXP pi, SEXP sigma2, SEXP V){
  double *yTx_p, *xTx_p, *beta_p, pi_c, sigma2_c;
  int v = INTEGER(V)[0];
  yTx_p  = NUMERIC_POINTER(yTx);
  xTx_p  = NUMERIC_POINTER(xTx);
  beta_p = NUMERIC_POINTER(beta);
  fvec yTx_c(yTx_p,  yTx_p  + v);
  fvec xTx_c(xTx_p,  xTx_p  + v*v);
  fvec beta_c(beta_p, beta_p + v);
  pi_c = REAL(pi)[0];
  sigma2_c = REAL(sigma2)[0];
  
  SEXP out = Ralloc_Real(1);
  REAL(out)[0] = pi_prime(yTx_c.begin(), xTx_c, beta_c.begin(), pi_c, sigma2_c, v);
  
  UNPROTECT(1);
  return out;
}

extern "C" SEXP compute_weightsR(SEXP yTx, SEXP xTx, SEXP beta, SEXP pi, SEXP G, SEXP K, SEXP n, SEXP V){
  double *yTx_p, *xTx_p, *beta_p, *pi_p;
  int GG, KK, nn, VV;
  yTx_p  = NUMERIC_POINTER(yTx);
  xTx_p  = NUMERIC_POINTER(xTx);
  beta_p = NUMERIC_POINTER(beta);
  pi_p   = NUMERIC_POINTER(pi);
  GG = INTEGER(G)[0];
  KK = INTEGER(K)[0];
  nn = INTEGER(n)[0];
  VV = INTEGER(V)[0];
  fvec arg0(GG*KK); //WEIGHTS
  fvec arg1(yTx_p,  yTx_p  + GG*VV);
  fvec arg2(xTx_p,  xTx_p  + VV*VV);
  fvec arg3(beta_p, beta_p + KK*VV);
  fvec arg4(pi_p, pi_p + KK);
  
  compute_weights(arg0, arg1, arg2, arg3, arg4, GG, KK, nn, VV);
  
  SEXP out = Ralloc_Real(GG*KK);
  for(int i=0; i<GG*KK; i++)
    REAL(out)[i] = arg0[i];
  UNPROTECT(1);
  return out;
}

extern "C" SEXP quad_formR(SEXP nin, SEXP xin, SEXP Ain){
  int n = INTEGER(nin)[0];
  double *A, *x;
  A = NUMERIC_POINTER(Ain);
  x = NUMERIC_POINTER(xin);
  
  SEXP result = PROTECT(allocVector(REALSXP, 1));
  REAL(result)[0] = quad_form(n, x, A);
  UNPROTECT(1);
  return result;
}

extern "C" SEXP norm_exp_lpR(SEXP weights, SEXP K){
  int KK = INTEGER(K)[0];
  double *weight_ptr = NUMERIC_POINTER(weights);
  fvec weights_c(weight_ptr, weight_ptr + KK);
  norm_exp_lp(weights_c.begin(), KK);
  SEXP result = PROTECT(allocVector(REALSXP, KK));
  
  for(int i=0; i<KK; i++){
    REAL(result)[i] = weights_c[i];
  }
  UNPROTECT(1);
  return result;
}

extern "C" SEXP rcategoricalR(SEXP weights, SEXP K){
  int KK = INTEGER(K)[0];
  double *weight_ptr = NUMERIC_POINTER(weights);
  fvec weights_c(weight_ptr, weight_ptr + KK);
  GetRNGstate();
  unsigned int z = rcategorical(weights_c.begin(), KK);
  PutRNGstate();
  SEXP result = PROTECT(allocVector(INTSXP, 1));
  
  INTEGER(result)[0] = z;
  
  UNPROTECT(1);
  return result;
}

extern "C" SEXP draw_zR(SEXP weights, SEXP G, SEXP K){
  int GG = INTEGER(G)[0];
  int KK = INTEGER(K)[0];
  double *weight_ptr = NUMERIC_POINTER(weights);
  fvec weights_c(weight_ptr, weight_ptr + GG*KK);
  uvec z(GG);
  GetRNGstate();
  draw_z(weights_c, z, GG, KK);
  PutRNGstate();
  
  SEXP result = PROTECT(allocVector(INTSXP, GG));
  for(int i=0; i<GG; i++)
    INTEGER(result)[i] = z[i];
  UNPROTECT(1);
  return result;
}

extern "C" SEXP cluster_statsR(SEXP z, SEXP G, SEXP K){
  int GG = INTEGER(G)[0];
  int KK = INTEGER(K)[0];
  int *z_ptr = INTEGER_POINTER(z);
  
  SEXP result = PROTECT(allocVector(INTSXP, KK));
  uvec z_c(z_ptr, z_ptr + GG);
  ivec Gk(KK);
  
  cluster_stats(z_c, Gk, GG, KK);
  for(int i=0; i<KK; i++){
    INTEGER(result)[i] = Gk[i];
  }
  UNPROTECT(1);
  return result;
}
  
  