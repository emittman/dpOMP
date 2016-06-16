#include "includes.h"

#define Ralloc_Real(len)  PROTECT(allocVector(REALSXP, (len)))
#define Ralloc_List(len)  PROTECT(allocVector(VECSXP, (len)))

extern "C" SEXP dpgmmR(SEXP yTyR, SEXP xTyR, SEXP xTxR, SEXP G, SEXP K, SEXP N, SEXP V){
  int GG, KK, NN, VV, xTy_len, xTx_len, beta_len;
  GG = INTEGER(G)[0];
  KK = INTEGER(K)[0];
  NN = INTEGER(N)[0];
  VV = INTEGER(V)[0];
  
  xTy_len = GG*VV;
  xTx_len = VV*VV;
  beta_len = KK*VV;

  double yTy, *xTy_p, *xTx_p;
  yTy = REAL(yTyR)[0];
  xTy_p = NUMERIC_POINTER(xTyR);
  xTx_p = NUMERIC_POINTER(xTxR);
  
  fvec xTy(xTy_p, xTy_p + xTy_len);
  fvec xTx(xTx_p, xTx_p + xTx_len);
  
  fvec beta(KK*VV);
  fvec pi(KK);
  uvec z(GG);
  fvec weights(GG*KK);
  fvec Gk(KK);
  double sigma2 = 1;
  
  SEXP beta_out = Ralloc_Real(beta_len);
  SEXP pi_out = Ralloc_Real(KK);
  SEXP list_out = Ralloc_List(2);
  
  draw_z(weights, z, GG, KK);
  draw_theta(beta, &sigma2, z, &yTy, xTy, xTx, Gk, GG, KK, VV, NN);
  draw_pi(pi, Gk, KK, 1.0);
  
  for(int i=0; i<beta_len; i++)
    REAL(beta_out)[i] = beta[i];
  
  for(int i=0; i<KK; i++)
    REAL(pi_out)[i] = pi[i];
  
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

extern "C" SEXP ddotR(SEXP n, SEXP x, SEXP y){
  int N = INTEGER(n)[0];
  double *X = NUMERIC_POINTER(x), *Y = NUMERIC_POINTER(y);
  SEXP result = Ralloc_Real(1);
  REAL(result)[0] = ddot(N, X, Y);
  UNPROTECT(1);
  return result;
}


//A is symmetric matrix, x is vector, compute xTAx
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

  
extern "C" SEXP cluster_statsR(SEXP k, SEXP xTy, SEXP xTx, SEXP G, SEXP V, SEXP n, SEXP z){
  int kk = INTEGER(k)[0];
  int GG = INTEGER(G)[0];
  int VV = INTEGER(V)[0];
  int nn = INTEGER(n)[0];
  int *z_ptr = INTEGER_POINTER(z);
  double *xTy_ptr = NUMERIC_POINTER(xTy);
  double *xTx_ptr = NUMERIC_POINTER(xTx);

  fvec xTx_c(xTx_ptr, xTx_ptr + VV*VV);  
  fvec xTy_c(xTy_ptr, xTy_ptr + GG*VV);
  uvec z_c(z_ptr, z_ptr + GG);
  fvec Gkk(1);
  fvec beta_hat(VV);
  fvec chol_S(VV*VV);
  fvec IGscale(1);
  cluster_stats(kk, xTy_c, xTx_c, GG, VV, nn, z_c, Gkk.begin(), beta_hat, chol_S, IGscale.begin());
  Rprintf("%d", z_c[0]);
  
  SEXP result = Ralloc_List(4);
  SEXP GkR = Ralloc_Real(1);
  SEXP beta_hatR = Ralloc_Real(VV);
  SEXP chol_SR = Ralloc_Real(VV*VV);
  SEXP IGscaleR = Ralloc_Real(1);
  
  REAL(GkR)[0] = Gkk[0];
  
  for(int i=0; i<VV; i++){
    REAL(beta_hatR)[i] = beta_hat[i];
  }
  for(int i=0; i<VV*VV; i++){
    REAL(chol_SR)[i] = chol_S[i];
  }
  
  REAL(IGscaleR)[0] = IGscale[0];
  
  SET_VECTOR_ELT(result, 0, GkR);
  SET_VECTOR_ELT(result, 1, beta_hatR);
  SET_VECTOR_ELT(result, 2, chol_SR);
  SET_VECTOR_ELT(result, 3, IGscaleR);
  
  
  UNPROTECT(5);
  return result;
}
  
extern "C" SEXP transposeR(SEXP matin, SEXP M, SEXP N){
    int MM = INTEGER(M)[0];
    int NN = INTEGER(N)[0];
    double *in_ptr = NUMERIC_POINTER(matin);
    
    SEXP result = PROTECT(allocVector(REALSXP, MM*NN));
    double *out_ptr = NUMERIC_POINTER(result);
    transpose(in_ptr, out_ptr, MM, NN);
    
    UNPROTECT(1);
    return result;
}

extern "C" SEXP solve_normaleq_symm_matR(SEXP n, SEXP A, SEXP B){
  int N = INTEGER(n)[0];
  int info;
  double *AA = NUMERIC_POINTER(A);
  double *BB = NUMERIC_POINTER(B);
  
  fvec a(AA, AA + N*N);
  fvec b(BB, BB + N);
  
  SEXP out  = Ralloc_List(3);
  SEXP chol = Ralloc_Real(N*N);
  SEXP sol  = Ralloc_Real(N);
  SEXP status = PROTECT(allocVector(INTSXP, 1));
  
  info = solve_normaleq_symm_mat(N, &(a[0]), &(b[0]));
  if(info == 0){
    info = dpotri(N, &(a[0]));
  }
  
  INTEGER(status)[0] = info;
  
  for(int i=0; i<N; i++)
    REAL(sol)[i] = b[i];
  
  for(int j=0; j<N*N; j++)
    REAL(chol)[j] = a[j];
  
  SET_VECTOR_ELT(out, 0, chol);
  SET_VECTOR_ELT(out, 1, sol);
  SET_VECTOR_ELT(out, 2, status);
  
  UNPROTECT(4);
  return out;
}

extern "C" SEXP get_k_indicesR(SEXP z, SEXP G, SEXP k){
  int kk = INTEGER(k)[0];
  int GG = INTEGER(G)[0];
  int *z_ptr = INTEGER_POINTER(z);
  uvec z_c(z_ptr, z_ptr+GG);
  fvec indicC(GG);
  get_k_indices(z_c, indicC, kk);
  
  SEXP result = Ralloc_Real(GG);
  for(int i=0; i<GG; i++) REAL(result)[i] = indicC[i];
  UNPROTECT(1);
  return result;
}

extern "C" SEXP rinvgammaR(SEXP n, SEXP shape, SEXP scale){
  int N = INTEGER(n)[0];
  double Sh = REAL(shape)[0];
  double Sc = REAL(scale)[0];
  SEXP result = Ralloc_Real(N);
  GetRNGstate();
  for(int i=0; i<N; i++)
    REAL(result)[i] = rinvgamma(Sh, Sc);
  PutRNGstate();
  UNPROTECT(1);
  return result;
}

extern "C" SEXP tail_sums(SEXP counts, SEXP len){
  int K = INTEGER(len)[0];
  double *Counts = NUMERIC_POINTER(counts);
  fvec C(Counts, Counts + K);
  fvec Tail(K);
  get_tail_sums(C, Tail);
  SEXP result = Ralloc_Real(K);
  for(int i=0; i<K; i++)
    REAL(result)[i] = Tail[i];
  UNPROTECT(1);
  return result;
}

extern "C" SEXP draw_piR(SEXP GK, SEXP K){
  int KK = INTEGER(K)[0];
  double *Gk_ptr = NUMERIC_POINTER(GK);
  fvec Gk(Gk_ptr, Gk_ptr + KK);
  fvec pi(KK);
  GetRNGstate();
  draw_pi(pi, Gk, KK, 1.0);
  PutRNGstate();
  SEXP out = Ralloc_Real(KK);
  for(int i=0; i<KK; i++)
    REAL(out)[i] = pi[i];
  UNPROTECT(1);
  return out;
}

extern "C" SEXP draw_thetaR(SEXP z, SEXP yTy, SEXP xTy, SEXP xTx, SEXP G, SEXP K, SEXP V, SEXP N){
  int GG, KK, VV, NN;
  GG = INTEGER(G)[0];
  KK = INTEGER(K)[0];
  VV = INTEGER(V)[0];
  NN = INTEGER(N)[0];
  int *z_p = INTEGER_POINTER(z);
  double yTyC = REAL(yTy)[0];
  double *xTy_p = NUMERIC_POINTER(xTy);
  double *xTx_p = NUMERIC_POINTER(xTx);
  uvec zC(z_p, z_p + GG);
  fvec xTyC(xTy_p, xTy_p + GG*VV);
  fvec xTxC(xTx_p, xTx_p + VV*VV);
  
  fvec betaC(KK*VV);
  double sigma2C = 1;
  fvec Gk(KK);
  
  GetRNGstate();
  draw_theta(betaC, &sigma2C, zC, &yTyC, xTyC, xTxC, Gk, GG, KK, VV, NN);
  PutRNGstate();
  
  SEXP out = Ralloc_List(2);
  SEXP beta = Ralloc_Real(KK*VV);
  SEXP sigma2 = Ralloc_Real(1);
  
  for(int i=0; i<(KK*VV); i++)
    REAL(beta)[i] = betaC[i];
  REAL(sigma2)[0] = sigma2C;
  
  SET_VECTOR_ELT(out, 0, beta);
  SET_VECTOR_ELT(out, 1, sigma2);
  UNPROTECT(3);
  return out;
}

extern "C" SEXP linear_comb_vecR(SEXP n, SEXP alpha, SEXP x, SEXP y){
  int N = INTEGER(n)[0];
  double a = REAL(alpha)[0];
  double *X = NUMERIC_POINTER(x);
  double *Y = NUMERIC_POINTER(y);
  fvec xC(X, X + N);
  fvec yC(Y, Y + N);
  linear_comb_vec(N, a, xC, yC);
  SEXP yout = Ralloc_Real(N);
  for(int i=0; i<N; i++)
    REAL(yout)[i] = yC[i];
  UNPROTECT(1);
  return yout;
}