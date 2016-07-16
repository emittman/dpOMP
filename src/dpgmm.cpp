#include "includes.h"

#define Ralloc_Real(len)  PROTECT(allocVector(REALSXP, (len)))
#define Ralloc_Int(len)   PROTECT(allocVector(INTSXP, (len)))
#define Ralloc_List(len)  PROTECT(allocVector(VECSXP, (len)))

extern "C" SEXP dpgmm_initR(SEXP yTyR, SEXP xTyR, SEXP xTxR, SEXP betaR, SEXP piR,
                           SEXP sigma2R, SEXP lambda2R, SEXP alphaR, SEXP G, SEXP V, SEXP K, SEXP N, SEXP iter){
  int GG, KK, NN, VV, beta_len, I;
  GG = INTEGER(G)[0];
  KK = INTEGER(K)[0];
  NN = INTEGER(N)[0];
  VV = INTEGER(V)[0];
  I = INTEGER(iter)[0];
  beta_len = KK*VV;
  
  double *yTy_p, *xTy_p, *xTx_p, *beta_p, *pi_p, lambda2, alpha;
  yTy_p = NUMERIC_POINTER(yTyR);
  xTy_p = NUMERIC_POINTER(xTyR);
  xTx_p = NUMERIC_POINTER(xTxR);
  beta_p = NUMERIC_POINTER(betaR);
  pi_p = NUMERIC_POINTER(piR);
  lambda2 = REAL(lambda2R)[0];
  alpha = REAL(alphaR)[0];
  
  chain_t chain = construct_chain(yTy_p, xTy_p, xTx_p, lambda2, alpha, GG, KK, VV, NN);
  
  //initialize_chain(chain);
  std::copy(beta_p, beta_p + beta_len, chain.beta.begin());
  std::copy(pi_p, pi_p + KK, chain.pi.begin());
  chain.sigma2 = REAL(sigma2R)[0];
  
//   SEXP beta_out = Ralloc_Real(beta_len*I);
  SEXP beta_g_out = Ralloc_Real(GG*VV*I);
//   SEXP sigma2_out = Ralloc_Real(KK*I);
  SEXP sigma2_out = Ralloc_Real(I);
//   SEXP pi_out = Ralloc_Real(KK*I);
  SEXP max_index = Ralloc_Int(I);
  SEXP list_out = Ralloc_List(3);
  
  for(int i=0; i<I; i++){
    // print_mat(weights, GG, KK);
    compute_weights(chain);
    draw_z(chain);
    draw_theta(chain);
    draw_pi(chain, chain.alpha);

    int offset = i*beta_len;
//     for(int j=0; j<beta_len; j++)
//       REAL(beta_out)[offset + j] = chain.beta[j];
//     
//     offset = i*KK;
//     for(int j=0; j<KK; j++){
//       REAL(pi_out)[offset + j] = chain.pi[j];
//       REAL(sigma2_out)[offset + j] = chain.sigma2[j];
//     }
    
    offset = i*GG*VV;
    for(int j=0; j<GG; j++)
      for(int k=0; k<VV; k++)
        REAL(beta_g_out)[offset + j*VV + k] = chain.beta[chain.z[j]*VV + k];
    
    REAL(sigma2_out)[i] = chain.sigma2;
    
    INTEGER(max_index)[i] = (int) *max_element(chain.z.begin(), chain.z.end());
  }
  
  // SET_VECTOR_ELT(list_out, 0, beta_out);
  // SET_VECTOR_ELT(list_out, 1, pi_out);
  SET_VECTOR_ELT(list_out, 0, beta_g_out);
  // SET_VECTOR_ELT(list_out, 3, sigma2_out);
  SET_VECTOR_ELT(list_out, 1, sigma2_out);
  SET_VECTOR_ELT(list_out, 2, max_index);
  
  UNPROTECT(4);
  
  return list_out;
}

extern "C" SEXP dpgmmR(SEXP yTyR, SEXP xTyR, SEXP xTxR, SEXP lambda2R, SEXP alphaR, SEXP G, SEXP V, SEXP K, SEXP N, SEXP iter){
  int GG, KK, NN, VV, beta_len, I;
  GG = INTEGER(G)[0];
  KK = INTEGER(K)[0];
  NN = INTEGER(N)[0];
  VV = INTEGER(V)[0];
  I = INTEGER(iter)[0];
  beta_len = KK*VV;
  
  double *yTy_p, *xTy_p, *xTx_p, lambda2, alpha;
  yTy_p = NUMERIC_POINTER(yTyR);
  xTy_p = NUMERIC_POINTER(xTyR);
  xTx_p = NUMERIC_POINTER(xTxR);
  lambda2 = REAL(lambda2R)[0];
  alpha = REAL(alphaR)[0];
  
  chain_t chain = construct_chain(yTy_p, xTy_p, xTx_p, lambda2, alpha, GG, KK, VV, NN);
  
  initialize_chain(chain);
  
  SEXP beta_out = Ralloc_Real(beta_len*I);
//   SEXP beta_g_out = Ralloc_Real(GG*VV*I);
  SEXP sigma2_out = Ralloc_Real(I);
  // SEXP sigma2_g_out = Ralloc_Real(GG*I);
  SEXP pi_out = Ralloc_Real(KK*I);
  SEXP list_out = Ralloc_List(3);
  
  for(int i=0; i<I; i++){
    // print_mat(weights, GG, KK);
    compute_weights(chain);
    draw_z(chain);
    //     Rprintf("allocation iter %d:\n",i);
    //     print_mat(z, 1, GG);
    // print_mat(weights, GG, KK);
    // print_mat(z, 1, GG);
    draw_theta(chain);
    draw_pi(chain, chain.alpha); //second arg is alpha of DP(alpha P_0) fame
    // Rprintf("iter %d: sigma2 = %lf\n", i, chain.sigma2);
    
    int offset = i*beta_len;
    for(int j=0; j<beta_len; j++){
      for(int k=0; k<KK; k++){
        REAL(beta_out)[offset + j] = chain.beta[j];
      }
    }
    
    offset = i*KK;
    for(int j=0; j<KK; j++){
      REAL(pi_out)[offset + j] = chain.pi[j];
    }

    REAL(sigma2_out)[i] = chain.sigma2;
    
        
//     offset = i*GG*VV;
//     for(int j=0; j<GG; j++)
//       for(int k=0; k<VV; k++)
//         REAL(beta_g_out)[offset + j*VV + k] = chain.beta[chain.z[j]*VV + k];
//     
//     offset = i*GG;
//     for(int j=0; j<GG; j++)
//       REAL(sigma2_g_out)[offset + j] = chain.sigma2[chain.z[j]];
    
  }
  SET_VECTOR_ELT(list_out, 0, beta_out);
  // SET_VECTOR_ELT(list_out, 2, beta_g_out);
  SET_VECTOR_ELT(list_out, 1, sigma2_out);
  SET_VECTOR_ELT(list_out, 2, pi_out);
  // SET_VECTOR_ELT(list_out, 4, sigma2_g_out);
  
  UNPROTECT(4);
  
  return list_out;
}

extern "C" SEXP pi_primeR(SEXP yTx, SEXP xTx, SEXP beta, SEXP pi, SEXP sigma2, SEXP V){
  double *yTx_p, *xTx_p, *beta_p, pi_c, sigma2_c;
  int v = INTEGER(V)[0];
  // yTy_c = REAL(yTy)[0];
  yTx_p  = NUMERIC_POINTER(yTx);
  xTx_p  = NUMERIC_POINTER(xTx);
  beta_p = NUMERIC_POINTER(beta);
  fvec xTx_c(xTx_p, xTx_p + v*v);
  pi_c = REAL(pi)[0];
  sigma2_c = REAL(sigma2)[0];
  
  SEXP out = Ralloc_Real(1);
  REAL(out)[0] = pi_prime(yTx_p, xTx_c, beta_p, pi_c, sigma2_c, v);
  
  UNPROTECT(1);
  return out;
}

extern "C" SEXP compute_weightsR(SEXP yTy, SEXP xTy, SEXP xTx, SEXP beta, SEXP pi, SEXP sigma2, SEXP G, SEXP K, SEXP V){
  double *yTy_p, *xTy_p, *xTx_p, *beta_p, *pi_p;
  int GG, KK, VV;
  GG = INTEGER(G)[0];
  KK = INTEGER(K)[0];
  VV = INTEGER(V)[0];
  int NN = 1; //arbitrary
  yTy_p  = NUMERIC_POINTER(yTy);
  xTy_p  = NUMERIC_POINTER(xTy);
  xTx_p  = NUMERIC_POINTER(xTx);
  chain_t chain = construct_chain(yTy_p, xTy_p, xTx_p, 1, 1, GG, KK, VV, NN);
  
  beta_p = NUMERIC_POINTER(beta);
  pi_p   = NUMERIC_POINTER(pi);
  std::copy(beta_p, beta_p + KK*VV, chain.beta.begin());
  std::copy(pi_p, pi_p + KK, chain.pi.begin());
  chain.sigma2 = REAL(sigma2)[0];
  compute_weights(chain);
  
  SEXP out = Ralloc_Real(GG*KK);
  for(int i=0; i<GG*KK; i++)
    REAL(out)[i] = chain.weights[i];
  UNPROTECT(1);
  return out;
}

extern "C" SEXP ddotR(SEXP n, SEXP x, SEXP y){
  int N = INTEGER(n)[0];
  double *X = NUMERIC_POINTER(x), *Y = NUMERIC_POINTER(y);
  SEXP result = Ralloc_Real(1);
  REAL(result)[0] = inner_prod_vec(N, X, Y);
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
  double yTy = 0;
  fvec xTy(GG);
  fvec xTx(1);
  chain_t chain = construct_chain(&yTy, &(xTy[0]), &(xTx[0]), 1, 1, GG, KK, 1, 1);
  std::copy(weight_ptr, weight_ptr + GG*KK, chain.weights.begin());
  
  GetRNGstate();
  draw_z(chain);
  PutRNGstate();
  
  SEXP result = PROTECT(allocVector(INTSXP, GG));
  for(int i=0; i<GG; i++)
    INTEGER(result)[i] = chain.z[i];
  UNPROTECT(1);
  return result;
}


extern "C" SEXP cluster_sumsR(SEXP k, SEXP xTy, SEXP G, SEXP V, SEXP z){
  int kk = INTEGER(k)[0];
  int GG = INTEGER(G)[0];
  int VV = INTEGER(V)[0];
  int *z_ptr = INTEGER_POINTER(z);
  double *xTy_ptr = NUMERIC_POINTER(xTy);

  fvec xTy_c(xTy_ptr, xTy_ptr + GG*VV);
  uvec z_c(z_ptr, z_ptr + GG);
  fvec Gkk(1);
  fvec xTyk(VV);
  cluster_sums(kk, xTy_c, GG, VV, z_c, Gkk.begin(), xTyk);
  
  SEXP result = Ralloc_List(2);
  SEXP GkR = Ralloc_Real(1);
  SEXP xTykR = Ralloc_Real(VV);
  
  REAL(GkR)[0] = Gkk[0];
  for(int i=0; i<VV; i++){
    REAL(xTykR)[i] = xTyk[i];
  }
  
  SET_VECTOR_ELT(result, 0, GkR);
  SET_VECTOR_ELT(result, 1, xTykR);
  
  
  UNPROTECT(3);
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
  double yTy = 0;
  fvec xTy(1);
  fvec xTx(1);
  chain_t chain = construct_chain(&yTy, &(xTy[0]), &(xTx[0]), 1, 1, 1, KK, 1, 1);
  std::copy(Gk_ptr, Gk_ptr + KK, chain.Gk.begin());
  GetRNGstate();
  draw_pi(chain, 1.0);
  PutRNGstate();
  SEXP out = Ralloc_Real(KK);
  for(int i=0; i<KK; i++)
    REAL(out)[i] = chain.pi[i];
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
  double yTy_c = REAL(yTy)[0];
  double *xTy_p = NUMERIC_POINTER(xTy);
  double *xTx_p = NUMERIC_POINTER(xTx);
  // assume that lambda2 = 1.0
  chain_t chain = construct_chain(&yTy_c, xTy_p, xTx_p, 1.0, 1.0, GG, KK, VV, NN);
  initialize_chain(chain);
  std::copy(z_p, z_p+GG, chain.z.begin());
  GetRNGstate();
  draw_theta(chain);
  PutRNGstate();
  SEXP out = Ralloc_List(2);
  SEXP beta = Ralloc_Real(KK*VV);
  SEXP sigma2 = Ralloc_Real(1);

  for(int i=0; i<(KK*VV); i++)
    REAL(beta)[i] = chain.beta[i];

  REAL(sigma2)[0] = chain.sigma2;
  
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

extern "C" SEXP send_chain(SEXP yTy, SEXP xTy, SEXP xTx, SEXP lambda2R, SEXP alphaR, SEXP G, SEXP K, SEXP V, SEXP N){
  int GG = INTEGER(G)[0];
  int KK = INTEGER(K)[0];
  int VV = INTEGER(V)[0];
  int NN = INTEGER(N)[0];
  double *yTy_p = NUMERIC_POINTER(yTy);
  double *xTy_p = NUMERIC_POINTER(xTy);
  double *xTx_p = NUMERIC_POINTER(xTx);
  double lambda2 = REAL(lambda2R)[0];
  double alpha = REAL(alphaR)[0];
  chain_t chain = construct_chain(yTy_p, xTy_p, xTx_p, lambda2, alpha, GG, KK, VV, NN);
  initialize_chain(chain);
  print_chain_state(chain);
  SEXP out = Ralloc_Int(1);
  INTEGER(out)[0] = 0;
  UNPROTECT(1);
  return out;
}                                                