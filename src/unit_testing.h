#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include "types.h"
#include "cluster_stats.h"

#define Ralloc_Real(len)  PROTECT(allocVector(REALSXP, (len)))
#define Ralloc_Int(len)   PROTECT(allocVector(INTSXP, (len)))
#define Ralloc_List(len)  PROTECT(allocVector(VECSXP, (len)))


extern "C"{
  
  SEXP construct_precision_matR(SEXP xTx, SEXP Gk, SEXP sigma2, SEXP lambda2, SEXP V){
    
    double *xTx_p = NUMERIC_POINTER(xTx),
           Gkk   = REAL(Gk)[0],
           sigma2C = REAL(sigma2)[0],
           lambda2C = REAL(lambda2)[0];
    int VV = INTEGER(V)[0];
    fvec xTxC(xTx_p, xTx_p + VV*VV);
    fvec Sinv(VV*VV);
    construct_precision_mat(xTxC, Sinv, Gkk, sigma2C, lambda2C, VV);
    SEXP out = Ralloc_Real(VV*VV);
    for(int vv=0; vv<VV*VV; vv++){
      REAL(out)[vv] = Sinv[vv];
    }
    UNPROTECT(1);
    return out;
  }
  
  SEXP calculate_IGscaleR( SEXP yTy, SEXP xTy, SEXP xTx, SEXP beta, SEXP Gk, SEXP V){
    double yTyk = REAL(yTy)[0],
           *xTy_p = NUMERIC_POINTER(xTy),
           *xTx_p = NUMERIC_POINTER(xTx),
           *beta_p = NUMERIC_POINTER(beta),
           GkC = REAL(Gk)[0];
    int VV = INTEGER(V)[0];
    double IGscale;
    fvec xTyk(xTy_p, xTy_p + VV);
    fvec xTxC(xTx_p, xTx_p + VV*VV);
    fvec betak(beta_p, beta_p + VV);
    calculate_IGscale(&IGscale, yTyk, xTyk, xTxC, betak, GkC, VV);
    SEXP out = Ralloc_Real(1);
    REAL(out)[0] = IGscale;
    UNPROTECT(1);
    return out;
  }

}