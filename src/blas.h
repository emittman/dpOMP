#ifndef BLAS_H
#define BLAS_H
//#include <R_ext/BLAS.h>
#include "types.h"
// extern void daxpy(int*, double*, double*, int*, double*, int*); 
extern "C"{
  double ddot_(int *n, double *x, int *incx, double *y, int *incy);
  void daxpy_(int *n, double *alpha, const double *x, int *incx, double *y, int *incy);
  void dsymv_(char *uplo, int *n, double *alpha, double *A, int *lda, double *x,
            int *incx, double *beta, double *y, int *incy);
  void dgemv_(char *trans, int *m, int *n, double *alpha, const double *A, int *lda, const double *x,
                      int *incx, double *beta, double *y, int *incy);
  void dposv_(char *uplo, int *n, int *nrhs, double *A, int *lda, double *B, int *ldb, int *info);
  
  void dpotri_(char *uplo, int *n, double *A, int *lda, int *info);
  void dtrtri_(char *uplo, char *diag, int *n, double *A, int *lda, int *info);
  void dtrmv_(char *uplo, char *trans, char *diag, int *N, const double *A, int *lda, double *x, int *incx);
}

void multiply_lowertri_vec(int n, const fvec &L, fvec &x){
  char uplo[] = "L", trans[] = "N", diag[] = "N";
  int lda = n, incx = 1;
  dtrmv_(uplo, trans, diag, &n, &(L[0]), &lda, &(x[0]), &incx);
}

int solve_normaleq_symm_mat(int n, double *A, double *B){
  /*solve Ax = B; store solution in B, chol(A) in A */
  char uplo[] = "L";
  int nrhs = 1;
  int lda = n;
  int ldb = n;
  int info = 0;
  dposv_(uplo, &n, &nrhs, A, &lda, B, &ldb, &info);
  return info;
}

/* Find inverse of symmetric matrix given cholesky factor*/
int dpotri(int n, double *A){
  char uplo[] = "L";
  int lda = n;
  int info = 0;
  dpotri_(uplo, &n, A, &lda, &info);
  return info;
}

/* Find inverse of lower triangular matrix*/
int invert_lower_tri(int n, double *A){
  char uplo[] = "L", diag[] = "N";
  int lda = n;
  int info = 0;
  dtrtri_(uplo, diag, &n, A, &lda, &info);
  return info;
}

double inner_prod_vec(int n, double *x, double *y){
  int incx = 1, incy = 1;
  double out;
  out = ddot_(&n, x, &incx, y, &incy);
  return out;
}

void linear_comb_vec(int n, double alpha, const fvec &x, fvec &y){
  /* y <- alpha * x + y; scalars{alpha}, vectors{x, y}*/
  int incx=1, incy=1;
  daxpy_(&n, &alpha, &(x[0]), &incx, &(y[0]), &incy);
}

void linear_comb_mat_vec(int n, double alpha, double *A, double *x, double beta, double *y){
  /* y <- alpha * A %*% x + beta * y; scalars{alpha, beta}, vectors{x, y}, matrices{A}*/
 char uplo[] = "U";
 int lda = n, incx = 1, incy = 1;
 dsymv_(uplo, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}

void multiply_mat_vec(int m, int n, const fvec &A, const fvec &invec, fvec &outvec, int tr){
  /* y <- alpha * A %*% x + beta * y; scalars{alpha, beta}, vectors{x, y}, matrices{A}*/
  char trans = tr > 0 ? 'T' : 'N';
  int lda = m, incx =1, incy = 1;
  double alpha = 1.0, beta = 1.0;
  dgemv_(&trans, &m, &n, &alpha, &(A[0]), &lda, &(invec[0]), &incx, &beta, &(outvec[0]), &incy);
}

/* Write custom function calling BLAS routines */
double quad_form(int n, double *x, double *A){
  double out;
  fvec y(n);
  linear_comb_mat_vec(n, 1, A, x, 1, &(y[0])); //y <- A%*%x
  out = inner_prod_vec(n, x, &(y[0]));
  return out;
}

#endif//BLAS_H