#ifndef BLAS_H
#define BLAS_H
//#include <R_ext/BLAS.h>
// extern void daxpy(int*, double*, double*, int*, double*, int*); 
extern "C"{
  double ddot_(int *n, double *x, int *incx, double *y, int *incy);
  void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
  void dsymv_(char *uplo, int *n, double *alpha, double *A, int *lda, double *x,
            int *incx, double *beta, double *y, int *incy);
  void dgemv_(char *trans, int *m, int *n, double *alpha, double *A, int *lda, double *x,
                      int *incx, double *beta, double *y, int *incy);
  void dposv_(char *uplo, int *n, int *nrhs, double *A, int *lda, double *B, int *ldb, int *info);
  
  void dpotri_(char *uplo, int *n, double *A, int *lda, int *info);
}

int dposv(int n, double *A, double *B){
  char uplo[] = "L";
  int nrhs = 1;
  int lda = n;
  int ldb = n;
  int info = 0;
  dposv_(uplo, &n, &nrhs, A, &lda, B, &ldb, &info);
  return info;
}

int dpotri(int n, double *A){
  char uplo[] = "L";
  int lda = n;
  int info = 0;
  dpotri_(uplo, &n, A, &lda, &info);
  return info;
}

double ddot(int n, double *x, double *y){
  int incx = 1, incy = 1;
  double out;
  out = ddot_(&n, x, &incx, y, &incy);
  return out;
}

void daxpy(int n, double alpha, double *x, double *y){
  /* y <- alpha * x + y; scalars{alpha}, vectors{x, y}*/
  int incx=1, incy=1;
  daxpy_(&n, &alpha, x, &incx, y, &incy);
}

void dsymv(int n, double alpha, double *A, double *x, double beta, double *y){
  /* y <- alpha * A %*% x + beta * y; scalars{alpha, beta}, vectors{x, y}, matrices{A}*/
 char uplo[] = "U";
 int lda = n, incx = 1, incy = 1;
 dsymv_(uplo, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}

void dgemv(int m, int n, double alpha, double *A, double *x, double beta, double *y){
  /* y <- alpha * A %*% x + beta * y; scalars{alpha, beta}, vectors{x, y}, matrices{A}*/
char trans[] = "N";
  int lda = m, incx =1, incy = 1;
  dgemv_(trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}

/* Write custom function calling BLAS routines */
double quad_form(int n, double *x, double *A){
  double *y, out;
  y = (double *) calloc(n, sizeof(double));
  dsymv(n, 1, A, x, 1, y); //y <- A%*%x
  out = ddot(n, x, y);
  return out;
}

#endif//BLAS_H