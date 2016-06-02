#ifndef BLAS_H
#define BLAS_H
//#include <R_ext/BLAS.h>
// extern void daxpy(int*, double*, double*, int*, double*, int*); 
extern "C" double ddot_(int *n, double *x, int *incx, double *y, int *incy);
extern "C" void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
extern "C" void dsymv_(char *uplo, int *n, double *alpha, double *A, int *lda, double *x,
            int *incx, double *beta, double *y, int *incy);

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


/* Write custom function calling BLAS routines */
double quad_form(int n, double *x, double *A){
  double *y, out;
  y = (double *) calloc(n, sizeof(double));
  dsymv(n, 1, A, x, 1, y); //y <- A%*%x
  out = ddot(n, x, y);
  return out;
}

#endif//BLAS_H