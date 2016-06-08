#ifndef TRANSPOSE_H
#define TRANSPOSE_H

// "matin" is (column major) M by N matrix
// "matout is (column major) N by M matrix
// reference: http://stackoverflow.com/questions/16737298/what-is-the-fastest-way-to-transpose-a-matrix-in-c
void transpose(double *matin, double *matout, int M, int N){
int n, i, j;
#pragma omp parallel for
  for(n = 0; n < M*N; n++){
    i = n/N;
    j = n%N;
    matout[n] = matin[M*j + i];
  }
}

#endif