#include <armadillo>
#include <iostream>

int main(){
  int N = 10;
  
  arma::uvec vec(N);
  
#pragma omp parallel for
  for(int i=0; i<N; i++)
    vec[i] = i * i;
  
  #pragma omp parallel for
  for(int i=0; i<6; i++){
    std::cout << "In thread " << i << " the value is: " << vec[i] << ".\n\n";
  }
  
  return 0;
}