#include <vector>
#include <iostream>

typedef std::vector<float> fvec;
typedef std::vector<unsigned int> uvec;
typedef std::vector<int> ivec;


int main(){
  int N = 10;
  
  uvec vec(N);
  
  #pragma omp parallel for
  for(int i=0; i<N; i++){
    vec[i] = i * i;
    std::cout <<"tid: " << i <<"\n";
  }
  
  for(int i=0; i<6; i++){
    std::cout << "In thread " << i << " the value is: " << vec[i] << ".\n\n";
  }
  
  return 0;
}
