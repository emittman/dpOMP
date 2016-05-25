#include <vector>
#include <iostream>

int main(){
  std::vector<float> vec = {0,1,2,3,4,5};
  
  #pragma omp parallel
  for(int i=0; i<6; i++){
    std::cout << "In thread " << i << " the value is: " << vec[i] << ".\n\n";
  }
  
  return 0;
}