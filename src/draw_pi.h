#include "types.h"
#include <numeric>
#include <algorithm>
#include <Rmath.h>
#include "sample_wrappers.h"
#include "debug_print.h"

void get_tail_sums(const fvec &counts, fvec &tail_sums){
  
  // example:
  // in = (1, 2, 3, 4, 5) -> out: (2+3+4+5, 3+4+5, 4+5, 5, 0)
  // .rend() - 1   -> "don't include the first element in any sums"
  // .rbegin() + 1 -> "start writing at the second to last position"
  //print_fmat(tail_sum, 1, K);
  std::partial_sum(counts.rbegin(), counts.rend()-1, tail_sums.rbegin()+1);
}



void draw_pi(fvec &pi, const fvec &Gk, int K, double alpha){
  fvec tail_sums(K);
  double V;
  double log_1minusV = 0;
  //print_fmat(Gk, 1, K);
  get_tail_sums(Gk, tail_sums);
  
  for(int i=0; i<K; i++){
    if(i == 0){
      V = rbeta(Gk[i] + 1, tail_sums[i] + alpha);
      pi[i] = V;
      log_1minusV = log(1 - V);
    }else if(i == K-1){
      pi[i] = exp(log_1minusV);
    } else {
      V = rbeta(Gk[i] + 1, tail_sums[i] + alpha);
      pi[i] = V * exp(log_1minusV);
      log_1minusV = log(1 - V) + log_1minusV;
    }
  }
  // print_fmat(pi, 1, K);
}