#ifndef RCATEGORICAL_H
#define RCATEGORICAL_H

#include <numeric> //accumulate
#include <algorithm> //transform
#include <functional>
#include <Rmath.h> //r_unif()
// #include <random> //discrete distribution

double logsumexp(double sum, double addend){
  double newsum = log(exp(sum) + exp(addend));
  return newsum;
}

void norm_exp_lp(fveci lp, int K){
  double max = *std::max_element(lp, lp + K);
  std::transform(lp, lp + K, lp, std::bind2nd(std::minus<double>(), max));
  //   double sum = std::accumulate(lp, lp + K, 0, logsumexp);
  //   std::transform(lp, lp + K, lp, std::bind2nd(std::minus<double>(), sum));
  std::transform(lp, lp + K, lp, exp);
}

unsigned int rcategorical(fveci lp, int K){
  norm_exp_lp(lp, K);
  std::partial_sum(lp, lp + K, lp, std::plus<double>());
  
  double u = runif(0,1) * lp[K-1]; //rand_r(&seed)/ (double) RAND_MAX * lp[K-1];
  // printf("u is %lf",u);
  unsigned int z=0;
  while(u>=lp[z]){
    z++;
  }
  return z;
}

#endif