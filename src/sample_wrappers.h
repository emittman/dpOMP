#ifndef SAMPLE_WRAPPERS_H
#define SAMPLE_WRAPPERS_H
#include <Rmath.h>

double rinvgamma(double shape, double scale){
  double out = 1/rgamma(shape, 1/scale);
  return out;
}

#endif