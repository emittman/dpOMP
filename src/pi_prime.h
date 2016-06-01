#ifndef PI_PRIME_H
#define PI_PRIME_H

double quad_form(fvec &A, fvec &x, int dim);

double pi_prime(fvec &yTx, fvec &xTx, fvec &beta, double pi, double sigma2, int V){

  double out = log(pi)  - 1.0 / (2.0 * sigma2) * (quad_form(xTx, beta, V) -
                   2.0 * std::inner_product(yTx.begin(), yTx.end(), beta.begin(), 0.0));

  return out;
}

double quad_form(fvec &A, fvec &x, int dim){
  double out = 0;

  //diagonal elements:
  for(int i=0; i<dim; i++)
    out += x[i]*x[i]*A[i*dim + i];

  for(int i=1; i<dim; i++)
    for(int j=0; j<i; j++)
      out += 2 * x[i]*x[j]*A[i*dim + j];

  return out;
}

#endif // PI_PRIME_H
