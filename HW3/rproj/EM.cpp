#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector exponential_mixture_EM(NumericVector y,
                                     NumericVector inits,
                                     double tol = 1e-4,
                                     int max_iter = 1000){
  double n = (double)y.size();
  double err = 0.0;

  NumericVector curr(3);
  NumericVector new_val(3);
  
  NumericVector a(n, 0);
  NumericVector b(n, 0);
  NumericVector w(n, 0);
  
  curr = inits;
  for(int i = 0; i < max_iter; i++){
    a = curr(0) * dexp(y, 1.0 / curr(1));
    b = (1.0 - curr(0)) * dexp(y, 1.0 / curr(2));
    w = a / (a + b);

    new_val(0) = sum(w) / n;
    new_val(1) = sum(w) / sum(w * y);
    new_val(2) = sum(1.0 - w) / sum((1.0 - w) * y);

    err = sqrt(sum(pow(new_val - curr, 2.0)));

    curr = new_val;

    if(err < tol){
      break;
    }
  }

  return(curr);
}
