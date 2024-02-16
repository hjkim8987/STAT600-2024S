#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Calculate the first derivative of the log-density of Cauchy distribution
//
// @param theta double
// @param x double
// @return first derivative of log-density of Cauchy for the given theta and x
// @export
// [[Rcpp::export]]
double get_ell_prime(double theta, double x){
  return 2.0 * (x - theta) / (1.0 + pow(x - theta, 2.0));
}

// Calculate the second derivative of the log-density of Cauchy distribution
//
// @param theta double
// @param x double
// @return second derivative of log-density of Cauchy for the given theta and x
// @export
// [[Rcpp::export]]
double get_ell_double_prime(double theta, double x){
  double deviation_squared = pow(x - theta, 2.0);
  return 2.0 * (deviation_squared - 1.0) /
    pow(1.0 + deviation_squared, 2.0);
}

// Calculate the first derivative of the log-likelihood
// using \ell = \sum log(density)
//
// @param theta double
// @param x arma::vec
// @param ftn Function
// @return the first derivative of the log-likelihood
// @export
// [[Rcpp::export]]
double get_log_likelihood(double theta, arma::vec x, Function ftn){
  int n = x.n_elem;
  double res = 0.0;
  
  for(int i = 0; i < n; i++){
    res += as<double> (ftn(theta, x(i)));
  }
  
  return res;
}

// Implement Bisection method to find the root of given `ftn`
//
// @param ftn Function
// @param x arma::vec
// @param a double
// @param b double
// @param tol double
// @return the root of given function, ftn
// @export
// [[Rcpp::export]]
List do_bisection(Function ftn,
                  arma::vec x,
                  double a, double b, double tol){
  double c;
  double f_a;
  double f_c;
  double eps = 1.0;
  
  int niter = 1;
  while(eps >= tol){
    c = (a + b) / 2.0;
    f_a = get_log_likelihood(a, x, ftn);
    f_c = get_log_likelihood(c, x, ftn);
    if(f_c == 0.0){
      break;
    }else if(f_c * f_a < 0.0){
      eps = abs(b - c) / abs(c);
      b = c;
    }else{
      eps = abs(a - c) / abs(c);
      a = c;
    }
    niter += 1;
  }
  return(List::create(
    Named("root") = c,
    Named("niter") = niter,
    Named("convergence") = eps
  ));
}

// Implement Newton-Raphson to find the root of given `ftn`
//
// @param ell_prime Function
// @param ell_double_prime Function
// @param x arma::vec
// @param init double
// @param tol double
// @return the root of given function, ftn
// @export
// [[Rcpp::export]]
List do_newton_raphson(Function ell_prime,
                       Function ell_double_prime,
                       arma::vec x,
                       double init, double tol){
  double new_val;
  double eps = 1.0;
  
  int niter = 1;
  while(eps >= tol){
    double l1 = get_log_likelihood(init, x, ell_prime);
    double l2 = get_log_likelihood(init, x, ell_double_prime);
    new_val = init - l1 / l2;
    eps = abs(new_val - init) / abs(init);
    init = new_val;
    niter += 1;
  }
  return(List::create(
      Named("root") = init,
      Named("niter") = niter,
      Named("convergence") = eps
  ));
}

// Implement Fisher Scoring to find the root of given `ftn`
//
// @param ell_prime Function
// @param ell_double_prime Function
// @param x arma::vec
// @param init double
// @param tol double
// @return the root of given function, ftn
// @export
// [[Rcpp::export]]
List do_fisher_scoring(Function ell_prime,
                       Function ell_double_prime,
                       arma::vec x,
                       double init, double tol){
  double new_val;
  double eps = 1.0;
  
  int niter = 1;
  while(eps >= tol){
    double l1 = get_log_likelihood(init, x, ell_prime);
    double I = -get_log_likelihood(init, x, ell_double_prime);
    new_val = init + l1 / I;
    eps = abs(new_val - init) / abs(init);
    init = new_val;
    niter += 1;
  }
  return(List::create(
      Named("root") = init,
      Named("niter") = niter,
      Named("convergence") = eps
  ));
}

// Implement Secant method to find the root of given `ftn`
//
// @param ftn Function
// @param x arma::vec
// @param x0 double
// @param x1 double
// @param tol double
// @return the root of given function, ftn
// @export
// [[Rcpp::export]]
List do_secant(Function ftn,
               arma::vec x,
               double x0, double x1, double tol){
  double x2;
  double eps = 1.0;
  
  int niter = 1;
  while(eps >= tol){
    double f_x0 = get_log_likelihood(x0, x, ftn);
    double f_x1 = get_log_likelihood(x1, x, ftn);
    x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);
    eps = abs(x2 - x1) / abs(x1);
    x0 = x1;
    x1 = x2;
    niter += 1;
  }
  return(List::create(
      Named("root") = x2,
      Named("niter") = niter,
      Named("convergence") = eps
  ));
}

// [[Rcpp::export]]
List logit(arma::vec y,
           arma::mat X,
           arma::vec beta_init,
           double tol){
  int n = y.n_elem;
  int p = X.n_cols;
  double eps = 1.0;
  
  arma::vec beta_new(p, arma::fill::zeros);
  arma::mat hessian(p, p, arma::fill::zeros);
  arma::vec yhat(n, arma::fill::zeros);
  arma::vec Xbeta(n, arma::fill::zeros);
  
  int niter = 1;
  while(eps >= tol){
    Xbeta = X * beta_init;
    for(int i = 0; i < n; i++){
      yhat(i) = exp(Xbeta(i)) / (1 + exp(Xbeta(i)));
    }
    hessian = X.t() * diagmat(yhat % (1.0 - yhat)) * X;
    beta_new = beta_init + inv(hessian) * X.t() * (y - yhat);
    
    double d1 = 0.0;
    double d0 = 0.0;
    for(int j = 0; j < p; j++){
      d1 += pow(beta_new(j) - beta_init(j), 2.0);
      d0 += pow(beta_init(j), 2.0);
    }
    eps = d1 / d0;
    beta_init = beta_new;
    niter += 1;
  }

  return(List::create(
    Named("coefs") = beta_new,
    Named("hessian") = hessian,
    Named("yhat") = yhat,
    Named("niter") = niter
  ));
}